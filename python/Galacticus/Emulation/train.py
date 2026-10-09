"""Collect training sets from a campaign of Galacticus models, and train emulators of them.

This is the reference trainer for emulator-assisted calibration. Given a design file (written by the ``emulatorDesign``
task) and the outputs of the runs of a campaign over it, it:

1. *collects* a training set for each observable from the ``analyses/<label>`` group of each run's output---the
   standard group written by ``outputAnalysis`` classes, with attributes ``xDataset``, ``yDataset``, ``yCovariance``,
   ``yDatasetTarget`` and ``yCovarianceTarget`` naming its datasets---averaging repeated realizations of a design
   point, and excluding runs which failed;
2. *transforms* each observable (optionally to log10, with a floor for empty bins, and replacing bins which are
   undefined in some runs);
3. *trains* an emulator of each observable: the bins are standardized, compressed by principal components analysis
   (PCA), and one Gaussian process (GP) is fitted to each (standardized) principal component coefficient, using
   scikit-learn, with a Matérn-5/2 kernel with a separate length scale for each input, and the finite-sampling noise of
   each run as heteroscedastic noise;
4. optionally *cross-validates* each emulator; and
5. writes everything to an emulator file (see :mod:`Galacticus.Emulation.emulatorFile`).

The emulator's inputs are the prior quantiles of the design parameters. scikit-learn is needed only for training
(:func:`train`); collection does not need it.

Andrew Benson, Claude (2026)
"""

from __future__ import annotations

import os
import warnings
from dataclasses import dataclass

import h5py
import numpy as np

from Galacticus.Emulation import emulatorFile as ef
from Galacticus.Emulation.design import read_design

__all__ = [
    'Observable',
    'CollectionReport',
    'collect',
    'fit_emulator',
    'cross_validate',
    'train',
]

# Bounds on the hyperparameters of the Gaussian processes (amplitude, and length scales in units of prior quantiles).
AMPLITUDE_BOUNDS = (1.0e-3, 1.0e3)
LENGTH_SCALE_BOUNDS = (1.0e-2, 1.0e1)
JITTER = 1.0e-10


@dataclass
class Observable:
    """How to collect and transform one observable.

    ``label`` is that of the ``analyses/<label>`` group. ``transform`` is ``identity`` or ``log10``. For ``log10``, values
    at or below ``10**floor`` (including empty bins) are set to ``floor``, with root variance ``rootVarianceFloored``
    (in dex). If ``undefined`` is ``median``, a bin whose value is undefined in a run (non-finite, or zero with zero
    variance, as for a mean relation with no galaxies in that bin) is set to the median of the bin over the runs in which
    it is defined, with root variance ``rootVarianceUndefined``.
    """

    label: str
    transform: str = 'identity'
    floor: float | None = None
    root_variance_floored: float | None = None
    undefined: str | None = None
    root_variance_undefined: float | None = None

    def validate(self):
        if self.transform not in ('identity', 'log10'):
            raise ValueError(f"observable '{self.label}': unknown transform '{self.transform}'")
        if self.transform == 'log10' and (self.floor is None) != (self.root_variance_floored is None):
            raise ValueError(f"observable '{self.label}': `floor` and `rootVarianceFloored` must be given together")
        if self.undefined not in (None, 'median'):
            raise ValueError(f"observable '{self.label}': unknown treatment of undefined bins '{self.undefined}'")
        if self.undefined == 'median' and self.root_variance_undefined is None:
            raise ValueError(f"observable '{self.label}': `rootVarianceUndefined` is required with `undefined=median`")

    def attributes(self):
        """The attributes recording this transform in a training set."""
        attributes = {'transform': self.transform}
        if self.floor is not None:
            attributes.update({'floor': self.floor, 'rootVarianceFloored': self.root_variance_floored})
        if self.undefined is not None:
            attributes.update({'undefined': self.undefined, 'rootVarianceUndefined': self.root_variance_undefined})
        return attributes


@dataclass
class CollectionReport:
    """A summary of the collection of a training set: which design points were used, and which runs failed."""

    points_used: list
    points_missing: list
    runs_failed: dict

    def __str__(self):
        text = f'{len(self.points_used)} design points used, {len(self.points_missing)} missing'
        if self.points_missing:
            text += ' (points ' + ', '.join(str(point) for point in self.points_missing) + ')'
        return text


def _read_analysis(path, label):
    """Read the ``analyses/<label>`` group of a model output: (x, y, covariance, yTarget, covarianceTarget, attributes)."""
    with h5py.File(path, 'r') as file:
        if 'analyses' not in file or label not in file['analyses']:
            raise KeyError(f"no analysis '{label}' in '{path}'")
        group = file['analyses'][label]
        attributes = {name: (value.decode() if isinstance(value, bytes) else value) for name, value in group.attrs.items()}
        if attributes.get('type') != 'function1D':
            raise ValueError(f"analysis '{label}' in '{path}' is of type '{attributes.get('type')}', not function1D")

        def dataset(name):
            key = attributes.get(name)
            return group[key][()] if key is not None and key in group else None

        return (dataset('xDataset'), dataset('yDataset'), dataset('yCovariance'), dataset('yDatasetTarget'),
                dataset('yCovarianceTarget'), attributes)


def _transform(observable, y, variance, mask):
    """Transform values and variances (each of shape (N, B)) for one observable, updating the mask in place."""
    rootVariance = np.sqrt(np.maximum(variance, 0.0))
    if observable.transform == 'log10':
        floor = observable.floor if observable.floor is not None else -np.inf
        with np.errstate(divide='ignore', invalid='ignore'):
            positive = y > 10.0 ** floor
            yTransformed = np.where(positive, np.log10(np.where(positive, y, 1.0)), floor)
            rootVarianceTransformed = np.where(positive, rootVariance / (np.where(positive, y, 1.0) * np.log(10.0)),
                                               observable.root_variance_floored if observable.floor is not None else np.inf)
        mask |= ~positive
        return yTransformed, rootVarianceTransformed
    return y.copy(), rootVariance


def _transform_target(observable, yTarget, covarianceTarget):
    """Transform the target data and its covariance (to first order, for a log10 transform)."""
    if yTarget is None:
        return None, None
    if observable.transform == 'log10':
        with np.errstate(divide='ignore', invalid='ignore'):
            jacobian = np.where(yTarget > 0.0, 1.0 / (yTarget * np.log(10.0)), 0.0)
            yTransformed = np.where(yTarget > 0.0, np.log10(np.where(yTarget > 0.0, yTarget, 1.0)), -np.inf)
        covariance = None if covarianceTarget is None else covarianceTarget * np.outer(jacobian, jacobian)
        return yTransformed, covariance
    return yTarget, covarianceTarget


def collect(design_file, observables, model_files=None, directory='.'):
    """Collect a training set for each observable from the outputs of a campaign.

    ``model_files`` lists the output file of each run of the design (in run order); by default these are the output files
    named in the design file, resolved relative to ``directory``. A run whose output is missing or lacks an observable is
    excluded from that observable's training set. Repeated realizations of a design point are averaged, and the variance
    of the average is that of a single realization divided by the number averaged.

    Returns ``(design, training_sets, reports)``: the :class:`Galacticus.Emulation.emulatorFile.Design`, and dictionaries,
    keyed by label, of :class:`Galacticus.Emulation.emulatorFile.TrainingSet` and :class:`CollectionReport`.
    """
    for observable in observables:
        observable.validate()
    design_file_content = read_design(design_file)
    design = design_file_content.design
    runs = design_file_content.runs
    countPoints = design.quantiles.shape[0]
    if model_files is None:
        model_files = [os.path.join(directory, name) for name in runs['outputFileName']]
    if len(model_files) != design_file_content.count_runs:
        raise ValueError(f'{len(model_files)} model files given for {design_file_content.count_runs} runs')
    training_sets = {}
    reports = {}
    for observable in observables:
        sums = {}
        failed = {}
        x = yTarget = covarianceTarget = None
        for index, path in enumerate(model_files):
            point = int(runs['pointIndex'][index])
            try:
                xRun, y, covariance, yTargetRun, covarianceTargetRun, _ = _read_analysis(path, observable.label)
            except (OSError, KeyError, ValueError) as error:
                failed[index] = str(error)
                continue
            if x is None:
                x, yTarget, covarianceTarget = xRun, yTargetRun, covarianceTargetRun
            elif xRun.shape != x.shape or not np.allclose(xRun, x, rtol=1.0e-10, atol=0.0):
                raise ValueError(f"observable '{observable.label}': bins of run {index} differ from those of other runs")
            elif yTargetRun is not None and not np.allclose(yTargetRun, yTarget, rtol=1.0e-10, atol=0.0):
                raise ValueError(f"observable '{observable.label}': target data of run {index} differ from those of other runs")
            variance = np.diag(covariance) if covariance is not None else np.zeros_like(y)
            entry = sums.setdefault(point, [np.zeros_like(y), np.zeros_like(y), 0, np.zeros_like(y, dtype=bool)])
            undefined = ~np.isfinite(y) | ((y == 0.0) & (variance == 0.0)) if observable.undefined is not None else np.zeros_like(y, dtype=bool)
            entry[0] += np.where(undefined, 0.0, y)
            entry[1] += np.where(undefined, 0.0, variance)
            entry[2] += 1
            entry[3] |= undefined
        if x is None:
            raise ValueError(f"observable '{observable.label}': not found in any run")
        points = sorted(sums)
        yRaw = np.array([sums[point][0] / sums[point][2] for point in points])
        varianceRaw = np.array([sums[point][1] / sums[point][2] ** 2 for point in points])
        undefinedBins = np.array([sums[point][3] for point in points])
        mask = np.zeros_like(yRaw, dtype=bool)
        y, rootVariance = _transform(observable, yRaw, varianceRaw, mask)
        # Replace bins undefined in a run by the median of that bin over the runs in which it is defined.
        if observable.undefined == 'median' and undefinedBins.any():
            for bin_ in range(y.shape[1]):
                undefined = undefinedBins[:, bin_]
                if not undefined.any():
                    continue
                defined = ~undefined
                if not defined.any():
                    raise ValueError(f"observable '{observable.label}': bin {bin_} is undefined in every run")
                y[undefined, bin_] = np.median(y[defined, bin_])
                rootVariance[undefined, bin_] = observable.root_variance_undefined
                mask[undefined, bin_] = True
        yTargetTransformed, covarianceTargetTransformed = _transform_target(observable, yTarget, covarianceTarget)
        if yTargetTransformed is None:
            yTargetTransformed = np.full(x.shape, np.nan)
            covarianceTargetTransformed = np.full((x.size, x.size), np.nan)
        elif covarianceTargetTransformed is None:
            covarianceTargetTransformed = np.full((x.size, x.size), np.nan)
        training_sets[observable.label] = ef.TrainingSet(
            x=np.asarray(x, dtype=float),
            y=y,
            root_variance=rootVariance,
            mask=mask.astype(np.int8),
            point_index=np.array(points, dtype=np.int64),
            y_target=np.asarray(yTargetTransformed, dtype=float),
            covariance_target=np.asarray(covarianceTargetTransformed, dtype=float),
            y_raw=yRaw,
            attributes=observable.attributes(),
        )
        reports[observable.label] = CollectionReport(
            points_used=points,
            points_missing=sorted(set(range(countPoints)) - set(points)),
            runs_failed=failed,
        )
    return design, training_sets, reports


def fit_emulator(inputs, input_names, y, root_variance, pca_variance_retained=0.99, restarts=4, seed=0):
    """Fit an emulator to training data.

    ``inputs`` (N, d) are prior quantiles; ``y`` and ``root_variance`` (N, B) are the (transformed) training values and
    their uncertainties. Returns a :class:`Galacticus.Emulation.emulatorFile.Emulator`.
    """
    from sklearn.exceptions import ConvergenceWarning
    from sklearn.gaussian_process import GaussianProcessRegressor
    from sklearn.gaussian_process.kernels import ConstantKernel, Matern

    countPoints, countInputs = inputs.shape
    # Standardize the bins. A bin which does not vary is given unit scale.
    binMean = y.mean(axis=0)
    binScale = y.std(axis=0)
    binScale = np.where(binScale > 0.0, binScale, 1.0)
    standardized = (y - binMean) / binScale
    # Principal components analysis: keep the fewest components which retain the requested fraction of the variance.
    _, singularValues, components = np.linalg.svd(standardized, full_matrices=False)
    varianceFraction = singularValues ** 2 / np.sum(singularValues ** 2) if np.sum(singularValues ** 2) > 0.0 else np.ones_like(singularValues)
    countComponents = int(np.searchsorted(np.cumsum(varianceFraction), pca_variance_retained - 1.0e-12) + 1)
    countComponents = min(countComponents, components.shape[0])
    components = components[:countComponents]
    coefficients = standardized @ components.T
    # Propagate the noise in each bin to the coefficients (ignoring covariance between bins).
    coefficientVariance = (root_variance / binScale) ** 2 @ (components ** 2).T
    coefficientMean = coefficients.mean(axis=0)
    coefficientScale = coefficients.std(axis=0)
    coefficientScale = np.where(coefficientScale > 0.0, coefficientScale, 1.0)
    gpComponents = []
    for k in range(countComponents):
        targets = (coefficients[:, k] - coefficientMean[k]) / coefficientScale[k]
        noiseVariance = coefficientVariance[:, k] / coefficientScale[k] ** 2
        kernel = (ConstantKernel(1.0, AMPLITUDE_BOUNDS)
                  * Matern(length_scale=np.ones(countInputs), length_scale_bounds=LENGTH_SCALE_BOUNDS, nu=2.5))
        regressor = GaussianProcessRegressor(kernel=kernel, alpha=noiseVariance + JITTER, normalize_y=False,
                                             n_restarts_optimizer=restarts, random_state=seed + k)
        with warnings.catch_warnings():
            # A hyperparameter at a bound is reported, not fatal; the fit is checked by cross-validation.
            warnings.simplefilter('ignore', ConvergenceWarning)
            regressor.fit(inputs, targets)
        amplitude = regressor.kernel_.k1.constant_value
        lengthScales = np.atleast_1d(regressor.kernel_.k2.length_scale) * np.ones(countInputs)
        gpComponents.append(ef.Component(
            log_amplitude=float(np.log(amplitude)),
            log_length_scales=np.log(lengthScales),
            targets=targets,
            noise_variance=noiseVariance,
            alpha=np.asarray(regressor.alpha_).ravel(),
            log_marginal_likelihood=float(regressor.log_marginal_likelihood_value_),
            cholesky_factor=np.asarray(regressor.L_),
        ))
    return ef.Emulator(
        kernel='matern52ARD',
        jitter=JITTER,
        input_names=list(input_names),
        inputs=inputs,
        bin_mean=binMean,
        bin_scale=binScale,
        pca_components=components,
        coefficient_mean=coefficientMean,
        coefficient_scale=coefficientScale,
        components=gpComponents,
        pca_variance_retained=float(pca_variance_retained),
    )


def cross_validate(inputs, input_names, y, root_variance, folds=5, seed=0, **fit_options):
    """k-fold cross-validation of an emulator: refit on each training fold and predict the held-out fold.

    Returns a :class:`Galacticus.Emulation.emulatorFile.Validation`. Standardized residuals and coverage use the combined
    variance of the emulator and of the finite-sampling noise of the held-out runs.
    """
    countPoints = inputs.shape[0]
    if folds < 2 or folds > countPoints:
        raise ValueError(f'the number of folds must be between 2 and the number of training points ({countPoints})')
    order = np.random.default_rng(seed).permutation(countPoints)
    foldOf = np.empty(countPoints, dtype=np.int64)
    for fold, indices in enumerate(np.array_split(order, folds)):
        foldOf[indices] = fold
    prediction = np.empty_like(y)
    variance = np.empty_like(y)
    for fold in range(folds):
        heldOut = foldOf == fold
        emulator = fit_emulator(inputs[~heldOut], input_names, y[~heldOut], root_variance[~heldOut], seed=seed, **fit_options)
        prediction[heldOut], variance[heldOut] = ef.predict(emulator, inputs[heldOut])
    residual = prediction - y
    varianceTotal = variance + root_variance ** 2
    standardized = residual / np.sqrt(varianceTotal)
    varianceBins = y.var(axis=0)
    with np.errstate(divide='ignore', invalid='ignore'):
        r2 = np.where(varianceBins > 0.0, 1.0 - np.mean(residual ** 2, axis=0) / varianceBins, np.nan)
    return ef.Validation(
        folds=foldOf,
        held_out_prediction=prediction,
        held_out_variance=variance,
        rmse=np.sqrt(np.mean(residual ** 2, axis=0)),
        r2=r2,
        rms_standardized_residual=np.sqrt(np.mean(standardized ** 2, axis=0)),
        coverage_1_sigma=np.mean(np.abs(standardized) <= 1.0, axis=0),
        coverage_2_sigma=np.mean(np.abs(standardized) <= 2.0, axis=0),
    )


def train(design, training_sets, pca_variance_retained=0.99, restarts=4, folds=0, seed=0, store_cholesky=True):
    """Train an emulator for each training set, optionally cross-validating it.

    Returns an :class:`Galacticus.Emulation.emulatorFile.EmulatorFile`, ready to be written. ``folds`` of zero disables
    cross-validation. The Cholesky factors of the Gaussian processes (8N^2 bytes each) are stored only if
    ``store_cholesky`` is true; otherwise a reader recomputes them.
    """
    emulators = {}
    validation = {}
    for label, training_set in training_sets.items():
        inputs = design.quantiles[training_set.point_index]
        emulator = fit_emulator(inputs, design.names, training_set.y, training_set.root_variance,
                                pca_variance_retained=pca_variance_retained, restarts=restarts, seed=seed)
        if not store_cholesky:
            for component in emulator.components:
                component.cholesky_factor = None
        emulators[label] = emulator
        if folds:
            validation[label] = cross_validate(inputs, design.names, training_set.y, training_set.root_variance, folds=folds,
                                               seed=seed, pca_variance_retained=pca_variance_retained, restarts=restarts)
    return ef.EmulatorFile(design=design, training_sets=training_sets, emulators=emulators, validation=validation,
                           attributes={'creator': 'Galacticus.Emulation.train'})
