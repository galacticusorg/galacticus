"""Fit parameters which change only the likelihood of a model's analyses, not its merger trees, by direct runs.

Some parameters of a calibration enter only through the analyses of the model - for example, the parameters of the model
of N-body halo mass errors (``nbodyHaloMassError``), which change a progenitor mass function only through its convolution
with those errors. At a fixed random number seed, runs which differ only in such parameters build identical merger trees,
so the likelihood is an exact, noise-free function of them. Their effects are often too small for an emulator to resolve
(about a percent per bin, below the emulator's accuracy, yet worth tens in the log-likelihood over many bins), so they are
better fit directly, and then fixed while the remaining parameters are emulated.

The procedure (see the "Emulator-Assisted Calibration" chapter of the user guide) has two stages:

1. A **star** design about a center, at one seed: the center, one step either side along each parameter, and the four
   sign combinations of each pair of parameters (``1 + 2n + 2n(n-1)`` points for ``n`` parameters, enough to fit a full
   quadratic). The direct log-likelihood at each point is fit by a quadratic, which predicts the optimum (within a
   bounded number of steps of the center) and the gain in log-likelihood there.
2. A **line** search from the center toward that optimum, at several seeds. The optimum of a single seed is fit to that
   seed's noise, so the gain is checked at independent seeds: at each point along the line, the log-likelihood of the
   predictions averaged over ``N`` seeds is extrapolated linearly in ``1/N`` to infinitely many seeds - the noise-free
   model which an emulator approximates. (A single run's own noise always lowers its log-likelihood, because the misfit
   is quadratic in it.)

Each stage is run as one campaign (:mod:`Galacticus.Emulation.campaign`) per *component* - a base parameter file and the
analyses whose likelihoods it contributes - over a design written here in the layout of the ``emulatorDesign`` task. The
log-likelihood of a point is the sum over components and analyses of the direct likelihood
(:mod:`Galacticus.Emulation.likelihood`) of the analyses, in the forms used by the emulated likelihood. The fractional root
variance of counts of the ``poissonCounts`` form may be fixed, or maximized at each point (profiled), as for a parameter
of the likelihood sampled in the calibration.

The command-line interface is ``scripts/emulation/emulatorProfile.py``.

Andrew Benson, Claude (28-October-2026).
"""

from __future__ import annotations

import itertools
import json
import os
import re
from dataclasses import dataclass, field, replace

import h5py
import numpy as np
from scipy.optimize import minimize

from Galacticus.Emulation import emulatorFile as ef
from Galacticus.Emulation.campaign import Campaign, output_file_path
from Galacticus.Emulation.design import DesignFile, write_design
from Galacticus.Emulation.likelihood import FORMS, LOG_IMPROBABLE, log_likelihood, profile_root_variance
from Galacticus.Emulation.train import Observable, collect

__all__ = [
    'STAGES',
    'ProfileParameter',
    'ProfileObservable',
    'ProfileComponent',
    'Profile',
    'star_offsets',
    'fit_quadratic',
]

STAGES = ('star', 'line')
# The path given, in place of its output, for a run excluded from an average over seeds.
_EXCLUDED = '/nonexistent/excluded-from-average.hdf5'


@dataclass
class ProfileParameter:
    """A parameter to profile: its path in the parameter file, the center of the profile, and the size of one step.

    With ``scale`` ``log`` the step is in the natural logarithm of the parameter (which must then be positive), otherwise
    it is linear.
    """

    name: str
    center: float
    step: float
    scale: str = 'linear'

    def validate(self):
        if self.scale not in ('linear', 'log'):
            raise ValueError(f"parameter '{self.name}': unknown scale '{self.scale}'")
        if not self.step > 0.0:
            raise ValueError(f"parameter '{self.name}': the step must be positive")
        if self.scale == 'log' and not self.center > 0.0:
            raise ValueError(f"parameter '{self.name}': a log-scaled parameter must have a positive center")

    def value(self, offset):
        """The value of the parameter ``offset`` steps from the center."""
        return self.center * np.exp(offset * self.step) if self.scale == 'log' else self.center + offset * self.step


@dataclass
class ProfileObservable:
    """An analysis contributing to the likelihood, with the form of its likelihood.

    ``observable`` gives the label and transform exactly as for training an emulator. If ``pattern`` is true, its label is
    a regular expression, matched in full against the labels of the analyses in the component's output.
    """

    observable: Observable
    form: str
    pattern: bool = False

    def validate(self):
        if self.form not in FORMS:
            raise ValueError(f"observable '{self.observable.label}': unknown likelihood form '{self.form}'")
        if not self.pattern:
            self.observable.validate()


@dataclass
class ProfileComponent:
    """A base parameter file, and the analyses of its output which contribute to the likelihood."""

    name: str
    base_parameters: str
    observables: list = field(default_factory=list)


def star_offsets(count):
    """The offsets (in steps) of the points of a star design in ``count`` parameters: the center, ±1 along each axis, and
    the four sign combinations of each pair of axes."""
    offsets = [np.zeros(count)]
    for k in range(count):
        for sign in (+1.0, -1.0):
            offset = np.zeros(count)
            offset[k] = sign
            offsets.append(offset)
    for k, l in itertools.combinations(range(count), 2):
        for signK, signL in itertools.product((+1.0, -1.0), repeat=2):
            offset = np.zeros(count)
            offset[k], offset[l] = signK, signL
            offsets.append(offset)
    return np.array(offsets)


def _features(offset):
    count = offset.size
    return np.concatenate([[1.0], offset, [offset[k] * offset[l] * (0.5 if k == l else 1.0) for k in range(count) for l in range(k, count)]])


def fit_quadratic(offsets, values, offset_maximum=2.0):
    """Fit ``values = c + g.u + u.H.u/2`` to the log-likelihoods ``values`` at ``offsets`` (in steps), and maximize the fit
    within ``offset_maximum`` steps of the center in each parameter.

    Returns a dictionary of the fit (``constant``, ``gradient``, ``hessian``, and the root-mean-square and maximum
    ``residual``), the ``optimum`` (in steps), the log-likelihood ``predicted`` there, and, if the Hessian is negative
    definite, the ``covariance`` of the parameters (in steps squared) from it.
    """
    offsets = np.asarray(offsets, dtype=float)
    values = np.asarray(values, dtype=float)
    count = offsets.shape[1]
    design = np.array([_features(offset) for offset in offsets])
    if design.shape[0] < design.shape[1]:
        raise ValueError(f'{design.shape[0]} points are too few to fit a quadratic in {count} parameters')
    coefficients, *_ = np.linalg.lstsq(design, values, rcond=None)
    residual = values - design @ coefficients
    gradient = coefficients[1:1 + count]
    hessian = np.zeros((count, count))
    m = 1 + count
    for k in range(count):
        for l in range(k, count):
            hessian[k, l] = hessian[l, k] = coefficients[m]
            m += 1
    # The Hessian need not be negative definite (a profile may be flat and then fall), so maximize within bounds.
    negative = lambda u: -(gradient @ u + 0.5 * u @ hessian @ u)  # noqa: E731
    optimum = minimize(negative, np.zeros(count), jac=lambda u: -(gradient + hessian @ u), method='L-BFGS-B',
                       bounds=[(-offset_maximum, offset_maximum)] * count).x
    result = {
        'constant': float(coefficients[0]),
        'gradient': gradient.tolist(),
        'hessian': hessian.tolist(),
        'residualRMS': float(np.sqrt(np.mean(residual ** 2))),
        'residualMaximum': float(np.max(np.abs(residual))),
        'optimum': optimum.tolist(),
        'predicted': float(coefficients[0] - negative(optimum)),
        'covariance': None,
    }
    if np.all(np.linalg.eigvalsh(hessian) < 0.0):
        result['covariance'] = np.linalg.inv(-hessian).tolist()
    return result


@dataclass
class Profile:
    """A profile of parameters by direct runs: the parameters, the components, and the options of each stage.

    ``directory`` holds every file written (designs, change files, outputs, logs, manifests and results), in a
    subdirectory for each stage. ``root_variance_target_fractional`` is the fractional root variance of counts of the
    ``poissonCounts`` form: a number, or None to maximize the likelihood over it (within
    ``root_variance_target_fractional_bounds``) at each point. If ``model_variance``, the variance of the model at each point
    (the variance of its predictions averaged over seeds, from the analyses' own covariances) is included in the over-dispersion
    of the ``poissonCounts`` form, as the progenitor mass function analysis does with ``[covarianceTargetOnly]`` false (see
    :func:`Galacticus.Emulation.likelihood.log_likelihood`). This reduces the bias of the likelihood of a model built from one,
    or a few, realizations of its merger trees, whose noise otherwise adds to its misfit. Bins whose value was floored (or
    replaced) in collection are given no model variance.
    """

    directory: str
    parameters: list
    components: list
    executable: str | None = None
    seed: int = 1
    line_seeds: list = field(default_factory=lambda: [101, 202, 303, 404])
    line_fractions: list = field(default_factory=lambda: [0.0, 0.5, 1.0, 1.5])
    offset_maximum: float = 2.0
    root_variance_target_fractional: float | None = 0.0
    root_variance_target_fractional_bounds: tuple = (1.0e-3, 3.0)
    model_variance: bool = False

    def validate(self):
        if not self.parameters:
            raise ValueError('at least one parameter must be given')
        for parameter in self.parameters:
            parameter.validate()
        if not self.components:
            raise ValueError('at least one component must be given')
        names = [component.name for component in self.components]
        if len(set(names)) != len(names):
            raise ValueError('component names must be unique')
        for component in self.components:
            if not component.observables:
                raise ValueError(f"component '{component.name}': at least one observable must be given")
            for observable in component.observables:
                observable.validate()
        if len(self.line_seeds) < 2:
            raise ValueError('the line search needs at least two seeds')
        if self.model_variance and any(observable.form != 'poissonCounts' for component in self.components
                                       for observable in component.observables):
            raise ValueError('the variance of the model can be included only for observables of the `poissonCounts` form')

    # Paths.

    def path(self, stage, *names):
        if stage not in STAGES:
            raise ValueError(f"unknown stage '{stage}'")
        return os.path.join(os.path.abspath(self.directory), stage, *names)

    def manifest(self, stage, component):
        return self.path(stage, component.name, 'campaign.json')

    def result_file(self, stage):
        return self.path(stage, 'result.json')

    def campaigns(self, stage):
        """The campaign of each component for a stage (in the order of the components)."""
        return [Campaign.load(self.manifest(stage, component)) for component in self.components]

    # Designs.

    def stage_points(self, stage):
        """The offsets (in steps from the center) of the points of a stage, and the seeds of each point."""
        if stage == 'star':
            return star_offsets(len(self.parameters)), [self.seed]
        star = self.read_result('star')
        if star is None:
            raise RuntimeError('the star stage must be analyzed before the line search is prepared')
        optimum = np.array(star['fit']['optimum'])
        return np.array([fraction * optimum for fraction in self.line_fractions]), list(self.line_seeds)

    def prepare(self, stage, overwrite=False):
        """Write the design, change files and campaign manifest of each component for a stage. Returns the campaigns."""
        offsets, seeds = self.stage_points(stage)
        values = np.array([[parameter.value(offset[k]) for k, parameter in enumerate(self.parameters)] for offset in offsets])
        campaigns = []
        for component in self.components:
            manifest = self.manifest(stage, component)
            if os.path.exists(manifest) and not overwrite:
                raise FileExistsError(f"stage '{stage}' is already prepared ('{manifest}' exists)")
            directory = self.path(stage, component.name)
            os.makedirs(os.path.join(directory, 'changes'), exist_ok=True)
            os.makedirs(os.path.join(directory, 'models'), exist_ok=True)
            runs = {'pointIndex': [], 'realizationIndex': [], 'changeFileName': [], 'outputFileName': [], 'seed': []}
            for point, row in enumerate(values):
                for realization, seed in enumerate(seeds):
                    label = f'point{point:04d}_{realization}'
                    change_file = os.path.join(directory, 'changes', f'{label}.xml')
                    output_file = os.path.join(directory, 'models', f'{label}.hdf5')
                    _write_change_file(change_file, [(parameter.name, value) for parameter, value in zip(self.parameters, row)],
                                       output_file=output_file, seed=seed)
                    for name, value in (('pointIndex', point), ('realizationIndex', realization), ('changeFileName', change_file),
                                        ('outputFileName', output_file), ('seed', seed)):
                        runs[name].append(value)
            design_file = os.path.join(directory, 'design.hdf5')
            write_design(design_file, _design(self.parameters, values, len(seeds), stage, runs), overwrite=overwrite)
            campaigns.append(Campaign.create(manifest, design_file, os.path.abspath(component.base_parameters),
                                             executable=self.executable, working_directory=directory, overwrite=overwrite))
        return campaigns

    # Likelihoods.

    def _observables(self, component, design_file):
        """The observables of a component, with any label patterns expanded against the first output found which has
        analyses (the output of a run which failed part way may have none)."""
        if not any(observable.pattern for observable in component.observables):
            return [(observable.observable, observable.form) for observable in component.observables]
        with h5py.File(design_file, 'r') as file:
            outputs = [name.decode() for name in file['runs']['outputFileName'][()]]
        labels = None
        for output in outputs:
            try:
                with h5py.File(output_file_path(output), 'r') as file:
                    if 'analyses' in file and len(file['analyses']) > 0:
                        labels = sorted(file['analyses'].keys())
                        break
            except OSError:
                continue
        if labels is None:
            raise RuntimeError(f"component '{component.name}': no output with analyses found in which to match observable labels")
        observables = []
        for observable in component.observables:
            if not observable.pattern:
                observables.append((observable.observable, observable.form))
                continue
            matched = [label for label in labels if re.fullmatch(observable.observable.label, label)]
            if not matched:
                raise RuntimeError(f"component '{component.name}': no analysis matches '{observable.observable.label}'")
            for label in matched:
                expanded = replace(observable.observable, label=label)
                expanded.validate()
                observables.append((expanded, observable.form))
        return observables

    def point_likelihoods(self, stage, realizations=None):
        """Return the direct log-likelihood at each point of a stage, and the fractional root variance of counts at which
        it is evaluated (None if fixed).

        The predictions of each point are averaged over its realizations (seeds), or over those listed in
        ``realizations`` only. A point missing from any component (its runs failed) has a log-likelihood of not-a-number.
        """
        terms = None
        for component in self.components:
            design_file = self.path(stage, component.name, 'design.hdf5')
            with h5py.File(design_file, 'r') as file:
                outputs = [name.decode() for name in file['runs']['outputFileName'][()]]
                realization = file['runs']['realizationIndex'][()]
                countPoints = int(file['design'].attrs['countPoints'])
            if terms is None:
                terms = [[] for _ in range(countPoints)]
                present = np.ones(countPoints, dtype=bool)
            if realizations is not None:
                outputs = [output if index in realizations else _EXCLUDED for output, index in zip(outputs, realization)]
            observables = self._observables(component, design_file)
            _, training_sets, reports = collect(design_file, [observable for observable, _ in observables], model_files=outputs)
            for observable, form in observables:
                training_set = training_sets[observable.label]
                missing = np.ones(countPoints, dtype=bool)
                for row, point in enumerate(training_set.point_index):
                    missing[point] = False
                    varianceModel = np.where(training_set.mask[row] != 0, 0.0, training_set.root_variance[row] ** 2) if self.model_variance else None
                    terms[point].append((training_set, form, training_set.y[row], varianceModel))
                present &= ~missing
        rootVariance = self.root_variance_target_fractional
        logLikelihoods, rootVariances = np.full(len(terms), np.nan), np.full(len(terms), np.nan)
        for point, pointTerms in enumerate(terms):
            if not present[point]:
                continue

            def evaluate(C, pointTerms=pointTerms):
                total = 0.0
                for training_set, form, mean, varianceModel in pointTerms:
                    value = log_likelihood(training_set, form, mean, root_variance_target_fractional=C, variance_model=varianceModel)
                    if value == LOG_IMPROBABLE:
                        return LOG_IMPROBABLE
                    total += value
                return total

            if rootVariance is None and any(form == 'poissonCounts' for _, form, _, _ in pointTerms):
                logLikelihoods[point], rootVariances[point] = profile_root_variance(evaluate, self.root_variance_target_fractional_bounds)
            else:
                logLikelihoods[point] = evaluate(rootVariance or 0.0)
        return logLikelihoods, (None if rootVariance is not None else rootVariances)

    # Analysis.

    def analyze(self, stage):
        """Analyze a stage whose runs are complete, write its result file, and return the result."""
        offsets, seeds = self.stage_points(stage)
        result = {'stage': stage, 'parameters': [parameter.name for parameter in self.parameters],
                  'offsets': offsets.tolist(), 'seeds': seeds}
        if stage == 'star':
            logLikelihoods, rootVariances = self.point_likelihoods(stage)
            used = np.isfinite(logLikelihoods) & (logLikelihoods > LOG_IMPROBABLE)
            result['logLikelihood'] = _list(logLikelihoods)
            result['rootVarianceTargetFractional'] = None if rootVariances is None else _list(rootVariances)
            if not used[0]:
                raise RuntimeError('the likelihood at the center of the star design could not be evaluated')
            result['fit'] = fit_quadratic(offsets[used], logLikelihoods[used], self.offset_maximum)
            result['pointsUsed'] = int(used.sum())
            optimum = np.array(result['fit']['optimum'])
            result['optimum'] = self.values(optimum)
            result['gainPredicted'] = result['fit']['predicted'] - float(logLikelihoods[0])
        else:
            # Log-likelihoods of the predictions averaged over N seeds, for N = 1, 2, 4, ..., averaged over the subsets of seeds
            # of that size, and extrapolated linearly in 1/N.
            counts, means = [], []
            count = 1
            while count <= len(seeds):
                values = [self.point_likelihoods(stage, realizations=set(subset))[0]
                          for subset in itertools.combinations(range(len(seeds)), count)]
                counts.append(count)
                means.append(np.mean(values, axis=0))
                count *= 2
            means = np.array(means)
            extrapolated = np.full(offsets.shape[0], np.nan)
            slopes = np.full(offsets.shape[0], np.nan)
            for point in range(offsets.shape[0]):
                if np.all(np.isfinite(means[:, point])):
                    slopes[point], extrapolated[point] = np.polyfit(1.0 / np.array(counts), means[:, point], 1)
            result['fractions'] = list(self.line_fractions)
            result['countsSeeds'] = counts
            result['logLikelihoodMean'] = [_list(row) for row in means]
            result['logLikelihoodExtrapolated'] = _list(extrapolated)
            result['noisePenaltySingleSeed'] = _list(slopes)
            if not np.isfinite(extrapolated).any():
                raise RuntimeError('no point of the line search could be evaluated')
            best = int(np.nanargmax(extrapolated))
            result['best'] = best
            result['optimum'] = self.values(offsets[best])
            result['gain'] = float(extrapolated[best] - extrapolated[0]) if np.isfinite(extrapolated[0]) else None
        with open(self.result_file(stage), 'w') as handle:
            json.dump(result, handle, indent=1)
        return result

    def read_result(self, stage):
        """The result of an analyzed stage, or None."""
        path = self.result_file(stage)
        if not os.path.exists(path):
            return None
        with open(path) as handle:
            return json.load(handle)

    def values(self, offset):
        """The values of the parameters at ``offset`` steps from the center, as a dictionary."""
        return {parameter.name: float(parameter.value(offset[k])) for k, parameter in enumerate(self.parameters)}

    def write_fixed(self, path, values):
        """Write a change file setting the parameters to ``values`` (a dictionary)."""
        _write_change_file(path, list(values.items()))


def _list(array):
    return [None if not np.isfinite(value) else float(value) for value in array]


def _design(parameters, values, realizations, stage, runs):
    """A design (in the layout of the ``emulatorDesign`` task) for the points of a profile. Each parameter is given a
    uniform prior spanning its values at the points, so that its quantiles lie in [0, 1]; they are not used."""
    priors, quantiles = [], np.zeros_like(values)
    for k in range(values.shape[1]):
        lower, upper = float(values[:, k].min()), float(values[:, k].max())
        if upper <= lower:
            upper = lower + 1.0
        priors.append(ef.Prior('uniform', {'limitLower': lower, 'limitUpper': upper}))
        quantiles[:, k] = (values[:, k] - lower) / (upper - lower)
    design = ef.Design(names=[parameter.name for parameter in parameters], priors=priors, mappers=['identity'] * len(parameters),
                       quantiles=np.clip(quantiles, 0.0, 1.0), values=values,
                       attributes={'type': f'profile:{stage}', 'realizationsPerPoint': np.int32(realizations), 'seedPerPoint': np.int32(1)})
    runs = {name: (np.array(value, dtype=np.int32) if name in ('pointIndex', 'realizationIndex') else
                   np.array(value, dtype=np.int64) if name == 'seed' else value) for name, value in runs.items()}
    return DesignFile(design=design, runs=runs)


def _write_change_file(path, changes, output_file=None, seed=None):
    """Write a Galacticus change file setting parameters (a list of (path, value)), and optionally the output file and seed."""
    lines = ['<?xml version="1.0" encoding="UTF-8"?>', '<!-- Written by Galacticus.Emulation.profile. -->', '<changes>']
    for name, value in changes:
        lines.append(f'  <change type="update" path="{_escape(name)}" value="{float(value):.16e}"/>')
    if output_file is not None:
        lines += ['  <change type="replaceOrAppend" path="outputFileName">',
                  f'    <outputFileName value="{_escape(output_file)}"/>',
                  '  </change>']
    if seed is not None:
        lines += ['  <change type="replaceOrAppend" path="randomNumberGenerator/seed">',
                  f'    <seed value="{int(seed)}"/>',
                  '  </change>']
    lines.append('</changes>')
    with open(path, 'w') as handle:
        handle.write('\n'.join(lines) + '\n')


def _escape(text):
    return text.replace('&', '&amp;').replace('"', '&quot;').replace('<', '&lt;').replace('>', '&gt;')
