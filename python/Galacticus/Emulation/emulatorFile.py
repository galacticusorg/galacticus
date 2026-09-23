"""Read, write and evaluate Galacticus emulator files.

An emulator file is an HDF5 file holding everything needed to evaluate, audit and reproduce a trained
emulator of Galacticus model predictions: the design (the parameters varied, their priors, and the
points at which the model was run), the training set for each emulated observable, the trained
emulator for each observable, and (optionally) cross-validation results. The format is specified in
the user guide (``docs/manuals/user-guide/emulation.rst``), and this module is its reference
implementation: :func:`predict` defines exactly how an emulator is evaluated, and the Fortran
``emulatorGaussianProcess`` class must reproduce it.

Conventions:

* Array shapes are given in C (row-major) order, as seen by ``h5py`` and NumPy. Fortran, being
  column-major, sees each array with its dimensions reversed.
* All strings, in attributes and datasets, are fixed-length, null-padded ASCII, which is what
  Galacticus' own HDF5 layer writes and reads.
* Indices stored in the file (``pointIndex``, ``folds``) are zero-based.

Andrew Benson, Claude (2026)
"""

from __future__ import annotations

import datetime
from dataclasses import dataclass, field, fields

import h5py
import numpy as np

__all__ = [
    'FORMAT_NAME',
    'FORMAT_VERSION',
    'KERNELS',
    'Prior',
    'Design',
    'TrainingSet',
    'Component',
    'Emulator',
    'Validation',
    'EmulatorFile',
    'kernel_matern52',
    'predict',
    'read',
    'write',
]

FORMAT_NAME = 'galacticusEmulator'
FORMAT_VERSION = 1
KERNELS = ('matern52ARD',)


def _camel(name):
    """Convert a snake_case field name to the camelCase name used in the file."""
    head, *tail = name.split('_')
    return head + ''.join(part[:1].upper() + part[1:] for part in tail)


def _encode(value):
    """Convert a value for storage, turning strings into fixed-length ASCII."""
    if isinstance(value, str):
        return np.bytes_(value.encode('ascii'))
    if isinstance(value, bool):
        return np.int8(value)
    if isinstance(value, (list, tuple)) and value and all(isinstance(v, str) for v in value):
        return np.array([v.encode('ascii') for v in value], dtype=f'S{max(1, max(len(v) for v in value))}')
    return value


def _decode(value):
    """Convert a stored value back, turning fixed-length strings into ``str``."""
    if isinstance(value, bytes):
        return value.decode('ascii')
    if isinstance(value, np.ndarray) and value.dtype.kind == 'S':
        return [v.decode('ascii') for v in value.tolist()]
    if isinstance(value, np.generic):
        return value.item()
    return value


@dataclass
class Prior:
    """A prior, described by its Galacticus ``distributionFunction1D`` class name and parameters.

    ``parameters`` maps each parameter name of that class (e.g. ``limitLower``, ``x0``, ``sigma``) to
    its value. It must hold every parameter needed to reconstruct the prior exactly, since the
    emulated likelihood checks that the priors of an MCMC agree with those the emulator was trained
    under.
    """

    class_name: str
    parameters: dict = field(default_factory=dict)


@dataclass
class Design:
    """The design: the varied parameters, their priors, and the points at which the model was run.

    ``quantiles`` has shape ``(N, d)`` and holds the prior cumulative probability of each parameter
    at each point; these are the emulator's inputs. ``values`` has the same shape and holds the
    physical parameter values. ``names`` are the parameter paths, exactly as given in each
    ``modelParameter``'s ``name``. ``mappers`` are the ``operatorUnaryMapper`` class names.
    ``attributes`` records how the design was generated (e.g. ``type='sobol'``, ``seed``,
    ``randomShift``).
    """

    names: list
    priors: list
    mappers: list
    quantiles: np.ndarray
    values: np.ndarray
    attributes: dict = field(default_factory=dict)

    def validate(self):
        count_parameters = len(self.names)
        if len(self.priors) != count_parameters or len(self.mappers) != count_parameters:
            raise ValueError('design: names, priors and mappers must have the same length')
        for name in ('quantiles', 'values'):
            array = getattr(self, name)
            if array.ndim != 2 or array.shape[1] != count_parameters:
                raise ValueError(f'design: {name} must have shape (N, {count_parameters})')
        if self.values.shape != self.quantiles.shape:
            raise ValueError('design: quantiles and values must have the same shape')
        if np.any(self.quantiles < 0.0) or np.any(self.quantiles > 1.0):
            raise ValueError('design: quantiles must lie in [0, 1]')


@dataclass
class TrainingSet:
    """The training set for one observable, after transforms, floors and masking.

    ``y`` and ``root_variance`` have shape ``(N, B)``: the transformed model prediction in each of B
    bins at each of N used design points, and its finite-sampling uncertainty. ``point_index`` (N)
    gives the design point (zero-based) of each row, so points which failed or were excluded are
    simply absent. ``mask`` (N, B) is nonzero where a value was floored or replaced. The raw values
    (``y_raw``, ``covariance_raw``) are optional. ``attributes`` records the transform
    (``transform``, ``floor``, ``rootVarianceFloored``, ``undefined``, ``rootVarianceUndefined``,
    ``xAxisIsLog``, ``yAxisIsLog``, ...).
    """

    x: np.ndarray
    y: np.ndarray
    root_variance: np.ndarray
    mask: np.ndarray
    point_index: np.ndarray
    y_target: np.ndarray
    covariance_target: np.ndarray
    y_raw: np.ndarray | None = None
    covariance_raw: np.ndarray | None = None
    attributes: dict = field(default_factory=dict)

    def validate(self, count_points):
        count_bins = self.x.shape[0]
        count_used = self.point_index.shape[0]
        for name in ('y', 'root_variance', 'mask'):
            if getattr(self, name).shape != (count_used, count_bins):
                raise ValueError(f'training set: {name} must have shape ({count_used}, {count_bins})')
        if self.y_target.shape != (count_bins,) or self.covariance_target.shape != (count_bins, count_bins):
            raise ValueError('training set: target shapes do not match the number of bins')
        if np.any(self.point_index < 0) or np.any(self.point_index >= count_points):
            raise ValueError('training set: point index out of range of the design')
        if len(np.unique(self.point_index)) != count_used:
            raise ValueError('training set: point indices must be unique')
        if self.y_raw is not None and self.y_raw.shape != (count_used, count_bins):
            raise ValueError('training set: yRaw has the wrong shape')
        if self.covariance_raw is not None and self.covariance_raw.shape != (count_used, count_bins, count_bins):
            raise ValueError('training set: covarianceRaw has the wrong shape')


@dataclass
class Component:
    """One Gaussian-process component: the GP emulating one (standardized) PCA coefficient.

    ``targets`` (N) are the standardized coefficient values at the training inputs, ``noise_variance``
    (N) their variance (the heteroscedastic noise added to the diagonal of the covariance matrix),
    and ``alpha`` (N) the solution of ``K alpha = targets``. ``cholesky_factor`` (N, N), the lower
    triangular factor of K, is optional: it is fully determined by the other data, and costs 8 N^2
    bytes, so a reader may recompute it instead.
    """

    log_amplitude: float
    log_length_scales: np.ndarray
    targets: np.ndarray
    noise_variance: np.ndarray
    alpha: np.ndarray
    log_marginal_likelihood: float = float('nan')
    cholesky_factor: np.ndarray | None = None


@dataclass
class Emulator:
    """The emulator for one observable.

    ``input_names`` are the design parameters the emulator uses, in its own input order (a subset of,
    or reordering of, the design's ``names``). ``inputs`` (N, d_obs) are the training inputs, i.e. the
    design quantiles of those parameters at the training set's points. ``bin_mean`` and ``bin_scale``
    (B) standardize the bins; ``pca_components`` (K, B) are the retained principal components; and
    ``coefficient_mean`` and ``coefficient_scale`` (K) standardize the coefficients.
    """

    kernel: str
    jitter: float
    input_names: list
    inputs: np.ndarray
    bin_mean: np.ndarray
    bin_scale: np.ndarray
    pca_components: np.ndarray
    coefficient_mean: np.ndarray
    coefficient_scale: np.ndarray
    components: list
    pca_variance_retained: float = float('nan')

    def validate(self, count_bins):
        if self.kernel not in KERNELS:
            raise ValueError(f"emulator: unknown kernel '{self.kernel}'")
        count_inputs = len(self.input_names)
        count_points = self.inputs.shape[0]
        count_components = len(self.components)
        if self.inputs.shape != (count_points, count_inputs):
            raise ValueError('emulator: inputs have the wrong shape')
        if self.bin_mean.shape != (count_bins,) or self.bin_scale.shape != (count_bins,):
            raise ValueError('emulator: bin standardization has the wrong shape')
        if self.pca_components.shape != (count_components, count_bins):
            raise ValueError('emulator: PCA components have the wrong shape')
        if self.coefficient_mean.shape != (count_components,) or self.coefficient_scale.shape != (count_components,):
            raise ValueError('emulator: coefficient standardization has the wrong shape')
        for component in self.components:
            if component.log_length_scales.shape != (count_inputs,):
                raise ValueError('emulator: length scales have the wrong shape')
            for name in ('targets', 'noise_variance', 'alpha'):
                if getattr(component, name).shape != (count_points,):
                    raise ValueError(f'emulator: component {name} has the wrong shape')
            if component.cholesky_factor is not None and component.cholesky_factor.shape != (count_points, count_points):
                raise ValueError('emulator: Cholesky factor has the wrong shape')


@dataclass
class Validation:
    """Cross-validation results for one observable (optional).

    ``folds`` (N) gives the fold (zero-based) in which each training row was held out;
    ``held_out_prediction`` and ``held_out_variance`` (N, B) are the predictions for each row when it
    was held out. The remaining arrays (B) are per-bin summary statistics.
    """

    folds: np.ndarray
    held_out_prediction: np.ndarray
    held_out_variance: np.ndarray
    rmse: np.ndarray
    r2: np.ndarray
    rms_standardized_residual: np.ndarray
    coverage_1_sigma: np.ndarray
    coverage_2_sigma: np.ndarray


@dataclass
class EmulatorFile:
    """The complete content of an emulator file.

    ``training_sets``, ``emulators`` and ``validation`` are keyed by observable label (normally the
    ``outputAnalysis`` label). Every emulator must have a training set of the same label.
    ``attributes`` holds provenance (``created``, ``creator``, ``gitHash``, ...).
    """

    design: Design
    training_sets: dict
    emulators: dict
    validation: dict = field(default_factory=dict)
    attributes: dict = field(default_factory=dict)

    def validate(self):
        self.design.validate()
        count_points = self.design.quantiles.shape[0]
        for label, training_set in self.training_sets.items():
            try:
                training_set.validate(count_points)
            except ValueError as error:
                raise ValueError(f"'{label}': {error}") from None
        for label, emulator in self.emulators.items():
            if label not in self.training_sets:
                raise ValueError(f"emulator '{label}' has no training set")
            missing = [name for name in emulator.input_names if name not in self.design.names]
            if missing:
                raise ValueError(f"emulator '{label}' uses inputs not in the design: {missing}")
            try:
                emulator.validate(self.training_sets[label].x.shape[0])
            except ValueError as error:
                raise ValueError(f"'{label}': {error}") from None
        for label in self.validation:
            if label not in self.training_sets:
                raise ValueError(f"validation '{label}' has no training set")


def kernel_matern52(inputs1, inputs2, log_amplitude, log_length_scales):
    """Evaluate the Matern-5/2 ARD kernel.

    k(u, u') = A (1 + sqrt(5) r + 5 r^2 / 3) exp(-sqrt(5) r), with r^2 = sum_i ((u_i - u'_i) / l_i)^2,
    A = exp(log_amplitude) and l_i = exp(log_length_scales_i). Returns an array of shape
    (len(inputs1), len(inputs2)).
    """
    length_scales = np.exp(log_length_scales)
    difference = (np.atleast_2d(inputs1)[:, None, :] - np.atleast_2d(inputs2)[None, :, :]) / length_scales
    radius = np.sqrt(np.sum(difference**2, axis=-1))
    return np.exp(log_amplitude) * (1.0 + np.sqrt(5.0) * radius + 5.0 * radius**2 / 3.0) * np.exp(-np.sqrt(5.0) * radius)


def _cholesky_factor(emulator, component):
    """Return the lower Cholesky factor of a component's covariance matrix, computing it if not stored."""
    if component.cholesky_factor is not None:
        return component.cholesky_factor
    covariance = kernel_matern52(emulator.inputs, emulator.inputs, component.log_amplitude, component.log_length_scales)
    covariance[np.diag_indices_from(covariance)] += component.noise_variance + emulator.jitter
    return np.linalg.cholesky(covariance)


def predict(emulator, quantiles):
    """Evaluate an emulator, returning its predicted mean and variance in each bin.

    ``quantiles`` is a vector of the emulator's inputs (prior quantiles, in the order of
    ``emulator.input_names``), or an array of shape (M, d_obs) of such vectors. For each component k
    the GP gives a mean m_k = k_*^T alpha_k and variance v_k = k_** - |L_k^-1 k_*|^2, in standardized
    coefficient space. These are unstandardized to c_k = mean_k + scale_k m_k and
    sigma_k^2 = scale_k^2 v_k, and projected onto the bins:

        y_b = binMean_b + binScale_b sum_k c_k E_kb
        sigma_b^2 = binScale_b^2 sum_k E_kb^2 sigma_k^2

    The variance ignores covariance between components and PCA truncation error. Returns arrays of
    shape (B,) for a single input vector, or (M, B) otherwise.
    """
    quantiles = np.asarray(quantiles, dtype=float)
    single = quantiles.ndim == 1
    quantiles = np.atleast_2d(quantiles)
    count_components = len(emulator.components)
    means = np.empty((quantiles.shape[0], count_components))
    variances = np.empty((quantiles.shape[0], count_components))
    for k, component in enumerate(emulator.components):
        cross = kernel_matern52(quantiles, emulator.inputs, component.log_amplitude, component.log_length_scales)
        factor = _cholesky_factor(emulator, component)
        solved = np.linalg.solve(factor, cross.T)
        means[:, k] = cross @ component.alpha
        variances[:, k] = np.maximum(np.exp(component.log_amplitude) - np.sum(solved**2, axis=0), 0.0)
    coefficients = emulator.coefficient_mean + emulator.coefficient_scale * means
    coefficient_variances = emulator.coefficient_scale**2 * variances
    mean = emulator.bin_mean + emulator.bin_scale * (coefficients @ emulator.pca_components)
    variance = emulator.bin_scale**2 * (coefficient_variances @ emulator.pca_components**2)
    if single:
        return mean[0], variance[0]
    return mean, variance


def _write_attributes(obj, attributes):
    for name, value in attributes.items():
        if value is None:
            continue
        obj.attrs[name] = _encode(value)


def _read_attributes(obj):
    return {name: _decode(value) for name, value in obj.attrs.items()}


def _write_datasets(group, record, exclude=()):
    for item in fields(record):
        if item.name in exclude:
            continue
        value = getattr(record, item.name)
        if value is None:
            continue
        group.create_dataset(_camel(item.name), data=np.asarray(value))


def _read_datasets(group, cls, exclude=()):
    values = {}
    for item in fields(cls):
        if item.name in exclude:
            continue
        name = _camel(item.name)
        if name in group:
            values[item.name] = group[name][()]
    return values


def write(path, content, overwrite=False):
    """Write an :class:`EmulatorFile` to ``path``, after validating it."""
    content.validate()
    with h5py.File(path, 'w' if overwrite else 'w-') as file:
        attributes = {'created': datetime.datetime.now(datetime.timezone.utc).isoformat(timespec='seconds')}
        attributes.update(content.attributes)
        attributes.update({'format': FORMAT_NAME, 'formatVersion': FORMAT_VERSION})
        _write_attributes(file, attributes)
        # Design.
        design = content.design
        group = file.create_group('design')
        _write_attributes(group, design.attributes)
        group.create_dataset('parameterNames', data=_encode(list(design.names)))
        group.create_dataset('mappers', data=_encode(list(design.mappers)))
        group.create_dataset('quantiles', data=design.quantiles)
        group.create_dataset('values', data=design.values)
        priors = group.create_group('priors')
        for i, prior in enumerate(design.priors):
            prior_group = priors.create_group(f'prior{i + 1}')
            _write_attributes(prior_group, {'class': prior.class_name, **prior.parameters})
        # Training sets.
        group = file.create_group('trainingSets')
        for label, training_set in content.training_sets.items():
            subgroup = group.create_group(label)
            _write_attributes(subgroup, training_set.attributes)
            _write_datasets(subgroup, training_set, exclude=('attributes',))
        # Emulators.
        group = file.create_group('emulators')
        for label, emulator in content.emulators.items():
            subgroup = group.create_group(label)
            _write_attributes(subgroup, {
                'kernel': emulator.kernel,
                'jitter': emulator.jitter,
                'pcaVarianceRetained': emulator.pca_variance_retained,
                'countComponents': len(emulator.components),
            })
            subgroup.create_dataset('inputNames', data=_encode(list(emulator.input_names)))
            _write_datasets(subgroup, emulator, exclude=('kernel', 'jitter', 'input_names', 'components', 'pca_variance_retained'))
            for k, component in enumerate(emulator.components):
                component_group = subgroup.create_group(f'component{k + 1}')
                _write_attributes(component_group, {
                    'logAmplitude': component.log_amplitude,
                    'logMarginalLikelihood': component.log_marginal_likelihood,
                })
                _write_datasets(component_group, component, exclude=('log_amplitude', 'log_marginal_likelihood'))
        # Validation.
        if content.validation:
            group = file.create_group('validation')
            for label, validation in content.validation.items():
                _write_datasets(group.create_group(label), validation)


def read(path):
    """Read an emulator file, returning an :class:`EmulatorFile`."""
    with h5py.File(path, 'r') as file:
        attributes = _read_attributes(file)
        if attributes.get('format') != FORMAT_NAME:
            raise ValueError(f"'{path}' is not a Galacticus emulator file")
        if attributes.get('formatVersion') != FORMAT_VERSION:
            raise ValueError(f"'{path}' has format version {attributes.get('formatVersion')}; this reader supports version {FORMAT_VERSION}")
        # Design.
        group = file['design']
        count_parameters = len(group['parameterNames'])
        priors = []
        for i in range(count_parameters):
            prior_attributes = _read_attributes(group['priors'][f'prior{i + 1}'])
            class_name = prior_attributes.pop('class')
            priors.append(Prior(class_name, prior_attributes))
        design = Design(
            names=_decode(group['parameterNames'][()]),
            priors=priors,
            mappers=_decode(group['mappers'][()]),
            quantiles=group['quantiles'][()],
            values=group['values'][()],
            attributes=_read_attributes(group),
        )
        # Training sets.
        training_sets = {}
        for label, subgroup in file['trainingSets'].items():
            training_sets[label] = TrainingSet(attributes=_read_attributes(subgroup), **_read_datasets(subgroup, TrainingSet, exclude=('attributes',)))
        # Emulators.
        emulators = {}
        for label, subgroup in file['emulators'].items():
            emulator_attributes = _read_attributes(subgroup)
            components = []
            for k in range(emulator_attributes['countComponents']):
                component_group = subgroup[f'component{k + 1}']
                component_attributes = _read_attributes(component_group)
                components.append(Component(
                    log_amplitude=component_attributes['logAmplitude'],
                    log_marginal_likelihood=component_attributes['logMarginalLikelihood'],
                    **_read_datasets(component_group, Component, exclude=('log_amplitude', 'log_marginal_likelihood')),
                ))
            emulators[label] = Emulator(
                kernel=emulator_attributes['kernel'],
                jitter=emulator_attributes['jitter'],
                pca_variance_retained=emulator_attributes['pcaVarianceRetained'],
                input_names=_decode(subgroup['inputNames'][()]),
                components=components,
                **_read_datasets(subgroup, Emulator, exclude=('kernel', 'jitter', 'input_names', 'components', 'pca_variance_retained')),
            )
        # Validation.
        validation = {}
        if 'validation' in file:
            for label, subgroup in file['validation'].items():
                validation[label] = Validation(**_read_datasets(subgroup, Validation))
    content = EmulatorFile(design=design, training_sets=training_sets, emulators=emulators, validation=validation, attributes=attributes)
    content.validate()
    return content
