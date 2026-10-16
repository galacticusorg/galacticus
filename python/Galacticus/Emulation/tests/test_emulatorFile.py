"""Tests for `Galacticus.Emulation.emulatorFile` (the emulator file format)."""

import h5py
import numpy as np
import pytest

from Galacticus.Emulation import emulatorFile as ef

COUNT_POINTS = 16
COUNT_BINS = 3
NOISE_VARIANCE = 1.0e-10


def _design(rng):
    quantiles = rng.uniform(size=(COUNT_POINTS, 2))
    return ef.Design(
        names=['coolingRate/multiplier', "nodeOperator/nodeOperator[@value='blackHolesSeed']/blackHoleSeeds/mass"],
        priors=[
            ef.Prior('uniform', {'limitLower': 0.05, 'limitUpper': 1.2}),
            ef.Prior('logNormal', {'x0': 3000.0, 'sigma': 4.0, 'limitLower': 100.0, 'limitUpper': 1.0e5}),
        ],
        mappers=['identity', 'logarithm'],
        quantiles=quantiles,
        values=quantiles * 2.0,  # Values need only have the right shape here.
        attributes={'type': 'sobol', 'seed': 42, 'randomShift': True},
    )


def _emulator(quantiles, y, store_factor):
    """Build an emulator of `y` by hand: standardize the bins, keep every principal component, and fit one GP per
    (standardized) coefficient with fixed hyperparameters and negligible noise, so that it interpolates the training data."""
    bin_mean = y.mean(axis=0)
    bin_scale = y.std(axis=0)
    standardized = (y - bin_mean) / bin_scale
    _, _, components = np.linalg.svd(standardized, full_matrices=False)
    coefficients = standardized @ components.T
    coefficient_mean = coefficients.mean(axis=0)
    coefficient_scale = coefficients.std(axis=0)
    log_length_scales = np.log([0.3, 0.5])
    gp_components = []
    for k in range(components.shape[0]):
        targets = (coefficients[:, k] - coefficient_mean[k]) / coefficient_scale[k]
        noise_variance = np.full(COUNT_POINTS, NOISE_VARIANCE)
        covariance = ef.kernel_matern52(quantiles, quantiles, 0.0, log_length_scales)
        covariance[np.diag_indices_from(covariance)] += noise_variance
        factor = np.linalg.cholesky(covariance)
        gp_components.append(ef.Component(
            log_amplitude=0.0,
            log_length_scales=log_length_scales,
            targets=targets,
            noise_variance=noise_variance,
            alpha=np.linalg.solve(covariance, targets),
            log_marginal_likelihood=-1.5,
            cholesky_factor=factor if store_factor else None,
        ))
    return ef.Emulator(
        kernel='matern52ARD',
        jitter=0.0,
        input_names=['coolingRate/multiplier', "nodeOperator/nodeOperator[@value='blackHolesSeed']/blackHoleSeeds/mass"],
        inputs=quantiles,
        bin_mean=bin_mean,
        bin_scale=bin_scale,
        pca_components=components,
        coefficient_mean=coefficient_mean,
        coefficient_scale=coefficient_scale,
        components=gp_components,
        pca_variance_retained=1.0,
    )


def _content(store_factor=True):
    rng = np.random.default_rng(1234)
    design = _design(rng)
    # A smooth function of the inputs, in each of three bins.
    u = design.quantiles
    y = np.stack([np.sin(3.0 * u[:, 0]) + u[:, 1], np.cos(2.0 * u[:, 1]) * u[:, 0], u[:, 0] ** 2 - u[:, 1]], axis=1)
    training_set = ef.TrainingSet(
        x=np.array([9.0, 10.0, 11.0]),
        y=y,
        root_variance=np.full(y.shape, 0.1),
        mask=np.zeros(y.shape, dtype=np.int8),
        point_index=np.arange(COUNT_POINTS),
        y_target=np.array([0.1, 0.2, 0.3]),
        covariance_target=np.diag([0.01, 0.02, 0.03]),
        attributes={'transform': 'log10', 'floor': -6.64, 'rootVarianceFloored': 0.5, 'xAxisIsLog': True},
    )
    validation = ef.Validation(
        folds=np.arange(COUNT_POINTS) % 4,
        held_out_prediction=y + 0.01,
        held_out_variance=np.full(y.shape, 0.01),
        rmse=np.full(COUNT_BINS, 0.01),
        r2=np.full(COUNT_BINS, 0.99),
        rms_standardized_residual=np.full(COUNT_BINS, 1.0),
        coverage_1_sigma=np.full(COUNT_BINS, 0.68),
        coverage_2_sigma=np.full(COUNT_BINS, 0.95),
    )
    label = 'massFunctionStellarTomczak2014ZFOURGEz0'
    return ef.EmulatorFile(
        design=design,
        training_sets={label: training_set},
        emulators={label: _emulator(u, y, store_factor)},
        validation={label: validation},
        attributes={'creator': 'test_emulatorFile'},
    ), label


def _assert_equal(a, b):
    """Recursively compare dataclasses, dicts, lists, arrays and scalars."""
    if hasattr(a, '__dataclass_fields__'):
        assert type(a) is type(b)
        for name in a.__dataclass_fields__:
            _assert_equal(getattr(a, name), getattr(b, name))
    elif isinstance(a, dict):
        assert set(a) == set(b)
        for key in a:
            _assert_equal(a[key], b[key])
    elif isinstance(a, (list, tuple)):
        assert len(a) == len(b)
        for x, y in zip(a, b):
            _assert_equal(x, y)
    elif isinstance(a, np.ndarray):
        np.testing.assert_array_equal(a, b)
    elif isinstance(a, float) and np.isnan(a):
        assert np.isnan(b)
    else:
        assert a == b


def test_round_trip(tmp_path):
    content, _ = _content()
    path = tmp_path / 'emulator.hdf5'
    ef.write(path, content)
    restored = ef.read(path)
    assert restored.attributes['format'] == ef.FORMAT_NAME
    assert restored.attributes['formatVersion'] == ef.FORMAT_VERSION
    assert 'created' in restored.attributes
    for name in ('format', 'formatVersion', 'created'):
        restored.attributes.pop(name)
    _assert_equal(content, restored)


def test_optional_datasets_are_omitted(tmp_path):
    content, label = _content(store_factor=False)
    path = tmp_path / 'emulator.hdf5'
    ef.write(path, content)
    with h5py.File(path, 'r') as file:
        assert 'choleskyFactor' not in file['emulators'][label]['component1']
        assert 'yRaw' not in file['trainingSets'][label]
    assert ef.read(path).emulators[label].components[0].cholesky_factor is None


def test_strings_are_fixed_length(tmp_path):
    """Galacticus reads fixed-length strings, so no string anywhere in the file may be variable-length."""
    content, _ = _content()
    path = tmp_path / 'emulator.hdf5'
    ef.write(path, content)
    found = []

    def check(name, obj):
        items = [(f'{name}@{key}', obj.attrs.get_id(key).get_type()) for key in obj.attrs]
        if isinstance(obj, h5py.Dataset):
            items.append((name, obj.id.get_type()))
        for item, datatype in items:
            if isinstance(datatype, h5py.h5t.TypeStringID):
                found.append(item)
                assert not datatype.is_variable_str(), f'{item} is a variable-length string'

    with h5py.File(path, 'r') as file:
        check('/', file)
        file.visititems(check)
    # Confirm that the check did examine strings, so that a pass is meaningful.
    assert '/@format' in found
    assert 'design/parameterNames' in found


def test_predict_interpolates_training_data():
    """With every principal component retained and negligible noise, the emulator must reproduce its training data,
    with negligible variance, at the training inputs."""
    content, label = _content()
    emulator = content.emulators[label]
    training_set = content.training_sets[label]
    mean, variance = ef.predict(emulator, emulator.inputs)
    np.testing.assert_allclose(mean, training_set.y, rtol=0.0, atol=1.0e-6)
    assert np.all(variance < 1.0e-6)
    # Away from the training inputs the variance must be substantial.
    _, variance = ef.predict(emulator, np.array([0.999, 0.001]))
    assert np.all(variance > 1.0e-3)


def test_predict_without_stored_factor():
    """Predictions must not depend on whether the Cholesky factor is stored or recomputed."""
    content_stored, label = _content(store_factor=True)
    content_computed, _ = _content(store_factor=False)
    points = np.random.default_rng(99).uniform(size=(5, 2))
    mean_stored, variance_stored = ef.predict(content_stored.emulators[label], points)
    mean_computed, variance_computed = ef.predict(content_computed.emulators[label], points)
    np.testing.assert_allclose(mean_computed, mean_stored, rtol=1.0e-12)
    np.testing.assert_allclose(variance_computed, variance_stored, rtol=1.0e-8, atol=1.0e-14)


def test_predict_single_point_shape():
    content, label = _content()
    mean, variance = ef.predict(content.emulators[label], np.array([0.5, 0.5]))
    assert mean.shape == (COUNT_BINS,)
    assert variance.shape == (COUNT_BINS,)


def test_kernel_matern52():
    """Check the kernel against its closed form at a known separation."""
    log_length_scales = np.log([2.0, 0.5])
    value = ef.kernel_matern52(np.array([0.0, 0.0]), np.array([1.0, 0.25]), np.log(3.0), log_length_scales)
    radius = np.sqrt((1.0 / 2.0) ** 2 + (0.25 / 0.5) ** 2)
    expected = 3.0 * (1.0 + np.sqrt(5.0) * radius + 5.0 * radius**2 / 3.0) * np.exp(-np.sqrt(5.0) * radius)
    assert value.shape == (1, 1)
    assert value[0, 0] == pytest.approx(expected, rel=1.0e-14)


@pytest.mark.parametrize('mutate, message', [
    (lambda c, l: setattr(c.design, 'quantiles', c.design.quantiles * 2.0), 'quantiles must lie in'),
    (lambda c, l: setattr(c.training_sets[l], 'point_index', np.zeros(COUNT_POINTS, dtype=int)), 'must be unique'),
    (lambda c, l: setattr(c.training_sets[l], 'point_index', np.arange(COUNT_POINTS) + 1), 'out of range'),
    (lambda c, l: setattr(c.emulators[l], 'input_names', ['coolingRate/multiplier', 'notInDesign']), 'not in the design'),
    (lambda c, l: setattr(c.emulators[l], 'kernel', 'squaredExponential'), 'unknown kernel'),
    (lambda c, l: setattr(c.emulators[l], 'pca_components', c.emulators[l].pca_components[:2]), 'PCA components'),
])
def test_validation_rejects_malformed_content(tmp_path, mutate, message):
    content, label = _content()
    mutate(content, label)
    with pytest.raises(ValueError, match=message):
        ef.write(tmp_path / 'emulator.hdf5', content)
    assert not (tmp_path / 'emulator.hdf5').exists()


def test_read_rejects_other_files(tmp_path):
    path = tmp_path / 'other.hdf5'
    with h5py.File(path, 'w') as file:
        file.attrs['format'] = np.bytes_(b'somethingElse')
    with pytest.raises(ValueError, match='not a Galacticus emulator file'):
        ef.read(path)
    content, _ = _content()
    path = tmp_path / 'future.hdf5'
    ef.write(path, content)
    with h5py.File(path, 'r+') as file:
        file.attrs['formatVersion'] = ef.FORMAT_VERSION + 1
    with pytest.raises(ValueError, match='format version'):
        ef.read(path)


def test_write_refuses_to_overwrite(tmp_path):
    content, _ = _content()
    path = tmp_path / 'emulator.hdf5'
    ef.write(path, content)
    with pytest.raises(FileExistsError):
        ef.write(path, content)
    ef.write(path, content, overwrite=True)
