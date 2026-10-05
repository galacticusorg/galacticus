"""Tests for `Galacticus.Emulation.train` (collecting training sets, and training emulators)."""

import h5py
import numpy as np
import pytest

from Galacticus.Emulation import emulatorFile as ef
from Galacticus.Emulation.train import Observable, collect, cross_validate, fit_emulator, train

COUNT_BINS = 8
NOISE_FRACTIONAL = 0.05


def truth_number_density(u, b):
    """A smooth number density, in each bin b, as a function of the inputs (prior quantiles) u."""
    return 10.0 ** (-2.0 - 0.4 * b + 0.8 * np.sin(2.0 * u[0]) + 0.5 * u[1] * b / COUNT_BINS)


def truth_relation(u, b):
    """A smooth mean relation."""
    return 1.0 + u[0] * b / 4.0 + u[1] ** 2


def write_design(path, quantiles, point_index, realization_index):
    """Write a design file in the layout of the Fortran `emulatorDesign` task."""
    countRuns = len(point_index)
    with h5py.File(path, 'w') as file:
        file.attrs['format'] = np.bytes_(b'galacticusDesign')
        file.attrs['formatVersion'] = np.int32(1)
        group = file.create_group('design')
        group.attrs['countPoints'] = np.int32(quantiles.shape[0])
        group.attrs['realizationsPerPoint'] = np.int32(1)
        group.attrs['seedPerPoint'] = np.int32(0)
        group.create_dataset('parameterNames', data=np.array([b'a/b', b'c/d']))
        group.create_dataset('mappers', data=np.array([b'identity', b'identity']))
        group.create_dataset('quantiles', data=quantiles)
        group.create_dataset('values', data=quantiles)
        priors = group.create_group('priors')
        for i in range(2):
            prior = priors.create_group(f'prior{i + 1}')
            prior.attrs['class'] = np.bytes_(b'uniform')
            prior.attrs['descriptor'] = np.bytes_(b'limitLower:0.0_limitUpper:1.0')
        runs = file.create_group('runs')
        runs.create_dataset('pointIndex', data=np.asarray(point_index, dtype=np.int32))
        runs.create_dataset('realizationIndex', data=np.asarray(realization_index, dtype=np.int32))
        runs.create_dataset('changeFileName', data=np.array([f'changes/point{i}.xml'.encode() for i in range(countRuns)]))
        runs.create_dataset('outputFileName', data=np.array([f'models/run{i}.hdf5'.encode() for i in range(countRuns)]))


def write_analysis(file, label, x, y, covariance, y_target=None, covariance_target=None):
    """Write an `analyses/<label>` group following the standard function1D attribute convention."""
    group = file.require_group('analyses').create_group(label)
    group.attrs['type'] = np.bytes_(b'function1D')
    group.attrs['xDataset'] = np.bytes_(b'x')
    group.attrs['yDataset'] = np.bytes_(b'y')
    group.attrs['yCovariance'] = np.bytes_(b'yCovariance')
    group.create_dataset('x', data=x)
    group.create_dataset('y', data=y)
    group.create_dataset('yCovariance', data=covariance)
    if y_target is not None:
        group.attrs['yDatasetTarget'] = np.bytes_(b'yTarget')
        group.attrs['yCovarianceTarget'] = np.bytes_(b'yCovarianceTarget')
        group.create_dataset('yTarget', data=y_target)
        group.create_dataset('yCovarianceTarget', data=covariance_target)


@pytest.fixture
def campaign(tmp_path):
    """A campaign of 40 runs (one per design point), with two observables. Run 3 failed (no output); in run 5 the last bin
    of the number density is empty; in run 7 the last bin of the relation is undefined."""
    rng = np.random.default_rng(7)
    countPoints = 40
    quantiles = rng.uniform(size=(countPoints, 2))
    write_design(tmp_path / 'design.hdf5', quantiles, np.arange(countPoints), np.zeros(countPoints))
    (tmp_path / 'models').mkdir()
    x = np.arange(COUNT_BINS, dtype=float)
    yTarget = np.array([truth_number_density([0.5, 0.5], b) for b in range(COUNT_BINS)])
    for i in range(countPoints):
        if i == 3:
            continue
        u = quantiles[i]
        density = np.array([truth_number_density(u, b) for b in range(COUNT_BINS)])
        relation = np.array([truth_relation(u, b) for b in range(COUNT_BINS)])
        densityCovariance = np.diag((NOISE_FRACTIONAL * density) ** 2)
        relationCovariance = np.diag(np.full(COUNT_BINS, 0.01 ** 2))
        if i == 5:
            density[-1] = 0.0
            densityCovariance[-1, -1] = 0.0
        if i == 7:
            relation[-1] = 0.0
            relationCovariance[-1, -1] = 0.0
        with h5py.File(tmp_path / 'models' / f'run{i}.hdf5', 'w') as file:
            write_analysis(file, 'numberDensity', x, density, densityCovariance, yTarget, np.diag((0.1 * yTarget) ** 2))
            write_analysis(file, 'relation', x, relation, relationCovariance)
    observables = [
        Observable('numberDensity', transform='log10', floor=-6.0, root_variance_floored=0.5),
        Observable('relation', undefined='median', root_variance_undefined=5.0),
    ]
    return tmp_path, quantiles, observables


def test_collect(campaign):
    path, quantiles, observables = campaign
    design, training_sets, reports = collect(path / 'design.hdf5', observables, directory=path)
    density = training_sets['numberDensity']
    relation = training_sets['relation']
    # The failed run is excluded, and reported.
    assert reports['numberDensity'].points_missing == [3]
    assert 3 not in density.point_index and density.y.shape == (39, COUNT_BINS)
    # Values are transformed to log10, with uncertainties converted to dex.
    row = int(np.where(density.point_index == 0)[0][0])
    expected = np.log10([truth_number_density(quantiles[0], b) for b in range(COUNT_BINS)])
    np.testing.assert_allclose(density.y[row], expected, rtol=1.0e-12)
    np.testing.assert_allclose(density.root_variance[row], NOISE_FRACTIONAL / np.log(10.0), rtol=1.0e-12)
    # The empty bin is floored, with the inflated uncertainty, and masked.
    row = int(np.where(density.point_index == 5)[0][0])
    assert density.y[row, -1] == -6.0 and density.root_variance[row, -1] == 0.5 and density.mask[row, -1] == 1
    assert density.mask.sum() == 1
    # The undefined bin is replaced by the median of that bin over the other runs, and masked.
    row = int(np.where(relation.point_index == 7)[0][0])
    others = relation.point_index != 7
    assert relation.y[row, -1] == pytest.approx(np.median(relation.y[others, -1]))
    assert relation.root_variance[row, -1] == 5.0 and relation.mask[row, -1] == 1 and relation.mask.sum() == 1
    # The target is transformed too, and the transform is recorded.
    np.testing.assert_allclose(density.y_target, np.log10([truth_number_density([0.5, 0.5], b) for b in range(COUNT_BINS)]))
    np.testing.assert_allclose(np.sqrt(np.diag(density.covariance_target)), 0.1 / np.log(10.0), rtol=1.0e-12)
    assert density.attributes == {'transform': 'log10', 'floor': -6.0, 'rootVarianceFloored': 0.5}
    assert design.names == ['a/b', 'c/d']


def test_collect_averages_realizations(tmp_path):
    """Two realizations of each of two design points are averaged, and the variance of the average halved."""
    quantiles = np.array([[0.2, 0.3], [0.7, 0.6]])
    write_design(tmp_path / 'design.hdf5', quantiles, [0, 0, 1, 1], [0, 1, 0, 1])
    (tmp_path / 'models').mkdir()
    values = [[1.0, 2.0], [3.0, 6.0], [5.0, 5.0], [7.0, 9.0]]
    for i, y in enumerate(values):
        with h5py.File(tmp_path / 'models' / f'run{i}.hdf5', 'w') as file:
            write_analysis(file, 'relation', np.array([0.0, 1.0]), np.array(y), np.diag([0.04, 0.16]))
    _, training_sets, _ = collect(tmp_path / 'design.hdf5', [Observable('relation')], directory=tmp_path)
    training_set = training_sets['relation']
    np.testing.assert_array_equal(training_set.point_index, [0, 1])
    np.testing.assert_allclose(training_set.y, [[2.0, 4.0], [6.0, 7.0]])
    np.testing.assert_allclose(training_set.root_variance, np.sqrt([[0.02, 0.08], [0.02, 0.08]]))


def test_fit_matches_scikit_learn():
    """The emulator file evaluates exactly what scikit-learn fitted: for a single bin (so that PCA is trivial), predictions
    by `emulatorFile.predict` must equal those of a scikit-learn regressor with the fitted hyperparameters."""
    sklearn = pytest.importorskip('sklearn')  # noqa: F841
    from sklearn.gaussian_process import GaussianProcessRegressor
    from sklearn.gaussian_process.kernels import ConstantKernel, Matern

    rng = np.random.default_rng(3)
    inputs = rng.uniform(size=(30, 2))
    y = (np.sin(3.0 * inputs[:, 0]) + inputs[:, 1] ** 2)[:, None]
    rootVariance = np.full_like(y, 0.02)
    emulator = fit_emulator(inputs, ['a', 'b'], y, rootVariance, restarts=1)
    component = emulator.components[0]
    kernel = (ConstantKernel(np.exp(component.log_amplitude), 'fixed')
              * Matern(length_scale=np.exp(component.log_length_scales), length_scale_bounds='fixed', nu=2.5))
    regressor = GaussianProcessRegressor(kernel=kernel, alpha=component.noise_variance + emulator.jitter, optimizer=None)
    regressor.fit(inputs, component.targets)
    points = rng.uniform(size=(10, 2))
    meanStandardized, sigmaStandardized = regressor.predict(points, return_std=True)
    scale = emulator.bin_scale[0] * emulator.pca_components[0, 0] * emulator.coefficient_scale[0]
    meanExpected = emulator.bin_mean[0] + emulator.bin_scale[0] * emulator.pca_components[0, 0] * (emulator.coefficient_mean[0] + emulator.coefficient_scale[0] * meanStandardized)
    mean, variance = ef.predict(emulator, points)
    np.testing.assert_allclose(mean[:, 0], meanExpected, rtol=1.0e-9, atol=1.0e-12)
    np.testing.assert_allclose(np.sqrt(variance[:, 0]), np.abs(scale) * sigmaStandardized, rtol=1.0e-7, atol=1.0e-12)


def test_train_and_predict(campaign, tmp_path):
    """Train on the campaign, write and read the emulator file, and check predictions at points not in the design."""
    pytest.importorskip('sklearn')
    path, _, observables = campaign
    design, training_sets, _ = collect(path / 'design.hdf5', observables, directory=path)
    content = train(design, training_sets, restarts=2, folds=0)
    ef.write(tmp_path / 'emulator.hdf5', content)
    restored = ef.read(tmp_path / 'emulator.hdf5')
    rng = np.random.default_rng(11)
    points = rng.uniform(0.1, 0.9, size=(20, 2))
    mean, variance = ef.predict(restored.emulators['numberDensity'], points)
    truth = np.log10([[truth_number_density(u, b) for b in range(COUNT_BINS)] for u in points])
    rmse = np.sqrt(np.mean((mean - truth) ** 2))
    assert rmse < 0.03, f'emulator RMSE {rmse} dex'
    assert np.all(variance > 0.0)
    # The fraction of retained PCA variance is respected, and fewer components than bins are needed for a smooth function.
    assert 1 <= len(restored.emulators['numberDensity'].components) < COUNT_BINS
    # Predictions do not depend on whether the Cholesky factors are stored.
    contentNoFactor = train(design, training_sets, restarts=2, folds=0, store_cholesky=False)
    meanNoFactor, varianceNoFactor = ef.predict(contentNoFactor.emulators['numberDensity'], points)
    np.testing.assert_allclose(meanNoFactor, mean, rtol=1.0e-8)
    np.testing.assert_allclose(varianceNoFactor, variance, rtol=1.0e-5, atol=1.0e-12)


def test_cross_validation(campaign):
    """Cross-validation predicts each held-out point, with sensible coverage."""
    pytest.importorskip('sklearn')
    path, quantiles, observables = campaign
    _, training_sets, _ = collect(path / 'design.hdf5', observables[:1], directory=path)
    training_set = training_sets['numberDensity']
    validation = cross_validate(quantiles[training_set.point_index], ['a/b', 'c/d'], training_set.y,
                                training_set.root_variance, folds=5, restarts=1)
    assert validation.held_out_prediction.shape == training_set.y.shape
    assert sorted(set(validation.folds)) == [0, 1, 2, 3, 4]
    # Away from the floored bin the emulator should predict held-out points well, and its uncertainties be roughly calibrated.
    assert np.all(validation.r2[:-1] > 0.9)
    assert np.all(validation.coverage_2_sigma[:-1] > 0.75)
    assert np.all(np.isfinite(validation.rms_standardized_residual))


def test_observable_validation():
    with pytest.raises(ValueError, match='unknown transform'):
        Observable('a', transform='sqrt').validate()
    with pytest.raises(ValueError, match='given together'):
        Observable('a', transform='log10', floor=-6.0).validate()
    with pytest.raises(ValueError, match='rootVarianceUndefined'):
        Observable('a', undefined='median').validate()


def test_command_line(campaign, tmp_path, capsys):
    """The command-line trainer reads its configuration, trains, and writes an emulator file."""
    pytest.importorskip('sklearn')
    import importlib.util
    import pathlib
    path, _, _ = campaign
    script = pathlib.Path(__file__).resolve().parents[4] / 'scripts' / 'emulation' / 'emulatorTrainReference.py'
    specification = importlib.util.spec_from_file_location('emulatorTrainReference', script)
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    configuration = tmp_path / 'train.xml'
    configuration.write_text(f'''<parameters>
  <task value="emulatorTrain">
    <designFileName       value="{path / 'design.hdf5'}"/>
    <emulatorFileName     value="{tmp_path / 'emulator.hdf5'}"/>
    <restartsOptimizer    value="1"/>
    <foldsCrossValidation value="3"/>
    <observable label="numberDensity" transform="log10" floor="-6.0" rootVarianceFloored="0.5"/>
    <observable label="relation" undefined="median" rootVarianceUndefined="5.0"/>
  </task>
</parameters>
''')
    assert module.main([str(configuration), '--directory', str(path)]) == 0
    output = capsys.readouterr().out
    assert 'numberDensity: 39 design points used, 1 missing (points 3)' in output
    assert 'cross-validation' in output
    content = ef.read(tmp_path / 'emulator.hdf5')
    assert set(content.emulators) == {'numberDensity', 'relation'}
    assert set(content.validation) == {'numberDensity', 'relation'}
