"""Tests for `Galacticus.Emulation.likelihood`, against independent implementations of each distribution (`scipy.stats`)."""

import numpy as np
import pytest
from scipy import stats

from Galacticus.Emulation import emulatorFile as ef
from Galacticus.Emulation import likelihood as lk


def _training_set(y_target, covariance_target, transform='identity', count_target=None, count_per_value=None, **attributes):
    count = len(y_target)
    return ef.TrainingSet(x=np.arange(count, dtype=float), y=np.zeros((1, count)), root_variance=np.ones((1, count)),
                          mask=np.zeros((1, count), dtype=np.int8), point_index=np.array([0]),
                          y_target=np.asarray(y_target, dtype=float), covariance_target=np.asarray(covariance_target, dtype=float),
                          y_raw=np.zeros((1, count)), count_target=count_target, count_per_value=count_per_value,
                          attributes={'transform': transform, **attributes})


def test_gaussian_diagonal():
    rng = np.random.default_rng(1)
    target, sigma, mean, variance = rng.normal(size=6), rng.uniform(0.1, 1.0, 6), rng.normal(size=6), rng.uniform(0.0, 0.2, 6)
    target[2] = np.nan  # Excluded.
    training_set = _training_set(target, np.diag(sigma ** 2))
    used = np.isfinite(target)
    expected = np.sum(stats.norm.logpdf(target[used], mean[used], sigma[used]))
    assert lk.log_likelihood(training_set, 'gaussianDiagonal', mean) == pytest.approx(expected, rel=1e-12)
    expected = np.sum(stats.norm.logpdf(target[used], mean[used], np.sqrt(sigma[used] ** 2 + variance[used])))
    assert lk.log_likelihood(training_set, 'gaussianDiagonal', mean, variance=variance) == pytest.approx(expected, rel=1e-12)


def test_gaussian_covariance():
    rng = np.random.default_rng(2)
    factor = rng.normal(size=(5, 5))
    covariance = factor @ factor.T + 0.5 * np.eye(5)
    target, mean, variance = rng.normal(size=5), rng.normal(size=5), rng.uniform(0.0, 0.2, 5)
    training_set = _training_set(target, covariance)
    expected = stats.multivariate_normal.logpdf(target, mean, covariance)
    assert lk.log_likelihood(training_set, 'gaussianCovariance', mean) == pytest.approx(expected, rel=1e-12)
    expected = stats.multivariate_normal.logpdf(target, mean, covariance + np.diag(variance))
    assert lk.log_likelihood(training_set, 'gaussianCovariance', mean, variance=variance) == pytest.approx(expected, rel=1e-12)


def test_log_normal_floored():
    floor, rootVarianceFloored = -3.0, 0.5
    # Bins: both above the floor; target floored (empty, -infinity); model floored; both floored; no target datum.
    target = np.array([-1.0, -np.inf, -1.5, -4.0, np.nan])
    sigma = np.array([0.1, 0.1, 0.2, 0.1, 0.1])
    mean = np.array([-1.2, -2.0, -3.5, -3.2, -1.0])
    training_set = _training_set(target, np.diag(sigma ** 2), transform='log10', floor=floor, rootVarianceFloored=rootVarianceFloored)
    expected = (stats.norm.logpdf(-1.0, -1.2, 0.1)
                + stats.norm.logpdf(floor, -2.0, rootVarianceFloored)
                + stats.norm.logpdf(-1.5, floor, np.sqrt(0.2 ** 2 + rootVarianceFloored ** 2))
                + stats.norm.logpdf(floor, floor, np.sqrt(2.0) * rootVarianceFloored))
    assert lk.log_likelihood(training_set, 'logNormalFloored', mean) == pytest.approx(expected, rel=1e-12)
    with pytest.raises(ValueError):
        lk.log_likelihood(_training_set(target, np.diag(sigma ** 2)), 'logNormalFloored', mean)


@pytest.mark.parametrize('transform', ['identity', 'log10'])
def test_poisson_counts(transform):
    rng = np.random.default_rng(3)
    countTarget = rng.poisson(20.0, 8).astype(float)
    countTarget[5] = np.nan  # Outside the range of the likelihood.
    countPerValue = rng.uniform(50.0, 200.0, 8)
    value = rng.uniform(0.05, 0.3, 8)
    mean = np.log10(value) if transform == 'log10' else value
    training_set = _training_set(np.full(8, np.nan), np.full((8, 8), np.nan), transform=transform,
                                 count_target=countTarget, count_per_value=countPerValue)
    used = np.isfinite(countTarget)
    n, mu = countTarget[used], (value * countPerValue)[used]
    # Poisson, for C = 0.
    assert lk.log_likelihood(training_set, 'poissonCounts', mean) == pytest.approx(np.sum(stats.poisson.logpmf(n, mu)), rel=1e-12)
    # Negative binomial, with stopping-time parameter r = 1/C^2, for C > 0.
    C = 0.3
    r = 1.0 / C ** 2
    expected = np.sum(stats.nbinom.logpmf(n, r, r / (r + mu)))
    assert lk.log_likelihood(training_set, 'poissonCounts', mean, root_variance_target_fractional=C) == pytest.approx(expected, rel=1e-12)
    # A model with no counts in a bin with target counts is impossible.
    empty = mean.copy()
    empty[0] = -np.inf if transform == 'log10' else 0.0
    assert lk.log_likelihood(training_set, 'poissonCounts', empty, root_variance_target_fractional=C) == lk.LOG_IMPROBABLE


def test_poisson_counts_variance():
    """The variance of the prediction is combined with C as in the emulated likelihood."""
    countTarget, countPerValue = np.array([12.0, 30.0]), np.array([100.0, 100.0])
    C, value, variance = 0.2, np.array([0.1, 0.25]), np.array([1.0e-4, 4.0e-4])
    training_set = _training_set(np.full(2, np.nan), np.full((2, 2), np.nan), count_target=countTarget, count_per_value=countPerValue)
    dispersion = C ** 2 + variance / value ** 2 + C ** 2 * variance / value ** 2
    r = 1.0 / dispersion
    mu = value * countPerValue
    expected = np.sum(stats.nbinom.logpmf(countTarget, r, r / (r + mu)))
    assert lk.log_likelihood(training_set, 'poissonCounts', value, root_variance_target_fractional=C, variance=variance) == pytest.approx(expected, rel=1e-12)
    # For an observable emulated in its logarithm, the variance (in dex^2) is that of a lognormal.
    training_set.attributes['transform'] = 'log10'
    varianceFractional = np.expm1(np.log(10.0) ** 2 * variance)
    dispersion = C ** 2 + varianceFractional + C ** 2 * varianceFractional
    r = 1.0 / dispersion
    mu = value * np.sqrt(1.0 + varianceFractional) * countPerValue
    expected = np.sum(stats.nbinom.logpmf(countTarget, r, r / (r + mu)))
    assert lk.log_likelihood(training_set, 'poissonCounts', np.log10(value), root_variance_target_fractional=C,
                             variance=variance) == pytest.approx(expected, rel=1e-12)


def test_profile_root_variance():
    rng = np.random.default_rng(4)
    mu = np.full(4000, 50.0)
    C = 0.25
    r = 1.0 / C ** 2
    n = rng.negative_binomial(r, r / (r + mu)).astype(float)
    training_set = _training_set(np.full(n.size, np.nan), np.full((n.size, n.size), np.nan), count_target=n,
                                 count_per_value=np.ones(n.size))
    best, CBest = lk.profile_root_variance(lambda C: lk.log_likelihood(training_set, 'poissonCounts', mu, root_variance_target_fractional=C))
    assert CBest == pytest.approx(C, rel=0.05)
    grid = np.linspace(0.8 * CBest, 1.2 * CBest, 41)
    assert best >= max(lk.log_likelihood(training_set, 'poissonCounts', mu, root_variance_target_fractional=C) for C in grid) - 1e-6


def test_unknown_form():
    with pytest.raises(ValueError):
        lk.log_likelihood(_training_set([1.0], [[1.0]]), 'gaussian', np.array([1.0]))
