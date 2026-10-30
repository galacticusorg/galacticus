"""The likelihood forms of the emulated likelihood, evaluated for a given model prediction.

These are the forms of :galacticus-class:`posteriorSampleLikelihoodEmulated` (``gaussianDiagonal``, ``gaussianCovariance``,
``logNormalFloored`` and ``poissonCounts``), evaluated exactly as there: without the variance of the prediction (as with
``includeEmulatorVariance`` false), or with it, if one is given. Applied to the
prediction of a model run itself - a row of the ``y`` of a training set, collected from the run's output by
:func:`Galacticus.Emulation.train.collect` - they give the *direct* likelihood of that run: the likelihood which an MCMC with
an exact emulator would assign to its parameters. This is used to profile parameters by direct runs
(:mod:`Galacticus.Emulation.profile`), and to compare emulated with direct likelihoods.

Andrew Benson, Claude (28-October-2026).
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import gammaln

__all__ = [
    'FORMS',
    'LOG_IMPROBABLE',
    'log_likelihood',
    'profile_root_variance',
]

FORMS = ('gaussianDiagonal', 'gaussianCovariance', 'logNormalFloored', 'poissonCounts')
# The log-likelihood of an impossible model, as `logImprobable` in `source/models/likelihoods/constants.F90`.
LOG_IMPROBABLE = -1.0e-16 * np.finfo(float).max


def log_likelihood(training_set, form, mean, root_variance_target_fractional=0.0, variance=None):
    """Return the log-likelihood of the prediction ``mean`` of one observable.

    ``training_set`` is the observable's :class:`Galacticus.Emulation.emulatorFile.TrainingSet`, which supplies the target
    data (and target counts), and the transform of the observable. ``mean`` is the prediction in the transformed space in
    which the observable is emulated (as are the rows of ``training_set.y``). ``root_variance_target_fractional`` is the
    fractional root variance of counts (:math:`C`) of the ``poissonCounts`` form. ``variance`` is the variance of the
    prediction in each bin (as with ``includeEmulatorVariance`` true), or None for none. For that form, a bin in which no counts
    are expected but some are observed makes the model impossible, and the result is then ``LOG_IMPROBABLE``, as in the
    emulated likelihood.
    """
    if form not in FORMS:
        raise ValueError(f"unknown likelihood form '{form}'")
    mean = np.asarray(mean, dtype=float)
    variance = np.zeros_like(mean) if variance is None else np.asarray(variance, dtype=float)
    yTarget = np.asarray(training_set.y_target, dtype=float)
    varianceTarget = np.diag(np.asarray(training_set.covariance_target, dtype=float))
    transformLog = training_set.attributes.get('transform') == 'log10'
    if form == 'gaussianDiagonal':
        included = np.isfinite(yTarget) & np.isfinite(varianceTarget)
        included[included] &= varianceTarget[included] + variance[included] > 0.0
        residual, varianceTotal = yTarget[included] - mean[included], varianceTarget[included] + variance[included]
        return float(-0.5 * np.sum(residual ** 2 / varianceTotal + np.log(2.0 * np.pi * varianceTotal)))
    if form == 'gaussianCovariance':
        included = np.isfinite(yTarget) & np.isfinite(varianceTarget)
        included[included] &= varianceTarget[included] + variance[included] > 0.0
        if not included.any():
            return 0.0
        covariance = np.asarray(training_set.covariance_target, dtype=float)[np.ix_(included, included)] + np.diag(variance[included])
        if not np.all(np.isfinite(covariance)):
            raise ValueError('the target covariance is not finite between included bins')
        residual = yTarget[included] - mean[included]
        factor = np.linalg.cholesky(covariance)
        solution = np.linalg.solve(factor, residual)
        return float(-0.5 * (solution @ solution + 2.0 * np.sum(np.log(np.diag(factor))) + residual.size * np.log(2.0 * np.pi)))
    if form == 'logNormalFloored':
        if not transformLog or 'floor' not in training_set.attributes:
            raise ValueError('the `logNormalFloored` form requires an observable emulated in its logarithm with a floor')
        floor = float(training_set.attributes['floor'])
        rootVarianceFloored = float(training_set.attributes['rootVarianceFloored'])
        included = ~np.isnan(yTarget) & np.isfinite(varianceTarget)
        with np.errstate(invalid='ignore'):
            targetFloored = ~np.isfinite(yTarget) | (yTarget <= floor)
        valueTarget = np.where(targetFloored, floor, yTarget)
        varianceTotal = np.where(targetFloored, rootVarianceFloored ** 2, varianceTarget)
        modelFloored = ~(mean > floor)
        valueModel = np.where(modelFloored, floor, mean)
        varianceTotal = varianceTotal + np.where(modelFloored, rootVarianceFloored ** 2, 0.0) + variance
        included &= varianceTotal > 0.0
        residual, varianceTotal = valueTarget[included] - valueModel[included], varianceTotal[included]
        return float(-0.5 * np.sum(residual ** 2 / varianceTotal + np.log(2.0 * np.pi * varianceTotal)))
    # The `poissonCounts` form.
    if training_set.count_target is None or training_set.count_per_value is None:
        raise ValueError('the `poissonCounts` form requires an observable with target counts')
    countTarget = np.asarray(training_set.count_target, dtype=float)
    countPerValue = np.asarray(training_set.count_per_value, dtype=float)
    included = np.isfinite(countTarget) & np.isfinite(countPerValue)
    # The mean of the untransformed observable, and the fractional variance of the prediction.
    if transformLog:
        varianceFractional = np.expm1(np.log(10.0) ** 2 * variance[included])
        value = 10.0 ** mean[included] * np.sqrt(1.0 + varianceFractional)
    else:
        value = mean[included]
        with np.errstate(divide='ignore', invalid='ignore'):
            varianceFractional = np.where(value > 0.0, variance[included] / value ** 2, 0.0)
    countExpected = value * countPerValue[included]
    countTarget = countTarget[included]
    empty = countExpected <= 0.0
    if np.any(empty & (countTarget > 0.0)):
        return LOG_IMPROBABLE
    countExpected, countTarget, varianceFractional = countExpected[~empty], countTarget[~empty], varianceFractional[~empty]
    # Combine the fractional variances of the expected count and of the prediction, as in the emulated likelihood.
    C2 = root_variance_target_fractional ** 2
    dispersion = C2 + varianceFractional + C2 * varianceFractional
    logFactorial = gammaln(np.rint(countTarget) + 1.0)
    negativeBinomial = dispersion > 0.0
    terms = countTarget * np.log(countExpected) - logFactorial - countExpected
    if negativeBinomial.any():
        n, mu, d = countTarget[negativeBinomial], countExpected[negativeBinomial], dispersion[negativeBinomial]
        stopping = 1.0 / d
        terms[negativeBinomial] = (n * np.log(mu) - logFactorial[negativeBinomial] + gammaln(stopping + n) - gammaln(stopping)
                                   - n * np.log(stopping + mu) - stopping * np.log1p(mu * d))
    return float(np.sum(terms))


def profile_root_variance(evaluate, bounds=(1.0e-3, 3.0)):
    """Maximize a log-likelihood over the fractional root variance of counts, :math:`C`.

    ``evaluate(C)`` returns the log-likelihood at :math:`C`. The maximum is found over :math:`\\log C` within ``bounds``.
    Returns ``(logLikelihood, C)``.
    """
    if not 0.0 < bounds[0] < bounds[1]:
        raise ValueError('the bounds on C must satisfy 0 < lower < upper')
    result = minimize_scalar(lambda logC: -evaluate(np.exp(logC)), bounds=(np.log(bounds[0]), np.log(bounds[1])),
                             method='bounded', options={'xatol': 1.0e-6})
    return float(-result.fun), float(np.exp(result.x))
