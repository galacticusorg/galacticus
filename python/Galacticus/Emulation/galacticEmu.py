"""Import campaigns of the GalacticEmu package, and check Galacticus' emulation against it.

GalacticEmu (https://github.com/Andrew-Robertson/Galacticus-emulation) runs a Sobol design of Galacticus models, trains
PCA-compressed Gaussian-process emulators of their analysis outputs, and samples the posterior of the model parameters with
`emcee`. This module:

* imports the design of a GalacticEmu campaign (``samples.csv`` and ``campaign_design.json``) as a Galacticus design file
  (:func:`import_design`), so that the campaign's runs can be collected, and emulators trained, with
  :mod:`Galacticus.Emulation.train`;
* maps GalacticEmu's treatment of each observable onto a :class:`Galacticus.Emulation.train.Observable`
  (:func:`observable_from_galacticemu`);
* compares, observable by observable, the training sets, the emulators, and the emulated likelihood with those of
  GalacticEmu (:func:`compare_training`, :func:`compare_emulators`, :func:`compare_likelihood`).

Importing a design needs only the campaign's files. The comparisons, and :func:`observable_from_galacticemu`, use
GalacticEmu itself (including some of its internal functions), and were written against its tag
``roman-grs-calibration-arxiv-v1`` (the version used for arXiv:2609.10680).

Andrew Benson, Claude (2026)
"""

from __future__ import annotations

import json
import os

import h5py
import numpy as np

from Galacticus.Emulation import emulatorFile as ef
from Galacticus.Emulation.design import DesignFile, serialize_descriptor
from Galacticus.Emulation.emulatorFile import Design, Prior
from Galacticus.Emulation.train import Observable

__all__ = [
    'compare_emulators',
    'compare_likelihood',
    'compare_training',
    'import_design',
    'log_likelihood_gaussian_diagonal',
    'observable_from_galacticemu',
    'prior_quantiles',
]


def _prior(prior):
    """Convert a GalacticEmu prior to a Galacticus ``distributionFunction1D`` class and parameters."""
    parameters = prior['parameters']
    if prior['class'] == 'UniformPrior':
        return 'uniform', {'limitLower': float(parameters['lower']), 'limitUpper': float(parameters['upper'])}
    if prior['class'] == 'TruncatedLogNormalPrior':
        # GalacticEmu's truncated log-normal prior is a normal distribution in ln(x), with median x0 and width sigma,
        # truncated to [lower, upper] - as Galacticus' `logNormal` distribution.
        return 'logNormal', {'x0': float(parameters['x0']), 'sigma': float(parameters['sigma']),
                             'limitLower': float(parameters['lower']), 'limitUpper': float(parameters['upper'])}
    raise NotImplementedError(f"GalacticEmu prior class '{prior['class']}' is not supported")


def prior_quantiles(prior_class, parameters, values, inverse=False):
    """The prior cumulative probability of each of ``values`` (or, if ``inverse``, the values at cumulative probabilities
    ``values``) under a Galacticus ``uniform`` or ``logNormal`` prior."""
    from scipy import stats
    values = np.asarray(values, dtype=float)
    lower, upper = parameters['limitLower'], parameters['limitUpper']
    if prior_class == 'uniform':
        return lower + values * (upper - lower) if inverse else (values - lower) / (upper - lower)
    mu, sigma = np.log(parameters['x0']), parameters['sigma']
    a, b = (np.log(lower) - mu) / sigma, (np.log(upper) - mu) / sigma
    if inverse:
        return np.exp(stats.truncnorm.ppf(values, a, b, loc=mu, scale=sigma))
    return stats.truncnorm.cdf(np.log(values), a, b, loc=mu, scale=sigma)


def import_design(campaign_root, log_parameters=None, output_file_name='galacticus.hdf5', tolerance=1.0e-6):
    """Import the design of a GalacticEmu campaign, returning a :class:`Galacticus.Emulation.design.DesignFile` (to be
    written with :func:`Galacticus.Emulation.design.write_design`).

    The parameters are named by their Galacticus parameter paths, and their priors converted to the equivalent Galacticus
    distributions. The value of every parameter at every design point is checked against that implied by its quantile under
    the converted prior (to relative precision ``tolerance``). ``log_parameters`` lists the short names of the parameters
    to be sampled in their logarithm (``logarithm`` mapper); by default, those with log-normal priors. The output of each run
    is ``evaluations/<evaluation_id>/<output_file_name>`` under ``campaign_root``.
    """
    import pandas as pd
    specification = json.load(open(os.path.join(campaign_root, 'campaign_design.json')))
    samples = pd.read_csv(os.path.join(campaign_root, 'samples.csv'))
    free = specification['free_parameters']
    short_names = [parameter['short_name'] for parameter in free]
    quantiles = samples[[f'{name}_quantile' for name in short_names]].to_numpy(dtype=float)
    values = samples[short_names].to_numpy(dtype=float)
    priors, mappers = [], []
    for j, parameter in enumerate(free):
        prior_class, parameters = _prior(parameter['prior'])
        expected = prior_quantiles(prior_class, parameters, quantiles[:, j], inverse=True)
        scale = np.where(values[:, j] != 0.0, np.abs(values[:, j]), 1.0)
        error = np.max(np.abs(expected - values[:, j]) / scale)
        if not error <= tolerance:
            raise ValueError(f"parameter '{short_names[j]}': its values do not follow from its quantiles under the converted prior "
                             f"(relative error {error:.2e})")
        priors.append(Prior(prior_class, {**parameters, 'descriptor': serialize_descriptor(parameters)}))
        logarithmic = short_names[j] in log_parameters if log_parameters is not None else prior_class == 'logNormal'
        mappers.append('logarithm' if logarithmic else 'identity')
    design = Design(
        names=[parameter['path'] for parameter in free],
        priors=priors,
        mappers=mappers,
        quantiles=quantiles,
        values=values,
        attributes={'type': 'sobol', 'seed': int(specification.get('seed', -1)), 'importedFrom': 'GalacticEmu'},
    )
    identifiers = samples['evaluation_id'].tolist()
    count = len(identifiers)
    root = os.path.abspath(campaign_root)
    runs = {
        'pointIndex': np.arange(count, dtype=np.int32),
        'realizationIndex': np.zeros(count, dtype=np.int32),
        'changeFileName': [os.path.join(root, 'evaluations', identifier, 'model_changes.xml') for identifier in identifiers],
        'outputFileName': [os.path.join(root, 'evaluations', identifier, output_file_name) for identifier in identifiers],
    }
    return DesignFile(design=design, runs=runs, attributes={'campaign': root, 'shortNames': ' '.join(short_names)})


def observable_from_galacticemu(config, y_axis_is_log, min_log10_y=-7.0, conventions=True, root_variance_target_minimum=1.0e-3):
    """Return the :class:`Galacticus.Emulation.train.Observable` reproducing GalacticEmu's treatment of an observable.

    ``config`` is the GalacticEmu observable configuration (an entry of its ``OBSERVABLE_CONFIGS``), and ``y_axis_is_log``
    whether the analysis is emulated in its logarithm (its ``yAxisIsLog`` attribute). If ``conventions`` is true,
    GalacticEmu's treatment of zero training noise (``min_training_sigma``) and its floor on the target uncertainty (the
    ``target_sigma_floor`` of its MCMC, by default 1e-3) are also reproduced.

    Only the treatments needed for the standard observables of arXiv:2609.10680 are supported (for which the collected
    training sets were found to be identical to GalacticEmu's): values at or below a floor in log10 kept at the floor
    (``nonpositive`` with ``keep``), and undefined values of a relation replaced by the bin median (``bin_median``).
    """
    condition = config.get('bad_training_condition', 'nonfinite')
    fill = config.get('bad_training_value_fill', 'bin_median')
    sigma = float(config.get('bad_training_sigma', 5.0))
    options = {}
    if y_axis_is_log:
        if condition != 'nonpositive' or str(fill) not in ('keep', 'as_is'):
            raise NotImplementedError(f"log10 observable '{config.get('analysis')}': only bad_training_condition=nonpositive with "
                                      f"bad_training_value_fill=keep is supported (got {condition}, {fill})")
        options.update(transform='log10', floor=float(min_log10_y), root_variance_floored=sigma)
    else:
        if fill != 'bin_median' or condition not in ('zero_only', 'nonpositive', 'zero_and_zero_noise'):
            raise NotImplementedError(f"observable '{config.get('analysis')}': only bad_training_value_fill=bin_median is supported "
                                      f"(got {condition}, {fill})")
        options.update(transform='identity', undefined='median', root_variance_undefined=sigma)
    if conventions:
        options.update(root_variance_zero=float(config.get('min_training_sigma', 1.0e-3)),
                       root_variance_target_minimum=root_variance_target_minimum)
    observable = Observable(label=config['analysis'], **options)
    observable.validate()
    return observable


def galacticemu_observables(campaign_root, keys, conventions=True, root_variance_target_minimum=1.0e-3):
    """Return a :class:`Galacticus.Emulation.train.Observable` for each of the GalacticEmu observables ``keys``, reading
    whether each is emulated in its logarithm from the first run of the campaign. Requires GalacticEmu."""
    from galacticus_emu import interactive_observables as gi
    configs = gi._observable_configs(None)
    unknown = [key for key in keys if key not in configs]
    if unknown:
        raise ValueError(f'unknown GalacticEmu observable(s): {unknown}')
    observables = []
    for key in keys:
        attributes = gi._infer_analysis_attributes(_path(campaign_root), configs[key]['analysis'], 'galacticus.hdf5')
        observables.append(observable_from_galacticemu(configs[key], bool(attributes.get('yAxisIsLog', False)), conventions=conventions,
                                                       root_variance_target_minimum=root_variance_target_minimum))
    return observables


def _path(path):
    from pathlib import Path
    return Path(path)


def compare_training(campaign_root, design_file, keys, conventions=True):
    """Compare the training sets collected by :func:`Galacticus.Emulation.train.collect` with those GalacticEmu prepares,
    for each of the GalacticEmu observables ``keys``. Returns, for each, a dictionary of the maximum absolute differences in
    the training values (``y``), their relative differences in root variance (``rootVariance``), and the differences in the
    target data (``target``, ``rootVarianceTarget``; not-a-number bins excluded). Requires GalacticEmu."""
    from galacticus_emu import interactive_observables as gi
    from Galacticus.Emulation.train import collect
    configs = gi._observable_configs(None)
    observables = galacticemu_observables(campaign_root, keys, conventions=conventions)
    _, training_sets, _ = collect(design_file, observables)
    root = _path(campaign_root)
    columns = gi._input_columns_for_samples(gi.pd.read_csv(root / 'samples.csv'), None)
    specifications = gi._parameter_specs_for_columns(columns)
    results = {}
    for key, observable in zip(keys, observables):
        config = configs[key]
        (_, _, _, target, target_noise, y, y_noise, _, transform) = gi._load_observable_campaign(
            root, analysis=config['analysis'], hdf5_filename='galacticus.hdf5', min_log10_y=-7.0,
            input_columns=columns, input_parameter_specs=specifications)
        supported = np.ones(y.shape[1], dtype=bool)
        if config.get('drop_unsupported_bins', False):
            supported = gi._supported_bin_mask(y, y_noise=y_noise, bad_training_condition=config['bad_training_condition'])
        y_fit, sigma, _ = gi._prepare_training_targets(
            y[:, supported], None if y_noise is None else y_noise[:, supported], analysis=config['analysis'], use_training_alpha=True,
            bad_training_condition=config['bad_training_condition'], bad_training_value_fill=config['bad_training_value_fill'],
            bad_training_sigma=config['bad_training_sigma'], min_training_sigma=config['min_training_sigma'],
            bad_mask_override=gi.training_bad_mask_override_from_transform(
                transform, bad_training_condition=config['bad_training_condition'], bad_training_value_fill=config['bad_training_value_fill']))
        training_set = training_sets[observable.label]
        target_sigma = np.sqrt(np.diag(training_set.covariance_target))[supported]
        target_noise = np.maximum(target_noise, observable.root_variance_target_minimum or 0.0)
        with np.errstate(invalid='ignore'):
            results[key] = {
                'label': observable.label,
                'bins': int(supported.sum()),
                'binsUnsupported': int((~supported).sum()),
                'y': float(np.max(np.abs(training_set.y[:, supported] - y_fit))),
                'rootVariance': float(np.max(np.abs(training_set.root_variance[:, supported] / sigma - 1.0))),
                'target': float(np.nanmax(np.abs(training_set.y_target[supported] - target))),
                'rootVarianceTarget': float(np.nanmax(np.abs(target_sigma / target_noise - 1.0))),
            }
    return results


def _bundle_prediction(observable, quantiles):
    """Predict an observable with a GalacticEmu PCA emulator (without its display offset): mean and standard deviation."""
    from galacticus_emu.gp import predict_scaled_gp
    coefficient = np.zeros((quantiles.shape[0], len(observable['models'])))
    coefficient_std = np.zeros_like(coefficient)
    for k, model in enumerate(observable['models']):
        coefficient[:, k], coefficient_std[:, k] = predict_scaled_gp(model, float(observable['y_means'][k]), float(observable['y_stds'][k]), quantiles)
    components = observable['pca_components_matrix']
    mean = (coefficient @ components + observable['pca_mean']) * observable['preprocessor_scale'] + observable['preprocessor_mean']
    std = np.sqrt(np.maximum(coefficient_std ** 2 @ components ** 2, 0.0)) * observable['preprocessor_scale']
    return mean, std


def compare_emulators(bundle, content, keys, count_points=1000, seed=1234, quantiles=None):
    """Compare the emulators of a Galacticus emulator file (``content``, an :class:`Galacticus.Emulation.emulatorFile.EmulatorFile`)
    with those of a GalacticEmu bundle (as loaded by ``joblib``), for the GalacticEmu observables ``keys``.

    Returns, for each, the number of principal components in each, the absolute cosines between corresponding principal
    components, the log marginal likelihoods of the Gaussian process of each component, and statistics of the differences
    in the predicted mean (in units of the target uncertainty, and of the combined predictive uncertainty) and the ratio
    of predictive uncertainties, at ``count_points`` points drawn uniformly in the prior quantiles (or at ``quantiles``).
    Requires GalacticEmu."""
    if quantiles is None:
        quantiles = np.random.default_rng(seed).uniform(size=(count_points, len(bundle['input_columns'])))
    results = {}
    for key in keys:
        observable = bundle['observables'][key]
        emulator = content.emulators[observable['analysis']]
        training_set = content.training_sets[observable['analysis']]
        supported = np.asarray(observable['supported_bin_mask'], dtype=bool)
        mean_theirs, std_theirs = _bundle_prediction(observable, quantiles)
        mean_ours, variance_ours = ef.predict(emulator, quantiles)
        mean_ours, std_ours = mean_ours[:, supported], np.sqrt(variance_ours[:, supported])
        sigma_target = np.sqrt(np.diag(training_set.covariance_target))[supported]
        with np.errstate(divide='ignore', invalid='ignore'):
            difference_target = np.abs(mean_ours - mean_theirs) / sigma_target
        difference_emulator = np.abs(mean_ours - mean_theirs) / np.sqrt(std_ours ** 2 + std_theirs ** 2)
        ratio = std_ours / std_theirs
        count = min(len(emulator.components), len(observable['models']))
        results[key] = {
            'label': observable['analysis'],
            'components': (len(emulator.components), len(observable['models'])),
            'cosines': [abs(float(np.dot(emulator.pca_components[k], observable['pca_components_matrix'][k]))) for k in range(count)],
            'logMarginalLikelihood': ([component.log_marginal_likelihood for component in emulator.components],
                                      [float(model.log_marginal_likelihood_value_) for model in observable['models']]),
            'meanDifferenceTarget': (float(np.nanmedian(difference_target)), float(np.nanpercentile(difference_target, 95)), float(np.nanmax(difference_target))),
            'meanDifferenceEmulator': (float(np.median(difference_emulator)), float(np.percentile(difference_emulator, 95))),
            'stdRatio': (float(np.median(ratio)), float(np.percentile(ratio, 5)), float(np.percentile(ratio, 95))),
        }
    return results


def log_likelihood_gaussian_diagonal(content, labels, quantiles, include_emulator_variance=True):
    """The emulated log-likelihood of the ``gaussianDiagonal`` form of ``posteriorSampleLikelihoodEmulated``, summed over the
    observables ``labels`` of an emulator file, at each of ``quantiles`` (an array of shape (N, d) of prior quantiles). As in
    Galacticus, bins whose target datum or target variance is not finite are excluded, as is any bin whose total variance
    is not positive."""
    quantiles = np.atleast_2d(quantiles)
    total = np.zeros(quantiles.shape[0])
    for label in labels:
        mean, variance = ef.predict(content.emulators[label], quantiles)
        training_set = content.training_sets[label]
        varianceTarget = np.diag(training_set.covariance_target)
        included = np.isfinite(training_set.y_target) & np.isfinite(varianceTarget)
        varianceTotal = varianceTarget[None, included] + (variance[:, included] if include_emulator_variance else 0.0)
        residual = training_set.y_target[None, included] - mean[:, included]
        positive = varianceTotal > 0.0
        with np.errstate(divide='ignore', invalid='ignore'):
            terms = residual ** 2 / varianceTotal + np.log(2.0 * np.pi * varianceTotal)
        total += -0.5 * np.sum(np.where(positive, terms, 0.0), axis=1)
    return total


def compare_likelihood(bundle, inputs, content, quantiles):
    """Compare the emulated log-likelihood of an emulator file (``gaussianDiagonal`` with emulator variance, as used by
    arXiv:2609.10680) with that of a GalacticEmu MCMC, at each of ``quantiles``. ``bundle`` is the GalacticEmu emulator
    bundle and ``inputs`` the MCMC inputs (``*_mcmc_inputs.joblib``), each as loaded by ``joblib``. Returns the two
    log-likelihoods. Requires GalacticEmu."""
    from galacticus_emu.observable_plot_metadata import y_display_offset
    quantiles = np.atleast_2d(quantiles)
    means, variances = [], []
    for key in inputs['observable_keys']:
        observable = bundle['observables'][key]
        mean, std = _bundle_prediction(observable, quantiles)
        mask = np.asarray(inputs['masks'][key], dtype=bool)
        # GalacticEmu adds a display offset to both its predictions and its targets, so it cancels in the likelihood.
        means.append(mean[:, mask] + y_display_offset(observable))
        variances.append(std[:, mask] ** 2)
    mean, variance = np.concatenate(means, axis=1), np.concatenate(variances, axis=1)
    variance = np.maximum(np.asarray(inputs['target_variance'])[None, :] + variance, 1.0e-12)
    theirs = -0.5 * np.sum((np.asarray(inputs['target_vector'])[None, :] - mean) ** 2 / variance + np.log(2.0 * np.pi * variance), axis=1)
    labels = [bundle['observables'][key]['analysis'] for key in inputs['observable_keys']]
    ours = log_likelihood_gaussian_diagonal(content, labels, quantiles)
    return ours, theirs


def map_quantiles(results_file):
    """Return the prior quantiles, and log-probability, of the maximum a posteriori point of a GalacticEmu MCMC (from its
    ``*_mcmc_results.hdf5`` file)."""
    with h5py.File(results_file, 'r') as file:
        return file['map/quantile'][()], float(file['map'].attrs['log_probability'])
