"""Tests for `Galacticus.Emulation.galacticEmu` (importing GalacticEmu campaigns). These need neither GalacticEmu nor a
real campaign: a small campaign is written in GalacticEmu's layout."""

import csv
import json

import numpy as np
import pytest
from scipy import stats

from Galacticus.Emulation import galacticEmu as ge
from Galacticus.Emulation.design import parse_descriptor, read_design, write_design

PRIORS = [
    {'short_name': 'velocity', 'path': "nodeOperator/nodeOperator[@value='stellarFeedbackDisks']/stellarFeedbackOutflows/stellarFeedbackOutflows/velocityCharacteristic",
     'prior': {'class': 'TruncatedLogNormalPrior', 'parameters': {'lower': 25.0, 'upper': 300.0, 'x0': 150.0, 'sigma': 0.5}}},
    {'short_name': 'exponent', 'path': 'hotHaloOutflowReincorporation/delta1',
     'prior': {'class': 'UniformPrior', 'parameters': {'lower': -2.0, 'upper': 4.0}}},
]


def write_campaign(path, count=16, perturb=False):
    """Write a GalacticEmu campaign's design files: `campaign_design.json` and `samples.csv`."""
    rng = np.random.default_rng(11)
    quantiles = rng.uniform(size=(count, 2))
    p = PRIORS[0]['prior']['parameters']
    mu, sigma = np.log(p['x0']), p['sigma']
    a, b = (np.log(p['lower']) - mu) / sigma, (np.log(p['upper']) - mu) / sigma
    velocity = np.exp(stats.truncnorm.ppf(quantiles[:, 0], a, b, loc=mu, scale=sigma))
    exponent = -2.0 + 6.0 * quantiles[:, 1]
    if perturb:
        velocity[3] *= 1.01
    identifiers = [f'campaign-eval-{i:04d}' for i in range(count)]
    with open(path / 'samples.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['evaluation_id', 'velocity_quantile', 'exponent_quantile', 'velocity', 'exponent'])
        for row in zip(identifiers, quantiles[:, 0], quantiles[:, 1], velocity, exponent):
            writer.writerow([row[0]] + [repr(float(value)) for value in row[1:]])
    json.dump({'seed': 42, 'free_parameters': PRIORS}, open(path / 'campaign_design.json', 'w'))
    return quantiles, np.column_stack([velocity, exponent]), identifiers


def test_import_design(tmp_path):
    quantiles, values, identifiers = write_campaign(tmp_path)
    design_file = ge.import_design(tmp_path)
    design = design_file.design
    assert design.names == [prior['path'] for prior in PRIORS]
    assert np.array_equal(design.quantiles, quantiles) and np.array_equal(design.values, values)
    assert [prior.class_name for prior in design.priors] == ['logNormal', 'uniform']
    # Parameters with log-normal priors are sampled in their logarithm by default.
    assert design.mappers == ['logarithm', 'identity']
    assert ge.import_design(tmp_path, log_parameters=[]).design.mappers == ['identity', 'identity']
    assert design_file.runs['outputFileName'][5].endswith(f'evaluations/{identifiers[5]}/galacticus.hdf5')
    # The converted priors reproduce the quantiles of the design.
    for j, prior in enumerate(design.priors):
        assert np.allclose(ge.prior_quantiles(prior.class_name, prior.parameters, values[:, j]), quantiles[:, j], rtol=0.0, atol=1.0e-10)
    # The design file round-trips, with descriptors that parse to the prior parameters.
    write_design(tmp_path / 'design.hdf5', design_file)
    with pytest.raises(FileExistsError):
        write_design(tmp_path / 'design.hdf5', design_file)
    read = read_design(tmp_path / 'design.hdf5')
    assert read.design.names == design.names and read.design.mappers == design.mappers
    assert np.array_equal(read.design.quantiles, design.quantiles) and read.count_runs == len(identifiers)
    assert parse_descriptor(read.design.priors[0].parameters['descriptor']) == {'x0': 150.0, 'sigma': 0.5, 'limitLower': 25.0, 'limitUpper': 300.0}
    assert read.runs['outputFileName'] == design_file.runs['outputFileName']


def test_import_design_inconsistent(tmp_path):
    write_campaign(tmp_path, perturb=True)
    with pytest.raises(ValueError, match="parameter 'velocity'"):
        ge.import_design(tmp_path)


def test_observable_from_galacticemu():
    logarithmic = ge.observable_from_galacticemu(
        {'analysis': 'massFunction', 'bad_training_condition': 'nonpositive', 'bad_training_value_fill': 'keep',
         'bad_training_sigma': 0.5, 'min_training_sigma': 1.0e-3}, y_axis_is_log=True)
    assert (logarithmic.transform, logarithmic.floor, logarithmic.root_variance_floored) == ('log10', -7.0, 0.5)
    assert (logarithmic.root_variance_zero, logarithmic.root_variance_target_minimum) == (1.0e-3, 1.0e-3)
    relation = ge.observable_from_galacticemu(
        {'analysis': 'relation', 'bad_training_condition': 'zero_only', 'bad_training_value_fill': 'bin_median',
         'bad_training_sigma': 5.0}, y_axis_is_log=False, conventions=False)
    assert (relation.transform, relation.undefined, relation.root_variance_undefined) == ('identity', 'median', 5.0)
    assert relation.root_variance_zero is None and relation.root_variance_target_minimum is None
    with pytest.raises(NotImplementedError):
        ge.observable_from_galacticemu({'analysis': 'x', 'bad_training_condition': 'nonfinite', 'bad_training_value_fill': 'bin_median'},
                                       y_axis_is_log=True)
