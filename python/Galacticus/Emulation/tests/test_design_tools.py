"""Tests for the splitting and pooling of designs (`Galacticus.Emulation.design`, `scripts/emulation/emulatorDesignTools.py`),
and the comparison of posteriors (`Galacticus.Emulation.chains`, `scripts/emulation/emulatorChainCompare.py`)."""

import os
import sys
import warnings

import h5py
import numpy as np
import pytest
from scipy import stats

from Galacticus.Emulation import chains as ch
from Galacticus.Emulation import design as ds
from Galacticus.Emulation.emulatorFile import Design, Prior
from Galacticus.Emulation.train import Observable, collect

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(ds.__file__)), '..', '..', '..', 'scripts', 'emulation'))


def test_prior_cumulative():
    x = np.linspace(0.5, 4.0, 15)
    cases = [
        (Prior('uniform', {'limitLower': 0.2, 'limitUpper': 5.0}), stats.uniform(0.2, 4.8).cdf(x)),
        (Prior('logUniform', {'limitLower': 0.2, 'limitUpper': 5.0}), stats.loguniform(0.2, 5.0).cdf(x)),
        (Prior('normal', {'mean': 2.0, 'variance': 0.64}), stats.norm(2.0, 0.8).cdf(x)),
        (Prior('normal', {'mean': 2.0, 'variance': 0.64, 'limitLower': 0.3, 'limitUpper': 4.5}),
         stats.truncnorm((0.3 - 2.0) / 0.8, (4.5 - 2.0) / 0.8, 2.0, 0.8).cdf(x)),
        (Prior('logNormal', {'x0': 1.5, 'sigma': 0.6, 'limitLower': 0.3, 'limitUpper': 6.0}),
         stats.truncnorm((np.log(0.3) - np.log(1.5)) / 0.6, (np.log(6.0) - np.log(1.5)) / 0.6, np.log(1.5), 0.6).cdf(np.log(x))),
    ]
    # A log-normal given by its mean and variance (of x): the variance of ln x is ln(1 + variance/mean^2).
    mean, variance = 2.0, 1.5
    sigma = np.sqrt(np.log(1.0 + variance / mean ** 2))
    cases.append((Prior('logNormal', {'mean': mean, 'variance': variance}), stats.lognorm(sigma, scale=mean * np.exp(-0.5 * sigma ** 2)).cdf(x)))
    for prior, expected in cases:
        assert np.allclose(ds.prior_cumulative(prior, x), expected, rtol=0.0, atol=1e-13), prior
    with pytest.raises(NotImplementedError):
        ds.prior_cumulative(Prior('beta', {}), x)


def _design_file(values, prior, seeds, directory, label):
    """A design of the given values (one parameter), each run once, with relative run file names."""
    count = values.shape[0]
    design = Design(names=['a'], priors=[prior], mappers=['identity'], quantiles=ds.prior_cumulative(prior, values)[:, None],
                    values=values[:, None], attributes={'countPoints': np.int32(count), 'seedPerPoint': np.int32(1)})
    runs = {'pointIndex': np.arange(count, dtype=np.int32), 'realizationIndex': np.zeros(count, dtype=np.int32),
            'changeFileName': [f'{label}/changes/point{i:04d}.xml' for i in range(count)],
            'outputFileName': [f'{label}/models/point{i:04d}.hdf5' for i in range(count)],
            'seed': np.array(seeds, dtype=np.int64)}
    return ds.DesignFile(design=design, runs=runs)


def _write_output(path, value):
    """A model output holding one analysis, whose value is linear in the parameter."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with h5py.File(path, 'w') as file:
        group = file.create_group('analyses/linear')
        group.attrs['type'] = np.bytes_(b'function1D')
        group.attrs['xDataset'] = np.bytes_(b'x')
        group.attrs['yDataset'] = np.bytes_(b'y')
        group['x'] = np.array([0.0, 1.0])
        group['y'] = np.array([value, 2.0 * value])


def test_split_and_pool(tmp_path):
    priorFirst = Prior('uniform', {'limitLower': 0.0, 'limitUpper': 1.0})
    first = _design_file(np.linspace(0.05, 0.95, 10), priorFirst, np.arange(10), tmp_path, 'first')
    # Parts: even and odd points, with their runs, renumbered.
    parts = ds.split_design(first, 2)
    assert [part.design.values[:, 0].tolist() for part in parts] == [first.design.values[0::2, 0].tolist(), first.design.values[1::2, 0].tolist()]
    assert parts[1].runs['pointIndex'].tolist() == list(range(5))
    assert parts[1].runs['outputFileName'] == first.runs['outputFileName'][1::2]
    assert int(parts[1].design.attributes['splitPart']) == 1
    with pytest.raises(ValueError):
        ds.split_design(first, 1)
    # Pool with a second design over a wider box: its points lie outside the prior of the first design.
    second = _design_file(np.linspace(0.6, 1.9, 6), Prior('uniform', {'limitLower': 0.5, 'limitUpper': 2.0}), np.arange(100, 106),
                          tmp_path, 'second')
    with pytest.raises(ValueError, match='outside the support'):
        ds.pool_designs([first, second])
    union = Prior('uniform', {'limitLower': 0.0, 'limitUpper': 2.0})
    pooled = ds.pool_designs([first, second], [union])
    assert pooled.design.quantiles.shape == (16, 1)
    assert np.allclose(pooled.design.quantiles[:, 0], pooled.design.values[:, 0] / 2.0)
    assert pooled.runs['pointIndex'].tolist() == list(range(16))
    # Shared seeds make the noise of the pooled runs correlated: warned.
    with pytest.warns(UserWarning, match='share a random number seed'):
        ds.pool_designs([first, _design_file(np.linspace(0.6, 1.9, 6), union, np.arange(6), tmp_path, 'third')], [union])
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        ds.pool_designs([first, second], [union])
    # A value at or below zero has no quantile under a log prior: it is outside the support.
    negative = _design_file(np.array([-0.5, 0.5]), Prior('uniform', {'limitLower': -1.0, 'limitUpper': 1.0}), [200, 201], tmp_path, 'negative')
    with pytest.raises(ValueError, match='outside the support'):
        ds.pool_designs([first, negative], [Prior('logUniform', {'limitLower': 1.0e-3, 'limitUpper': 2.0})])
    # Designs of different parameters can not be pooled.
    other = _design_file(np.linspace(0.1, 0.9, 4), priorFirst, np.arange(4), tmp_path, 'other')
    other.design.names = ['b']
    with pytest.raises(ValueError, match='same parameters'):
        ds.pool_designs([first, other])


def test_command_line(tmp_path, monkeypatch):
    import emulatorDesignTools
    monkeypatch.chdir(tmp_path)
    prior = Prior('uniform', {'limitLower': 0.0, 'limitUpper': 1.0, 'descriptor': 'limitLower:0.0_limitUpper:1.0'})
    union = Prior('uniform', {'limitLower': 0.0, 'limitUpper': 2.0, 'descriptor': 'limitLower:0.0_limitUpper:2.0'})
    designs = {'first': _design_file(np.linspace(0.05, 0.95, 8), prior, np.arange(8), tmp_path, 'first'),
               'second': _design_file(np.linspace(1.05, 1.95, 8), union, np.arange(100, 108), tmp_path, 'second')}
    for label, content in designs.items():
        # The runs of each design were made from its own directory.
        os.makedirs(label)
        ds.write_design(f'{label}/design.hdf5', content)
        for path, value in zip(content.runs['outputFileName'], content.design.values[:, 0]):
            _write_output(os.path.join(label, path), value)
    # Split the first design; the run file names are made absolute.
    assert emulatorDesignTools.main(['split', 'first/design.hdf5', '--output', 'half', '--directory', 'first']) == 0
    half = ds.read_design('half1.hdf5')
    assert half.design.values[:, 0].tolist() == designs['first'].design.values[1::2, 0].tolist()
    assert all(os.path.isabs(path) and os.path.exists(path) for path in half.runs['outputFileName'])
    # A training set can be collected from a part, from any directory.
    _, training_sets, _ = collect('half1.hdf5', [Observable('linear')])
    assert np.allclose(training_sets['linear'].y[:, 0], half.design.values[:, 0])
    # Pool the two designs under the union prior, read from the second design.
    assert emulatorDesignTools.main(['pool', 'pooled.hdf5', 'first/design.hdf5', 'second/design.hdf5', '--directories', 'first', 'second',
                                     '--priors', 'second/design.hdf5']) == 0
    pooled = ds.read_design('pooled.hdf5')
    assert pooled.design.quantiles.shape == (16, 1)
    assert pooled.design.priors[0].parameters['limitUpper'] == 2.0
    _, training_sets, reports = collect('pooled.hdf5', [Observable('linear')])
    assert reports['linear'].points_missing == [] and np.allclose(training_sets['linear'].y[:, 0], pooled.design.values[:, 0])
    # Priors may instead come from the `modelParameter` blocks of a parameter file.
    (tmp_path / 'priors.xml').write_text('''<parameters>
  <posteriorSampleSimulation value="differentialEvolution">
    <modelParameter value="active">
      <name value="a"/>
      <distributionFunction1DPrior value="uniform">
        <limitLower value="0.0"/>
        <limitUpper value="2.0"/>
      </distributionFunction1DPrior>
      <operatorUnaryMapper value="identity"/>
      <distributionFunction1DPerturber value="cauchy">
        <median value="0.0"/>
        <scale value="1.0e-3"/>
      </distributionFunction1DPerturber>
    </modelParameter>
  </posteriorSampleSimulation>
</parameters>
''')
    assert emulatorDesignTools.main(['pool', 'pooled2.hdf5', 'first/design.hdf5', 'second/design.hdf5', '--directories', 'first', 'second',
                                     '--priors', 'priors.xml']) == 0
    assert np.array_equal(ds.read_design('pooled2.hdf5').design.quantiles, pooled.design.quantiles)


def test_priors_from_parameter_file():
    """Priors are read from the `modelParameter` blocks of a real parameter file (that of the design test)."""
    path = os.path.join(os.path.dirname(os.path.abspath(ds.__file__)), '..', '..', '..', 'testSuite', 'parameters', 'emulatorDesign.xml')
    names = ['hotHaloOutflowReincorporation/multiplier', 'blackHoleWind/efficiencyWind']
    priors = ds.priors_from_file(path, names)
    assert [prior.class_name for prior in priors] == ['uniform', 'logUniform']
    assert priors[1].parameters['limitLower'] == 2.4e-4 and priors[1].parameters['limitUpper'] == 2.4e-2
    with pytest.raises(ValueError, match='gives no prior'):
        ds.priors_from_file(path, ['noSuchParameter'])


def _write_chains(root, samples, converged_from):
    """Write chain logs in the format of Galacticus posterior sampling simulations: one file per chain."""
    for chain, states in enumerate(samples):
        with open(f'{root}_{chain:04d}.log', 'w') as handle:
            handle.write('# Simulation state chain file\n# Columns:\n#    1 = Simulation step\n#    2 = Chain index\n'
                         '#    3 = Evaluation time (s)\n#    4 = Chain is converged? [T/F]\n#    5 = log posterior\n#    6 = log likelihood\n'
                         '#    7 = Parameter `a`\n#    8 = Parameter `b`\n')
            for step, state in enumerate(states):
                flag = 'T' if step >= converged_from else 'F'
                handle.write(f'{step + 1:12d} {chain:11d}   1.0E-03 {flag}  -1.0  -1.0  {state[0]:.17g}  {state[1]:.17g}\n')


def test_compare_posteriors(tmp_path, capsys):
    import emulatorChainCompare
    rng = np.random.default_rng(6)
    # The reference posterior, and two posteriors shifted in `a` by +0.3 and -0.1 of its width. Unconverged states lie far away.
    roots = {}
    for label, shift in (('full', 0.0), ('half0', 0.3), ('half1', -0.1)):
        samples = []
        for _ in range(4):
            states = np.column_stack([rng.normal(1.0 + shift, 1.0, 20000), rng.normal(5.0, 0.1, 20000)])
            states[:500] = 100.0
            samples.append(states)
        roots[label] = str(tmp_path / label)
        _write_chains(roots[label], samples, converged_from=500)
    names, states = ch.read_chains(roots['full'])
    assert names == ['a', 'b'] and states.shape == (4 * 19500, 2)
    result = ch.compare_posteriors(ch.read_chains(roots['full']), [ch.read_chains(roots['half0']), ch.read_chains(roots['half1'])])
    assert result[0]['name'] == 'a' and result[0]['sigma'] == pytest.approx(1.0, rel=0.03)
    assert result[0]['shifts'] == pytest.approx([0.3, -0.1], abs=0.03)
    assert result[0]['ratio'] == pytest.approx(0.3, abs=0.03)
    assert result[1]['ratio'] < 0.1
    # Ignoring the convergence flags (with no burn-in), the unconverged states are included.
    assert ch.read_chains(roots['full'], converged_only=False, burn_fraction=0.0)[1].shape == (4 * 20000, 2)
    assert emulatorChainCompare.main([roots['full'], roots['half0'], roots['half1']]) == 0
    assert 'systematic' in capsys.readouterr().out
    with pytest.raises(FileNotFoundError):
        ch.read_chains(str(tmp_path / 'missing'))
    # A root which is a prefix of another (`full` and `full_extra`) reads only its own chains.
    _write_chains(str(tmp_path / 'full_extra'), [np.full((10, 2), 50.0)], converged_from=0)
    assert ch.read_chains(roots['full'])[1].shape == (4 * 19500, 2)
