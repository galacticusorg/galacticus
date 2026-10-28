"""Tests for `Galacticus.Emulation.profile` and `scripts/emulation/emulatorProfile.py`.

A stub executable stands in for Galacticus. It writes one counts-based analysis, whose model value in each bin depends on two
parameters (`a`, an amplitude, and `b`, a slope) exactly, times a noise factor fixed by the random number seed (as merger
trees fixed by the seed would give). Its target counts are drawn once from the noise-free model at `a=1`, `b=0.5`.
"""

import os
import stat
import sys

import h5py
import numpy as np
import pytest

from Galacticus.Emulation import profile as pf
from Galacticus.Emulation.train import Observable

STUB = '''#!{python}
import re, sys
import numpy as np, h5py
text = ''.join(open(f).read() for f in sys.argv[1:])
def last(pattern):
    return re.findall(pattern, text)[-1]
values = dict(re.findall(r'path="(\\w+)" value="([^"]+)"', text))
base = dict(re.findall(r'<(\\w+) value="([^"]+)"/>', open(sys.argv[1]).read()))
a = float(values.get('a', base['a']))
b = float(values.get('b', base['b']))
seed = int(last(r'<seed value="(\\d+)"'))
output = last(r'<outputFileName value="([^"]+)"')
x = np.linspace(0.0, 2.0, 12)
countPerValue = np.full(x.size, 400.0)
truth = np.exp(-0.5 * x)
countTarget = np.random.default_rng(0).poisson(truth * countPerValue).astype(float)
noise = 1.0 + 0.03 * np.random.default_rng(seed).normal(size=x.size)
y = a * np.exp(-b * x) * noise
with h5py.File(output, 'w') as file:
    group = file.create_group('analyses/counts')
    group.attrs['type'] = np.bytes_(b'function1D')
    group.attrs['xDataset'] = np.bytes_(b'x')
    group.attrs['yDataset'] = np.bytes_(b'y')
    group['x'] = x
    group['y'] = y
    group['countTarget'] = countTarget
    group['countPerValue'] = countPerValue
print('stub Galacticus: done')
'''


@pytest.fixture
def workspace(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    executable = tmp_path / 'stub.py'
    executable.write_text(STUB.format(python=sys.executable))
    executable.chmod(executable.stat().st_mode | stat.S_IEXEC)
    (tmp_path / 'base.xml').write_text('<parameters>\n  <a value="1.0"/>\n  <b value="0.5"/>\n'
                                       '  <randomNumberGenerator value="GSL">\n    <seed value="1"/>\n  </randomNumberGenerator>\n'
                                       '  <outputFileName value="galacticus.hdf5"/>\n</parameters>\n')
    (tmp_path / 'profile.xml').write_text(f'''<emulatorProfile>
  <directory     value="profile"/>
  <executable    value="{executable}"/>
  <lineSeeds     value="11 12 13 14"/>
  <lineFractions value="0 0.5 1"/>
  <rootVarianceTargetFractional value="profile" minimum="1.0e-3" maximum="1.0"/>
  <parameter name="a" center="1.05" step="0.05" scale="log"/>
  <parameter name="b" center="0.48" step="0.02"/>
  <component name="stub" baseParameters="base.xml">
    <observable label="co.*" pattern="true" likelihoodForm="poissonCounts" transform="log10"/>
  </component>
</emulatorProfile>
''')
    return tmp_path


def test_star_offsets():
    offsets = pf.star_offsets(3)
    assert offsets.shape == (1 + 2 * 3 + 2 * 3 * 2, 3)
    assert np.all(offsets[0] == 0.0)
    assert len({tuple(offset) for offset in offsets}) == offsets.shape[0]


def test_fit_quadratic():
    rng = np.random.default_rng(5)
    hessian = np.array([[-4.0, 1.0, 0.0], [1.0, -3.0, 0.5], [0.0, 0.5, -2.0]])
    gradient = np.array([2.0, -1.0, 0.5])
    offsets = pf.star_offsets(3)
    values = 7.0 + offsets @ gradient + 0.5 * np.einsum('ij,jk,ik->i', offsets, hessian, offsets)
    fit = pf.fit_quadratic(offsets, values)
    assert np.allclose(fit['hessian'], hessian) and np.allclose(fit['gradient'], gradient)
    assert np.allclose(fit['optimum'], np.linalg.solve(-hessian, gradient), atol=1e-6)
    assert np.allclose(fit['covariance'], np.linalg.inv(-hessian))
    # The optimum is bounded.
    fit = pf.fit_quadratic(offsets, 10.0 * offsets[:, 0] - 0.1 * offsets[:, 0] ** 2 + rng.normal(0.0, 1e-9, offsets.shape[0]), offset_maximum=2.0)
    assert fit['optimum'][0] == pytest.approx(2.0)


def test_profile_command_line(workspace, capsys):
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(pf.__file__)), '..', '..', '..', 'scripts', 'emulation'))
    import emulatorProfile
    configuration = str(workspace / 'profile.xml')
    profile = emulatorProfile.read_configuration(configuration)
    # The line search can not be prepared before the star stage is analyzed.
    with pytest.raises(RuntimeError):
        profile.prepare('line')
    assert emulatorProfile.main(['prepare', configuration, 'star']) == 0
    assert emulatorProfile.main(['analyze', configuration, 'star']) == 1  # Runs not yet made.
    assert emulatorProfile.main(['run', configuration, 'star', '--parallel', '4']) == 0
    assert emulatorProfile.main(['analyze', configuration, 'star']) == 0
    star = profile.read_result('star')
    assert star['pointsUsed'] == 9
    # At one seed the likelihood is a smooth function of the parameters (the trees, here the noise, are fixed), so a quadratic
    # fits it well, with a negative-definite curvature near the truth (a=1, b=0.5). Its optimum is shifted by the noise of that
    # seed: here it overshoots in `a` to the bound of two steps.
    assert star['fit']['residualRMS'] < 0.1 * star['gainPredicted']
    assert star['fit']['covariance'] is not None
    assert star['fit']['optimum'][0] == pytest.approx(-2.0)
    # The direct log-likelihood at the center is that of the model output there, with C maximized.
    from Galacticus.Emulation import likelihood as lk
    from Galacticus.Emulation.train import collect
    _, training_sets, _ = collect(profile.path('star', 'stub', 'design.hdf5'), [Observable('counts', transform='log10')])
    C = star['rootVarianceTargetFractional'][0]
    assert lk.log_likelihood(training_sets['counts'], 'poissonCounts', training_sets['counts'].y[0],
                             root_variance_target_fractional=C) == pytest.approx(star['logLikelihood'][0], rel=1e-10)
    # The line search, at four other seeds.
    assert emulatorProfile.main(['prepare', configuration, 'line']) == 0
    assert emulatorProfile.main(['run', configuration, 'line', '--parallel', '4']) == 0
    assert emulatorProfile.main(['analyze', configuration, 'line']) == 0
    line = profile.read_result('line')
    assert line['countsSeeds'] == [1, 2, 4]
    # Free of the noise of the star's seed, the line search stops half way, much nearer the truth than the center.
    assert line['best'] == 1 and line['gain'] > 0.0
    assert line['optimum']['a'] == pytest.approx(1.0, abs=0.01)
    assert line['optimum']['b'] == pytest.approx(0.5, abs=0.035)
    # Averaging over more seeds reduces the noise, which raises the log-likelihood.
    means = np.array(line['logLikelihoodMean'])
    assert np.all(means[2] > means[0])
    assert all(penalty < 0.0 for penalty in line['noisePenaltySingleSeed'])
    output = capsys.readouterr().out
    assert 'N=infinity' in output
    # Fix the parameters at the optimum of the line search, in a change file and in a copy of the base parameter file.
    assert emulatorProfile.main(['fix', configuration, 'line', '--output', 'fixed.xml', '--baseDirectory', 'fixed']) == 0
    import lxml.etree as ET
    fixed = ET.parse(str(workspace / 'fixed' / 'base.xml')).getroot()
    for name in ('a', 'b'):
        assert float(fixed.find(name).get('value')) == pytest.approx(line['optimum'][name], rel=1e-14)


def test_observable_patterns(tmp_path):
    """Label patterns are matched in the first output which has analyses: a run which failed part way leaves none."""
    from types import SimpleNamespace
    outputs = [tmp_path / 'failed.hdf5', tmp_path / 'complete.hdf5']
    with h5py.File(outputs[0], 'w') as file:
        file.create_group('Parameters')
    with h5py.File(outputs[1], 'w') as file:
        file.create_group('analyses/countsA')
        file.create_group('analyses/other')
    with h5py.File(tmp_path / 'design.hdf5', 'w') as file:
        file['runs/outputFileName'] = np.array([str(output).encode() for output in outputs])
    component = SimpleNamespace(name='stub', observables=[SimpleNamespace(pattern=True, observable=Observable('co.*'), form='poissonCounts')])
    observables = pf.Profile._observables(None, component, str(tmp_path / 'design.hdf5'))
    assert [(observable.label, form) for observable, form in observables] == [('countsA', 'poissonCounts')]
