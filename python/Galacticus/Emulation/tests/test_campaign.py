"""Tests for `Galacticus.Emulation.design` and `Galacticus.Emulation.campaign`.

A stub executable stands in for Galacticus. It reads the run's change file and behaves according to a keyword placed
in it: `FAILEXIT` (exit with status 3), `MARKER` (print a failure marker but exit 0), `NOOUTPUT` (exit 0 without
writing the output file), or `BADBOUND` (fail, but only in a dry run). Otherwise it writes the output file named in the
change file and exits 0.
"""

import json
import os
import stat
import subprocess
import sys

import h5py
import numpy as np
import pytest

from Galacticus.Emulation import campaign as cm
from Galacticus.Emulation.design import parse_descriptor, read_design

COUNT_POINTS = 6

STUB = '''#!{python}
import re, sys
arguments = [a for a in sys.argv[1:] if a != '--dry-run']
dry_run = '--dry-run' in sys.argv
text = ''.join(open(f).read() for f in arguments[1:])
output = re.findall(r'<outputFileName value="([^"]+)"', text)[-1]
print('stub Galacticus: running', arguments)
if 'BADBOUND' in text and dry_run:
    print('Fatal error:')
    print('parameter [x] has value -1, but must be greater than 0.0')
    print(' Occurred at: somewhere')
    sys.exit(1)
if dry_run:
    sys.exit(0)
if 'FAILEXIT' in text:
    sys.exit(3)
if 'MARKER' in text:
    print('an aborted evolution')
if 'NOOUTPUT' not in text:
    open(output, 'w').write('model')
sys.exit(0)
'''


@pytest.fixture
def workspace(tmp_path, monkeypatch):
    """A working directory holding a design file, change files, a base parameter file and a stub executable."""
    monkeypatch.chdir(tmp_path)
    (tmp_path / 'changes').mkdir()
    (tmp_path / 'models').mkdir()
    (tmp_path / 'base.xml').write_text('<parameters/>\n')
    executable = tmp_path / 'stub.py'
    executable.write_text(STUB.format(python=sys.executable))
    executable.chmod(executable.stat().st_mode | stat.S_IEXEC)
    rng = np.random.default_rng(3)
    quantiles = rng.uniform(size=(COUNT_POINTS, 2))
    values = np.column_stack([1.0 + 9.0 * quantiles[:, 0], 10.0 ** (-3.0 + 2.0 * quantiles[:, 1])])
    change_files, output_files = [], []
    for i in range(COUNT_POINTS):
        change_file = f'changes/point{i:04d}.xml'
        output_file = f'models/point{i:04d}.hdf5'
        (tmp_path / change_file).write_text(
            '<changes>\n'
            f'  <change type="update" path="a/b" value="{values[i, 0]}"/>\n'
            f'  <change type="replaceOrAppend" path="outputFileName"><outputFileName value="{output_file}"/></change>\n'
            '</changes>\n')
        change_files.append(change_file)
        output_files.append(output_file)
    # Write a design file in the layout of the Fortran `emulatorDesign` task.
    with h5py.File(tmp_path / 'design.hdf5', 'w') as file:
        file.attrs['format'] = np.bytes_(b'galacticusDesign')
        file.attrs['formatVersion'] = np.int32(1)
        group = file.create_group('design')
        group.attrs['countPoints'] = np.int32(COUNT_POINTS)
        group.attrs['realizationsPerPoint'] = np.int32(1)
        group.attrs['seedPerPoint'] = np.int32(0)
        group.create_dataset('parameterNames', data=np.array([b'a/b', b"c/d[@value='e']/f"]))
        group.create_dataset('mappers', data=np.array([b'identity', b'logarithm']))
        group.create_dataset('quantiles', data=quantiles)
        group.create_dataset('values', data=values)
        priors = group.create_group('priors')
        prior = priors.create_group('prior1')
        prior.attrs['class'] = np.bytes_(b'uniform')
        prior.attrs['descriptor'] = np.bytes_(b'limitLower:1.0_limitUpper:10.0')
        prior = priors.create_group('prior2')
        prior.attrs['class'] = np.bytes_(b'logUniform')
        prior.attrs['descriptor'] = np.bytes_(b'limitLower:1.0e-3_limitUpper:1.0e-1')
        runs = file.create_group('runs')
        runs.create_dataset('pointIndex', data=np.arange(COUNT_POINTS, dtype=np.int32))
        runs.create_dataset('realizationIndex', data=np.zeros(COUNT_POINTS, dtype=np.int32))
        runs.create_dataset('changeFileName', data=np.array([name.encode() for name in change_files]))
        runs.create_dataset('outputFileName', data=np.array([name.encode() for name in output_files]))
    campaign = cm.Campaign.create('campaign.json', 'design.hdf5', 'base.xml', executable=str(executable))
    return tmp_path, campaign, values


def _mark(workspace_path, index, keyword):
    """Add a keyword to a run's change file, to direct the stub executable."""
    path = workspace_path / f'changes/point{index:04d}.xml'
    path.write_text(path.read_text().replace('</changes>', f'<!-- {keyword} -->\n</changes>'))


def test_parse_descriptor():
    parsed = parse_descriptor('x0:250.0_sigma:0.5_name:abc_randomNumberGenerator:GSL{seed:42_offset:false}_last:1e-3')
    assert parsed == {'x0': 250.0, 'sigma': 0.5, 'name': 'abc',
                      'randomNumberGenerator': ('GSL', {'seed': 42.0, 'offset': 'false'}), 'last': 1.0e-3}
    assert parse_descriptor('a:1{b:2{c:3}}_d:4') == {'a': (1.0, {'b': (2.0, {'c': 3.0})}), 'd': 4.0}


def test_read_design(workspace):
    path, _, values = workspace
    design = read_design(path / 'design.hdf5')
    assert design.count_runs == COUNT_POINTS
    assert design.design.names[1] == "c/d[@value='e']/f"
    assert design.design.priors[1].class_name == 'logUniform'
    assert design.design.priors[1].parameters['limitUpper'] == pytest.approx(0.1)
    np.testing.assert_array_equal(design.design.values, values)
    assert design.runs['changeFileName'][2] == 'changes/point0002.xml'


def test_manifest_round_trip(workspace):
    _, campaign, _ = workspace
    loaded = cm.Campaign.load('campaign.json')
    assert loaded == campaign
    assert [run.status for run in loaded.runs] == ['pending'] * COUNT_POINTS
    assert loaded.runs[3].log_file == os.path.join('logs', 'point0003.log')
    with pytest.raises(FileExistsError):
        cm.Campaign.create('campaign.json', 'design.hdf5', 'base.xml')


def test_run_local_statuses(workspace):
    path, campaign, _ = workspace
    _mark(path, 1, 'FAILEXIT')
    _mark(path, 2, 'MARKER')
    _mark(path, 3, 'NOOUTPUT')
    summary = campaign.run_local(parallel=3)
    assert summary == {'complete': 3, 'failed': 3}
    runs = cm.Campaign.load('campaign.json').runs
    assert runs[0].status == 'complete' and runs[0].exit_status == 0
    assert runs[1].status == 'failed' and runs[1].exit_status == 3 and 'status 3' in runs[1].message
    assert runs[2].status == 'failed' and 'aborted' in runs[2].message
    assert runs[3].status == 'failed' and 'output file is missing' in runs[3].message
    assert all(run.attempts == 1 for run in runs)
    # Running again retries only the failed runs.
    campaign.run_local()
    runs = cm.Campaign.load('campaign.json').runs
    assert [run.attempts for run in runs] == [1, 2, 2, 2, 1, 1]


def test_refresh_in_progress(workspace):
    path, campaign, _ = workspace
    run = campaign.runs[0]
    run.status = 'submitted'
    (path / 'logs').mkdir()
    (path / run.log_file).write_text('evolving trees...\n')
    campaign.refresh()
    assert campaign.runs[0].status == 'running'
    (path / run.output_file).write_text('model')
    with open(path / run.log_file, 'a') as handle:
        handle.write('EXIT_STATUS=0\n')
    campaign.refresh()
    assert campaign.runs[0].status == 'complete'


def test_array_specification():
    assert cm._array_specification([0, 1, 2, 3]) == '0-3'
    assert cm._array_specification([7, 0, 1, 9, 10, 3]) == '0-1,3,7,9-10'
    assert cm._array_specification([5]) == '5'


class FakeManager:
    """A stand-in for a queueManager manager, recording the job it is given."""

    def __init__(self):
        self.jobs = []

    def submit_detached(self, job):
        self.jobs.append(job)
        return '12345'


def test_submit_array_job(workspace):
    path, campaign, _ = workspace
    _mark(path, 4, 'FAILEXIT')
    campaign.run_local(indices=[4, 5])
    manager = FakeManager()
    # Submitting by default takes the pending and failed runs: all but run 5.
    job_id = campaign.submit(manager, job_options={'walltime': '1:00:00'})
    assert job_id == '12345'
    job = manager.jobs[0]
    assert job['array'] == '0-4'
    assert job['walltime'] == '1:00:00'
    runs = cm.Campaign.load('campaign.json').runs
    assert [run.status for run in runs] == ['submitted'] * 5 + ['complete']
    assert runs[4].attempts == 2 and runs[4].job_id == '12345'
    # The log of the earlier failed attempt is removed on resubmission.
    assert not (path / runs[4].log_file).exists()
    # Execute the array task for run 2 as Slurm would, then check that its status is picked up.
    environment = dict(os.environ, SLURM_ARRAY_TASK_ID='2')
    completed = subprocess.run(['bash', '-c', job['command']], env=environment, capture_output=True, text=True)
    assert completed.returncode == 0, completed.stderr
    summary = campaign.refresh()
    assert campaign.runs[2].status == 'complete'
    assert summary['submitted'] == 4


def test_validate(workspace):
    path, campaign, values = workspace
    extremes = campaign.extreme_runs()
    expected = {0} | {int(np.argmin(values[:, j])) for j in range(2)} | {int(np.argmax(values[:, j])) for j in range(2)}
    assert set(extremes) == expected
    results = campaign.validate(indices=[0, 1, 2])
    assert all(result.ok for result in results)
    _mark(path, 1, 'BADBOUND')
    results = {result.run: result for result in campaign.validate(indices=[0, 1, 2], parallel=2)}
    assert results[0].ok and results[2].ok
    assert not results[1].ok
    assert results[1].message == 'parameter [x] has value -1, but must be greater than 0.0'
    # Validation writes no model output and leaves the manifest unchanged.
    assert not any((path / 'models').iterdir())
    assert json.load(open('campaign.json'))['runs'][1]['status'] == 'pending'
