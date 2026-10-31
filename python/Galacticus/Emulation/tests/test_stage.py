"""Tests for `Galacticus.Emulation.stage` and `scripts/emulation/emulatorStageWorker.py`.

A stub executable stands in for Galacticus. It writes one analysis (`obs`, a smooth function of the design parameter `a`),
and reads keywords placed in a run's change file: `FAILEXIT` (exit with status 3) and `SLEEP` (sleep for fifteen seconds
first). Emulators are trained from the stub's outputs by the real training script.
"""

import os
import signal
import stat
import subprocess
import sys
import time

import numpy as np
import pytest

from Galacticus.Emulation import design as ds
from Galacticus.Emulation import stage as sg
from Galacticus.Emulation.campaign import Campaign
from Galacticus.Emulation.emulatorFile import Design, Prior

STUB = '''#!{python}
import re, sys, time
import numpy as np, h5py
text = ''.join(open(f).read() for f in sys.argv[2:])
a = float(re.findall(r'path="a" value="([^"]+)"', text)[-1])
output = re.findall(r'<outputFileName value="([^"]+)"', text)[-1]
if 'SLEEP' in text:
    time.sleep(15)
if 'FAILEXIT' in text:
    sys.exit(3)
x = np.linspace(0.0, 1.0, 5)
with h5py.File(output, 'w') as file:
    group = file.create_group('analyses/obs')
    group.attrs['type'] = np.bytes_(b'function1D')
    group.attrs['xDataset'] = np.bytes_(b'x')
    group.attrs['yDataset'] = np.bytes_(b'y')
    group.attrs['yCovariance'] = np.bytes_(b'yCovariance')
    group['x'] = x
    group['y'] = np.sin(3.0 * a + x)
    group['yCovariance'] = 1.0e-4 * np.eye(x.size)
print('stub Galacticus: done')
'''

COUNT_POINTS = 8


def _component(root, name, executable, keywords=None):
    """Write the design, change files, campaign manifest and training configurations of one component."""
    directory = root / name
    (directory / 'changes').mkdir(parents=True)
    (directory / 'models').mkdir()
    values = np.linspace(0.05, 0.95, COUNT_POINTS)
    runs = {'pointIndex': np.arange(COUNT_POINTS, dtype=np.int32), 'realizationIndex': np.zeros(COUNT_POINTS, dtype=np.int32),
            'changeFileName': [], 'outputFileName': []}
    for i, value in enumerate(values):
        change = directory / 'changes' / f'point{i:04d}.xml'
        output = directory / 'models' / f'point{i:04d}.hdf5'
        keyword = (keywords or {}).get(i, '')
        change.write_text(f'<changes>\n  <change type="update" path="a" value="{float(value)!r}"/>\n'
                          f'  <change type="replaceOrAppend" path="outputFileName"><outputFileName value="{output}"/></change>\n'
                          f'  <!-- {keyword} -->\n</changes>\n')
        runs['changeFileName'].append(str(change))
        runs['outputFileName'].append(str(output))
    prior = Prior('uniform', {'limitLower': 0.0, 'limitUpper': 1.0, 'descriptor': 'limitLower:0.0_limitUpper:1.0'})
    design = ds.DesignFile(design=Design(names=['a'], priors=[prior], mappers=['identity'], quantiles=values[:, None],
                                         values=values[:, None], attributes={'countPoints': np.int32(COUNT_POINTS)}), runs=runs)
    ds.write_design(str(directory / 'design.hdf5'), design)
    for part, content in enumerate(ds.split_design(design)):
        ds.write_design(str(directory / f'half{part}.hdf5'), content)
    (root / 'base.xml').write_text('<parameters/>\n')
    Campaign.create(str(directory / 'campaign.json'), str(directory / 'design.hdf5'), str(root / 'base.xml'), executable=str(executable),
                    working_directory=str(directory))
    trainings = []
    for label, design_file in (('full', 'design.hdf5'), ('half0', 'half0.hdf5'), ('half1', 'half1.hdf5')):
        configuration = directory / f'train_{label}.xml'
        configuration.write_text(f'''<emulatorTrain>
  <designFileName       value="{directory / design_file}"/>
  <emulatorFileName     value="{directory / f'emulator_{label}.hdf5'}"/>
  <restartsOptimizer    value="1"/>
  <foldsCrossValidation value="0"/>
  <observable label="obs"/>
</emulatorTrain>
''')
        trainings.append(sg.Training(label, str(configuration), str(directory / f'emulator_{label}.hdf5')))
    return sg.StageComponent(name=name, manifest=str(directory / 'campaign.json'), trainings=trainings, cost_mean=10.0, cost_maximum=20.0)


@pytest.fixture
def workspace(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv('GALACTICUS_EXEC_PATH', os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(sg.__file__)), '..', '..', '..')))
    for name in ('SLURM_JOB_ID', 'SLURM_ARRAY_JOB_ID', 'SLURM_ARRAY_TASK_ID'):
        monkeypatch.delenv(name, raising=False)
    executable = tmp_path / 'stub.py'
    executable.write_text(STUB.format(python=sys.executable))
    executable.chmod(executable.stat().st_mode | stat.S_IEXEC)
    return tmp_path, executable


def _messages():
    messages = []
    return messages, messages.append


def test_stage_round_trip(workspace):
    root, executable = workspace
    components = [_component(root, 'first', executable)]
    stage = sg.Stage.create(str(root / 'stage.json'), components, training_processes=2)
    loaded = sg.Stage.load(str(root / 'stage.json'))
    assert loaded.components == stage.components and loaded.training_processes == 2
    with pytest.raises(FileExistsError):
        sg.Stage.create(str(root / 'stage.json'), components)
    with pytest.raises(ValueError):
        sg.Stage.create(str(root / 'other.json'), components + components)


def test_work(workspace):
    root, executable = workspace
    components = [_component(root, 'first', executable, keywords={3: 'FAILEXIT'}), _component(root, 'second', executable)]
    stage = sg.Stage.create(str(root / 'stage.json'), components, training_processes=1)
    assert len(stage.check()) == 6
    messages, log = _messages()
    stage.work(slots=4, log=log)
    # Every run ended, and every emulator was trained (the failed run is left out of training).
    assert stage.check() == []
    assert any('first run 3 failed' in message for message in messages)
    first = Campaign.load(components[0].manifest)
    assert [run.status for run in first.runs].count('complete') == COUNT_POINTS - 1 and first.runs[3].status == 'failed'
    assert Campaign.load(components[1].manifest).summary() == {'complete': COUNT_POINTS}
    from Galacticus.Emulation import emulatorFile as ef
    content = ef.read(components[0].trainings[0].emulator)
    assert content.training_sets['obs'].point_index.tolist() == [0, 1, 2, 4, 5, 6, 7]
    # Planning reports the failed run, and leaves the emulators trained.
    summary = stage.plan()
    assert summary['runs'] == 0 and summary['trainings'] == 0 and [entry[:2] for entry in summary['failed']] == [('first', 3)]
    assert stage.check() == []
    # Retrying failed runs sets the emulators of their component aside, for retraining.
    summary = stage.plan(retry_failed=True)
    assert summary['runs'] == 1 and summary['trainings'] == 3
    assert os.path.exists(components[0].trainings[0].emulator + '.previous')
    assert stage.check() == [('first', 'full'), ('first', 'half0'), ('first', 'half1')]


def test_killed_worker(workspace):
    root, executable = workspace
    components = [_component(root, 'first', executable, keywords={0: 'SLEEP'})]
    stage = sg.Stage.create(str(root / 'stage.json'), components)
    # Run a worker in another process, and kill it while its first (slow) run is in progress.
    worker = subprocess.Popen([sys.executable, os.path.join(os.environ['GALACTICUS_EXEC_PATH'], 'scripts', 'emulation', 'emulatorStageWorker.py'),
                               'work', str(root / 'stage.json'), '--slots', '1'], env=dict(os.environ, PYTHONPATH=os.pathsep.join(sys.path)),
                              stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    process = stage.state(components[0], 'run000000.process')
    for _ in range(100):
        if os.path.exists(process):
            break
        time.sleep(0.1)
    worker.send_signal(signal.SIGKILL)
    worker.wait()
    # Its model run is orphaned, and still running: a new worker does the other runs, but does not take that one over.
    messages, log = _messages()
    stage.work(slots=2, log=log)
    assert not stage.run_ended(components[0], 0)
    # No emulator of the component is trained until all its runs have ended.
    assert stage.check() == [('first', 'full'), ('first', 'half0'), ('first', 'half1')]
    assert all(stage.run_ended(components[0], index) for index in range(1, COUNT_POINTS))
    # Once the orphaned run has ended (without recording its outcome), the next worker takes it over.
    pid = int(open(process).read().split(':')[1])
    for _ in range(300):
        try:
            os.kill(pid, 0)
        except ProcessLookupError:
            break
        time.sleep(0.1)
    messages, log = _messages()
    stage.work(slots=2, log=log)
    assert any('took over run000000.claim' in message for message in messages)
    assert stage.check() == []
    assert Campaign.load(components[0].manifest).summary() == {'complete': COUNT_POINTS}


def test_command_line(workspace, capsys):
    sys.path.insert(0, os.path.join(os.environ['GALACTICUS_EXEC_PATH'], 'scripts', 'emulation'))
    import emulatorStageWorker
    root, executable = workspace
    sg.Stage.create(str(root / 'stage.json'), [_component(root, 'first', executable)])
    assert emulatorStageWorker.main(['check', str(root / 'stage.json')]) == 1
    assert emulatorStageWorker.main(['plan', str(root / 'stage.json')]) == 0
    assert f'remaining: {COUNT_POINTS} model runs' in capsys.readouterr().out
    assert emulatorStageWorker.main(['work', str(root / 'stage.json'), '--slots', '4']) == 0
    assert emulatorStageWorker.main(['check', str(root / 'stage.json')]) == 0
    assert 'all 3 emulators are trained' in capsys.readouterr().out


def test_holder_ended():
    import socket
    host = socket.gethostname()
    # A Slurm array task is ended once it is no longer listed; nothing is known to have ended if the listing failed.
    assert sg._holder_ended(f'123_4:{host}:1', {'123', '123_5'})
    assert not sg._holder_ended(f'123_4:{host}:1', {'123', '123_4'})
    assert not sg._holder_ended(f'123_4:{host}:1', None)
    # A local worker is ended once its process has ended, which can only be checked on its own host.
    process = subprocess.Popen([sys.executable, '-c', 'pass'])
    process.wait()
    assert sg._holder_ended(f'local:{host}:{process.pid}', None)
    assert not sg._holder_ended(f'local:{host}:{os.getpid()}', None)
    assert not sg._holder_ended(f'local:another-host:{process.pid}', None)
    assert not sg._holder_ended('malformed', set())


def test_plan_keeps_orphans(workspace):
    import socket
    root, executable = workspace
    components = [_component(root, 'first', executable)]
    stage = sg.Stage.create(str(root / 'stage.json'), components)
    os.makedirs(stage.state(components[0]))
    # Records of model processes left by killed workers: one still running (this process), one ended.
    ended = subprocess.Popen([sys.executable, '-c', 'pass'])
    ended.wait()
    for index, pid in ((0, os.getpid()), (1, ended.pid)):
        sg._write(stage.state(components[0], f'run{index:06d}.process'), f'{socket.gethostname()}:{pid}')
        os.symlink(f'local:{socket.gethostname()}:{pid}', stage.state(components[0], f'run{index:06d}.claim'))
    # Planning clears the queue, but keeps the record of the process still running, so that no worker starts its run.
    summary = stage.plan()
    assert summary['orphaned'] == 1 and summary['runs'] == COUNT_POINTS
    assert sorted(os.listdir(stage.state(components[0]))) == ['run000000.process']
    messages, log = _messages()
    stage.work(slots=2, log=log)
    assert not stage.run_ended(components[0], 0)
    assert all(stage.run_ended(components[0], index) for index in range(1, COUNT_POINTS))


def test_live_jobs(monkeypatch):
    # A listing of the live jobs is trusted only if it includes this worker's own job.
    monkeypatch.setenv('SLURM_JOB_ID', '200')
    monkeypatch.setenv('SLURM_ARRAY_JOB_ID', '123')
    monkeypatch.setenv('SLURM_ARRAY_TASK_ID', '4')
    for listing, expected in ((set(), None), ({'123', '123_5'}, None), (None, None), ({'123', '123_4', '99'}, {'123', '123_4', '99'})):
        monkeypatch.setattr(sg, 'slurm_live_jobs', lambda listing=listing: listing)
        assert sg._live_jobs() == expected
    monkeypatch.delenv('SLURM_JOB_ID')
    monkeypatch.setattr(sg, 'slurm_live_jobs', lambda: {'1'})
    assert sg._live_jobs() is None
