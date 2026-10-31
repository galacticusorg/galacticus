"""Run the model runs and emulator trainings of a stage from a shared work queue.

A *stage* is a set of *components*, each a campaign of model runs (:mod:`Galacticus.Emulation.campaign`) with one or
more emulators to train from its runs (e.g. one on the whole design, and one on each half; see
:func:`Galacticus.Emulation.design.split_design`). A calibration against many simulations (the progenitor mass functions
of tens of N-body simulations, say) is such a stage: tens of thousands of single-threaded runs of very different costs.
Rather than a Slurm array task per run, any number of *worker* jobs (on any nodes, starting at any time) share the work:

* Each worker takes work by creating a *claim* - an atomic symbolic link naming the worker - in
  ``<stage>/state/<component>/``: ``runNNNNNN.claim`` for a model run, ``train-<name>.claim`` for a training. When the
  work ends, the worker writes ``runNNNNNN.complete`` or ``runNNNNNN.failed`` (holding the reason), or
  ``train-<name>.complete`` or ``train-<name>.failed``.
* Runs are taken slowest component first (by the costs given for the components), single-threaded, filling the worker's
  slots. A component's emulators are trained as soon as all its runs have ended (training takes precedence over new
  runs, since it is on the critical path to sampling the posterior); runs that failed are left out of training.
* A worker in a Slurm job does not start a run that could not finish before its time limit, and takes over the claims of
  workers whose jobs Slurm no longer lists (killed, or out of time). A worker outside Slurm takes over the claims of
  workers on the same host whose processes have ended. A worker exits when nothing is left that it can take.
* The campaign manifests are read, but not written, while workers run, except by the worker training a component, which
  records that component's runs in its manifest.

Before (re)starting workers, with none running, :meth:`Stage.plan` brings the manifests up to date, clears the work queue,
and reports the work remaining; runs that failed are retried only if asked. :meth:`Stage.check` reports whether every
emulator is trained. The command-line interface is ``scripts/emulation/emulatorStageWorker.py``.

Andrew Benson, Claude (28-October-2026).
"""

from __future__ import annotations

import json
import os
import re
import socket
import subprocess
import sys
import time
from dataclasses import asdict, dataclass, field

from Galacticus.Emulation.campaign import (COMPLETE, EXIT_SENTINEL, FAILED, PENDING, SUBMITTED, Campaign, Run,
                                           run_status_from_log, slurm_live_jobs)

__all__ = [
    'STAGE_FORMAT',
    'STAGE_VERSION',
    'Training',
    'StageComponent',
    'Stage',
    'owner_identity',
]

STAGE_FORMAT = 'galacticusStage'
STAGE_VERSION = 1


@dataclass
class Training:
    """An emulator to train from a component's runs: a name (unique within the component), and the training configuration
    (see ``scripts/emulation/emulatorTrainReference.py``) and the emulator file it writes."""

    name: str
    configuration: str
    emulator: str


@dataclass
class StageComponent:
    """A component of a stage: its campaign manifest, the emulators to train from its runs, and (optionally) the mean and
    maximum wall time, in seconds, of one of its runs."""

    name: str
    manifest: str
    trainings: list = field(default_factory=list)
    cost_mean: float | None = None
    cost_maximum: float | None = None


@dataclass
class Stage:
    """A stage, persisted as a JSON file. Its work queue is the ``state`` directory alongside that file.

    ``training_processes`` is the number of processes used by each training (and so the number of a worker's slots it
    occupies); ``training_seconds`` is the wall time reserved for one training when deciding whether one can finish before
    a job's time limit. ``python`` is the interpreter used to run training.
    """

    stage_file: str
    components: list
    training_processes: int = 1
    training_seconds: float = 3600.0
    python: str = sys.executable

    # Construction and persistence.

    @classmethod
    def create(cls, stage_file, components, overwrite=False, **options):
        if os.path.exists(stage_file) and not overwrite:
            raise FileExistsError(f"stage '{stage_file}' already exists")
        names = [component.name for component in components]
        if len(set(names)) != len(names) or any(not re.fullmatch(r'[A-Za-z0-9_.:+-]+', name) for name in names):
            raise ValueError('component names must be unique, and made of letters, digits and `_.:+-`')
        for component in components:
            trainings = [training.name for training in component.trainings]
            if len(set(trainings)) != len(trainings):
                raise ValueError(f"component '{component.name}': training names must be unique")
        stage = cls(stage_file=os.path.abspath(stage_file), components=components, **options)
        stage.save()
        return stage

    @classmethod
    def load(cls, stage_file):
        with open(stage_file) as handle:
            data = json.load(handle)
        if data.get('format') != STAGE_FORMAT or data.get('formatVersion') != STAGE_VERSION:
            raise ValueError(f"'{stage_file}' is not a Galacticus stage file (version {STAGE_VERSION})")
        components = [StageComponent(**{**component, 'trainings': [Training(**training) for training in component['trainings']]})
                      for component in data['components']]
        options = {key: value for key, value in data.items() if key not in ('format', 'formatVersion', 'components')}
        return cls(stage_file=os.path.abspath(stage_file), components=components, **options)

    def save(self):
        data = {'format': STAGE_FORMAT, 'formatVersion': STAGE_VERSION}
        data.update({key: value for key, value in asdict(self).items() if key != 'stage_file'})
        temporary = self.stage_file + '.tmp'
        with open(temporary, 'w') as handle:
            json.dump(data, handle, indent=1)
        os.replace(temporary, self.stage_file)

    # Paths and state.

    @property
    def directory(self):
        return os.path.dirname(self.stage_file)

    def state(self, component, name=''):
        return os.path.join(self.directory, 'state', component.name, name)

    def campaign(self, component):
        return Campaign.load(component.manifest)

    def run_ended(self, component, index, names=None):
        """Whether a run has ended (completed or failed) in this pass of the workers."""
        names = names if names is not None else set(os.listdir(self.state(component)))
        return f'run{index:06d}.complete' in names or f'run{index:06d}.failed' in names

    def trained(self, component, training):
        """Whether an emulator is trained: its training was recorded as complete in this pass, or (from an earlier pass,
        cleared by :meth:`plan`) its emulator exists and no worker has claimed its training."""
        return (os.path.exists(self.state(component, f'train-{training.name}.complete'))
                or (os.path.exists(training.emulator) and not os.path.lexists(self.state(component, f'train-{training.name}.claim'))))

    def runs_ended(self, component, campaign=None):
        """Whether every run of a component has ended (in the manifest, or in this pass of the workers)."""
        campaign = campaign or self.campaign(component)
        names = set(os.listdir(self.state(component))) if os.path.isdir(self.state(component)) else set()
        return all(run.status in (COMPLETE, FAILED) or self.run_ended(component, run.index, names) for run in campaign.runs)

    # Planning and checking (with no worker running).

    def plan(self, retry_failed=False):
        """Bring the manifests up to date from the outputs of the runs, clear the work queue, and return a summary of the
        work remaining. Call only when no worker of the stage is running.

        Runs which completed are recorded as complete; runs which were started but did not end, or were never started, are
        made pending. Runs which failed stay failed (and are left out of training) unless ``retry_failed``, when they are
        made pending, and the emulators of their components are set aside (renamed ``<emulator>.previous``) for retraining.
        An emulator whose training was claimed but did not complete (possibly written by a killed training) is removed.
        The record of a model process started by an earlier worker and still running (a local worker's model runs outlive
        it if it is killed) is kept, so that no worker starts that run again until the process has ended.
        """
        summary = {'runs': 0, 'seconds': 0.0, 'failed': [], 'trainings': 0, 'orphaned': 0}
        for component in self.components:
            state = self.state(component)
            if os.path.isdir(state):
                for training in component.trainings:
                    claimed = os.path.lexists(os.path.join(state, f'train-{training.name}.claim'))
                    complete = os.path.exists(os.path.join(state, f'train-{training.name}.complete'))
                    if claimed and not complete and os.path.exists(training.emulator):
                        os.remove(training.emulator)
                for name in os.listdir(state):
                    if name.endswith('.process') and _process_alive(os.path.join(state, name)):
                        summary['orphaned'] += 1
                        continue
                    os.remove(os.path.join(state, name))
                if not os.listdir(state):
                    os.rmdir(state)
            campaign = self.campaign(component)
            for run in campaign.runs:
                if run.status != COMPLETE and os.path.exists(campaign.path(run.log_file)):
                    run.status = SUBMITTED
                    run.job_id = None
            campaign.refresh()
            retried = False
            for run in campaign.runs:
                if run.status == FAILED:
                    if retry_failed:
                        run.status, retried = PENDING, True
                    else:
                        summary['failed'].append((component.name, run.index, run.message))
                elif run.status != COMPLETE:
                    run.status = PENDING
            campaign.save()
            if retried:
                for training in component.trainings:
                    if os.path.exists(training.emulator):
                        os.replace(training.emulator, training.emulator + '.previous')
            pending = [run for run in campaign.runs if run.status == PENDING]
            summary['runs'] += len(pending)
            summary['seconds'] += len(pending) * (component.cost_mean or 0.0)
            summary['trainings'] += sum(1 for training in component.trainings if not os.path.exists(training.emulator))
        return summary

    def check(self):
        """Return the trainings (as ``(component, training)`` names) whose emulators are not trained."""
        return [(component.name, training.name) for component in self.components for training in component.trainings
                if not self.trained(component, training)]

    # Working.

    def work(self, slots=None, log=print):
        """Work on the stage until nothing is left that this worker can take. ``slots`` is the number of runs to make at
        once (default: the number of CPUs)."""
        slots = slots or os.cpu_count()
        trainSlots = max(1, min(self.training_processes, slots))
        owner = owner_identity()
        left = time_left()
        deadline = None if left is None else time.time() + left
        environment = dict(os.environ, OMP_NUM_THREADS='1')
        log(f'worker {owner}: {slots} slots, time left {"unlimited" if left is None else f"{left} s"}, {len(self.components)} components')
        campaigns = {component.name: self.campaign(component) for component in self.components}
        for component in self.components:
            os.makedirs(self.state(component), exist_ok=True)
        # Runs to do, slowest component first.
        order = sorted(self.components, key=lambda component: -(component.cost_mean or 0.0))
        todo = [(component, run) for component in order for run in campaigns[component.name].runs if run.status == PENDING]
        deferred = []
        running = {}
        free = slots
        live, liveTime = None, 0.0

        def remaining():
            return float('inf') if deadline is None else deadline - time.time()

        def trainable():
            return [(component, training) for component in self.components for training in component.trainings
                    if not self.trained(component, training)
                    and not os.path.exists(self.state(component, f'train-{training.name}.failed'))
                    and self.runs_ended(component, campaigns[component.name])]

        def start_run(component, run):
            campaign = campaigns[component.name]
            logFile = campaign.path(run.log_file)
            os.makedirs(os.path.dirname(logFile) or '.', exist_ok=True)
            handle = open(logFile, 'w')
            try:
                process = subprocess.Popen(campaign.command(run), cwd=campaign.working_directory, env=environment,
                                           stdout=handle, stderr=subprocess.STDOUT)
            except OSError as error:
                handle.write(f'could not launch {campaign.executable}: {error}\n')
                handle.close()
                _write(self.state(component, f'run{run.index:06d}.failed'), f'could not launch: {error}')
                log(f'{component.name} run {run.index} could not be launched: {error}')
                return False
            # Record the process, so that the run is not taken over while it may still be running (a local worker's model
            # runs outlive it if it is killed).
            _write(self.state(component, f'run{run.index:06d}.process'), f'{socket.gethostname()}:{process.pid}')
            running[process] = ('run', component, run, 1, handle)
            return True

        def start_training(component, training):
            campaign = Campaign.load(component.manifest)
            # Record the runs of this component in its manifest (only the worker training it writes it now).
            for run in campaign.runs:
                if run.status != COMPLETE:
                    run.status, run.job_id = SUBMITTED, None
            campaign.refresh()
            if os.path.exists(training.emulator):
                os.remove(training.emulator)
            handle = open(self.state(component, f'train-{training.name}.log'), 'w')
            command = [self.python, os.path.join(os.environ.get('GALACTICUS_EXEC_PATH', '.'), 'scripts', 'emulation', 'emulatorTrainReference.py'),
                       training.configuration, '--overwrite', '--processes', str(trainSlots)]
            process = subprocess.Popen(command, cwd=os.path.dirname(os.path.abspath(training.configuration)), stdout=handle,
                                       stderr=subprocess.STDOUT)
            running[process] = ('train', component, training, trainSlots, handle)
            log(f'{component.name}: training {training.name} started ({trainSlots} processes)')

        while True:
            # Reap work which has ended.
            for process in [process for process in running if process.poll() is not None]:
                kind, component, item, used, handle = running.pop(process)
                free += used
                handle.close()
                if kind == 'run':
                    campaign = campaigns[component.name]
                    with open(campaign.path(item.log_file), 'a') as logHandle:
                        logHandle.write(f'\n{EXIT_SENTINEL}{process.returncode}\n')
                    resolved = Run(**{**asdict(item), 'log_file': campaign.path(item.log_file), 'output_file': campaign.path(item.output_file)})
                    status, _, message = run_status_from_log(resolved)
                    outcome = 'complete' if status == COMPLETE else 'failed'
                    _write(self.state(component, f'run{item.index:06d}.{outcome}'), message)
                    if outcome == 'failed':
                        log(f'{component.name} run {item.index} failed: {message}')
                else:
                    outcome = 'complete' if process.returncode == 0 and os.path.exists(item.emulator) else 'failed'
                    _write(self.state(component, f'train-{item.name}.{outcome}'), '')
                    log(f'{component.name}: training {item.name} {outcome}')
            # The live jobs are needed only to take over the claims of ended jobs: refresh them at most every five minutes.
            if time.time() - liveTime > 300.0:
                live, liveTime = _live_jobs(), time.time()
            # Training first.
            waiting = False
            for component, training in trainable():
                claimPath = self.state(component, f'train-{training.name}.claim')
                if remaining() < self.training_seconds or not _claimable(claimPath, live):
                    continue
                if free < trainSlots:
                    waiting = True
                    break
                if _claim(claimPath, owner, live, log):
                    start_training(component, training)
                    free -= trainSlots
            # Then model runs.
            reserve = trainSlots if waiting else 0
            while free > reserve and todo:
                component, run = todo.pop(0)
                if self.run_ended(component, run.index):
                    continue
                if remaining() < 1.5 * (component.cost_maximum or 0.0) + 60.0:
                    continue
                if _process_alive(self.state(component, f'run{run.index:06d}.process')):
                    # Started by a worker which has ended, but its model process is still running.
                    deferred.append((component, run))
                    continue
                if _claim(self.state(component, f'run{run.index:06d}.claim'), owner, live, log):
                    if start_run(component, run):
                        free -= 1
                else:
                    deferred.append((component, run))
            if running or waiting:
                time.sleep(1.0)
                continue
            # Nothing of this worker's is running, and no unclaimed run is left. Take over the runs of workers which have
            # ended, if any; then any training which can be started; otherwise stop.
            live, liveTime = _live_jobs(), time.time()
            stale = [(component, run) for component, run in deferred
                     if not self.run_ended(component, run.index) and _claimable(self.state(component, f'run{run.index:06d}.claim'), live)
                     and not _process_alive(self.state(component, f'run{run.index:06d}.process'))]
            if stale:
                deferred = [entry for entry in deferred if entry not in stale]
                todo = stale
                continue
            if any(remaining() >= self.training_seconds and _claimable(self.state(component, f'train-{training.name}.claim'), live)
                   for component, training in trainable()):
                continue
            break
        log(f'worker {owner} finished: no work left that it can take')


def owner_identity():
    """``<job>:<host>:<process ID>``, unique to this worker, where ``<job>`` names the Slurm job (``<array>_<task>`` for an
    array task) whose end makes the worker's claims stale, or is ``local`` outside Slurm."""
    if 'SLURM_ARRAY_TASK_ID' in os.environ and 'SLURM_ARRAY_JOB_ID' in os.environ:
        job = f"{os.environ['SLURM_ARRAY_JOB_ID']}_{os.environ['SLURM_ARRAY_TASK_ID']}"
    else:
        job = os.environ.get('SLURM_JOB_ID', 'local')
    return f'{job}:{socket.gethostname()}:{os.getpid()}'


def time_left():
    """Seconds until this Slurm job's time limit (None outside Slurm, or if unlimited or unknown)."""
    if 'SLURM_JOB_ID' not in os.environ:
        return None
    job = (f"{os.environ['SLURM_ARRAY_JOB_ID']}_{os.environ['SLURM_ARRAY_TASK_ID']}"
           if 'SLURM_ARRAY_TASK_ID' in os.environ and 'SLURM_ARRAY_JOB_ID' in os.environ else os.environ['SLURM_JOB_ID'])
    try:
        result = subprocess.run(['squeue', '-h', '-j', job, '-o', '%L'], capture_output=True, text=True, timeout=120)
    except (OSError, subprocess.TimeoutExpired):
        return None
    match = re.fullmatch(r'(?:(\d+)-)?(?:(\d+):)?(\d+):(\d+)', result.stdout.strip())
    if result.returncode != 0 or not match:
        return None
    days, hours, minutes, seconds = (int(value) if value else 0 for value in match.groups())
    return ((days * 24 + hours) * 60 + minutes) * 60 + seconds


def _live_jobs():
    """The live Slurm jobs (see :func:`slurm_live_jobs`), or None if they are unknown: outside Slurm, if Slurm can not be
    queried, or if the listing omits this worker's own job (which is running, so the listing can not be trusted, and
    taking it as complete would make every other worker look ended)."""
    if 'SLURM_JOB_ID' not in os.environ:
        return None
    live = slurm_live_jobs()
    if live is None or owner_identity().split(':')[0] not in live:
        return None
    return live


def _holder_ended(holder, live):
    """Whether the worker holding a claim has ended: a local worker on this host whose process has ended, or a Slurm job
    not in ``live`` (None if the live jobs are unknown, when no Slurm holder is taken to have ended)."""
    fields = holder.split(':')
    if len(fields) != 3:
        return False
    job, host, pid = fields
    if job == 'local':
        if host != socket.gethostname():
            return False
        try:
            os.kill(int(pid), 0)
        except ProcessLookupError:
            return True
        except (PermissionError, ValueError):
            return False
        return False
    return live is not None and job not in live


def _process_alive(path):
    """Whether the model process recorded in a ``.process`` file may still be running: it is running on this host, or
    is on another host (where it can not be checked, but a Slurm job's processes end with the job)."""
    try:
        with open(path) as handle:
            host, pid = handle.read().strip().split(':')
    except (OSError, ValueError):
        return False
    if host != socket.gethostname():
        return False
    try:
        os.kill(int(pid), 0)
    except ProcessLookupError:
        return False
    except (PermissionError, ValueError):
        return True
    return True


def _claimable(path, live):
    """Whether a claim could be taken: it does not exist, or its holder has ended."""
    try:
        holder = os.readlink(path)
    except FileNotFoundError:
        return True
    except OSError:
        return False
    return _holder_ended(holder, live)


def _claim(path, owner, live, log):
    """Claim a piece of work: create the claim, or take it over if its holder has ended."""
    try:
        os.symlink(owner, path)
        return True
    except FileExistsError:
        pass
    try:
        holder = os.readlink(path)
    except OSError:
        return False
    if not _holder_ended(holder, live):
        return False
    # Whichever worker renames the stale claim away first takes it over.
    try:
        os.rename(path, f'{path}.stale.{socket.gethostname()}.{os.getpid()}')
    except OSError:
        return False
    try:
        os.symlink(owner, path)
    except FileExistsError:
        return False
    log(f'took over {os.path.basename(path)} from {holder}')
    return True


def _write(path, text):
    with open(path, 'w') as handle:
        handle.write(text)
