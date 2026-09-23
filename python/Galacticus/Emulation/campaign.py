"""Run and track a campaign of Galacticus models over an emulator design.

A campaign is the set of model runs listed in a design file (written by the ``emulatorDesign`` task): one run per
realization per design point, each applying that run's change file to a common base parameter file. This module keeps a
*manifest* of the campaign (a JSON file recording each run's status), and can validate runs before submission, run them
locally, submit them to Slurm as an array job, track their status, and resubmit failures. The command-line interface is
``scripts/emulation/emulatorCalibration.py``; the functions here are importable so that other pipelines (e.g. the
dark-matter constraint pipeline) can drive a campaign from within their own state machines.

A run is *complete* when its command exited with status zero, its log contains none of the failure markers used by the
test suite, and its output file exists. Every run's command appends an ``EXIT_STATUS=<n>`` line to its log, so that
completion is judged from the exit status of Galacticus itself, not inferred from the output.

Andrew Benson, Claude (2026)
"""

from __future__ import annotations

import json
import os
import shlex
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict, dataclass, field

import numpy as np

from Galacticus.Emulation.design import read_design

__all__ = [
    'MANIFEST_FORMAT',
    'MANIFEST_VERSION',
    'FAILURE_MARKERS',
    'Run',
    'Campaign',
    'ValidationResult',
    'run_status_from_log',
]

MANIFEST_FORMAT = 'galacticusCampaign'
MANIFEST_VERSION = 1
# The markers which the test suite and CI treat as failure, matched as substrings anywhere in a log.
FAILURE_MARKERS = ('FAIL', 'fatal', 'aborted', 'ODE integration failed', 'unrecognized parameter')
EXIT_SENTINEL = 'EXIT_STATUS='
# Run states.
PENDING = 'pending'
SUBMITTED = 'submitted'
RUNNING = 'running'
COMPLETE = 'complete'
FAILED = 'failed'


@dataclass
class Run:
    """One run of the model: a design point (and realization) with its change, output and log files."""

    index: int
    point: int
    realization: int
    change_file: str
    output_file: str
    log_file: str
    status: str = PENDING
    attempts: int = 0
    job_id: str | None = None
    exit_status: int | None = None
    message: str = ''


@dataclass
class ValidationResult:
    """The result of validating one run: whether Galacticus accepted it, and why not if it did not."""

    run: int
    ok: bool
    message: str = ''


def run_status_from_log(run):
    """Determine a run's status (and exit status, and a message) from its log and output files.

    Returns ``(status, exit_status, message)``. A run with no log is left in its current state (it has not started); a
    log without an exit sentinel means the run is in progress (or was killed without the chance to write one, which
    only the queue manager can tell).
    """
    if not os.path.exists(run.log_file):
        return run.status, run.exit_status, run.message
    with open(run.log_file, errors='replace') as handle:
        log = handle.read()
    exit_status = None
    for line in log.splitlines():
        if line.startswith(EXIT_SENTINEL):
            try:
                exit_status = int(line[len(EXIT_SENTINEL):].strip())
            except ValueError:
                pass
    if exit_status is None:
        return RUNNING, None, ''
    markers = [marker for marker in FAILURE_MARKERS if marker in log]
    if exit_status != 0:
        return FAILED, exit_status, f'exited with status {exit_status}'
    if markers:
        return FAILED, exit_status, 'log contains failure markers: ' + ', '.join(markers)
    if not os.path.exists(run.output_file):
        return FAILED, exit_status, 'output file is missing'
    return COMPLETE, exit_status, ''


@dataclass
class Campaign:
    """A campaign of model runs, persisted as a JSON manifest.

    Paths in the manifest (the design file, base parameter file, change, output and log files) are stored as given, and
    are interpreted relative to the directory from which Galacticus is run (``working_directory``), exactly as the
    design task wrote them.
    """

    manifest_file: str
    design_file: str
    base_parameters: str
    executable: str
    working_directory: str
    runs: list = field(default_factory=list)
    seed_per_point: bool = False

    # Construction and persistence.

    @classmethod
    def create(cls, manifest_file, design_file, base_parameters, executable=None, working_directory=None,
               log_directory=None, overwrite=False):
        """Create a campaign from a design file, and write its manifest.

        Each run's log is written to ``log_directory`` (by default a ``logs`` directory alongside the change files),
        named after its change file.
        """
        if os.path.exists(manifest_file) and not overwrite:
            raise FileExistsError(f"manifest '{manifest_file}' already exists")
        working_directory = os.path.abspath(working_directory or os.getcwd())
        if executable is None:
            executable = os.path.join(os.environ.get('GALACTICUS_EXEC_PATH', '.'), 'Galacticus.exe')
        design = read_design(_resolve(working_directory, design_file))
        runs = []
        for index in range(design.count_runs):
            change_file = design.runs['changeFileName'][index]
            label = os.path.splitext(os.path.basename(change_file))[0]
            directory = log_directory if log_directory is not None else os.path.join(os.path.dirname(os.path.dirname(change_file)) or '.', 'logs')
            runs.append(Run(
                index=index,
                point=int(design.runs['pointIndex'][index]),
                realization=int(design.runs['realizationIndex'][index]),
                change_file=change_file,
                output_file=design.runs['outputFileName'][index],
                log_file=os.path.normpath(os.path.join(directory, label + '.log')),
            ))
        campaign = cls(
            manifest_file=manifest_file,
            design_file=design_file,
            base_parameters=base_parameters,
            executable=executable,
            working_directory=working_directory,
            runs=runs,
            seed_per_point=bool(design.design.attributes.get('seedPerPoint', 0)),
        )
        campaign.save()
        return campaign

    @classmethod
    def load(cls, manifest_file):
        """Load a campaign from its manifest."""
        with open(manifest_file) as handle:
            data = json.load(handle)
        if data.get('format') != MANIFEST_FORMAT:
            raise ValueError(f"'{manifest_file}' is not a Galacticus campaign manifest")
        if data.get('formatVersion') != MANIFEST_VERSION:
            raise ValueError(f"'{manifest_file}' has format version {data.get('formatVersion')}; this reader supports version {MANIFEST_VERSION}")
        runs = [Run(**run) for run in data.pop('runs')]
        for key in ('format', 'formatVersion'):
            data.pop(key)
        return cls(manifest_file=manifest_file, runs=runs, **data)

    def save(self):
        """Write the manifest (atomically, so that an interrupted write never leaves it truncated)."""
        data = {'format': MANIFEST_FORMAT, 'formatVersion': MANIFEST_VERSION}
        data.update({key: value for key, value in asdict(self).items() if key not in ('manifest_file', 'runs')})
        data['runs'] = [asdict(run) for run in self.runs]
        directory = os.path.dirname(os.path.abspath(self.manifest_file))
        os.makedirs(directory, exist_ok=True)
        handle, temporary = tempfile.mkstemp(dir=directory, suffix='.tmp')
        with os.fdopen(handle, 'w') as file:
            json.dump(data, file, indent=1)
        os.replace(temporary, self.manifest_file)

    # Helpers.

    def path(self, path):
        """Resolve a path from the manifest against the working directory."""
        return _resolve(self.working_directory, path)

    def select(self, statuses=None, indices=None):
        """Return the runs with the given statuses and/or indices."""
        runs = self.runs
        if indices is not None:
            wanted = set(int(index) for index in indices)
            runs = [run for run in runs if run.index in wanted]
        if statuses is not None:
            runs = [run for run in runs if run.status in statuses]
        return runs

    def command(self, run, extra_change_files=(), dry_run=False):
        """Return the command (as a list) which runs Galacticus for a run."""
        command = [self.executable, self.base_parameters, run.change_file, *extra_change_files]
        if dry_run:
            command.append('--dry-run')
        return command

    def summary(self):
        """Return a count of runs in each state."""
        counts = {}
        for run in self.runs:
            counts[run.status] = counts.get(run.status, 0) + 1
        return counts

    # Validation.

    def extreme_runs(self):
        """Return the indices of the runs holding the smallest and largest value of each parameter.

        Parameter bounds are limits on each parameter separately, so a design which violates any bound does so at one of
        these runs. Validating them (plus the first run) checks the whole design against the bounds at a cost of at most
        2d+1 runs.
        """
        design = read_design(self.path(self.design_file))
        points = set([0])
        for column in design.design.values.T:
            points.add(int(np.argmin(column)))
            points.add(int(np.argmax(column)))
        return sorted({run.index for run in self.runs if run.point in points and run.realization == 0})

    def validate(self, indices=None, parallel=1):
        """Validate runs by a Galacticus dry run: construct the model at each run's parameters without evolving it.

        This catches parameter values outside declared bounds, paths in the change files which do not exist in the base
        parameter file, and (with per-run seeds) a base parameter file which lacks a ``randomNumberGenerator``. The
        output file of each dry run is redirected to a temporary directory. By default the runs at the extremes of each
        parameter are validated (see :meth:`extreme_runs`). Returns a list of :class:`ValidationResult`.
        """
        if indices is None:
            indices = self.extreme_runs()
        runs = self.select(indices=indices)
        with tempfile.TemporaryDirectory(prefix='galacticusValidate') as directory:

            def check(run):
                redirect = os.path.join(directory, f'run{run.index}.xml')
                with open(redirect, 'w') as handle:
                    handle.write('<changes>\n')
                    handle.write('  <change type="replaceOrAppend" path="outputFileName">\n')
                    handle.write(f'    <outputFileName value="{os.path.join(directory, f"run{run.index}.hdf5")}"/>\n')
                    handle.write('  </change>\n')
                    handle.write('</changes>\n')
                completed = subprocess.run(self.command(run, extra_change_files=(redirect,), dry_run=True),
                                           cwd=self.working_directory, capture_output=True, text=True)
                if completed.returncode == 0:
                    return ValidationResult(run.index, True)
                return ValidationResult(run.index, False, _fatal_message(completed.stdout + completed.stderr))

            with ThreadPoolExecutor(max_workers=max(1, parallel)) as pool:
                return list(pool.map(check, runs))

    # Running.

    def run_local(self, indices=None, parallel=1, threads=None):
        """Run the given runs (by default all pending or failed runs) locally, ``parallel`` at a time."""
        runs = self.select(indices=indices) if indices is not None else self.select(statuses=(PENDING, FAILED))
        environment = os.environ.copy()
        if threads is not None:
            environment['OMP_NUM_THREADS'] = str(threads)

        def execute(run):
            log_file = self.path(run.log_file)
            os.makedirs(os.path.dirname(log_file) or '.', exist_ok=True)
            with open(log_file, 'w') as log:
                completed = subprocess.run(self.command(run), cwd=self.working_directory, env=environment,
                                           stdout=log, stderr=subprocess.STDOUT)
            with open(log_file, 'a') as log:
                log.write(f'\n{EXIT_SENTINEL}{completed.returncode}\n')

        for run in runs:
            run.attempts += 1
            run.status = SUBMITTED
            run.job_id = 'local'
        self.save()
        with ThreadPoolExecutor(max_workers=max(1, parallel)) as pool:
            list(pool.map(execute, runs))
        return self.refresh()

    def array_script(self, runs, job_options=None):
        """Return a Slurm job description (for :mod:`queueManager`) which runs the given runs as one array job.

        The array index is the run index. Each array task reads its change and log file from a task list written
        alongside the manifest (one line per run, in run order), runs Galacticus, and appends the exit status to its log.
        """
        task_list = os.path.splitext(os.path.abspath(self.manifest_file))[0] + '.tasks'
        with open(task_list, 'w') as handle:
            for run in self.runs:
                handle.write(f'{run.change_file} {run.log_file}\n')
        base = shlex.quote(self.base_parameters)
        executable = shlex.quote(self.executable)
        command = '\n'.join([
            f'cd {shlex.quote(self.working_directory)}',
            f'read -r changeFile logFile <<< "$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" {shlex.quote(task_list)})"',
            'mkdir -p "$(dirname "$logFile")"',
            f'{executable} {base} "$changeFile" > "$logFile" 2>&1',
            'status=$?',
            f'echo "{EXIT_SENTINEL}$status" >> "$logFile"',
            'exit $status',
        ])
        label = os.path.splitext(os.path.basename(self.manifest_file))[0]
        job = {
            'label': label,
            'launchFile': os.path.splitext(os.path.abspath(self.manifest_file))[0] + '.slurm',
            'command': command,
            'array': _array_specification([run.index for run in runs]),
            'logOutput': os.path.splitext(os.path.abspath(self.manifest_file))[0] + '_%A_%a.out',
        }
        if job_options:
            job.update(job_options)
        return job

    def submit(self, manager, indices=None, job_options=None):
        """Submit the given runs (by default all pending or failed runs) as a Slurm array job, via a queue manager.

        ``manager`` is a :mod:`queueManager` manager (e.g. from ``queueManager.factory(args)``); its
        ``submit_detached(job)`` method is called, and must return the job ID.
        """
        runs = self.select(indices=indices) if indices is not None else self.select(statuses=(PENDING, FAILED))
        if not runs:
            return None
        job = self.array_script(runs, job_options)
        job_id = manager.submit_detached(job)
        for run in runs:
            run.attempts += 1
            run.status = SUBMITTED
            run.job_id = str(job_id)
            run.exit_status = None
            run.message = ''
            # Remove any log from an earlier attempt, so that its exit sentinel is not mistaken for this attempt's.
            log_file = self.path(run.log_file)
            if os.path.exists(log_file):
                os.remove(log_file)
        self.save()
        return job_id

    def refresh(self):
        """Update the status of every submitted or running run from its log and output files, and save the manifest."""
        for run in self.runs:
            if run.status in (SUBMITTED, RUNNING):
                resolved = Run(**{**asdict(run), 'log_file': self.path(run.log_file), 'output_file': self.path(run.output_file)})
                run.status, run.exit_status, run.message = run_status_from_log(resolved)
                if run.status == RUNNING and not os.path.exists(resolved.log_file):
                    run.status = SUBMITTED
        self.save()
        return self.summary()


def _resolve(directory, path):
    return path if os.path.isabs(path) else os.path.join(directory, path)


def _array_specification(indices):
    """Compress a list of indices into a Slurm ``--array`` specification, e.g. ``0-3,7,9-10``."""
    indices = sorted(set(indices))
    ranges = []
    start = previous = indices[0]
    for index in indices[1:]:
        if index == previous + 1:
            previous = index
            continue
        ranges.append((start, previous))
        start = previous = index
    ranges.append((start, previous))
    return ','.join(str(a) if a == b else f'{a}-{b}' for a, b in ranges)


def _fatal_message(output):
    """Extract the message of a Galacticus fatal error from its output (or the output's end, if there is none)."""
    lines = output.splitlines()
    for i, line in enumerate(lines):
        if line.startswith('Fatal error'):
            message = []
            for following in lines[i + 1:]:
                if following.strip().startswith('Occurred at') or following.startswith('#'):
                    break
                message.append(following.strip())
            return ' '.join(message)
    return ' '.join(line.strip() for line in lines[-3:])
