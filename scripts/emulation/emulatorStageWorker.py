#!/usr/bin/env python3
"""Run the model runs and emulator trainings of a stage from a shared work queue.

A stage (a set of campaigns, each with emulators to train from its runs) is described by a stage file, written with
`Galacticus.Emulation.stage.Stage.create`. The commands, each taking the stage file, are:

  plan    with no worker running: bring the campaign manifests up to date, clear the work queue, and report the work
          remaining (--retryFailed to run failed runs again)
  work    work on the stage until nothing is left that this worker can take (--slots runs at once; default: all CPUs)
  submit  plan, then submit --workers worker jobs to Slurm as an array job (each with --cpus CPUs)
  check   exit 0 if every emulator is trained, else list those which are not

Workers may run on any nodes and start at any time: the first to start takes all the work it can, and one which starts
after all the work is taken exits at once. On a busy queue, ask for more, shorter worker jobs, which start sooner. See
the "Emulator-Assisted Calibration" chapter of the user guide.

Andrew Benson, Claude (28-October-2026).
"""

import argparse
import os
import shlex
import sys
import time

sys.path.insert(0, os.path.join(os.environ.get('GALACTICUS_EXEC_PATH', '.'), 'python'))

from Galacticus.Emulation.stage import Stage  # noqa: E402


def _plan(stage, retry_failed):
    summary = stage.plan(retry_failed=retry_failed)
    for component, index, message in summary['failed']:
        print(f'{component} run {index} failed (left out of training; --retryFailed to run it again): {message}')
    if summary['orphaned'] > 0:
        print(f"{summary['orphaned']} model runs started by an earlier worker are still running: no worker will start them again until they end")
    hours = summary['seconds'] / 3600.0
    print(f"remaining: {summary['runs']} model runs" + (f' (about {hours:.1f} CPU-hours)' if hours > 0.0 else '')
          + f", {summary['trainings']} emulators to train")
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    subparsers = parser.add_subparsers(dest='command', required=True)
    commands = {}
    for name in ('plan', 'work', 'submit', 'check'):
        commands[name] = subparsers.add_parser(name)
        commands[name].add_argument('stage', help='the stage file')
    for name in ('plan', 'submit'):
        commands[name].add_argument('--retryFailed', action='store_true', help='run failed runs again')
    commands['work'].add_argument('--slots', type=int, help='the number of runs to make at once (default: the number of CPUs)')
    commands['submit'].add_argument('--workers', type=int, required=True, help='the number of worker jobs')
    commands['submit'].add_argument('--cpus', type=int, required=True, help='CPUs per worker job (the runs each makes at once)')
    commands['submit'].add_argument('--walltime', help='the walltime of each worker job, e.g. `24:00:00`')
    commands['submit'].add_argument('--memory', type=int, help='memory per worker job, in MB')
    commands['submit'].add_argument('--partition', help='the Slurm partition (default: from galacticusConfig.xml)')
    args = parser.parse_args(argv)
    stage = Stage.load(args.stage)

    if args.command == 'plan':
        _plan(stage, args.retryFailed)
        return 0
    if args.command == 'check':
        missing = stage.check()
        if missing:
            print(f'{len(missing)} emulators are not trained: ' + ', '.join(f'{component}/{training}' for component, training in missing))
            return 1
        print(f'all {sum(len(component.trainings) for component in stage.components)} emulators are trained')
        return 0
    if args.command == 'work':
        stage.work(slots=args.slots, log=lambda message: print(time.strftime('%Y-%m-%d %H:%M:%S'), message, flush=True))
        return 0
    # Submit.
    import queueManager
    summary = _plan(stage, args.retryFailed)
    if summary['runs'] == 0 and summary['trainings'] == 0:
        print('nothing to do')
        return 0
    root = os.path.join(stage.directory, 'workers')
    os.makedirs(root, exist_ok=True)
    script = os.path.abspath(__file__)
    job = {
        'label': 'emulatorStageWorker',
        'launchFile': os.path.join(root, 'workers.slurm'),
        'command': f'{shlex.quote(stage.python)} {shlex.quote(script)} work {shlex.quote(stage.stage_file)} --slots {args.cpus}',
        'array': f'0-{args.workers - 1}',
        'logOutput': os.path.join(root, 'worker_%A_%a.out'),
        'cpusPerTask': args.cpus,
        'countOpenMPThreads': 1,
    }
    if args.walltime is not None:
        job['walltime'] = args.walltime
    if args.memory is not None:
        job['memory'] = args.memory
    job_id = queueManager.factory(args).submit_detached(job)
    print(f'submitted {args.workers} workers as array job {job_id}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
