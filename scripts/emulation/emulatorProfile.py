#!/usr/bin/env python3
"""Fit parameters which change only the analyses of a model (not its merger trees) by direct runs, before emulating the rest.

The configuration is an XML file:

  <emulatorProfile>
    <directory     value="profile"/>                    <!-- where every file is written -->
    <executable    value="Galacticus.exe"/>             <!-- optional; default: $GALACTICUS_EXEC_PATH/Galacticus.exe -->
    <seed          value="1"/>                          <!-- optional; the seed of the star design -->
    <lineSeeds     value="101 202 303 404"/>            <!-- optional; the seeds of the line search -->
    <lineFractions value="0 0.5 1 1.5"/>                <!-- optional; points of the line search, as fractions of the way to the star optimum -->
    <offsetMaximum value="2"/>                          <!-- optional; the star optimum is sought within this many steps of the center -->
    <rootVarianceTargetFractional value="profile" minimum="1.0e-3" maximum="3.0"/>   <!-- optional; or a fixed value (default 0) -->
    <parameter name="progenitorMassFunctionParameters/errorAmplitude" center="0.16" step="0.15" scale="log"/>
    <parameter name="progenitorMassFunctionParameters/errorExponent"  center="0.50" step="0.05"/>
    <component name="MDPL2" baseParameters="progenitorMassFunctionBaseMDPL2.xml">
      <observable label="progenitorMassFunction.*" pattern="true" likelihoodForm="poissonCounts"
                  transform="log10" floor="-13.6" rootVarianceFloored="1.0" likelihoodRangeOnly="true"/>
    </component>
  </emulatorProfile>

Each parameter is stepped linearly, or (with scale="log") in its natural logarithm. Each component is a base parameter file
and the analyses of its output which contribute to the likelihood: observables are given as for training an emulator
(see emulatorTrainReference.py), plus the likelihood form of each (as in the emulated likelihood). With pattern="true" the
label is a regular expression matched against the labels of the analyses. With <rootVarianceTargetFractional
value="profile"/> the fractional root variance of counts of the poissonCounts form is maximized at each point.

The commands, each taking the configuration and a stage (`star` or `line`), are:

  prepare   write the designs, change files and campaign manifests of the stage
  run       run the stage's pending and failed runs locally (all components sharing --parallel slots)
  submit    submit the stage's pending and failed runs to Slurm (one array job per component)
  status    report the status of the stage's runs
  analyze   compute the direct log-likelihoods, and report (and write to <directory>/<stage>/result.json) the result
  fix       write a change file setting the parameters to the stage's optimum (and, with --baseDirectory, copies of
            the base parameter files with that change applied)

A profile is a `star` stage (prepare, run, analyze), which predicts an optimum from a quadratic fit at one seed, then a
`line` stage toward that optimum at several seeds, which checks the gain free of seed noise. Iterate by centering a new
profile on the result. See the "Emulator-Assisted Calibration" chapter of the user guide.

Andrew Benson, Claude (28-October-2026).
"""

import argparse
import os
import sys

import lxml.etree as ET
import numpy as np

sys.path.insert(0, os.path.join(os.environ.get('GALACTICUS_EXEC_PATH', '.'), 'python'))

from Galacticus.Emulation import campaign as cm  # noqa: E402
from Galacticus.Emulation.profile import STAGES, Profile, ProfileComponent, ProfileObservable, ProfileParameter  # noqa: E402
from Galacticus.Emulation.train import Observable  # noqa: E402


def _true(text):
    return text is not None and text.strip().lower() in ('true', 't', '1', 'yes')


def _value(root, name, default, kind=str):
    element = root.find(name)
    return default if element is None else kind(element.get('value'))


def read_configuration(path):
    """Read a profile configuration, returning a :class:`Galacticus.Emulation.profile.Profile`."""
    root = ET.parse(path).getroot()
    if root.tag != 'emulatorProfile':
        raise ValueError(f"'{path}' is not a profile configuration (its root element must be <emulatorProfile>)")
    base = os.path.dirname(os.path.abspath(path))

    def resolve(name):
        return name if os.path.isabs(name) else os.path.join(base, name)

    parameters = [ProfileParameter(name=element.get('name'), center=float(element.get('center')), step=float(element.get('step')),
                                   scale=element.get('scale', 'linear'))
                  for element in root.findall('parameter')]
    components = []
    for element in root.findall('component'):
        observables = []
        for child in element.findall('observable'):

            def number(name, child=child):
                text = child.get(name)
                return None if text is None else float(text)

            observable = Observable(
                label=child.get('label'),
                transform=child.get('transform', 'identity'),
                floor=number('floor'),
                root_variance_floored=number('rootVarianceFloored'),
                undefined=child.get('undefined'),
                root_variance_undefined=number('rootVarianceUndefined'),
                root_variance_zero=number('rootVarianceZero'),
                root_variance_target_minimum=number('rootVarianceTargetMinimum'),
                likelihood_range_only=_true(child.get('likelihoodRangeOnly')),
            )
            if child.get('likelihoodForm') is None:
                raise ValueError(f"observable '{observable.label}': `likelihoodForm` is required")
            observables.append(ProfileObservable(observable, child.get('likelihoodForm'), pattern=_true(child.get('pattern'))))
        components.append(ProfileComponent(name=element.get('name'), base_parameters=resolve(element.get('baseParameters')),
                                           observables=observables))
    rootVariance, bounds = 0.0, (1.0e-3, 3.0)
    element = root.find('rootVarianceTargetFractional')
    if element is not None:
        if element.get('value') == 'profile':
            rootVariance = None
            bounds = (float(element.get('minimum', bounds[0])), float(element.get('maximum', bounds[1])))
        else:
            rootVariance = float(element.get('value'))
    executable = _value(root, 'executable', None)
    profile = Profile(
        directory=resolve(_value(root, 'directory', 'profile')),
        parameters=parameters,
        components=components,
        executable=None if executable is None else resolve(executable),
        seed=_value(root, 'seed', 1, int),
        line_seeds=[int(seed) for seed in _value(root, 'lineSeeds', '101 202 303 404').split()],
        line_fractions=[float(fraction) for fraction in _value(root, 'lineFractions', '0 0.5 1 1.5').split()],
        offset_maximum=_value(root, 'offsetMaximum', 2.0, float),
        root_variance_target_fractional=rootVariance,
        root_variance_target_fractional_bounds=bounds,
    )
    profile.validate()
    return profile


def _report_status(profile, stage, campaigns):
    total = {}
    for component, campaign in zip(profile.components, campaigns):
        summary = campaign.summary()
        for status, count in summary.items():
            total[status] = total.get(status, 0) + count
        for run in campaign.select(statuses=('failed',)):
            print(f'  {component.name} run {run.index} failed: {run.message} [log: {run.log_file}]')
    print(f'{stage}: ' + '  '.join(f'{status}: {count}' for status, count in sorted(total.items())))
    return total


def _report_result(profile, result):
    names = [parameter.name for parameter in profile.parameters]
    if result['stage'] == 'star':
        fit = result['fit']
        print(f"star: {result['pointsUsed']} of {len(result['offsets'])} points; log-likelihood at the center "
              f"{result['logLikelihood'][0]:.2f}")
        if result['rootVarianceTargetFractional'] is not None:
            print(f"  fractional root variance of counts at the center: {result['rootVarianceTargetFractional'][0]:.4g}")
        print(f"  quadratic fit: residual RMS {fit['residualRMS']:.2f}, maximum {fit['residualMaximum']:.2f}")
        print('  gradient (per step):  ' + '  '.join(f'{name} {value:+.2f}' for name, value in zip(names, fit['gradient'])))
        print('  curvature (diagonal): ' + '  '.join(f'{name} {value:+.2f}' for name, value in zip(names, np.diag(fit['hessian']))))
        print(f"  optimum (steps): {np.round(fit['optimum'], 3).tolist()}; predicted gain {result['gainPredicted']:+.2f}")
        if any(abs(offset) >= profile.offset_maximum - 1.0e-6 for offset in fit['optimum']):
            print(f'  the optimum lies at the bound of {profile.offset_maximum} steps: recenter and repeat')
        if fit['covariance'] is not None:
            sigma = np.sqrt(np.diag(fit['covariance'])) * np.array([parameter.step for parameter in profile.parameters])
            print('  1-sigma (from the curvature; in the log for log-scaled parameters): '
                  + '  '.join(f'{name} {value:.3g}' for name, value in zip(names, sigma)))
        else:
            print('  the curvature is not negative definite: no uncertainties')
    else:
        print('line: fraction  ' + '  '.join(f'N={count}' for count in result['countsSeeds']) + '  N=infinity')
        for i, fraction in enumerate(result['fractions']):
            means = '  '.join('-' if row[i] is None else f'{row[i]:.2f}' for row in result['logLikelihoodMean'])
            extrapolated = result['logLikelihoodExtrapolated'][i]
            print(f"  {fraction:6.3f}  {means}  {'-' if extrapolated is None else f'{extrapolated:.2f}'}")
        gain = result['gain']
        print(f"  best: fraction {result['fractions'][result['best']]}" + ('' if gain is None else f'; gain {gain:+.2f}'))
    print('  values: ' + '  '.join(f'{name} {value:.8g}' for name, value in result['optimum'].items()))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    subparsers = parser.add_subparsers(dest='command', required=True)
    commands = {}
    for name, help_text in (('prepare', 'write the designs, change files and campaign manifests of a stage'),
                            ('run', 'run a stage locally'), ('submit', 'submit a stage to Slurm'),
                            ('status', 'report the status of the runs of a stage'), ('analyze', 'analyze a stage'),
                            ('fix', "write a change file setting the parameters to a stage's optimum")):
        commands[name] = subparsers.add_parser(name, help=help_text)
        commands[name].add_argument('configuration', help='the profile configuration (an XML file)')
        commands[name].add_argument('stage', choices=STAGES)
    commands['prepare'].add_argument('--overwrite', action='store_true', help='replace a stage already prepared')
    commands['run'].add_argument('--parallel', type=int, default=1, help='the number of runs to make at once')
    commands['run'].add_argument('--threads', type=int, help='OpenMP threads per run')
    commands['submit'].add_argument('--partition', help='the Slurm partition (default: from galacticusConfig.xml)')
    commands['submit'].add_argument('--walltime', help='the walltime of each run, e.g. `4:00:00`')
    commands['submit'].add_argument('--cpus', type=int, help='CPUs (and OpenMP threads) per run')
    commands['submit'].add_argument('--memory', type=int, help='memory per run, in MB')
    commands['fix'].add_argument('--output', required=True, help='the change file to write')
    commands['fix'].add_argument('--baseDirectory', help='write copies of the base parameter files, with the change applied, here')
    args = parser.parse_args(argv)
    profile = read_configuration(args.configuration)

    if args.command == 'prepare':
        campaigns = profile.prepare(args.stage, overwrite=args.overwrite)
        print(f'{args.stage}: prepared {len(campaigns)} components of {len(campaigns[0].runs)} runs each in {profile.path(args.stage)}')
        return 0

    if args.command == 'fix':
        result = profile.read_result(args.stage)
        if result is None:
            print(f'FAILED: stage `{args.stage}` has not been analyzed')
            return 1
        profile.write_fixed(args.output, result['optimum'])
        print(f'wrote {args.output}')
        if args.baseDirectory is not None:
            from Galacticus.Parameters.resolve import resolve_file
            os.makedirs(args.baseDirectory, exist_ok=True)
            for component in profile.components:
                output = os.path.join(args.baseDirectory, os.path.basename(component.base_parameters))
                resolve_file(component.base_parameters, change_files=[args.output], output=output)
                print(f'wrote {output}')
        return 0

    campaigns = profile.campaigns(args.stage)
    if args.command == 'run':
        cm.run_local(campaigns, parallel=args.parallel, threads=args.threads)
        _report_status(profile, args.stage, campaigns)
        return 0
    if args.command == 'submit':
        import queueManager
        options = {}
        if args.walltime is not None:
            options['walltime'] = args.walltime
        if args.cpus is not None:
            options['cpusPerTask'] = args.cpus
            options['countOpenMPThreads'] = args.cpus
        if args.memory is not None:
            options['memory'] = args.memory
        manager = queueManager.factory(args)
        live_jobs = cm.slurm_live_jobs()
        for component, campaign in zip(profile.components, campaigns):
            campaign.refresh(live_jobs=live_jobs)
            job_id = campaign.submit(manager, job_options=options)
            if job_id is not None:
                print(f'{component.name}: submitted array job {job_id}')
        return 0
    if args.command == 'status':
        live_jobs = cm.slurm_live_jobs() if any(run.job_id is not None and not run.job_id.startswith('local')
                                                for campaign in campaigns for run in campaign.runs) else None
        for campaign in campaigns:
            campaign.refresh(live_jobs=live_jobs)
        _report_status(profile, args.stage, campaigns)
        return 0
    if args.command == 'analyze':
        for campaign in campaigns:
            campaign.refresh()
        total = _report_status(profile, args.stage, campaigns)
        if set(total) - {'complete', 'failed'}:
            print(f'FAILED: stage `{args.stage}` still has runs in progress')
            return 1
        _report_result(profile, profile.analyze(args.stage))
        return 0
    return 1


if __name__ == '__main__':
    sys.exit(main())
