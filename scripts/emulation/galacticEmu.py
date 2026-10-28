#!/usr/bin/env python3
"""Import a GalacticEmu campaign, and check Galacticus' emulation against GalacticEmu.

GalacticEmu (https://github.com/Andrew-Robertson/Galacticus-emulation) calibrates Galacticus with PCA-compressed
Gaussian-process emulators trained on a Sobol design of model runs (arXiv:2609.10680). This script reproduces such a
calibration with Galacticus' own emulation tools. Its subcommands are:

  design CAMPAIGN DESIGN
      import the campaign's design (samples.csv, campaign_design.json) as a Galacticus design file;
  configure CAMPAIGN DESIGN --observables KEYS --emulator EMULATOR --training TRAINING [--mcmc MCMC]
      write a training configuration (for emulatorTrainReference.py) reproducing GalacticEmu's treatment of the given
      observables, and optionally a parameter file sampling the posterior with the emulated likelihood;
  parity-training CAMPAIGN DESIGN --observables KEYS
      compare the training sets collected by Galacticus with those prepared by GalacticEmu;
  parity-emulators BUNDLE EMULATOR --observables KEYS
      compare trained emulators with those of a GalacticEmu bundle;
  parity-likelihood BUNDLE INPUTS EMULATOR [--results RESULTS] [--grid GRID]
      compare the emulated likelihood with that of a GalacticEmu MCMC, at random points, at the MCMC's maximum a
      posteriori point, and (optionally) with the likelihoods computed by Galacticus in a `grid` simulation.

KEYS are GalacticEmu observable keys (e.g. smf_z0,mzr_blanc2019), BUNDLE a GalacticEmu emulator bundle
(*_bundle.joblib), and INPUTS and RESULTS the inputs and results of a GalacticEmu MCMC (*_mcmc_inputs.joblib,
*_mcmc_results.hdf5). Every subcommand other than `design` requires GalacticEmu to be importable (its tag
roman-grs-calibration-arxiv-v1 was used in developing this script). See the "Emulator-Assisted Calibration" chapter of the
user guide.

Andrew Benson, Claude (12-October-2026).
"""

import argparse
import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.environ.get('GALACTICUS_EXEC_PATH', '.'), 'python'))

from Galacticus.Emulation import emulatorFile as ef  # noqa: E402
from Galacticus.Emulation import galacticEmu as ge  # noqa: E402
from Galacticus.Emulation.design import parse_descriptor, read_design, write_design  # noqa: E402


def _keys(text):
    return [key.strip() for key in text.split(',') if key.strip()]


def command_design(args):
    log_parameters = _keys(args.log_parameters) if args.log_parameters is not None else None
    design_file = ge.import_design(args.campaign, log_parameters=log_parameters)
    write_design(args.design, design_file, overwrite=args.overwrite)
    design = design_file.design
    print(f'{design.quantiles.shape[0]} design points, {len(design.names)} parameters written to {args.design}')
    for name, prior, mapper in zip(design_file.attributes['shortNames'].split(), design.priors, design.mappers):
        print(f'  {name:45s} {prior.class_name:10s} {mapper}')
    return 0


def _model_parameters(design):
    """The active model parameters of a design, as XML."""
    text = ''
    for name, prior, mapper in zip(design.names, design.priors, design.mappers):
        parameters = parse_descriptor(prior.parameters['descriptor'])
        elements = ''.join(f'        <{key} value="{value!r}"/>\n' for key, value in parameters.items())
        text += (f'    <modelParameter value="active">\n      <name value="{name}"/>\n'
                 f'      <distributionFunction1DPrior value="{prior.class_name}">\n{elements}      </distributionFunction1DPrior>\n'
                 f'      <operatorUnaryMapper value="{mapper}"/>\n'
                 f'      <distributionFunction1DPerturber value="cauchy">\n        <median value="0.0"/>\n        <scale  value="1.0e-4"/>\n'
                 f'      </distributionFunction1DPerturber>\n    </modelParameter>\n')
    return text


def command_configure(args):
    observables = ge.galacticemu_observables(args.campaign, _keys(args.observables), conventions=not args.no_conventions,
                                             root_variance_target_minimum=args.target_sigma_minimum)
    # GalacticEmu's training: 99% of the variance retained by PCA, four restarts of the optimizer, no cross-validation, and
    # only the stated noise of each run.
    lines = ['<emulatorTrain>', f'  <designFileName       value="{args.design}"/>', f'  <emulatorFileName     value="{args.emulator}"/>',
             '  <pcaVarianceRetained  value="0.99"/>', '  <restartsOptimizer    value="4"/>', '  <foldsCrossValidation value="0"/>',
             '  <noiseLearned         value="false"/>']
    for observable in observables:
        attributes = {'label': observable.label}
        attributes.update({name: (f'{value!r}' if isinstance(value, float) else str(value)) for name, value in observable.attributes().items()})
        lines.append('  <observable ' + ' '.join(f'{name}="{value}"' for name, value in attributes.items()) + '/>')
    lines.append('</emulatorTrain>')
    with open(args.training, 'w') as file:
        file.write('\n'.join(lines) + '\n')
    print(f'wrote {args.training}')
    if args.mcmc:
        design = read_design(args.design).design
        labels = [observable.label for observable in observables]
        if args.likelihood_observables:
            keys = _keys(args.observables)
            labels = [labels[keys.index(key)] for key in _keys(args.likelihood_observables)]
        emulators = ''.join(f'    <emulator value="gaussianProcess">\n      <fileName value="{args.emulator}"/>\n'
                            f'      <label    value="{label}"/>\n    </emulator>\n' for label in labels)
        root = os.path.splitext(args.mcmc)[0]
        with open(args.mcmc, 'w') as file:
            file.write(f'''<parameters>
  <formatVersion>2</formatVersion>
  <outputFileName value="{root}.hdf5"/>
  <task value="posteriorSample">
    <initializeNodeClassHierarchy value="false"/>
  </task>
  <posteriorSampleLikelihood value="emulated">
{emulators}    <likelihoodForms         value="gaussianDiagonal"/>
    <includeEmulatorVariance value="true"/>
  </posteriorSampleLikelihood>
  <posteriorSampleSimulation value="differentialEvolution">
    <stepsMaximum           value="{args.steps}"/>
    <acceptanceAverageCount value="10"/>
    <stateSwapCount         value="100"/>
    <logFileRoot            value="{root}"/>
    <reportCount            value="100"/>
    <sampleOutliers         value="false"/>
    <logFlushCount          value="10"/>
    <proposalArchive        value="true"/>
    <archiveInitialCount    value="{10 * len(design.names)}"/>
{_model_parameters(design)}    <posteriorSampleState value="correlation">
      <acceptedStateCount value="100"/>
    </posteriorSampleState>
    <posteriorSampleStateInitialize value="latinHypercube">
      <maximinTrialCount value="100"/>
    </posteriorSampleStateInitialize>
    <posteriorSampleConvergence value="gelmanRubin">
      <thresholdHatR              value="1.2"/>
      <burnCount                  value="500"/>
      <testCount                  value="100"/>
      <outlierCountMaximum        value="1"/>
      <outlierSignificance        value="0.95"/>
      <outlierLogLikelihoodOffset value="60"/>
      <reportCount                value="10"/>
      <logFileName                value="{root}_convergence.log"/>
    </posteriorSampleConvergence>
    <posteriorSampleStoppingCriterion value="stepCount">
      <stopAfterCount value="{args.steps}"/>
    </posteriorSampleStoppingCriterion>
    <posteriorSampleDffrntlEvltnRandomJump value="adaptive"/>
    <posteriorSampleDffrntlEvltnProposalSize value="adaptive">
      <gammaInitial          value="0.5"/>
      <gammaAdjustFactor     value="1.1"/>
      <gammaMinimum          value="1.0e-4"/>
      <gammaMaximum          value="3.0"/>
      <acceptanceRateMinimum value="0.1"/>
      <acceptanceRateMaximum value="0.9"/>
      <updateCount           value="10"/>
      <logFileName           value="{root}_proposalSize.log"/>
    </posteriorSampleDffrntlEvltnProposalSize>
  </posteriorSampleSimulation>
  <randomNumberGenerator value="GSL">
    <seed          value="{args.seed}"/>
    <mpiRankOffset value="true"/>
  </randomNumberGenerator>
</parameters>
''')
        print(f'wrote {args.mcmc}')
    return 0


def command_parity_training(args):
    results = ge.compare_training(args.campaign, args.design, _keys(args.observables), conventions=not args.no_conventions)
    print(f'{"observable":36s} {"bins":>5s}  max|dy|    max|dσ/σ|  target max|d|  max|dσ_target/σ_target|')
    for key, result in results.items():
        print(f'{key:36s} {result["bins"]:5d}  {result["y"]:9.2e}  {result["rootVariance"]:9.2e}  {result["target"]:13.2e}  {result["rootVarianceTarget"]:9.2e}'
              + (f'  ({result["binsUnsupported"]} bins unsupported in GalacticEmu)' if result['binsUnsupported'] else ''))
    return 0


def _load(path):
    import joblib
    return joblib.load(path)


def command_parity_emulators(args):
    bundle = _load(args.bundle)
    content = ef.read(args.emulator)
    results = ge.compare_emulators(bundle, content, _keys(args.observables), count_points=args.points)
    for key, result in results.items():
        ours, theirs = result['logMarginalLikelihood']
        print(f'{key}: components {result["components"][0]} (GalacticEmu {result["components"][1]}); |cos| ' + ' '.join(f'{c:.4f}' for c in result['cosines']))
        print('    log marginal likelihood:            ' + ' '.join(f'{value:9.1f}' for value in ours))
        print('    log marginal likelihood GalacticEmu: ' + ' '.join(f'{value:9.1f}' for value in theirs))
        print('    mean difference / target sigma: median {:.3f}, 95% {:.3f}, max {:.3f}'.format(*result['meanDifferenceTarget']))
        print('    mean difference / emulator sigma: median {:.3f}, 95% {:.3f}'.format(*result['meanDifferenceEmulator']))
        print('    emulator sigma ratio: median {:.3f}, 5-95% {:.3f}-{:.3f}'.format(*result['stdRatio']))
    return 0


def command_parity_likelihood(args):
    bundle = _load(args.bundle)
    inputs = _load(args.inputs)
    content = ef.read(args.emulator)
    quantiles = np.random.default_rng(args.seed).uniform(size=(args.points, len(bundle['input_columns'])))
    ours, theirs = ge.compare_likelihood(bundle, inputs, content, quantiles)
    difference = ours - theirs
    print(f'{args.points} random points: log-likelihood {ours.min():.1f} to {ours.max():.1f}; difference from GalacticEmu: '
          f'median {np.median(difference):+.3f}, max |.| {np.max(np.abs(difference)):.3f}')
    if args.results:
        quantile_map, log_probability = ge.map_quantiles(args.results)
        ours_map, theirs_map = ge.compare_likelihood(bundle, inputs, content, quantile_map)
        print(f'maximum a posteriori point: log-likelihood {ours_map[0]:.4f} (GalacticEmu {theirs_map[0]:.4f}; its log-posterior {log_probability:.4f})')
    if args.grid:
        design = content.design
        rows = []
        for path in sorted(glob.glob(f'{args.grid}_[0-9][0-9][0-9][0-9].log')):
            for line in open(path):
                if not line.startswith('#'):
                    columns = line.split()
                    rows.append([float(columns[5])] + [float(value) for value in columns[6:]])
        if not rows:
            print(f'no grid states found in {args.grid}_*.log')
            return 1
        rows = np.array(rows)
        quantiles_grid = np.column_stack([ge.prior_quantiles(prior.class_name, prior.parameters, rows[:, 1 + j]) for j, prior in enumerate(design.priors)])
        labels = [bundle['observables'][key]['analysis'] for key in inputs['observable_keys']]
        reference = ge.log_likelihood_gaussian_diagonal(content, labels, quantiles_grid)
        error = np.max(np.abs(rows[:, 0] - reference) / np.maximum(np.abs(reference), 1.0))
        print(f'{len(rows)} Galacticus grid states: log-likelihood agrees with the Python reference to {error:.2e} (relative)')
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0], epilog=__doc__.split('\n\n', 1)[1],
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    subparsers = parser.add_subparsers(dest='command', required=True)
    sub = subparsers.add_parser('design', help='import the design of a GalacticEmu campaign')
    sub.add_argument('campaign', help='the campaign directory (containing samples.csv and campaign_design.json)')
    sub.add_argument('design', help='the design file to write')
    sub.add_argument('--log-parameters', help='short names of the parameters to sample in their logarithm (default: those with log-normal priors)')
    sub.add_argument('--overwrite', action='store_true', help='replace an existing design file')
    sub.set_defaults(function=command_design)
    sub = subparsers.add_parser('configure', help='write a training configuration, and optionally an MCMC parameter file')
    sub.add_argument('campaign')
    sub.add_argument('design')
    sub.add_argument('--observables', required=True, help='comma-separated GalacticEmu observable keys')
    sub.add_argument('--emulator', required=True, help='the emulator file to be trained')
    sub.add_argument('--training', required=True, help='the training configuration to write')
    sub.add_argument('--mcmc', help='the MCMC parameter file to write')
    sub.add_argument('--likelihood-observables', help='the observables (keys) to include in the likelihood (default: all)')
    sub.add_argument('--steps', type=int, default=40000, help='the number of MCMC steps')
    sub.add_argument('--seed', type=int, default=219, help='the random number seed of the MCMC')
    sub.add_argument('--target-sigma-minimum', type=float, default=1.0e-3, help='the floor on the target uncertainty (as GalacticEmu)')
    sub.add_argument('--no-conventions', action='store_true', help="do not reproduce GalacticEmu's treatment of zero training noise and its target floor")
    sub.set_defaults(function=command_configure)
    sub = subparsers.add_parser('parity-training', help='compare training sets with those of GalacticEmu')
    sub.add_argument('campaign')
    sub.add_argument('design')
    sub.add_argument('--observables', required=True)
    sub.add_argument('--no-conventions', action='store_true')
    sub.set_defaults(function=command_parity_training)
    sub = subparsers.add_parser('parity-emulators', help='compare emulators with those of a GalacticEmu bundle')
    sub.add_argument('bundle')
    sub.add_argument('emulator')
    sub.add_argument('--observables', required=True)
    sub.add_argument('--points', type=int, default=1000, help='the number of random points at which to compare predictions')
    sub.set_defaults(function=command_parity_emulators)
    sub = subparsers.add_parser('parity-likelihood', help='compare the emulated likelihood with that of a GalacticEmu MCMC')
    sub.add_argument('bundle')
    sub.add_argument('inputs')
    sub.add_argument('emulator')
    sub.add_argument('--results', help='the GalacticEmu MCMC results (to compare at its maximum a posteriori point)')
    sub.add_argument('--grid', help='the log file root of a Galacticus `grid` simulation using the emulated likelihood')
    sub.add_argument('--points', type=int, default=1000)
    sub.add_argument('--seed', type=int, default=1234)
    sub.set_defaults(function=command_parity_likelihood)
    args = parser.parse_args(argv)
    return args.function(args)


if __name__ == '__main__':
    sys.exit(main())
