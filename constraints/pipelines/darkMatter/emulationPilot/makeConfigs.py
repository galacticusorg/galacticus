#!/usr/bin/env python3
"""Write the training configuration, and the emulated-likelihood MCMC configuration, for an emulator campaign over the
progenitor mass function stage of the dark-matter pipeline.

* `train.xml`: one log10 observable per `progenitorMassFunction` analysis found in a completed run, collecting only the
  bins used in the likelihood. Each is floored at the value giving an expected count of `--countFloor` in its bin of
  largest count per value, so floored bins contribute negligible expected counts.
* `mcmc.xml`: a differential evolution simulation (one chain per MPI process) with the `emulated` likelihood, the
  `poissonCounts` form, and the stage's active parameters, unchanged: the emulator inputs, and
  `rootVarianceTargetFractional`, bound to the likelihood through `rootVarianceTargetFractionalParameter`.

Andrew Benson, Claude (13-October-2026)
"""
import argparse
import copy
import os

import h5py
import lxml.etree as ET
import numpy as np

parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
parser.add_argument('--pipelineOutput', required=True, help="the pipeline's output directory, holding the generated stage")
parser.add_argument('--campaignDirectory', required=True, help='the campaign directory (holding design.hdf5)')
parser.add_argument('--modelFile', required=True, help='the output of any completed run of the campaign')
parser.add_argument('--countFloor', type=float, default=0.01, help='the expected count at which to floor each observable')
parser.add_argument('--rootVarianceFloored', type=float, default=1.0, help='the root variance (dex) of floored values')
parser.add_argument('--folds', type=int, default=5, help='cross-validation folds (0 to disable)')
parser.add_argument('--steps', type=int, default=20000, help='the maximum number of MCMC steps')
parser.add_argument('--emulatorFile', default=None, help='the emulator file (default: <campaign>/emulator.hdf5)')
args = parser.parse_args()
pipeline = os.path.abspath(args.pipelineOutput) + '/'
campaign = os.path.abspath(args.campaignDirectory) + '/'
emulatorFile = os.path.abspath(args.emulatorFile) if args.emulatorFile else campaign + 'emulator.hdf5'

with h5py.File(args.modelFile, 'r') as file:
    labels = sorted(name for name in file['analyses'] if name.startswith('progenitorMassFunction'))
    floors = {label: np.log10(args.countFloor / file['analyses'][label]['countPerValue'][()].max()) for label in labels}
if not labels:
    raise RuntimeError(f'no progenitorMassFunction analyses in {args.modelFile}')

train = ET.Element('emulatorTrain')
ET.SubElement(train, 'designFileName', value=campaign + 'design.hdf5')
ET.SubElement(train, 'emulatorFileName', value=emulatorFile)
ET.SubElement(train, 'foldsCrossValidation', value=str(args.folds))
for label in labels:
    ET.SubElement(train, 'observable', label=label, transform='log10', floor=f'{floors[label]:.4f}',
                  rootVarianceFloored=str(args.rootVarianceFloored), likelihoodRangeOnly='true')
ET.ElementTree(train).write(campaign + 'train.xml', pretty_print=True, xml_declaration=True, encoding='UTF-8')

config = ET.parse(pipeline + 'progenitorMassFunctionConfig.xml').getroot()
mcmc = ET.Element('parameters')
ET.SubElement(mcmc, 'formatVersion').text = '2'
ET.SubElement(mcmc, 'verbosityLevel', value='standard')
ET.SubElement(mcmc, 'task', value='posteriorSample')
ET.SubElement(mcmc, 'outputFileName', value=campaign + 'mcmc.hdf5')
likelihood = ET.SubElement(mcmc, 'posteriorSampleLikelihood', value='emulated')
for label in labels:
    emulator = ET.SubElement(likelihood, 'emulator', value='gaussianProcess')
    ET.SubElement(emulator, 'fileName', value=emulatorFile)
    ET.SubElement(emulator, 'label', value=label)
ET.SubElement(likelihood, 'likelihoodForms', value='poissonCounts')
ET.SubElement(likelihood, 'includeEmulatorVariance', value='true')
ET.SubElement(likelihood, 'rootVarianceTargetFractionalParameter',
              value='progenitorMassFunctionParameters/rootVarianceTargetFractional')
ET.SubElement(mcmc, 'posteriorSampleStateInitialize', value='priorRandom')
proposal = ET.SubElement(mcmc, 'posteriorSampleDffrntlEvltnProposalSize', value='adaptive')
for name, value in (('gammaInitial', '0.5'), ('gammaAdjustFactor', '1.1'), ('gammaMinimum', '1.0e-4'), ('gammaMaximum', '3.0'),
                    ('acceptanceRateMinimum', '0.1'), ('acceptanceRateMaximum', '0.9'), ('updateCount', '10'),
                    ('logFileName', campaign + 'chains/gamma.log')):
    ET.SubElement(proposal, name, value=value)
simulation = ET.SubElement(mcmc, 'posteriorSampleSimulation', value='differentialEvolution')
for name, value in (('stepsMaximum', str(args.steps)), ('acceptanceAverageCount', '10'), ('stateSwapCount', '11'),
                    ('logFileRoot', campaign + 'chains/chain'), ('reportCount', '100'), ('sampleOutliers', 'false'),
                    ('logFlushCount', '1')):
    ET.SubElement(simulation, name, value=value)
state = ET.SubElement(simulation, 'posteriorSampleState', value='correlation')
ET.SubElement(state, 'acceptedStateCount', value='100')
ET.SubElement(simulation, 'posteriorSampleConvergence', value='never')
ET.SubElement(simulation, 'posteriorSampleStoppingCriterion', value='never')
ET.SubElement(simulation, 'posteriorSampleDffrntlEvltnRandomJump', value='adaptive')
for parameter in config.find('posteriorSampleSimulation').findall('modelParameter'):
    simulation.append(copy.deepcopy(parameter))
rng = ET.SubElement(mcmc, 'randomNumberGenerator', value='GSL')
ET.SubElement(rng, 'seed', value='9372')
ET.SubElement(rng, 'mpiRankOffset', value='true')
ET.ElementTree(mcmc).write(campaign + 'mcmc.xml', pretty_print=True, xml_declaration=True, encoding='UTF-8')
os.makedirs(campaign + 'chains', exist_ok=True)
print(f'wrote {campaign}train.xml and {campaign}mcmc.xml for {len(labels)} observables')
