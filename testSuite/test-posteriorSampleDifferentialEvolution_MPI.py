#!/usr/bin/env python3
"""Test that differential evolution MCMC samples a known posterior without bias.

The likelihood is a multivariate normal distribution in eight dimensions, with unequal widths and one strong correlation,
and the priors are uniform and much wider than it, so the posterior is that normal distribution. The posterior is sampled
twice (with one chain per MPI process; the likelihood is cheap, so the processes are oversubscribed as needed):

1. with sixteen chains, constructing proposals from the current states of other chains. When a chain could use its own state
   in constructing its proposal (an asymmetric proposal) the sampled widths were only about 0.6 of their true values with
   eight chains, and 0.75 with sixteen;
2. with only four chains (fewer than the number of parameters), constructing proposals from an archive of past states of all
   chains (``proposalArchive``).

Andrew Benson, Claude (12-October-2026).
"""

import argparse
import glob
import os
import shutil
import subprocess
import sys

import numpy as np

root = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
os.environ['GALACTICUS_EXEC_PATH'] = root
os.chdir(root)

directory = 'testSuite/outputs/posteriorSampleDifferentialEvolution'
countParameters = 8
countSteps = 20000
# The posterior: unequal widths, and a correlation of 0.8 between the first two parameters.
sigmaTrue = np.linspace(0.5, 2.0, countParameters)
correlationTrue = 0.8
covarianceTrue = np.diag(sigmaTrue ** 2)
covarianceTrue[0, 1] = covarianceTrue[1, 0] = correlationTrue * sigmaTrue[0] * sigmaTrue[1]
meanTrue = np.linspace(-1.0, 1.0, countParameters)


def values(array):
    return ' '.join(f'{value:.17g}' for value in np.ravel(array))


def writeParameters(name, archive):
    modelParameters = ''.join(f'''    <modelParameter value="active">
      <name value="p/x{i}"/>
      <distributionFunction1DPrior value="uniform">
        <limitLower value="-20.0"/>
        <limitUpper value="+20.0"/>
      </distributionFunction1DPrior>
      <operatorUnaryMapper value="identity"/>
      <distributionFunction1DPerturber value="cauchy">
        <median value="0.0"/>
        <scale  value="1.0e-4"/>
      </distributionFunction1DPerturber>
    </modelParameter>
''' for i in range(countParameters))
    archiveParameters = '''    <proposalArchive     value="true"/>
    <archiveThinning     value="10"/>
    <archiveInitialCount value="80"/>
''' if archive else ''
    path = f'{directory}/{name}.xml'
    with open(path, 'w') as file:
        file.write(f'''<parameters>
  <formatVersion>2</formatVersion>
  <outputFileName value="{directory}/{name}.hdf5"/>
  <task value="posteriorSample">
    <initializeNodeClassHierarchy value="false"/>
  </task>
  <posteriorSampleLikelihood value="multivariateNormal">
    <means      value="{values(meanTrue)}"/>
    <covariance value="{values(covarianceTrue)}"/>
  </posteriorSampleLikelihood>
  <posteriorSampleSimulation value="differentialEvolution">
    <stepsMaximum           value="{countSteps}"/>
    <acceptanceAverageCount value="10"/>
    <stateSwapCount         value="100"/>
    <logFileRoot            value="{directory}/{name}_chains"/>
    <reportCount            value="1000"/>
    <sampleOutliers         value="false"/>
    <logFlushCount          value="1000"/>
{archiveParameters}{modelParameters}    <posteriorSampleState value="correlation">
      <acceptedStateCount value="100"/>
    </posteriorSampleState>
    <posteriorSampleStateInitialize value="latinHypercube">
      <maximinTrialCount value="100"/>
    </posteriorSampleStateInitialize>
    <posteriorSampleConvergence value="gelmanRubin">
      <thresholdHatR              value="1.2"/>
      <burnCount                  value="1000"/>
      <testCount                  value="100"/>
      <outlierCountMaximum        value="1"/>
      <outlierSignificance        value="0.95"/>
      <outlierLogLikelihoodOffset value="60"/>
      <reportCount                value="10"/>
      <logFileName                value="{directory}/{name}_convergence.log"/>
    </posteriorSampleConvergence>
    <posteriorSampleStoppingCriterion value="stepCount">
      <stopAfterCount value="{countSteps}"/>
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
      <logFileName           value="{directory}/{name}_proposalSize.log"/>
    </posteriorSampleDffrntlEvltnProposalSize>
  </posteriorSampleSimulation>
  <randomNumberGenerator value="GSL">
    <seed          value="811"/>
    <mpiRankOffset value="true"/>
  </randomNumberGenerator>
</parameters>
''')
    return path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--processesPerNode', type=int, default=1)
    parser.add_argument('--allowRunAsRoot', type=str, default='no')
    args, _ = parser.parse_known_args()
    allowRunAsRoot = ['--allow-run-as-root'] if args.allowRunAsRoot == 'yes' else []
    shutil.rmtree(directory, ignore_errors=True)
    os.makedirs(directory)
    for name, countChains, archive in (('current', 16, False), ('archive', 4, True)):
        check(name, countChains, archive, allowRunAsRoot)


def check(name, countChains, archive, allowRunAsRoot):
    """Sample the posterior, and check the sampled mean, widths, and correlation against their true values."""
    label = f'{countChains} chains, proposals from {"an archive of past states" if archive else "current states"}'
    with open(f'{directory}/{name}.log', 'w') as log:
        status = subprocess.run(['mpirun', '--oversubscribe', '-np', str(countChains)] + allowRunAsRoot + ['./Galacticus.exe', writeParameters(name, archive)],
                                stdout=log, stderr=subprocess.STDOUT, env=dict(os.environ, OMP_NUM_THREADS='1')).returncode
    if status != 0:
        print(f'FAILED: posterior sampling ({label}) failed (see {directory}/{name}.log)')
        return
    # Gather the converged states from the chains.
    states = []
    for path in sorted(glob.glob(f'{directory}/{name}_chains_[0-9][0-9][0-9][0-9].log')):
        with open(path) as file:
            for line in file:
                if line.startswith('#'):
                    continue
                columns = line.split()
                if columns[3] == 'T':
                    states.append([float(value) for value in columns[6:6 + countParameters]])
    states = np.array(states)
    if states.shape[0] < 10000:
        print(f'FAILED: only {states.shape[0]} converged states ({label}) were found in the chains (see {directory}/{name}_convergence.log)')
        return
    meanMCMC = states.mean(axis=0)
    covarianceMCMC = np.cov(states.T)
    sigmaMCMC = np.sqrt(np.diag(covarianceMCMC))
    offset = np.abs(meanMCMC - meanTrue) / sigmaTrue
    ratio = sigmaMCMC / sigmaTrue
    correlationMCMC = covarianceMCMC[0, 1] / (sigmaMCMC[0] * sigmaMCMC[1])
    print(f'{label}: {states.shape[0]} converged states')
    message = f'{label}: posterior means are within {offset.max():.3f} sigma of the true values'
    print(f'SUCCESS: {message}' if offset.max() < 0.25 else f'FAILED: {message} (tolerance 0.25)')
    message = f'{label}: posterior widths are {ratio.min():.3f} to {ratio.max():.3f} times the true values'
    print(f'SUCCESS: {message}' if 0.9 < ratio.min() and ratio.max() < 1.1 else f'FAILED: {message} (tolerance 0.9-1.1)')
    message = f'{label}: posterior correlation is {correlationMCMC:.3f} (true value {correlationTrue:.3f})'
    print(f'SUCCESS: {message}' if abs(correlationMCMC - correlationTrue) < 0.05 else f'FAILED: {message} (tolerance 0.05)')

if __name__ == '__main__':
    try:
        main()
    except Exception as error:  # Report any unexpected error as a failure, rather than a bare traceback.
        print(f'FAILED: unexpected error: {type(error).__name__}: {error}')
    sys.exit(0)
