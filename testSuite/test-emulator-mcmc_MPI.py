#!/usr/bin/env python3
"""Test emulator-assisted calibration end to end, against a posterior known analytically.

The "model" is the synthetic analysis `outputAnalysisSyntheticLinear`, a straight line y(x) = a + b x, compared with target
data having independent Gaussian uncertainties, so that (with uniform priors much wider than the posterior) the posterior of
the intercept, a, and slope, b, is exactly Gaussian with known mean and covariance. The test runs every step of a
calibration:

1. generates a Sobol design over (a, b) with the `emulatorDesign` task;
2. runs the campaign of models over the design (`Galacticus.Emulation.campaign`);
3. collects the training set from the model outputs, and trains an emulator (`Galacticus.Emulation.train`);
4. samples the posterior with differential evolution MCMC using the emulated likelihood;

and then checks that the posterior mean and covariance found by the MCMC recover the known values.

Posterior sampling requires an MPI build of Galacticus; training requires scikit-learn.

Andrew Benson, Claude (7-October-2026).
"""

import argparse
import glob
import os
import shutil
import subprocess
import sys

import lxml.etree as ET
import numpy as np

root = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
os.environ['GALACTICUS_EXEC_PATH'] = root
sys.path.insert(0, os.path.join(root, 'python'))
os.chdir(root)

directory = 'testSuite/outputs/emulatorMCMC'
label = 'syntheticLinear'
# The synthetic data: a straight line, observed with independent Gaussian uncertainties.
x = np.linspace(0.0, 1.0, 8)
interceptTrue, slopeTrue, sigma = 0.5, 1.2, 0.1
yTarget = interceptTrue + slopeTrue * x + sigma * np.random.default_rng(3).normal(size=x.size)
# Uniform priors on the intercept and slope - wide compared to the posterior (whose widths are about 0.07 and 0.11, with a
# correlation of -0.84).
priors = {'outputAnalysis/intercept': (-0.5, 1.5), 'outputAnalysis/slope': (0.0, 2.5)}
countDesign = 64
countSteps = 2000
# The number of chains (one per MPI process). Differential evolution builds its proposals from differences between chains,
# and with few chains the population under-disperses: with 4 chains the posterior widths found here are only 0.6 of their
# true values, with 8 chains 0.8, and with 16 chains 0.95. Each evaluation of the emulated likelihood is cheap, so the
# processes are oversubscribed as needed.
countChains = 16


def posteriorAnalytic():
    """Return the mean and covariance of the posterior of (intercept, slope): for uniform priors, the generalized least squares
    solution."""
    design = np.column_stack([np.ones_like(x), x])
    precision = design.T @ design / sigma ** 2
    covariance = np.linalg.inv(precision)
    return covariance @ design.T @ yTarget / sigma ** 2, covariance


def values(array):
    return ' '.join(f'{value:.17g}' for value in array)


def modelParameters():
    """The active model parameters, as XML."""
    return ''.join(f'''    <modelParameter value="active">
      <name value="{name}"/>
      <distributionFunction1DPrior value="uniform">
        <limitLower value="{lower}"/>
        <limitUpper value="{upper}"/>
      </distributionFunction1DPrior>
      <operatorUnaryMapper value="identity"/>
      <distributionFunction1DPerturber value="cauchy">
        <median value="0.0"/>
        <scale  value="1.0e-4"/>
      </distributionFunction1DPerturber>
    </modelParameter>
''' for name, (lower, upper) in priors.items())


def writeBase():
    """Write the base parameter file: the quick test model, with its outputs replaced by the synthetic analysis."""
    tree = ET.parse('parameters/quickTest.xml')
    for outputter in tree.getroot().findall('mergerTreeOutputter'):
        tree.getroot().remove(outputter)
    ET.SubElement(tree.getroot(), 'mergerTreeOutputter', value='analyzer')
    analysis = ET.SubElement(tree.getroot(), 'outputAnalysis', value='syntheticLinear')
    for name, value in (('x', values(x)), ('yTarget', values(yTarget)), ('varianceTarget', values(np.full(x.size, sigma ** 2))),
                        ('intercept', '0.5'), ('slope', '1.0')):
        ET.SubElement(analysis, name, value=value)
    path = f'{directory}/base.xml'
    tree.write(path, xml_declaration=True, encoding='UTF-8')
    return path


def writeDesign():
    path = f'{directory}/design.xml'
    with open(path, 'w') as file:
        file.write(f'''<parameters>
  <formatVersion>2</formatVersion>
  <task value="emulatorDesign">
    <designFileName     value="{directory}/design.hdf5"/>
    <changeFilesRoot    value="{directory}/changes/point"/>
    <outputFileNameRoot value="{directory}/models/point"/>
{modelParameters()}    <posteriorSamples value="sobol">
      <countSamples value="{countDesign}"/>
      <randomShift  value="true"/>
    </posteriorSamples>
  </task>
  <randomNumberGenerator value="GSL">
    <seed value="11"/>
  </randomNumberGenerator>
</parameters>
''')
    return path


def writeSampler(emulatorFileName):
    path = f'{directory}/mcmc.xml'
    with open(path, 'w') as file:
        file.write(f'''<parameters>
  <formatVersion>2</formatVersion>
  <outputFileName value="{directory}/mcmc.hdf5"/>
  <task value="posteriorSample">
    <initializeNodeClassHierarchy value="false"/>
  </task>
  <posteriorSampleLikelihood value="emulated">
    <emulator value="gaussianProcess">
      <fileName value="{emulatorFileName}"/>
      <label    value="{label}"/>
    </emulator>
  </posteriorSampleLikelihood>
  <posteriorSampleSimulation value="differentialEvolution">
    <stepsMaximum           value="{countSteps}"/>
    <acceptanceAverageCount value="10"/>
    <stateSwapCount         value="100"/>
    <logFileRoot            value="{directory}/chains"/>
    <reportCount            value="100"/>
    <sampleOutliers         value="false"/>
    <logFlushCount          value="100"/>
{modelParameters()}    <posteriorSampleState value="correlation">
      <acceptedStateCount value="100"/>
    </posteriorSampleState>
    <posteriorSampleStateInitialize value="latinHypercube">
      <maximinTrialCount value="100"/>
    </posteriorSampleStateInitialize>
    <posteriorSampleConvergence value="gelmanRubin">
      <thresholdHatR              value="1.2"/>
      <burnCount                  value="100"/>
      <testCount                  value="50"/>
      <outlierCountMaximum        value="1"/>
      <outlierSignificance        value="0.95"/>
      <outlierLogLikelihoodOffset value="60"/>
      <reportCount                value="10"/>
      <logFileName                value="{directory}/convergence.log"/>
    </posteriorSampleConvergence>
    <posteriorSampleStoppingCriterion value="stepCount">
      <stopAfterCount value="{countSteps}"/>
    </posteriorSampleStoppingCriterion>
    <posteriorSampleDffrntlEvltnRandomJump value="adaptive"/>
    <posteriorSampleDffrntlEvltnProposalSize value="adaptive">
      <logFileName           value="{directory}/proposalSize.log"/>
      <gammaInitial          value="0.5"/>
      <gammaAdjustFactor     value="1.1"/>
      <gammaMinimum          value="1.0e-4"/>
      <gammaMaximum          value="3.0"/>
      <acceptanceRateMinimum value="0.1"/>
      <acceptanceRateMaximum value="0.9"/>
      <updateCount           value="10"/>
    </posteriorSampleDffrntlEvltnProposalSize>
  </posteriorSampleSimulation>
  <randomNumberGenerator value="GSL">
    <seed          value="219"/>
    <mpiRankOffset value="true"/>
  </randomNumberGenerator>
</parameters>
''')
    return path


def run(command, logName):
    """Run a command, logging its output, and return its exit status."""
    with open(f'{directory}/{logName}', 'w') as log:
        return subprocess.run(command, stdout=log, stderr=subprocess.STDOUT).returncode


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--processesPerNode', type=int, default=1)
    parser.add_argument('--allowRunAsRoot', type=str, default='no')
    args, _ = parser.parse_known_args()
    try:
        from Galacticus.Emulation import emulatorFile as ef
        from Galacticus.Emulation.campaign import Campaign
        from Galacticus.Emulation.train import Observable, collect, train
    except ImportError as error:
        print(f'FAILED: unable to import the emulation modules (scikit-learn is required): {error}')
        return
    processes = max(args.processesPerNode, 1)
    allowRunAsRoot = ['--allow-run-as-root'] if args.allowRunAsRoot == 'yes' else []
    if args.allowRunAsRoot == 'yes':
        # Allow the MPI build of Galacticus to run as root without `mpirun` (as for the design and the campaign).
        os.environ['OMPI_ALLOW_RUN_AS_ROOT'] = '1'
        os.environ['OMPI_ALLOW_RUN_AS_ROOT_CONFIRM'] = '1'
    shutil.rmtree(directory, ignore_errors=True)
    os.makedirs(directory)

    # Generate the design.
    base = writeBase()
    if run(['./Galacticus.exe', writeDesign()], 'design.log') != 0:
        print(f'FAILED: emulatorDesign task failed (see {directory}/design.log)')
        return
    # Run the campaign.
    campaign = Campaign.create(f'{directory}/campaign.json', f'{directory}/design.hdf5', base)
    summary = campaign.run_local(parallel=processes, threads=1)
    if summary.get('complete', 0) != countDesign:
        failed = [run_ for run_ in campaign.runs if run_.status != 'complete']
        print(f'FAILED: {len(failed)} of {countDesign} campaign runs did not complete (e.g. run {failed[0].index}: {failed[0].message}; '
              f'see {failed[0].log_file})')
        return
    print(f'SUCCESS: campaign of {countDesign} runs completed')
    # Collect the training set and train the emulator.
    design, trainingSets, reports = collect(f'{directory}/design.hdf5', [Observable(label)])
    content = train(design, trainingSets, restarts=2, folds=0, seed=0)
    emulatorFileName = f'{directory}/emulator.hdf5'
    ef.write(emulatorFileName, content, overwrite=True)
    # Check the emulator against the exact model at the posterior mean (as a diagnostic of any failure below).
    meanTrue, covarianceTrue = posteriorAnalytic()
    emulator = content.emulators[label]
    quantiles = [(value - priors[name][0]) / (priors[name][1] - priors[name][0]) for name, value in zip(priors, meanTrue)]
    mean, variance = ef.predict(emulator, np.array([quantiles[list(priors).index(name)] for name in emulator.input_names]))
    errorEmulator = np.max(np.abs(mean.ravel() - (meanTrue[0] + meanTrue[1] * x))) / sigma
    print(f'emulator error at the posterior mean: {errorEmulator:.1e} of the data uncertainty '
          f'(predicted uncertainty {np.sqrt(np.max(variance)) / sigma:.1e})')
    # Sample the posterior.
    status = run(['mpirun', '--oversubscribe', '-np', str(countChains)] + allowRunAsRoot + ['./Galacticus.exe', writeSampler(emulatorFileName)],
                 'mcmc.log')
    if status != 0:
        print(f'FAILED: posterior sampling failed (see {directory}/mcmc.log)')
        return
    # Gather the converged states from the chains.
    states, logLikelihoods = [], []
    for path in sorted(glob.glob(f'{directory}/chains_[0-9][0-9][0-9][0-9].log')):
        with open(path) as file:
            for line in file:
                if line.startswith('#'):
                    continue
                columns = line.split()
                if columns[3] == 'T':
                    states.append([float(value) for value in columns[6:8]])
                    logLikelihoods.append(float(columns[5]))
    states = np.array(states)
    if states.shape[0] < 1000:
        print(f'FAILED: only {states.shape[0]} converged states were found in the chains (see {directory}/convergence.log)')
        return
    if not np.all(np.isfinite(logLikelihoods)):
        print('FAILED: non-finite log-likelihoods were found in the chains')
        return
    # Compare with the analytic posterior. The posterior found by the MCMC is slightly wider than the analytic posterior, as
    # the emulated likelihood includes the uncertainty of the emulator.
    meanMCMC = states.mean(axis=0)
    covarianceMCMC = np.cov(states.T)
    sigmaTrue = np.sqrt(np.diag(covarianceTrue))
    sigmaMCMC = np.sqrt(np.diag(covarianceMCMC))
    offset = np.abs(meanMCMC - meanTrue) / sigmaTrue
    ratio = sigmaMCMC / sigmaTrue
    correlationTrue = covarianceTrue[0, 1] / np.prod(sigmaTrue)
    correlationMCMC = covarianceMCMC[0, 1] / np.prod(sigmaMCMC)
    print(f'{states.shape[0]} converged states; posterior mean {meanMCMC} (analytic {meanTrue}); '
          f'widths {sigmaMCMC} (analytic {sigmaTrue}); correlation {correlationMCMC:.3f} (analytic {correlationTrue:.3f})')
    for name, offset_, ratio_ in zip(('intercept', 'slope'), offset, ratio):
        message = f'posterior mean of the {name} is {offset_:.2f} sigma from the analytic value'
        print(f'SUCCESS: {message}' if offset_ < 0.25 else f'FAILED: {message} (tolerance 0.25)')
        message = f'posterior width of the {name} is {ratio_:.3f} times the analytic value'
        print(f'SUCCESS: {message}' if 0.85 < ratio_ < 1.20 else f'FAILED: {message} (tolerance 0.85-1.20)')
    message = f'posterior correlation is {correlationMCMC:.3f} (analytic {correlationTrue:.3f})'
    print(f'SUCCESS: {message}' if abs(correlationMCMC - correlationTrue) < 0.1 else f'FAILED: {message} (tolerance 0.1)')


if __name__ == '__main__':
    try:
        main()
    except Exception as error:  # Report any unexpected error as a failure, rather than a bare traceback.
        print(f'FAILED: unexpected error: {type(error).__name__}: {error}')
    sys.exit(0)
