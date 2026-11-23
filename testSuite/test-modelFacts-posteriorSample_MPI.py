#!/usr/bin/env python3
"""Check that the model fact registry is cleared, and correctly refilled, between models evaluated by a likelihood.

The `galaxyPopulation` likelihood builds a new model, from the same parameter document, for each state it evaluates, and
clears the model fact registry (`modelFactsReset()`) after each. The model used here (a single subhalo, with tidal heating
of its dark matter halo) requires the fact `satellite:tidalHeatingNormalized:evolved`, provided by the
`satelliteTidalHeating` node operator. A `grid` posterior sampling simulation, run on two MPI processes, evaluates two states, so that two models are
built and validated in turn on each process. It is checked that:

* with the node operator present, both models are evaluated - so the facts recorded by the second model satisfy its
  requirements after the registry was cleared following the first;
* with the node operator removed, evaluation stops with the error giving the reason for the requirement and naming the node
  operator that would satisfy it.

Andrew Benson, Claude
"""

import argparse
import os
import re
import shutil
import subprocess
import sys

root = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
os.environ['GALACTICUS_EXEC_PATH'] = root
os.chdir(root)

directory         = 'testSuite/outputs/modelFactsPosteriorSample'
parameterFileBase = 'testSuite/parameters/zeroVelocityDispersionDisruptingHalo.xml'
operatorElement   = '<nodeOperator value="satelliteTidalHeating"/>'
operatorPattern   = re.compile(r'<nodeOperator\s+value="satelliteTidalHeating"\s*/>')
countProcesses    = 2


def writeParameters(name, withOperator):
    """Write a base model (with or without the node operator), and a grid posterior sampling of it."""
    with open(parameterFileBase) as fileBase:
        parameters = fileBase.read()
    if len(operatorPattern.findall(parameters)) != 1:
        return None
    if not withOperator:
        parameters = operatorPattern.sub('', parameters)
    parameters = re.sub(r'<outputFileName\s+value="[^"]*"\s*/>', f'<outputFileName value="{directory}/{name}Model.hdf5"/>', parameters)
    with open(f'{directory}/{name}Base.xml', 'w') as file:
        file.write(parameters)
    # One grid step per process, so that each call to the likelihood evaluates (and resets between) two models.
    with open(f'{directory}/{name}.xml', 'w') as file:
        file.write(f'''<?xml version="1.0" encoding="UTF-8"?>
<parameters>
  <formatVersion>2</formatVersion>
  <verbosityLevel value="standard"/>
  <task value="posteriorSample">
    <initializeNodeClassHierarchy value="false"/>
  </task>
  <outputFileName value="{directory}/{name}.hdf5"/>
  <posteriorSampleLikelihood value="galaxyPopulation">
    <baseParametersFileName   value="{directory}/{name}Base.xml"/>
    <failedParametersFileName value="{directory}/{name}FailedParameters.xml"/>
    <evolveForestsVerbosity   value="standard"/>
  </posteriorSampleLikelihood>
  <posteriorSampleSimulation value="grid">
    <logFileRoot value="{directory}/{name}Chains"/>
    <posteriorSamples value="priorGrid">
      <countGrid value="{countProcesses}"/>
    </posteriorSamples>
    <modelParameter value="active">
      <name value="darkMatterProfileHeating/darkMatterProfileHeating/coefficientSecondOrder0"/>
      <distributionFunction1DPrior value="uniform">
        <limitLower value="0.20"/>
        <limitUpper value="0.22"/>
      </distributionFunction1DPrior>
      <operatorUnaryMapper value="identity"/>
      <distributionFunction1DPerturber value="cauchy">
        <median value="0.0"/>
        <scale  value="1.0e-3"/>
      </distributionFunction1DPerturber>
    </modelParameter>
  </posteriorSampleSimulation>
  <randomNumberGenerator value="GSL">
    <seed          value="219" />
    <mpiRankOffset value="true"/>
  </randomNumberGenerator>
</parameters>
''')
    return f'{directory}/{name}.xml'


def run(name, withOperator, allowRunAsRoot):
    """Run the posterior sampling, returning its exit status and log (or `None` if the parameters could not be written)."""
    path = writeParameters(name, withOperator)
    if path is None:
        return None, None
    with open(f'{directory}/{name}.log', 'w') as log:
        status = subprocess.run(['mpirun', '--oversubscribe', '-np', str(countProcesses)] + allowRunAsRoot + ['./Galacticus.exe', path],
                                stdout=log, stderr=subprocess.STDOUT, env=dict(os.environ, OMP_NUM_THREADS='1')).returncode
    with open(f'{directory}/{name}.log') as log:
        return status, log.read()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--processesPerNode', type=int, default=1)
    parser.add_argument('--allowRunAsRoot', type=str, default='no')
    args, _ = parser.parse_known_args()
    allowRunAsRoot = ['--allow-run-as-root'] if args.allowRunAsRoot == 'yes' else []
    shutil.rmtree(directory, ignore_errors=True)
    os.makedirs(directory)
    failures = []

    # With the provider, every model must be evaluated.
    status, log = run('provided', True, allowRunAsRoot)
    if status is None:
        print(f'FAILED: expected exactly one {operatorElement} in {parameterFileBase}')
        return
    countEvaluated = len(re.findall(r'Chain\s+\d+\s+has\s+log', log))
    if status != 0:
        failures.append(f'posterior sampling with the satelliteTidalHeating node operator failed (see {directory}/provided.log)')
    elif countEvaluated != countProcesses:
        failures.append(f'expected {countProcesses} models to be evaluated, but found {countEvaluated} (see {directory}/provided.log)')

    # Without the provider, evaluation must stop with the generated error.
    status, log = run('missing', False, allowRunAsRoot)
    if status == 0:
        failures.append('posterior sampling without the satelliteTidalHeating node operator ran to completion')
    else:
        for expected, description in (
                ('tidal heating of dark matter halos requires that the tidal heating of satellites is evolved', 'the reason for the requirement'),
                ("the fact 'satellite:tidalHeatingNormalized:evolved' is established by", 'the required fact'             ),
                (operatorElement                                                        , 'the XML selecting the provider'),
        ):
            if expected not in log:
                failures.append(f'the error report does not give {description} (see {directory}/missing.log)')

    if failures:
        for failure in failures:
            print(f'FAILED: {failure}')
    else:
        print('SUCCESS: model facts are cleared and refilled between models evaluated by a likelihood, and an unmet requirement is reported')


if __name__ == '__main__':
    main()
    sys.exit(0)
