#!/usr/bin/env python3
"""Test generation of an emulator design, and a campaign of model runs over it.

Runs the `emulatorDesign` task (testSuite/parameters/emulatorDesign.xml: three parameters of parameters/quickTest.xml
with uniform, log-uniform and truncated log-normal priors, and a 16-point Sobol design), then checks:

* the design file: the prior quantile recorded for each value agrees with an independent computation of each prior's
  cumulative distribution, and the design is balanced (one point in each 1/16 interval of every quantile);
* the change files: each sets its parameters to exactly the values in the design file;
* validation: a Galacticus dry run accepts the runs at the extremes of each parameter, and rejects a run which violates a
  declared parameter bound, reporting the bound;
* running: two runs of quickTest.xml complete, and each model was run with exactly the design's parameter values.

Andrew Benson, Claude (23-September-2026).
"""

import math
import os
import sys

import h5py
import numpy as np

root = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
os.environ['GALACTICUS_EXEC_PATH'] = root
sys.path.insert(0, os.path.join(root, 'python'))
os.chdir(root)

import subprocess  # noqa: E402

import lxml.etree as ET  # noqa: E402

from Galacticus.Emulation.campaign import Campaign  # noqa: E402
from Galacticus.Emulation.design import read_design  # noqa: E402

directory = 'testSuite/outputs/emulatorDesign'
countPoints = 16


def normalCDF(x):
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def priorCumulative(prior, x):
    """An independent implementation of the cumulative distribution of each prior used here."""
    p = prior.parameters
    if prior.class_name == 'uniform':
        return (x - p['limitLower']) / (p['limitUpper'] - p['limitLower'])
    if prior.class_name == 'logUniform':
        return math.log(x / p['limitLower']) / math.log(p['limitUpper'] / p['limitLower'])
    if prior.class_name == 'logNormal':
        # Truncated log-normal, parameterized by its median x0 and the standard deviation of ln x.
        z = lambda y: math.log(y / p['x0']) / p['sigma']  # noqa: E731
        low, high = normalCDF(z(p['limitLower'])), normalCDF(z(p['limitUpper']))
        return (normalCDF(z(x)) - low) / (high - low)
    # No reference implementation: return NaN, which fails the comparison.
    return float('nan')


def modelParameter(parameters, name):
    """Read a (possibly `[@value='...']`-selected) parameter from the `Parameters` group of a model output file."""
    group = parameters
    segments = name.split('/')
    for segment in segments[:-1]:
        if '[@value=' in segment:
            # Multiple copies of a parameter are written as groups `tag[N]`, with their values in attributes of the same names.
            tag, selector = segment.split('[@value=')
            value = selector.strip("']")
            key = [k for k, v in group.attrs.items() if k.startswith(tag + '[') and v.decode() == value][0]
            group = group[key]
        else:
            group = group[segment]
    return group.attrs[segments[-1]]


def check(condition, success, failure):
    print(f'SUCCESS: {success}' if condition else f'FAILED: {failure}')
    return condition


def main():
    # Generate the design.
    subprocess.run(['rm', '-rf', directory])
    with open('testSuite/outputs/test-emulator-design.log', 'w') as log:
        status = subprocess.run(['./Galacticus.exe', 'testSuite/parameters/emulatorDesign.xml'], stdout=log, stderr=subprocess.STDOUT).returncode
    if not check(status == 0, 'emulatorDesign task ran', f'emulatorDesign task exited with status {status} (see testSuite/outputs/test-emulator-design.log)'):
        return
    design = read_design(f'{directory}/design.hdf5')
    values = design.design.values
    quantiles = design.design.quantiles
    check(values.shape == (countPoints, 3) and design.count_runs == countPoints,
          'design has 16 points of 3 parameters, and 16 runs', f'design has shape {values.shape} and {design.count_runs} runs')
    check(np.all((quantiles > 0.0) & (quantiles < 1.0)), 'quantiles lie in (0,1)', 'quantiles lie outside (0,1)')
    # Quantiles against an independent computation of each prior's cumulative distribution.
    errors = [abs(priorCumulative(design.design.priors[j], values[i, j]) - quantiles[i, j]) for i in range(countPoints) for j in range(3)]
    errors = np.array(errors)
    check(np.max(errors) < 1.0e-8, f'quantiles agree with the prior distributions (maximum difference {np.max(errors):.1e})',
          f'quantiles disagree with the prior distributions (maximum difference {np.max(errors):.1e})')
    bins = np.minimum((quantiles * countPoints).astype(int), countPoints - 1)
    check(all(len(set(bins[:, j])) == countPoints for j in range(3)), 'design is balanced', 'design is not balanced')
    # Change files set exactly the design's values.
    exact = True
    for index in range(design.count_runs):
        changes = ET.parse(design.runs['changeFileName'][index]).getroot()
        updates = {change.get('path'): float(change.get('value')) for change in changes if change.get('type') == 'update'}
        point = design.runs['pointIndex'][index]
        exact = exact and all(updates[name] == values[point, j] for j, name in enumerate(design.design.names))
    check(exact, 'change files set exactly the design values', 'change files do not set exactly the design values')

    # Validation by dry run.
    campaign = Campaign.create(f'{directory}/campaign.json', f'{directory}/design.hdf5', 'parameters/quickTest.xml')
    results = campaign.validate(parallel=2)
    check(all(result.ok for result in results), f'dry runs accept the {len(results)} runs at the parameter extremes',
          'dry runs reject runs: ' + '; '.join(f'run {r.run}: {r.message}' for r in results if not r.ok))
    # A run which violates a declared bound must be rejected, with the bound reported.
    changeFile = design.runs['changeFileName'][1]
    with open(changeFile) as handle:
        original = handle.read()
    with open(changeFile, 'w') as handle:
        handle.write(original.replace('</changes>', '  <change type="update" path="cosmologyParameters/OmegaMatter" value="-0.3"/>\n</changes>'))
    try:
        result = campaign.validate(indices=[1])[0]
    finally:
        with open(changeFile, 'w') as handle:
            handle.write(original)
    check(not result.ok and 'must be greater than' in result.message, 'dry run rejects a run which violates a parameter bound',
          f'dry run did not reject a run violating a parameter bound as expected (ok={result.ok}, message="{result.message}")')

    # Run two points, and check that each model used exactly the design's values.
    summary = campaign.run_local(indices=[0, 1], parallel=2, threads=1)
    if not check(summary.get('complete', 0) == 2, 'two runs of quickTest.xml completed',
                 'runs did not complete: ' + '; '.join(f'run {run.index}: {run.message} [log: {run.log_file}]' for run in campaign.runs[:2] if run.status != 'complete')):
        return
    matched = True
    for run in campaign.runs[:2]:
        with h5py.File(run.output_file, 'r') as model:
            used = [modelParameter(model['Parameters'], name) for name in design.design.names]
        matched = matched and np.array_equal(np.array(used, dtype=float), values[run.point])
    check(matched, 'models ran with exactly the design values', 'models did not run with the design values')


if __name__ == '__main__':
    try:
        main()
    except Exception as error:  # Report any unexpected error as a failure, rather than a bare traceback.
        print(f'FAILED: unexpected error: {type(error).__name__}: {error}')
    sys.exit(0)
