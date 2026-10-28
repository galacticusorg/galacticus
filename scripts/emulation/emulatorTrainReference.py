#!/usr/bin/env python3
"""Train emulators of Galacticus model predictions from a campaign of runs over an emulator design.

The configuration is an XML file:

  <emulatorTrain>
    <designFileName       value="design.hdf5"  />
    <emulatorFileName     value="emulator.hdf5"/>
    <pcaVarianceRetained  value="0.99"/>      <!-- optional -->
    <restartsOptimizer    value="4"   />      <!-- optional -->
    <foldsCrossValidation value="5"   />      <!-- optional; 0 disables cross-validation -->
    <seed                 value="0"   />      <!-- optional -->
    <storeCholesky        value="true"/>      <!-- optional -->
    <collectOnly          value="false"/>     <!-- optional; write the training sets only -->
    <observable label="massFunctionStellarTomczak2014ZFOURGEz0" transform="log10" floor="-6.64" rootVarianceFloored="0.5"/>
    <observable label="massMetallicityBlanc2019" transform="identity" undefined="median" rootVarianceUndefined="5.0"/>
    <observable label="progenitorMassFunctionM1Z0" transform="log10" floor="-13.6" rootVarianceFloored="1.0" likelihoodRangeOnly="true"/>
  </emulatorTrain>

Each observable may also give ``rootVarianceZero`` (the root variance given to values whose root variance is zero) and
``rootVarianceTargetMinimum`` (a floor on the root variance of the target data). With ``pattern="true"`` the label is a
regular expression, matched in full against the labels of the analyses in the first run output found, and the observable
stands for every analysis it matches (each with the same options).

The outputs of the runs are those named in the design file, resolved relative to the working directory (or to
--directory). The emulator file format is described in the "Emulator-Assisted Calibration" chapter of the user guide.

Andrew Benson, Claude (5-October-2026).
"""

import argparse
import os
import sys

import lxml.etree as ET
import numpy as np

sys.path.insert(0, os.path.join(os.environ.get('GALACTICUS_EXEC_PATH', '.'), 'python'))

from Galacticus.Emulation import emulatorFile as ef  # noqa: E402
from Galacticus.Emulation.train import Observable, collect, train  # noqa: E402


def _value(task, name, default=None, kind=str):
    element = task.find(name)
    if element is None:
        if default is None:
            raise ValueError(f'`{name}` is required')
        return default
    text = element.get('value')
    if kind is bool:
        return text.strip().lower() in ('true', 't', '1', 'yes')
    return kind(text)


def read_configuration(path):
    """Read a training configuration, returning a dictionary of options and a list of observables."""
    task = ET.parse(path).getroot()
    if task.tag != 'emulatorTrain':
        raise ValueError(f"'{path}' is not an emulator training configuration (its root element must be <emulatorTrain>)")
    options = {
        'designFileName':       _value(task, 'designFileName'),
        'emulatorFileName':     _value(task, 'emulatorFileName'),
        'pcaVarianceRetained':  _value(task, 'pcaVarianceRetained' , 0.99 , float),
        'restartsOptimizer':    _value(task, 'restartsOptimizer'   , 4    , int  ),
        'foldsCrossValidation': _value(task, 'foldsCrossValidation', 5    , int  ),
        'seed':                 _value(task, 'seed'                , 0    , int  ),
        'storeCholesky':        _value(task, 'storeCholesky'       , True , bool ),
        'collectOnly':          _value(task, 'collectOnly'         , False, bool ),
    }
    observables = []
    for element in task.findall('observable'):

        def number(name):
            text = element.get(name)
            return None if text is None else float(text)

        observable = Observable(
            label=element.get('label'),
            transform=element.get('transform', 'identity'),
            floor=number('floor'),
            root_variance_floored=number('rootVarianceFloored'),
            undefined=element.get('undefined'),
            root_variance_undefined=number('rootVarianceUndefined'),
            root_variance_zero=number('rootVarianceZero'),
            root_variance_target_minimum=number('rootVarianceTargetMinimum'),
            likelihood_range_only=element.get('likelihoodRangeOnly', 'false').strip().lower() in ('true', 't', '1', 'yes'),
        )
        if element.get('pattern', 'false').strip().lower() in ('true', 't', '1', 'yes'):
            observable.pattern = True
        else:
            observable.validate()
        observables.append(observable)
    if not observables:
        raise ValueError('at least one <observable> must be given')
    return options, observables


def expand_patterns(observables, design_file, directory='.'):
    """Replace each observable whose label is a pattern by one observable per matching analysis label, found in the
    first run output (of the design's runs) which exists."""
    import re
    from dataclasses import replace
    import h5py
    from Galacticus.Emulation.campaign import output_file_path
    from Galacticus.Emulation.design import read_design
    if not any(getattr(observable, 'pattern', False) for observable in observables):
        return observables
    labels = None
    for name in read_design(design_file).runs['outputFileName']:
        path = output_file_path(os.path.join(directory, name))
        try:
            with h5py.File(path, 'r') as file:
                labels = sorted(file['analyses'].keys()) if 'analyses' in file else []
            break
        except OSError:
            continue
    if labels is None:
        raise ValueError('no run output exists in which to match observable label patterns')
    expanded = []
    for observable in observables:
        if not getattr(observable, 'pattern', False):
            expanded.append(observable)
            continue
        matched = [label for label in labels if re.fullmatch(observable.label, label)]
        if not matched:
            raise ValueError(f"no analysis label matches the pattern '{observable.label}'")
        for label in matched:
            copy = replace(observable, label=label)
            copy.validate()
            expanded.append(copy)
    return expanded


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('configuration', help='the training configuration (an XML file)')
    parser.add_argument('--directory', default='.', help='the directory relative to which the run outputs are resolved')
    parser.add_argument('--overwrite', action='store_true', help='replace an existing emulator file')
    parser.add_argument('--processes', type=int, default=1,
                        help='the number of processes in which to fit emulators (and cross-validation folds) in parallel')
    args = parser.parse_args(argv)
    options, observables = read_configuration(args.configuration)
    observables = expand_patterns(observables, options['designFileName'], args.directory)
    design, training_sets, reports = collect(options['designFileName'], observables, directory=args.directory)
    for label, report in reports.items():
        print(f'{label}: {report}')
        masked = training_sets[label].mask
        if masked.any():
            print(f'  {int(masked.sum())} of {masked.size} values floored or replaced')
    if options['collectOnly']:
        content = ef.EmulatorFile(design=design, training_sets=training_sets, emulators={},
                                  attributes={'creator': 'Galacticus.Emulation.train (collect only)'})
    else:
        content = train(design, training_sets, pca_variance_retained=options['pcaVarianceRetained'],
                        restarts=options['restartsOptimizer'], folds=options['foldsCrossValidation'], seed=options['seed'],
                        store_cholesky=options['storeCholesky'], processes=args.processes)
        for label, emulator in content.emulators.items():
            line = f'{label}: {len(emulator.components)} principal components'
            if label in content.validation:
                validation = content.validation[label]
                line += (f'; cross-validation: median RMSE {np.median(validation.rmse):.3g},'
                         f' median 1/2-sigma coverage {np.median(validation.coverage_1_sigma):.2f}/{np.median(validation.coverage_2_sigma):.2f}'
                         f' (nominal 0.68/0.95)')
            print(line)
    ef.write(options['emulatorFileName'], content, overwrite=args.overwrite)
    print(f"wrote {options['emulatorFileName']}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
