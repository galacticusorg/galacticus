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
  </emulatorTrain>

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
        )
        observable.validate()
        observables.append(observable)
    if not observables:
        raise ValueError('at least one <observable> must be given')
    return options, observables


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('configuration', help='the training configuration (an XML file)')
    parser.add_argument('--directory', default='.', help='the directory relative to which the run outputs are resolved')
    parser.add_argument('--overwrite', action='store_true', help='replace an existing emulator file')
    args = parser.parse_args(argv)
    options, observables = read_configuration(args.configuration)
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
                        store_cholesky=options['storeCholesky'])
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
