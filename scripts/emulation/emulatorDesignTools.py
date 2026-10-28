#!/usr/bin/env python3
"""Split or pool emulator designs whose runs are complete, reusing the outputs of their runs.

  split  write the parts of a design: part k holds the points whose index modulo the number of parts is k (by default,
         two parts: the even and the odd points). Emulators trained on each part give independent estimates of the
         posterior; their spread measures the systematic error of the emulators (see emulatorChainCompare.py).
  pool   write one design holding all the points and runs of several designs run with the same base parameter files and
         executable, with the quantiles of every point recomputed under the given priors (which must contain them all,
         e.g. the union of the boxes of the designs).

The run file names in the designs written are made absolute, resolving relative names against the directory from which
the runs were made (--directory, or --directories for pool; default: the current directory). An emulator is then trained
on a written design by naming it as the `designFileName` of a training configuration (see emulatorTrainReference.py).

Priors for pool are read from --priors: a design file, or a Galacticus parameter file holding the `modelParameter` blocks
of the parameters (e.g. that of the emulatorDesign task or of the posterior sampling). By default the priors of the first
design are used.

See the "Emulator-Assisted Calibration" chapter of the user guide.

Andrew Benson, Claude (28-October-2026).
"""

import argparse
import os
import sys

import lxml.etree as ET

sys.path.insert(0, os.path.join(os.environ.get('GALACTICUS_EXEC_PATH', '.'), 'python'))

from Galacticus.Emulation.design import pool_designs, read_design, resolve_paths, serialize_descriptor, split_design, write_design  # noqa: E402
from Galacticus.Emulation.emulatorFile import Prior  # noqa: E402


def priors_from_file(path, names):
    """Read the prior of each of ``names`` from a design file, or from the `modelParameter` blocks of a parameter file."""
    try:
        design = read_design(path).design
        priors = dict(zip(design.names, design.priors))
    except (OSError, ValueError):
        priors = {}
        for element in ET.parse(path).getroot().iter('modelParameter'):
            if element.get('value') != 'active' or element.find('name') is None:
                continue
            distribution = element.find('distributionFunction1D')
            if distribution is None:
                continue
            parameters = {child.tag: float(child.get('value')) for child in distribution if isinstance(child.tag, str)
                          and child.get('value') is not None and child.tag != 'randomNumberGenerator'}
            parameters['descriptor'] = serialize_descriptor(parameters)
            priors[element.find('name').get('value')] = Prior(distribution.get('value'), parameters)
    missing = [name for name in names if name not in priors]
    if missing:
        raise ValueError(f"'{path}' gives no prior for: {', '.join(missing)}")
    return [priors[name] for name in names]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    subparsers = parser.add_subparsers(dest='command', required=True)
    parser_split = subparsers.add_parser('split', help='write the parts of a design')
    parser_split.add_argument('design')
    parser_split.add_argument('--output', required=True, help='the root of the parts written: <output><k>.hdf5')
    parser_split.add_argument('--parts', type=int, default=2, help='the number of parts (default: 2)')
    parser_split.add_argument('--directory', default='.', help='the directory from which the runs were made')
    parser_split.add_argument('--overwrite', action='store_true')
    parser_pool = subparsers.add_parser('pool', help='pool several designs into one')
    parser_pool.add_argument('output')
    parser_pool.add_argument('designs', nargs='+')
    parser_pool.add_argument('--priors', help='a design or parameter file giving the priors of the pooled design')
    parser_pool.add_argument('--directories', nargs='+', help='the directory from which the runs of each design were made')
    parser_pool.add_argument('--overwrite', action='store_true')
    args = parser.parse_args(argv)

    if args.command == 'split':
        design = resolve_paths(read_design(args.design), args.directory)
        for part, content in enumerate(split_design(design, args.parts)):
            path = f'{args.output}{part}.hdf5'
            write_design(path, content, overwrite=args.overwrite)
            print(f'wrote part {part} ({content.design.quantiles.shape[0]} points, {len(content.runs["pointIndex"])} runs) to {path}')
        return 0

    directories = args.directories or ['.'] * len(args.designs)
    if len(directories) != len(args.designs):
        print('FAILED: give one directory per design')
        return 1
    designs = [resolve_paths(read_design(path), directory) for path, directory in zip(args.designs, directories)]
    priors = None if args.priors is None else priors_from_file(args.priors, designs[0].design.names)
    content = pool_designs(designs, priors)
    write_design(args.output, content, overwrite=args.overwrite)
    print(f'wrote pooled design ({content.design.quantiles.shape[0]} points, {len(content.runs["pointIndex"])} runs) to {args.output}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
