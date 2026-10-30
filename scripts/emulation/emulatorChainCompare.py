#!/usr/bin/env python3
"""Compare posteriors sampled with different emulators, to measure the systematic error of the emulators.

Train emulators on the whole of a design and on each of its parts (see emulatorDesignTools.py split), sample the posterior
with each, and compare:

  emulatorChainCompare.py <reference log file root> <log file root> [<log file root> ...]

For each parameter this reports the median and half-width (half the 16-84 percentile range) of the reference posterior,
the shift of the median of each other posterior from it, and the systematic error (the largest absolute shift) relative
to the half-width. A systematic comparable to the half-width means that the posterior is limited by the emulators: add
design points (e.g. by pooling a second design; see emulatorDesignTools.py pool), and quote the systematic with the result.

States flagged as converged are used, if any are; otherwise the second half of each chain (see --burnFraction).

Andrew Benson, Claude (28-October-2026).
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.environ.get('GALACTICUS_EXEC_PATH', '.'), 'python'))

from Galacticus.Emulation.chains import compare_posteriors, read_chains  # noqa: E402


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('reference', help='the log file root of the reference chains (e.g. emulators trained on the whole design)')
    parser.add_argument('others', nargs='+', help='the log file roots of the chains to compare (e.g. emulators trained on its parts)')
    parser.add_argument('--burnFraction', type=float, default=0.5,
                        help='the fraction of each chain discarded if no state is flagged as converged (default: 0.5)')
    parser.add_argument('--all', action='store_true', help='use all states after the burn fraction, ignoring convergence flags')
    args = parser.parse_args(argv)
    reference = read_chains(args.reference, converged_only=not args.all, burn_fraction=args.burnFraction)
    others = [read_chains(root, converged_only=not args.all, burn_fraction=args.burnFraction) for root in args.others]
    print(f'reference: {args.reference} ({reference[1].shape[0]} states)')
    for k, (root, other) in enumerate(zip(args.others, others)):
        print(f'posterior {k + 1}: {root} ({other[1].shape[0]} states)')
    width = max(len(name) for name in reference[0])
    print(f"{'parameter':{width}}  {'median':>12} {'sigma':>10}  " + ' '.join(f'{f"shift {k + 1}":>10}' for k in range(len(others)))
          + f"  {'systematic':>10} {'/sigma':>7}")
    for row in compare_posteriors(reference, others):
        print(f"{row['name']:{width}}  {row['median']:12.5g} {row['sigma']:10.3g}  " + ' '.join(f'{shift:+10.3g}' for shift in row['shifts'])
              + f"  {row['systematic']:10.3g} {row['ratio']:7.2f}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
