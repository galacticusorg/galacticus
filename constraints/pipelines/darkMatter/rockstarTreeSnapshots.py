#!/usr/bin/env python3
"""Extract the unique snapshots from a Rockstar/consistent-trees tree file.

Given a directory, this script locates a Rockstar tree file (e.g.
`tree_0_0_0.dat`) within it, scans the halo rows for the distinct
(snapshot number, expansion factor) pairs, and writes a `snapshots.txt`
listing the snapshot ID, redshift, and expansion factor back to that
directory.

Usage:
    rockstarTreeSnapshots.py <path>

Andrew Benson (2026)
"""

import argparse
import glob
import os
import sys


def find_tree_file(path):
    """Return the path to a Rockstar tree file found in `path`.

    Some tree files can be empty, so the largest file is chosen to ensure
    there is snapshot data to extract.
    """
    candidates = glob.glob(os.path.join(path, "tree_*.dat"))
    if not candidates:
        raise FileNotFoundError(f"no Rockstar tree file (tree_*.dat) found in '{path}'")
    return max(candidates, key=os.path.getsize)


def extract_snapshots(tree_file):
    """Scan `tree_file`, returning a dict of {snapshot_id: expansion_factor}.

    In consistent-trees output the column labels in the header carry their
    (zero-based) column index, so `scale(0)` is the expansion factor and
    `Snap_num(31)` is the snapshot number.
    """
    scale_column = 0
    snap_column = 31
    snapshots = {}
    with open(tree_file, "r") as tree:
        for line in tree:
            # Skip header/comment lines and the per-tree marker lines.
            if not line or line[0] == "#":
                continue
            fields = line.split()
            # Data rows have many columns; the tree-count and "#tree" marker
            # lines have too few, so guard against them.
            if len(fields) <= snap_column:
                continue
            try:
                scale = float(fields[scale_column])
                snap_id = int(fields[snap_column])
            except ValueError:
                continue
            # Record the expansion factor for this snapshot (identical across
            # all halos sharing the snapshot).
            snapshots.setdefault(snap_id, scale)
    return snapshots


def write_snapshots(snapshots, output_file):
    """Write the snapshot table (ID, redshift, expansion factor) to disk."""
    with open(output_file, "w") as out:
        out.write("# snapshotID    redshift    expansionFactor\n")
        for snap_id in sorted(snapshots):
            scale = snapshots[snap_id]
            redshift = 1.0 / scale - 1.0
            out.write(f"{snap_id:12d} {redshift:15.8f} {scale:15.8f}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Extract unique snapshots from a Rockstar tree file."
    )
    parser.add_argument(
        "path",
        help="directory containing a Rockstar tree file; snapshots.txt is written here",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="re-extract snapshots even if snapshots.txt already exists",
    )
    args = parser.parse_args()

    if not os.path.isdir(args.path):
        sys.exit(f"error: '{args.path}' is not a directory")

    output_file = os.path.join(args.path, "snapshots.txt")
    if os.path.exists(output_file) and not args.force:
        print(f"Snapshot file already exists (use --force to regenerate): {output_file}")
        return

    tree_file = find_tree_file(args.path)
    print(f"Reading tree file: {tree_file}")

    snapshots = extract_snapshots(tree_file)
    if not snapshots:
        sys.exit(f"error: no snapshot data found in '{tree_file}'")

    write_snapshots(snapshots, output_file)
    print(f"Wrote {len(snapshots)} snapshots to: {output_file}")


if __name__ == "__main__":
    main()
