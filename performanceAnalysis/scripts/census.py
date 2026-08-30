#!/usr/bin/env python3
"""Count live nodes / subhalos per output time in a standard-outputter file.

Usage: python3 census.py out/<name>.hdf5
"""
import sys

import h5py
import numpy as np


def main():
    path = sys.argv[1]
    with h5py.File(path, "r") as f:
        outputs = f["Outputs"]
        print(f"{'output':>8s} {'z':>7s} {'nodes':>9s} {'subhalos':>9s} {'centrals':>9s}")
        for name in sorted(outputs.keys(), key=lambda s: int(s.replace("Output", ""))):
            g = outputs[name]
            z = 1.0 / g.attrs["outputExpansionFactor"] - 1.0 if "outputExpansionFactor" in g.attrs else float("nan")
            nd = g["nodeData"]
            # nodeIsIsolated: 1 for centrals/isolated halos, 0 for subhalos.
            isolated = np.asarray(nd["nodeIsIsolated"]) if "nodeIsIsolated" in nd else None
            n = len(next(iter(nd.values()))) if len(nd) else 0
            nSub = int((isolated == 0).sum()) if isolated is not None else -1
            nCen = int((isolated == 1).sum()) if isolated is not None else -1
            print(f"{name:>8s} {z:7.2f} {n:9d} {nSub:9d} {nCen:9d}")


if __name__ == "__main__":
    main()
