#!/usr/bin/env python3
"""Compare evolver-profiler histograms across runs.

Usage: python3 profcompare.py diag_m1.0e9 diag_m1.0e8 [...]
"""
import sys

import h5py
import numpy as np

PERF = "/home/user/perf"


def main():
    print(f"{'run':>16s} {'steps':>10s} {'interrupted':>11s} {'evals':>10s} {'cpu[s]':>8s} "
          f"{'f(dt<1e-3)':>10s} {'f(dt<1e-2)':>10s} {'cpuf(dt<1e-2)':>13s}")
    for name in sys.argv[1:]:
        with h5py.File(f"{PERF}/out/{name}.hdf5", "r") as f:
            g = f["metaData/evolverProfiler"]
            ts = np.asarray(g["timeStep"])
            n = np.asarray(g["timeStepCount"]) + np.asarray(g["timeStepCountInterrupted"])
            ev = np.asarray(g["evaluationCount"]) + np.asarray(g["evaluationCountInterrupted"])
            cpu = np.asarray(g["timeCPU"]) + np.asarray(g["timeCPUInterrupted"])
            ni = np.asarray(g["timeStepCountInterrupted"])
            fSmall3 = n[ts < 1e-3].sum() / max(n.sum(), 1)
            fSmall2 = n[ts < 1e-2].sum() / max(n.sum(), 1)
            fCPU2 = cpu[ts < 1e-2].sum() / max(cpu.sum(), 1e-9)
            print(f"{name:>16s} {n.sum():10d} {ni.sum():11d} {ev.sum():10d} {cpu.sum():8.2f} "
                  f"{fSmall3:10.3f} {fSmall2:10.3f} {fCPU2:13.3f}")
            hits = {x.decode(): h for x, h in zip(np.asarray(g["propertyNames"]), np.asarray(g["propertyHitCount"]))}
            top = sorted(hits.items(), key=lambda t: -t[1])[:3]
            print(f"{'':>16s} top limiters: " + ", ".join(f"{k}={v}" for k, v in top))


if __name__ == "__main__":
    main()
