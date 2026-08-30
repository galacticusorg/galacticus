#!/usr/bin/env python3
"""Parse per-walk progress reports from a run log (verbosityLevel=warn).

Each throttled report carries: walk counter, live nodes in tree, nodes evolved
in that walk. Usage: python3 walkstats.py logs/<name>.log [...]
"""
import re
import sys

import numpy as np


def parse(path):
    """Return a list of per-tree (walks, live, evolved) arrays.

    The walk counter resets for each tree (and each evolve call), so a decrease
    in the counter marks a new segment.
    """
    segments = []
    walks, live, evolved = [], [], []
    pat_w = re.compile(r"Evolving tree \[\s*(\d+)\]")
    pat_n = re.compile(r"Nodes in tree:\s*(\d+)")
    pat_e = re.compile(r"Nodes evolved:\s*(\d+)")
    w = n = None
    for line in open(path, errors="replace"):
        m = pat_w.search(line)
        if m:
            w = int(m.group(1))
            if walks and w < walks[-1]:
                segments.append((np.array(walks), np.array(live), np.array(evolved)))
                walks, live, evolved = [], [], []
            continue
        m = pat_n.search(line)
        if m and w is not None:
            n = int(m.group(1))
            continue
        m = pat_e.search(line)
        if m and w is not None and n is not None:
            walks.append(w)
            live.append(n)
            evolved.append(int(m.group(1)))
            w = n = None
    if walks:
        segments.append((np.array(walks), np.array(live), np.array(evolved)))
    return segments


def main():
    for path in sys.argv[1:]:
        segments = parse(path)
        if not segments:
            print(f"{path}: no walk reports")
            continue
        print(f"{path}: {len(segments)} tree segments")
        totalVisits = totalEvolved = 0.0
        for i, (walks, live, evolved) in enumerate(segments):
            frac = evolved / np.maximum(live, 1)
            dw = np.diff(np.concatenate([[0], walks]))
            visits = np.sum(dw * live)
            evs = np.sum(dw * evolved)
            stalled = frac < 0.10
            share = np.sum(dw[stalled] * live[stalled]) / max(visits, 1)
            totalVisits += visits
            totalEvolved += evs
            print(f"  tree {i}: walks~{walks.max():5d} liveMed={np.median(live):6.0f} liveMax={live.max():6d} "
                  f"visits~{visits:.3g} evolved~{evs:.3g} stalledVisitShare={share:.0%}")
        print(f"  TOTAL visits~{totalVisits:.3g} evolved~{totalEvolved:.3g} ratio={totalVisits/max(totalEvolved,1):.1f}")


if __name__ == "__main__":
    main()
