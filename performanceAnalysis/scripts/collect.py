#!/usr/bin/env python3
"""Collect timing and instrumentation results from scaling-experiment runs.

Usage: python3 collect.py <runName> [<runName> ...]
Emits one JSON record per run to results.jsonl (appending, deduplicated by name)
and prints a human-readable summary table.
"""
import json
import os
import re
import sys

import h5py
import numpy as np

PERF = "/home/user/perf"


def parse_log(path):
    out = {}
    if not os.path.exists(path):
        return out
    text = open(path, errors="replace").read()
    m = re.search(r"Elapsed \(wall clock\) time.*: (.*)", text)
    if m:
        parts = m.group(1).strip().split(":")
        seconds = 0.0
        for p in parts:
            seconds = seconds * 60 + float(p)
        out["wallSeconds"] = seconds
    m = re.search(r"Maximum resident set size \(kbytes\): (\d+)", text)
    if m:
        out["maxRSSMB"] = int(m.group(1)) / 1024.0
    m = re.search(r"User time \(seconds\): ([\d.]+)", text)
    if m:
        out["userSeconds"] = float(m.group(1))
    # Tree walk counts: "Evolving tree [    N]" lines are throttled (printed when
    # count grew by >10%), so the maximum seen is a lower bound within ~10%.
    walks = [int(w) for w in re.findall(r"Evolving tree \[\s*(\d+)\]", text)]
    if walks:
        out["treeWalkCountMax"] = max(walks)
    out["failed"] = bool(re.search(r"FATAL|Aborted|exit=[1-9]", text))
    return out


def read_hdf5(path):
    out = {}
    if not os.path.exists(path):
        return out
    with h5py.File(path, "r") as f:
        if "metaData/treeTiming" in f:
            g = f["metaData/treeTiming"]
            for key, name in [("treeMass", "treeMasses"), ("timeConstruct", "timesConstruct"),
                              ("timeEvolve", "timesEvolve"), ("countNodes", "countsNodes"),
                              ("treeID", "treeIDs")]:
                if key in g:
                    out[name] = np.asarray(g[key]).tolist()
        if "metaData/evolverProfiler" in f:
            g = f["metaData/evolverProfiler"]
            prof = {}
            for key in g.keys():
                arr = np.asarray(g[key])
                if arr.dtype.kind in "SO":
                    arr = [x.decode() if isinstance(x, bytes) else str(x) for x in arr]
                    prof[key] = arr
                else:
                    prof[key] = arr.tolist()
            out["evolverProfiler"] = prof
    return out


def collect(name):
    rec = {"name": name}
    rec.update(parse_log(f"{PERF}/logs/{name}.log"))
    rec.update(read_hdf5(f"{PERF}/out/{name}.hdf5"))
    return rec


def main():
    names = sys.argv[1:]
    results_path = f"{PERF}/results.jsonl"
    existing = {}
    if os.path.exists(results_path):
        with open(results_path) as f:
            for line in f:
                r = json.loads(line)
                existing[r["name"]] = r
    for name in names:
        existing[name] = collect(name)
    with open(results_path, "w") as f:
        for r in existing.values():
            f.write(json.dumps(r) + "\n")
    # Summary table.
    print(f"{'run':32s} {'wall[s]':>9s} {'RSS[MB]':>8s} {'walks':>7s} {'nodes(mean)':>12s} {'tCon(mean)':>11s} {'tEvo(mean)':>11s}")
    for name in names:
        r = existing[name]
        nodes = np.mean(r["countsNodes"]) if "countsNodes" in r else float("nan")
        tcon = np.mean(r["timesConstruct"]) if "timesConstruct" in r else float("nan")
        tevo = np.mean(r["timesEvolve"]) if "timesEvolve" in r else float("nan")
        print(f"{name:32s} {r.get('wallSeconds', float('nan')):9.1f} {r.get('maxRSSMB', float('nan')):8.0f} "
              f"{r.get('treeWalkCountMax', 0):7d} {nodes:12.0f} {tcon:11.2f} {tevo:11.2f}")


if __name__ == "__main__":
    main()
