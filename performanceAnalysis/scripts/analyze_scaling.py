#!/usr/bin/env python3
"""Fit local power-law slopes of tree build/evolve CPU time vs mass resolution.

Usage: python3 analyze_scaling.py <prefix>
Reads results.jsonl records whose names start with <prefix>, expecting names of
the form <prefix>_m<massResolution>. Prints a table of per-resolution means and
local logarithmic slopes d log(time) / d log(1/m_res).
"""
import json
import re
import sys

import numpy as np

PERF = "/home/user/perf"


def main():
    prefix = sys.argv[1]
    recs = []
    with open(f"{PERF}/results.jsonl") as f:
        for line in f:
            r = json.loads(line)
            m = re.match(rf"^{re.escape(prefix)}_m([0-9.e+E]+)$", r["name"])
            if m and "timesEvolve" in r:
                r["massResolution"] = float(m.group(1))
                recs.append(r)
    recs.sort(key=lambda r: -r["massResolution"])
    if not recs:
        print("no records found")
        return
    rows = []
    for r in recs:
        te = np.array(r["timesEvolve"])
        tc = np.array(r["timesConstruct"])
        nn = np.array(r["countsNodes"], dtype=float)
        rows.append(dict(mres=r["massResolution"], nTree=len(te),
                         tEvolve=np.mean(te), tEvolveMed=np.median(te),
                         tConstruct=np.mean(tc), nodes=np.mean(nn),
                         walks=r.get("treeWalkCountMax", 0),
                         wall=r.get("wallSeconds", float("nan"))))
    print(f"{'m_res':>10s} {'nTree':>5s} {'nodes':>10s} {'tCon[s]':>9s} {'tEvo[s]':>9s} "
          f"{'walks':>7s} | {'a(nodes)':>8s} {'a(tCon)':>8s} {'a(tEvo)':>8s}")
    for i, row in enumerate(rows):
        slopes = ["", "", ""]
        if i > 0:
            p = rows[i - 1]
            dlm = np.log10(p["mres"] / row["mres"])  # decades of resolution increase
            for j, key in enumerate(["nodes", "tConstruct", "tEvolve"]):
                if row[key] > 0 and p[key] > 0:
                    slopes[j] = f"{np.log10(row[key] / p[key]) / dlm:8.2f}"
        print(f"{row['mres']:10.2e} {row['nTree']:5d} {row['nodes']:10.0f} {row['tConstruct']:9.2f} "
              f"{row['tEvolve']:9.2f} {row['walks']:7d} | {slopes[0]:>8s} {slopes[1]:>8s} {slopes[2]:>8s}")
    # Global fits over the top decade vs bottom decade.
    lm = np.log10([row["mres"] for row in rows])
    lt = np.log10([row["tEvolve"] for row in rows])
    if len(rows) >= 3:
        coeffs = np.polyfit(lm, lt, 1)
        print(f"\nGlobal power law: tEvolve ~ m_res^{coeffs[0]:.2f}  "
              f"(i.e. ~ (1/m_res)^{-coeffs[0]:.2f})")


if __name__ == "__main__":
    main()
