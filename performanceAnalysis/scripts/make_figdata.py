#!/usr/bin/env python3
"""Assemble figures.json for the report from results.jsonl, walk logs, and perf shares."""
import json
import re

import numpy as np

import walkstats

PERF = "/home/user/perf"


def series(prefix, massTree):
    recs = []
    for line in open(f"{PERF}/results.jsonl"):
        r = json.loads(line)
        m = re.match(rf"^{re.escape(prefix)}_m([0-9.e+E]+)$", r["name"])
        if m and "timesEvolve" in r:
            mres = float(m.group(1))
            recs.append(dict(
                massResolution=mres,
                ratio=massTree / mres,
                nTree=len(r["timesEvolve"]),
                tEvolve=float(np.mean(r["timesEvolve"])),
                tEvolveErr=float(np.std(r["timesEvolve"]) / max(np.sqrt(len(r["timesEvolve"])), 1)),
                tConstruct=float(np.mean(r["timesConstruct"])),
                nodes=float(np.mean(r["countsNodes"])),
                walks=r.get("treeWalkCountMax", 0),
            ))
    recs.sort(key=lambda x: x["ratio"])
    return recs


def stall_anatomy(log, segment=0):
    segs = walkstats.parse(f"{PERF}/logs/{log}")
    if not segs:
        return None
    walks, live, evolved = segs[segment]
    return dict(walks=walks.tolist(), live=live.tolist(), evolved=evolved.tolist())


def main():
    fig = {
        "series1e13": series("base", 1.0e13),
        "series1e12": series("base12", 1.0e12),
        "seriesDefaults": series("def", 1.0e13),
        "stall": stall_anatomy("base_m1.0e8.log", 0),
    }
    # Variant comparison at m_res=1e8.
    variants = {}
    for line in open(f"{PERF}/results.jsonl"):
        r = json.loads(line)
        if "timesEvolve" in r and (r["name"].startswith("var") or r["name"] == "base_m1.0e8"):
            variants[r["name"]] = dict(tEvolve=float(np.mean(r["timesEvolve"])),
                                       walks=r.get("treeWalkCountMax", 0))
    fig["variants"] = variants
    with open(f"{PERF}/figures.json", "w") as f:
        json.dump(fig, f, indent=1)
    for key, val in fig.items():
        if isinstance(val, list):
            print(key, [f"M/m={x['ratio']:.1e}: tEvo={x['tEvolve']:.2f}s" for x in val])
        elif key == "variants":
            print(key, {k: f"{v['tEvolve']:.2f}s" for k, v in val.items()})


if __name__ == "__main__":
    main()
