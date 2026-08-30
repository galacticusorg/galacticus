#!/usr/bin/env python3
"""Extract top function shares from `perf report --stdio` output.

Usage: python3 perfshares.py out/<name>.perf.flat.txt [more ...]
Prints aligned share tables and writes perfshares.json for the report figures.
"""
import json
import re
import sys


def parse(path):
    rows = []
    for line in open(path, errors="replace"):
        # Format: "    12.34%  Galacticus.exe  Galacticus.exe  [.] __module_MOD_proc"
        m = re.match(r"\s+([\d.]+)%\s+\S+\s+(\S+)\s+\[[.k]\]\s+(.+)$", line)
        if m:
            rows.append({"share": float(m.group(1)), "dso": m.group(2), "symbol": m.group(3).strip()})
    return rows


def main():
    out = {}
    for path in sys.argv[1:]:
        rows = parse(path)
        out[path] = rows
        print(f"== {path} (top 25) ==")
        for r in rows[:25]:
            print(f"  {r['share']:6.2f}%  {r['symbol'][:100]}")
    with open("/home/user/perf/perfshares.json", "w") as f:
        json.dump(out, f, indent=1)


if __name__ == "__main__":
    main()
