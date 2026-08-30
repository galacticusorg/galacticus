#!/usr/bin/env python3
"""Aggregate perf flat-profile shares into subsystem groups, compare runs.

Usage: python3 perfgroup.py prof_m1.0e9 prof_m1.0e8 [...]
"""
import re
import sys

PERF = "/home/user/perf"

GROUPS = [
    ("tree walk & evolvability", ["allnodesnext", "walker", "nodeisevolvable", "issatellite"]),
    ("timestep limits (timeEvolveTo)", ["timeevolveto", "evolve_to_time", "refusetoevolve"]),
    ("evolver-call overhead", ["serializ", "deserial", "serialization", "odestepinactive", "odestepanalytic",
                               "odestepscales", "calculations_reset", "calculationreset",
                               "differentialevolutionpre", "differentialevolutionpost", "standardevolve"]),
    ("ODE solver core", ["gsl_odeiv2", "rkf45", "odesolver", "solver_", "brent_iterate", "root_finder",
                         "rootfinder", "rootfunction"]),
    ("satellite/tidal physics RHS", ["tidal", "velocitydispersion", "kinematic", "chandrasekhar",
                                     "dynamicalfriction", "orbit", "heated", "heating"]),
    ("mass distributions / profiles", ["massdistribution", "nfw", "darkmatterprofile", "darkmatteronly",
                                       "kepler", "concentration"]),
    ("tree surgery & lists", ["promote", "lastsatellite", "removefromhost", "movecomponents", "destroy",
                              "mergee", "merge", "satellite_move"]),
    ("memory (alloc/copy)", ["malloc", "_int_free", "free", "memmove", "memcpy", "memset", "operator new",
                             "arena", "calloc"]),
    ("math library", ["ieee754", "log_fma", "exp_fma", "pow", "sincos", "atan", "sqrt"]),
]


def parse(path):
    rows = []
    for line in open(path, errors="replace"):
        m = re.match(r"\s+([\d.]+)%\s+\S+\s+(\S+)\s+\[[.k]\]\s+(.+)$", line)
        if m:
            rows.append((float(m.group(1)), m.group(3).strip().lower()))
    return rows


def group_shares(rows):
    shares = {name: 0.0 for name, _ in GROUPS}
    shares["other"] = 0.0
    for share, symbol in rows:
        for name, keys in GROUPS:
            if any(k in symbol for k in keys):
                shares[name] += share
                break
        else:
            shares["other"] += share
    return shares


def main():
    names = sys.argv[1:]
    all_shares = {}
    for name in names:
        rows = parse(f"{PERF}/out/{name}.perf.flat.txt")
        all_shares[name] = group_shares(rows)
    keys = [g for g, _ in GROUPS] + ["other"]
    print(f"{'group':36s} " + " ".join(f"{n[-12:]:>13s}" for n in names))
    for k in keys:
        print(f"{k:36s} " + " ".join(f"{all_shares[n][k]:12.1f}%" for n in names))
    import json
    json.dump(all_shares, open(f"{PERF}/perfgroups.json", "w"), indent=1)


if __name__ == "__main__":
    main()
