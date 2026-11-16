#!/usr/bin/env python3
"""Run the equilibrium iteration regression and coupled disk/spheroid reference tests.

Build tests.galactic_structure_equilibrium_iteration.exe and
tests.galactic_structure_equilibrium.exe first (optionally with SUFFIX).
The numerical test calls the same Fortran update used by the production solver.
The reference test exercises the full solver in contracting NFW haloes.

Contributions: Andrew Robertson, Codex.
"""
import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
TESTS = ("galactic_structure_equilibrium_iteration", "galactic_structure_equilibrium")


def main():
    """Run both compiled tests on one thread and reject assertions even when the executable exits zero."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--suffix", default="", help="Executable suffix passed to make as SUFFIX")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "testSuite/outputs/equilibriumOscillation")
    args = parser.parse_args()
    executables = [ROOT / ("tests." + name + ".exe" + args.suffix) for name in TESTS]
    for executable in executables:
        if not executable.is_file() or not os.access(executable, os.X_OK):
            parser.error(f"Build the test executable first: {executable}")
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    if any((output / (name + ".log")).exists() for name in TESTS):
        parser.error("Test logs already exist; choose a new --output-dir to preserve them")
    env = dict(os.environ, GALACTICUS_EXEC_PATH=str(ROOT), OMP_NUM_THREADS="1")
    records = []
    for name, executable in zip(TESTS, executables):
        log = output / (name + ".log")
        print(f"Running {name}; log: {log}", flush=True)
        with log.open("w") as stream:
            result = subprocess.run([str(executable)], cwd=ROOT, env=env, stdout=stream, stderr=subprocess.STDOUT)
        text = log.read_text()
        passed = re.findall(r"Tests passed:\s*(\d+)", text)
        failed = re.findall(r"Tests failed:\s*(\d+)", text)
        ok = (result.returncode == 0 and len(passed) == 1 and int(passed[0]) > 0
              and len(failed) == 1 and int(failed[0]) == 0 and "FAILED:" not in text)
        records.append(dict(test=name, executable=str(executable), threads=1, returncode=result.returncode,
                            passed=int(passed[0]) if len(passed) == 1 else None,
                            failed=int(failed[0]) if len(failed) == 1 else None, status="passed" if ok else "failed"))
        print(f"{'SUCCESS' if ok else 'FAILED'}: {name}", flush=True)
    (output / "summary.json").write_text(json.dumps(records, indent=2) + "\n")
    return 0 if all(row["status"] == "passed" for row in records) else 1


if __name__ == "__main__":
    sys.exit(main())
