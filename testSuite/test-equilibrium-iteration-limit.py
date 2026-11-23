#!/usr/bin/env python3
"""Test the equilibrium iteration limit using the physical-reference executable.

Code written by OpenAI Codex at Andrew Robertson's request.
"""

import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
TARGET = "tests.galactic_structure_equilibrium.exe"


def main():
    """Run successful, exhausted-limit and invalid-parameter cases separately."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path,
                        help="Use an already-built executable (including a suffixed build).")
    parser.add_argument("--output", type=Path,
                        default=ROOT / "testSuite/outputs/equilibriumIterationLimit")
    args = parser.parse_args()
    if args.executable is None:
        subprocess.run(["make", TARGET], cwd=ROOT, check=True)
        executable = ROOT / TARGET
    else:
        executable = args.executable.resolve()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    results = []
    reference = None
    cases = [("omitted", None, None), ("explicit100", 100, None),
             ("explicit1000", 1000, None),
             ("exhausted", 2, "failed to find converged solution"),
             ("invalid1", 1, "[countIterationsMaximum] must be at least 2"),
             ("invalid0", 0, "[countIterationsMaximum] must be at least 2"),
             ("invalidNegative", -1, "[countIterationsMaximum] must be at least 2")]
    environment = dict(os.environ, OMP_NUM_THREADS="1")
    for name, limit, error in cases:
        tree = ET.parse(ROOT / "testSuite/parameters/galacticStructureEquilibrium.xml")
        solver = tree.getroot().find("galacticStructureSolver")
        if limit is not None:
            ET.SubElement(solver, "countIterationsMaximum", value=str(limit))
        parameters = output / f"{name}.xml"
        tree.write(parameters, encoding="utf-8", xml_declaration=True)
        log = output / f"{name}.log"
        with log.open("w") as handle:
            result = subprocess.run([str(executable), str(parameters)], cwd=ROOT,
                                    env=environment, stdout=handle, stderr=subprocess.STDOUT)
        text = log.read_text(errors="replace")
        if error is None:
            # The program reports the two physical-reference comparisons for each of 14 models.
            comparisons = re.findall(r"model\s+\d+: fractional difference[^\n]*", text)
            passed = result.returncode == 0 and "FAILED" not in text and len(comparisons) == 14
            if name == "omitted":
                reference = comparisons
            else:
                # Raising the ceiling must not change the result of solves which already converge.
                passed = passed and comparisons == reference
        else:
            passed = result.returncode != 0 and error in text
        results.append(dict(case=name, passed=passed, returncode=result.returncode, log=str(log)))
        print(f"{'SUCCESS' if passed else 'FAILED'}: {name}; log: {log}")
    (output / "summary.json").write_text(json.dumps(results, indent=2) + "\n")
    # The regression harness determines failure from the FAILED markers above.


if __name__ == "__main__":
    main()
