#!/usr/bin/env python3
"""Regression test for the object-build recursion guards (issue #397).

A class that composites a member of its own class -- directly (e.g. the
`galacticFilterNot` filter, which wraps another `galacticFilter`) or via
another, mutually-compositing class (e.g. `criticalOverdensityEnvironmental`
<-> `haloEnvironmentNormal`) -- could, if no such object was provided
explicitly, search up the parameter tree, re-discover the object currently
being built, and attempt to build it again, leading to an unbounded recursion
that ultimately crashed the executable with a stack-overflow (SIGSEGV) and no
useful diagnostic.

The guards turn each of these into a clean, informative `Error_Report` abort.
This test runs two parameter files that each trigger one of the recursion
modes and checks that the model:

  * exits with a non-zero status (it must NOT silently "succeed"),
  * does NOT die from a segmentation fault (return code 139 / signal 11 ==
    the unguarded stack overflow we are protecting against), and
  * emits the expected recursion diagnostic.

Andrew Benson
"""

import subprocess
import sys

# Ensure output directory exists.
subprocess.run("mkdir -p outputs", shell=True)

# Each case: (label, parameter file, a string expected in the diagnostic output).
cases = [
    (
        "self-composition (galacticFilter 'not')",
        "testSuite/parameters/recursionGuardSelfComposition.xml",
        "composites a member of its own class",
    ),
    (
        "mutual recursion (criticalOverdensity <-> haloEnvironment)",
        "testSuite/parameters/recursionGuardMutual.xml",
        "recursive build of [criticalOverdensity] detected",
    ),
]

failed = False
for label, parameterFile, expected in cases:
    process = subprocess.run(
        f"cd ..; ./Galacticus.exe {parameterFile}",
        shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        universal_newlines=True,
    )
    output     = process.stdout or ""
    returnCode = process.returncode

    # A segmentation fault (the unguarded behaviour) reports 139 via the shell
    # (128 + SIGSEGV), or a negative signal number if not run through a shell.
    if returnCode == 139 or returnCode < 0:
        print(f"FAILED: {label}: model crashed (return code {returnCode}) "
              f"- recursion guard did not fire")
        failed = True
        continue
    # The recursive configuration is invalid, so the model must abort.
    if returnCode == 0:
        print(f"FAILED: {label}: model ran to completion but the "
              f"configuration is recursive and should have been rejected")
        failed = True
        continue
    # It aborted - confirm it did so with the expected recursion diagnostic.
    if expected not in output:
        print(f"FAILED: {label}: model aborted (return code {returnCode}) "
              f"but without the expected recursion diagnostic "
              f"['{expected}']")
        failed = True
        continue
    print(f"SUCCESS: {label}: recursion detected and reported cleanly "
          f"(return code {returnCode})")

if failed:
    print("FAILED: object-builder recursion guard test")
    sys.exit(1)

print("SUCCESS: object-builder recursion guard test")
