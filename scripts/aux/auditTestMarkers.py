#!/usr/bin/env python3
# Audit testSuite scripts against the outcome-marker convention.
# Andrew Benson with assistance from Claude (09-September-2026).

"""Audit ``testSuite`` scripts against the outcome-marker convention.

A test script signals its outcome by printing a marker line and then exiting
with status zero:

* ``FAILED: <what went wrong>`` - the test failed;
* ``SUCCESS: <what passed>``    - the test passed;
* ``SKIPPED: <why>``            - a prerequisite could not be met, so the test
  ran nothing at all.

The harnesses judge a test by scanning its log for ``FAIL`` (as a substring,
case sensitively), so a marker misspelled ``FAIL:`` or ``FAILURE`` is still
detected - but a *passing* run whose output happens to contain ``FAIL``
anywhere is reported as a failure. Both harnesses also fail a test whose script
exits non-zero, so an uncaught exception is detected, but as a traceback rather
than as a readable message.

This tool therefore reports:

* a script with no ``FAILED`` marker at all, which can never report a failure;
* markers spelled ``FAIL:``, ``FAILURE``, ``PASS:``, ``SKIP:`` or lower-case
  ``success:`` rather than the canonical spellings;
* ``FAILED`` appearing anywhere but the start of an output line, which risks
  turning a passing run red;
* a script with no ``SUCCESS`` marker;
* a ``raise`` or bare ``assert`` which would abort the script with a traceback
  instead of a marker;
* an explicit non-zero exit.

A string literal on a line ending in the pragma ``# markers: exempt`` is
ignored - use it for a literal which deliberately contains a marker word
without being one, such as a ``grep`` pattern.

Usage::

   ./scripts/aux/auditTestMarkers.py                 # report on testSuite/
   ./scripts/aux/auditTestMarkers.py --check         # exit 1 if anything is reported
   ./scripts/aux/auditTestMarkers.py --json out.json # also write a machine-readable report
"""

import argparse
import ast
import glob
import json
import os
import re
import sys

# Exceptions which are a legitimate way to leave a script: `SystemExit` is an exit, and argparse converts a
# `ArgumentTypeError` raised by a `type=` callable into its own usage message and exit.
raisesPermitted = ("SystemExit", "ArgumentTypeError")

# Scripts which are not tests, and so have no outcome to mark.
scriptsExempt = ("test-all.py",)


def docstringsOf(tree):
    """Return the set of string-literal nodes which are docstrings, and so are never printed."""
    docstrings = set()
    for node in ast.walk(tree):
        if not isinstance(node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            continue
        body = getattr(node, "body", None)
        if body and isinstance(body[0], ast.Expr) and isinstance(body[0].value, ast.Constant) \
           and isinstance(body[0].value.value, str):
            docstrings.add(id(body[0].value))
    return docstrings


def literalsOf(tree, source, lines):
    """Yield every string literal in the file as (lineNumber, value), less docstrings and those carrying the
    exemption pragma."""
    docstrings = docstringsOf(tree)
    for node in ast.walk(tree):
        if isinstance(node, ast.Constant) and isinstance(node.value, str):
            if id(node) in docstrings:
                continue
            lineNumber = node.lineno
            if lineNumber <= len(lines) and lines[lineNumber-1].rstrip().endswith("# markers: exempt"):
                continue
            yield lineNumber, node.value


def callName(call):
    """Return the dotted name of a called function, as far as it can be resolved statically."""
    function = call.func
    if isinstance(function, ast.Name):
        return function.id
    parts = []
    while isinstance(function, ast.Attribute):
        parts.append(function.attr)
        function = function.value
    if isinstance(function, ast.Name):
        parts.append(function.id)
    return ".".join(reversed(parts))


def auditScript(path):
    """Audit a single script, returning a list of issue strings."""
    source = open(path, errors="replace").read()
    lines  = source.splitlines()
    issues = []
    try:
        tree = ast.parse(source)
    except SyntaxError as error:
        return [f"syntax error: {error}"]

    literals = list(literalsOf(tree, source, lines))

    # Marker spellings. `FAILED` and `SUCCESS` may be assembled into a variable and interpolated later, so every literal
    # is considered, not only those appearing directly in a `print`.
    hasFailed  = any("FAILED"  in text for _, text in literals)
    hasSuccess = any("SUCCESS" in text for _, text in literals)
    if not hasFailed:
        issues.append("no FAILED marker anywhere: the script can never report a failure")
    if not hasSuccess:
        issues.append("no SUCCESS marker anywhere")

    for lineNumber, text in literals:
        stripped = text.strip()
        # A marker word which is not the canonical spelling.
        if "FAILURE" in text:
            issues.append(f"line {lineNumber}: 'FAILURE' rather than 'FAILED': {stripped[:60]!r}")
        elif re.search(r'\bFAIL(?!ED)', text):
            issues.append(f"line {lineNumber}: 'FAIL' rather than 'FAILED': {stripped[:60]!r}")
        elif "FAILED" in text and not stripped.startswith("FAILED"):
            # `FAILED` anywhere but at the start of the line will appear in the output of a passing run if the line is
            # printed unconditionally - and the harnesses match `FAIL` as a substring anywhere in the log.
            issues.append(f"line {lineNumber}: 'FAILED' is not at the start of the line: {stripped[:60]!r}")
        if stripped.startswith("PASS:"):
            issues.append(f"line {lineNumber}: 'PASS:' rather than 'SUCCESS:': {stripped[:60]!r}")
        if stripped.startswith("success:"):
            issues.append(f"line {lineNumber}: lower-case 'success:' rather than 'SUCCESS:': {stripped[:60]!r}")
        if re.match(r'SKIP(?!PED)\b', stripped):
            issues.append(f"line {lineNumber}: 'SKIP' rather than 'SKIPPED': {stripped[:60]!r}")

    # Exits and aborts. The convention is to exit zero always, so that the outcome is read from the marker; anything which
    # leaves the script by another route reports itself as a traceback rather than as a readable message.
    for node in ast.walk(tree):
        if isinstance(node, ast.Call) and callName(node) in ("sys.exit", "exit", "os._exit", "quit"):
            if node.args:
                code = node.args[0]
                if isinstance(code, ast.Constant) and code.value not in (0, None, False):
                    issues.append(f"line {node.lineno}: exits with non-zero status {code.value!r}")
        elif isinstance(node, ast.Raise):
            exception = node.exc
            name = ""
            if isinstance(exception, ast.Call):
                name = callName(exception).split(".")[-1]
            elif isinstance(exception, ast.Name):
                name = exception.id
            if name and name not in raisesPermitted:
                issues.append(f"line {node.lineno}: raises {name} - aborts with a traceback rather than a FAILED marker")
        elif isinstance(node, ast.Assert):
            issues.append(f"line {node.lineno}: bare assert - aborts with a traceback rather than a FAILED marker")

    # Early-exit guards. A guard which explains itself but prints no marker leaves the harnesses nothing to see: the script
    # exits zero, its log contains no `FAIL`, and the test is recorded as a pass while having tested nothing.
    for node in ast.walk(tree):
        if not isinstance(node, ast.If):
            continue
        for block in (node.body, node.orelse):
            if not block:
                continue
            last = block[-1]
            exits = isinstance(last, ast.Return)
            if isinstance(last, ast.Expr) and isinstance(last.value, ast.Call) \
               and callName(last.value) in ("sys.exit", "exit", "os._exit", "quit"):
                exits = True
            if not exits:
                continue
            printed = []
            for statement in block:
                for inner in ast.walk(statement):
                    if isinstance(inner, ast.Call) and callName(inner) == "print":
                        printed.extend(text for _, text in literalsOf(inner, source, lines))
            if printed and not any(marker in text for text in printed for marker in ("FAILED", "SUCCESS", "SKIPPED")):
                issues.append(f"line {last.lineno}: early exit which prints but sets no marker: {printed[0].strip()[:60]!r}")

    return issues


def main():
    parser = argparse.ArgumentParser(description="Audit testSuite scripts against the outcome-marker convention")
    parser.add_argument("root" , nargs="?", default=None, help="directory holding the test scripts (default: $GALACTICUS_EXEC_PATH/testSuite)")
    parser.add_argument("--check", action="store_true"  , help="exit with status 1 if any issue is reported")
    parser.add_argument("--json" , type=str, default=None, help="also write a machine-readable report to this file")
    arguments = parser.parse_args()

    root = arguments.root
    if root is None:
        root = os.path.join(os.environ.get("GALACTICUS_EXEC_PATH", "."), "testSuite")
    scripts = sorted(glob.glob(os.path.join(root, "test-*.py")) + glob.glob(os.path.join(root, "validate-*.py")))
    scripts = [script for script in scripts if os.path.basename(script) not in scriptsExempt]
    if not scripts:
        print(f"auditTestMarkers.py: no test scripts found under '{root}'")
        return 1

    report      = []
    countIssues = 0
    for script in scripts:
        issues = auditScript(script)
        countIssues += len(issues)
        report.append({"script": os.path.basename(script), "issues": issues})
        if issues:
            print(f"{os.path.relpath(script)}:")
            for issue in issues:
                print(f"   {issue}")

    countScripts = sum(1 for record in report if record["issues"])
    print(f"{len(scripts)} scripts audited; {len(scripts)-countScripts} clean; {countScripts} with issues ({countIssues} in total)")

    if arguments.json:
        json.dump(report, open(arguments.json, "w"), indent=1)

    return 1 if (arguments.check and countIssues > 0) else 0


if __name__ == "__main__":
    sys.exit(main())
