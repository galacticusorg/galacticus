#!/usr/bin/env python3
"""Differential test for the Python parameter-file resolver.

Compares `Galacticus.Parameters.resolve` against the authoritative Galacticus
oracle (`Galacticus.exe --output-processed-parameters`), which serializes the
fully-processed DOM. The oracle faithfully reflects the two transformations the
resolver implements at this stage -- XInclude expansion and change-file
application -- because Galacticus bakes both into the DOM before any value is
read (see the plan in ~/.claude/plans/galacticus-parameter-resolver.md).

Run from the `testSuite/` directory (as `test-all.py` does). Prints `FAILED` on
any divergence; `SKIP` (not a failure) if no built `Galacticus.exe` is found.

Andrew Benson (2026).
"""

import glob
import os
import subprocess
import sys
import tempfile
import time

REPO = os.path.abspath(os.environ.get("GALACTICUS_EXEC_PATH", ".."))
sys.path.insert(0, os.path.join(REPO, "python"))
from lxml import etree                                        # noqa: E402
from Galacticus.Parameters import resolve                     # noqa: E402

BINARY = os.path.join(REPO, "Galacticus.exe")
# Cap the (oracle-backed) XInclude corpus scan so the test stays quick; the cap
# is logged so partial coverage is never silent.
XINCLUDE_SAMPLE = int(os.environ.get("RESOLVER_XINCLUDE_SAMPLE", "40"))


def _norm(element):
    """Canonical form ignoring comments, whitespace, and attribute order."""
    if not isinstance(element.tag, str):
        return None
    children = tuple(c for c in (_norm(x) for x in element) if c is not None)
    return (element.tag, tuple(sorted(element.attrib.items())),
            (element.text or "").strip(), children)


def _first_diff(a, b, path="parameters"):
    if a is None or b is None:
        return f"{path}: node presence differs"
    if a[0] != b[0]:
        return f"{path}: tag {a[0]} vs {b[0]}"
    if a[1] != b[1]:
        return f"{path} <{a[0]}>: attribs {a[1]} vs {b[1]}"
    if a[2] != b[2]:
        return f"{path} <{a[0]}>: text {a[2]!r} vs {b[2]!r}"
    if len(a[3]) != len(b[3]):
        return (f"{path} <{a[0]}>: child count {len(a[3])} vs {len(b[3])} "
                f"({[c[0] for c in a[3]]} vs {[c[0] for c in b[3]]})")
    for ca, cb in zip(a[3], b[3]):
        deeper = _first_diff(ca, cb, f"{path}/{ca[0]}")
        if deeper:
            return deeper
    return None


def _run_oracle(param_file, change_files, out_path, timeout=40):
    """Run Galacticus.exe to serialize the processed tree; poll-and-kill once the
    output file is written (it is emitted before the heavy run begins)."""
    if os.path.exists(out_path):
        os.remove(out_path)
    cmd = [BINARY, param_file, *change_files,
           "--output-processed-parameters", out_path, "--dry-run"]
    proc = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    start, last, stable = time.time(), -1, 0
    try:
        while time.time() - start < timeout:
            if proc.poll() is not None:
                break
            if os.path.exists(out_path):
                size = os.path.getsize(out_path)
                stable = stable + 1 if size > 0 and size == last else 0
                last = size
                if stable >= 2:
                    break
            time.sleep(0.15)
    finally:
        if proc.poll() is None:
            proc.terminate()
            try:
                proc.wait(2)
            except subprocess.TimeoutExpired:
                proc.kill()
    return os.path.exists(out_path) and os.path.getsize(out_path) > 0


def _compare(param_file, change_files, label, out_path):
    # Compare only the oracle-faithful scope (XInclude + changes + reference
    # validation); the oracle keeps `active=` and does not prune, so conditional
    # evaluation is excluded here and covered by unit tests instead.
    if not _run_oracle(param_file, change_files, out_path):
        # Oracle errored before serialize; a faithful resolver must agree (error).
        try:
            resolve.resolve_file(param_file, change_files, conditionals=False)
        except Exception:
            return True, "both error (agreement)"
        return False, "oracle errored but resolver succeeded (too lenient)"
    try:
        oracle = etree.parse(out_path).getroot()
        mine = resolve.resolve_file(param_file, change_files,
                                    conditionals=False).getroot()
    except Exception as exc:                                  # pragma: no cover
        return False, f"comparison setup error: {exc}"
    diff = _first_diff(_norm(oracle), _norm(mine))
    return diff is None, diff or "ok"


# Synthetic base+change pairs exercising every change operation.
CHANGE_CASES = {
    "remove": ('<parameters><a value="1"/><b value="2"/></parameters>',
               '<changes><change type="remove" path="a"/></changes>'),
    "update": ('<parameters><a value="1"/></parameters>',
               '<changes><change type="update" path="a" value="9"/></changes>'),
    "update-append": ('<parameters><a value="1 2"/></parameters>',
               '<changes><change type="update" path="a" value=" 3" append="true"/></changes>'),
    "append": ('<parameters><a><x value="1"/></a></parameters>',
               '<changes><change type="append" path="a"><y value="2"/></change></changes>'),
    "insertBefore": ('<parameters><a value="1"/><b value="2"/></parameters>',
               '<changes><change type="insertBefore" path="b"><z value="0"/></change></changes>'),
    "insertAfter": ('<parameters><a value="1"/><b value="2"/></parameters>',
               '<changes><change type="insertAfter" path="a"><z value="0"/></change></changes>'),
    "replace": ('<parameters><a value="1"/></parameters>',
               '<changes><change type="replace" path="a"><b value="2"/><c value="3"/></change></changes>'),
    "replaceOrAppend-exists": ('<parameters><a value="1"/></parameters>',
               '<changes><change type="replaceOrAppend" path="a"><a value="2"/></change></changes>'),
    "replaceOrAppend-missing": ('<parameters><a><x value="1"/></a></parameters>',
               '<changes><change type="replaceOrAppend" path="a/y"><y value="2"/></change></changes>'),
    "replaceWith": ('<parameters><src value="k"><deep value="d"/></src><dst value="x"/></parameters>',
               '<changes><change type="replaceWith" path="dst" target="src"/></changes>'),
    "encapsulate": ('<parameters><t value="i"/></parameters>',
               '<changes><change type="encapsulate" path="t"><w value="o"/></change></changes>'),
    "predicate": ('<parameters><op value="a"/><op value="b"><inner value="1"/></op></parameters>',
               """<changes><change type="update" path="op[@value='b']/inner" value="9"/></changes>"""),
}


def main():
    if not os.path.isfile(BINARY):
        print(f"SKIP: no built Galacticus.exe at {BINARY}; "
              "differential resolver test requires the binary.")
        return 0

    failures = 0
    out_path = os.path.join(tempfile.gettempdir(), "resolver_oracle.xml")

    print("Change-operation parity vs oracle:")
    with tempfile.TemporaryDirectory() as work:
        for name, (base, changes) in CHANGE_CASES.items():
            base_file = os.path.join(work, "base.xml")
            change_file = os.path.join(work, "changes.xml")
            open(base_file, "w").write(base)
            open(change_file, "w").write(changes)
            ok, detail = _compare(base_file, [change_file], name, out_path)
            if ok:
                print(f"  ok: {name}")
            else:
                failures += 1
                print(f"  FAILED: {name} :: {detail}")

    print(f"XInclude parity vs oracle (first {XINCLUDE_SAMPLE} corpus files):")
    candidates = []
    for pattern in ("parameters/**/*.xml", "testSuite/parameters/**/*.xml"):
        for path in sorted(glob.glob(os.path.join(REPO, pattern), recursive=True)):
            try:
                root_tag = etree.parse(path).getroot().tag
            except Exception:
                continue
            if root_tag == "parameters" and b"xi:include" in open(path, "rb").read():
                candidates.append(path)
    scanned = candidates[:XINCLUDE_SAMPLE]
    print(f"  ({len(scanned)} of {len(candidates)} XInclude-using parameter files)")
    for path in scanned:
        ok, detail = _compare(path, [], os.path.relpath(path, REPO), out_path)
        if not ok:
            failures += 1
            print(f"  FAILED: {os.path.relpath(path, REPO)} :: {detail}")

    # Conditionals: the oracle keeps `active=` (cannot certify pruning), so verify
    # end-to-end instead -- fully resolve a real conditional file and confirm
    # Galacticus accepts the pruned result.
    tests_params = os.path.join(REPO, "testSuite", "parameters", "testsParameters.xml")
    if os.path.isfile(tests_params):
        print("Conditional end-to-end (testsParameters.xml):")
        resolved = os.path.join(tempfile.gettempdir(), "resolver_resolved.xml")
        try:
            root = resolve.resolve_file(tests_params, output=resolved).getroot()
        except Exception as exc:
            failures += 1
            print(f"  FAILED: resolve raised :: {exc}")
        else:
            active1 = [c.get("value") for c in root if c.tag == "active1"]
            if active1 != ["0.2"]:
                failures += 1
                print(f"  FAILED: conditional pruning -> active1 {active1}, expected ['0.2']")
            elif not _run_oracle(resolved, [], out_path):
                failures += 1
                print("  FAILED: Galacticus rejected the conditionally-pruned file")
            else:
                print("  ok: pruned to active1=0.2 and accepted by Galacticus")

    if failures:
        print(f"FAILED: {failures} resolver/oracle divergence(s)")
        return 1
    print("SUCCESS: resolver matches the Galacticus oracle "
          f"({len(CHANGE_CASES)} change ops, {len(scanned)} XInclude files)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
