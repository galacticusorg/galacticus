#!/usr/bin/env python3
# Audit the git revisions recorded by parameter files and parameter migrations.
# Andrew Benson with assistance from Claude (09-October-2026).

"""Audit the git revisions recorded by parameter files and parameter migrations.

Parameter migration is driven entirely by commit hashes:

* each parameter file records, in ``<lastModified revision="..."/>``, the
  revision to which it was last migrated;
* ``scripts/aux/migrations.xml`` registers each migration against the commit
  which introduced it.

``parametersMigrate.py`` applies every migration whose commit lies on the
ancestry path from a file's ``lastModified`` revision to ``HEAD``, and
Galacticus itself (when built with libgit2) performs the same ancestry test to
decide whether a parameter file is out of date. A hash which is not in the
history - typically because the branch on which it was recorded was rebased or
amended before it was merged - therefore breaks both: a missing ``lastModified``
revision makes migration fail outright, and a missing migration commit makes
that migration silently never fire.

This tool reports:

* a ``lastModified`` revision, in any tracked ``*.xml`` file, which does not
  resolve to a commit, or which is not an ancestor of ``HEAD``;
* a migration commit which is not a full, lower-case, 40 character hash (the
  migration is matched against ``git rev-list`` output by exact string
  comparison, so an abbreviated hash never matches), which does not resolve to a
  commit, or which is not an ancestor of ``HEAD``.

A full history is required: in a shallow clone every commit beyond the
shallow boundary would be reported as missing, so the tool refuses to run in
one (use ``fetch-depth: 0`` with ``actions/checkout``).

Usage::

   ./scripts/aux/auditRevisions.py           # report
   ./scripts/aux/auditRevisions.py --check   # exit 1 if anything is reported
"""

import argparse
import os
import re
import subprocess
import sys
import xml.etree.ElementTree as ElementTree

# A `lastModified` element may have its attributes split across lines, so match across newlines.
lastModifiedPattern = re.compile(r'<lastModified\b[^>]*?\brevision\s*=\s*"([^"]*)"', re.DOTALL)
fullHashPattern     = re.compile(r'^[0-9a-f]{40}$')
migrationsPath      = os.path.join("scripts", "aux", "migrations.xml")

def git(root, *arguments):
    """Run a git command in `root`, returning the completed process."""
    return subprocess.run(["git", "-C", root, *arguments], capture_output=True, text=True)

def resolves(root, revision):
    """Return true if `revision` resolves to a commit."""
    return git(root, "rev-parse", "--verify", "--quiet", f"{revision}^{{commit}}").returncode == 0

def isAncestor(root, revision):
    """Return true if `revision` is an ancestor of (or is) `HEAD`."""
    status = git(root, "merge-base", "--is-ancestor", revision, "HEAD").returncode
    if status not in (0, 1):
        raise RuntimeError(f"failed to determine whether '{revision}' is an ancestor of HEAD")
    return status == 0

def lastModifiedRevisions(root):
    """Return a dictionary mapping each `lastModified` revision to the tracked files which record it."""
    listing = git(root, "ls-files", "-z", "--", "*.xml")
    if listing.returncode != 0:
        raise RuntimeError(f"failed to list tracked files: {listing.stderr.strip()}")
    revisions = {}
    for fileName in filter(None, listing.stdout.split("\0")):
        try:
            with open(os.path.join(root, fileName), encoding="utf-8", errors="replace") as file:
                content = file.read()
        except FileNotFoundError:
            # Tracked but deleted in the working tree.
            continue
        for revision in set(lastModifiedPattern.findall(content)):
            revisions.setdefault(revision, []).append(fileName)
    return revisions

def migrationCommits(root):
    """Return the commit hashes registered in `migrations.xml`, in order."""
    tree = ElementTree.parse(os.path.join(root, migrationsPath))
    return [migration.get("commit") or "" for migration in tree.getroot().findall("migration")]

def audit(root):
    """Return a list of issues found in the repository at `root`."""
    issues = []
    # Last-modified revisions recorded by parameter files.
    for revision, fileNames in sorted(lastModifiedRevisions(root).items()):
        if not resolves(root, revision):
            problem = "does not resolve to a commit"
        elif not isAncestor(root, revision):
            problem = "is not an ancestor of HEAD"
        else:
            continue
        shown = ", ".join(fileNames[:3]) + (f", and {len(fileNames)-3} more" if len(fileNames) > 3 else "")
        issues.append(f"lastModified revision \"{revision}\" {problem}: recorded in {len(fileNames)} file(s): {shown}")
    # Commits registered for migrations.
    for commit in migrationCommits(root):
        if not fullHashPattern.match(commit):
            problem = "is not a full, lower-case, 40 character hash"
        elif not resolves(root, commit):
            problem = "does not resolve to a commit"
        elif not isAncestor(root, commit):
            problem = "is not an ancestor of HEAD"
        else:
            continue
        issues.append(f"{migrationsPath}: migration commit \"{commit}\" {problem}")
    return issues

def main():
    parser = argparse.ArgumentParser(description="Audit the git revisions recorded by parameter files and parameter migrations")
    parser.add_argument("root"   , nargs="?", default=None, help="repository root (default: $GALACTICUS_EXEC_PATH, else the current directory)")
    parser.add_argument("--check", action="store_true"    , help="exit with status 1 if any issue is reported")
    arguments = parser.parse_args()
    root = arguments.root or os.environ.get("GALACTICUS_EXEC_PATH", ".")
    shallow = git(root, "rev-parse", "--is-shallow-repository")
    if shallow.returncode != 0:
        print(f"error: '{root}' is not a git repository: {shallow.stderr.strip()}", file=sys.stderr)
        return 2
    if shallow.stdout.strip() == "true":
        print(f"error: '{root}' is a shallow clone; a full history is needed to audit revisions (try `git fetch --unshallow`)", file=sys.stderr)
        return 2
    issues = audit(root)
    for issue in issues:
        print(issue)
    if issues:
        print(f"\n{len(issues)} issue(s) found. A revision missing from the history usually means that the branch on which it was"
              " recorded was rebased or amended; replace it with the equivalent commit which was actually merged.")
    else:
        print("All lastModified revisions and migration commits are ancestors of HEAD.")
    return 1 if issues and arguments.check else 0

if __name__ == "__main__":
    sys.exit(main())
