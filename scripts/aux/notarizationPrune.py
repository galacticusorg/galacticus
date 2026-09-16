#!/usr/bin/env python3
"""Prune superseded macOS notarization bookkeeping assets from the `bleeding-edge` release.

Two small JSON assets track the notarization of each `master` deploy:

  * `notarization-requests-<sha>.json` -- written by `notarizationSubmit.py` and
    uploaded by the `Deploy` job; it records the `notarytool` submission ids.
  * `notarization-complete-<sha>.json` -- written by the `Notarize-MacOS` poller
    once the result has been reported to Slack, so that a later poll of the same
    commit does not report it a second time.

Both are keyed by commit, and `gh release upload --clobber` only ever replaces an
asset of the same name, so the pair accumulates on the rolling `bleeding-edge`
release: one per deploy, indefinitely.  They are tiny (tens of bytes each) but
they came to make up the bulk of that release's asset list, which is what a user
looking for a binary has to read past.

Only the pair for the *current* `master` HEAD is ever read -- `notarizationRetrieve.py`
resolves both names from `git rev-parse HEAD` -- so every older pair is history
with no reader.  This prunes them, subject to two guards which between them make
an unintended deletion of a live pair impossible:

  * the HEAD commit's pair is always kept, even if it is also the oldest; and
  * a pair is kept until it is `--min-age-days` old, whatever the retention count,
    so a burst of deploys cannot push a still-recent pair out of the window.

Assets are deleted in pairs, by commit.  Removing a `complete` marker while its
`requests` sibling remained would, for HEAD, invite a duplicate Slack
notification; keeping them together means the two are never out of step.
"""

import argparse
import json
import os
import re
import subprocess
import sys
from datetime import datetime, timedelta, timezone

RELEASE = "bleeding-edge"
DEFAULT_REPOSITORY = "galacticusorg/galacticus"

# `notarization-requests-<sha>.json` / `notarization-complete-<sha>.json`.  The
# commit is a full 40-character hash: both producers write `${GITHUB_SHA}` or
# `git rev-parse HEAD`, never an abbreviation.
_ASSET = re.compile(r"^notarization-(requests|complete)-([0-9a-f]{40})\.json$")

# How many commits' pairs to keep, newest first.  Deploys run on every push to
# `master`, so this is a few weeks of history at the usual rate -- enough to read
# back what happened to a recent release without the list growing without bound.
DEFAULT_KEEP = 10

# A pair younger than this is kept regardless of the retention count.
DEFAULT_MIN_AGE_DAYS = 7


def _run(command):
    """Run `command`, returning stdout; raises `CalledProcessError` on failure."""
    return subprocess.run(command, check=True, capture_output=True,
                          text=True).stdout


def _api(path):
    """Return the parsed JSON body of a `gh api` GET of `path`."""
    return json.loads(_run(["gh", "api", "-H", "Accept: application/vnd.github+json",
                            path]))


def repository():
    """The `owner/name` this is running against."""
    return os.environ.get("GITHUB_REPOSITORY") or DEFAULT_REPOSITORY


def head_commit():
    """The commit the checked-out tree is at, or None outside a repository.

    This matches how `notarizationRetrieve.py` names the assets it reads, which
    is what makes "the HEAD pair" the right thing to protect.
    """
    try:
        return _run(["git", "rev-parse", "HEAD"]).strip() or None
    except (subprocess.CalledProcessError, OSError):
        return None


def list_assets(repo):
    """Every asset on the release, as `{id, name, updated_at}` dictionaries.

    The assets are paged rather than read from the release object, which carries
    only the first page of them -- and this release has long since outgrown one
    page, precisely because of the assets being pruned here.
    """
    release = _api(f"repos/{repo}/releases/tags/{RELEASE}")
    assets, page = [], 1
    while True:
        batch = _api(f"repos/{repo}/releases/{release['id']}/assets"
                     f"?per_page=100&page={page}")
        assets.extend(batch)
        if len(batch) < 100:
            return assets
        page += 1


def group_by_commit(assets):
    """Group the notarization assets among `assets` by the commit they belong to.

    Returns `{sha: {"assets": [...], "updated": datetime}}`, where `updated` is
    the most recent of the pair -- the `complete` marker is written after the
    `requests` file, so the pair's age is that of its newest member.
    """
    groups = {}
    for asset in assets:
        matched = _ASSET.match(asset["name"])
        if not matched:
            continue
        sha = matched.group(2)
        updated = datetime.strptime(asset["updated_at"],
                                    "%Y-%m-%dT%H:%M:%SZ").replace(tzinfo=timezone.utc)
        group = groups.setdefault(sha, {"assets": [], "updated": updated})
        group["assets"].append(asset)
        group["updated"] = max(group["updated"], updated)
    return groups


def select_prunable(groups, head, keep, min_age_days, now=None):
    """The commits in `groups` whose assets may be deleted, oldest first.

    Kept: `head`, the `keep` most recently updated commits, and anything updated
    within `min_age_days`.  Everything else is returned.
    """
    now = now or datetime.now(timezone.utc)
    cutoff = now - timedelta(days=min_age_days)
    ordered = sorted(groups, key=lambda sha: groups[sha]["updated"], reverse=True)
    protected = set(ordered[:max(keep, 0)])
    if head is not None:
        protected.add(head)
    prunable = [sha for sha in ordered
                if sha not in protected and groups[sha]["updated"] < cutoff]
    prunable.reverse()               # oldest first, so the log reads chronologically
    return prunable


def delete_asset(repo, asset):
    """Delete one release asset by id."""
    _run(["gh", "api", "-X", "DELETE",
          f"repos/{repo}/releases/assets/{asset['id']}"])


def summarize(lines):
    """Append `lines` to the job summary, when running under Actions."""
    path = os.environ.get("GITHUB_STEP_SUMMARY")
    if not path:
        return
    with open(path, "a", encoding="utf-8") as summary:
        summary.write("## Notarization asset retention\n")
        for line in lines:
            summary.write(f"- {line}\n")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--keep", type=int, default=DEFAULT_KEEP,
                        help=f"commits' pairs to retain, newest first "
                             f"(default {DEFAULT_KEEP})")
    parser.add_argument("--min-age-days", type=float, default=DEFAULT_MIN_AGE_DAYS,
                        help=f"never delete a pair younger than this, whatever "
                             f"--keep says (default {DEFAULT_MIN_AGE_DAYS})")
    parser.add_argument("--dry-run", action="store_true",
                        help="report what would be deleted; delete nothing")
    arguments = parser.parse_args(argv)

    repo = repository()
    head = head_commit()
    groups = group_by_commit(list_assets(repo))
    prunable = select_prunable(groups, head, arguments.keep, arguments.min_age_days)

    if not prunable:
        message = (f"{len(groups)} notarization commit(s) on {RELEASE}; "
                   f"none older than the retention window.")
        print(message)
        summarize([message])
        return 0

    deleted, failed = 0, 0
    for sha in prunable:
        for asset in sorted(groups[sha]["assets"], key=lambda item: item["name"]):
            if arguments.dry_run:
                print(f"would delete {asset['name']}")
                deleted += 1
                continue
            try:
                delete_asset(repo, asset)
                print(f"deleted {asset['name']}")
                deleted += 1
            except subprocess.CalledProcessError as error:
                # Housekeeping must not fail the notarization report it runs
                # alongside, so a failure is annotated and the sweep continues.
                print(f"::warning::could not delete {asset['name']}: "
                      f"{error.stderr.strip()}")
                failed += 1

    verb = "would delete" if arguments.dry_run else "deleted"
    lines = [f"{verb} {deleted} asset(s) across {len(prunable)} superseded "
             f"commit(s); kept {len(groups) - len(prunable)}."]
    if failed:
        lines.append(f"{failed} asset(s) could not be deleted -- see the log.")
    print(lines[0])
    summarize(lines)
    return 0


if __name__ == "__main__":
    sys.exit(main())
