"""Tests for the retention selection in scripts/aux/notarizationPrune.py.

The prune runs unattended against the public `bleeding-edge` release, and
deleting a release asset cannot be undone, so what is exercised here is which
commits it is willing to delete. Two live pairs must survive every policy:

  * the `master` HEAD commit's pair, which `notarizationRetrieve.py` reads on
    every poll -- deleting its `requests` file would silently reduce the poll to
    "no-request" and the notarization result would never be reported; and
  * any recent pair, whatever the retention count, so that a burst of deploys
    (the count is per commit, not per day) cannot push a pair that a poll might
    still want out of the window.

Asset grouping is covered too, because the pair is keyed by a 40-character
commit hash embedded in the asset name: a pattern that also matched the other
`notarization-*` names, or that took the pair's age from the older member, would
retire a pair earlier than the policy says.
"""

import importlib.util
import os
from datetime import datetime, timedelta, timezone

import pytest

_SCRIPT = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), os.pardir, 'notarizationPrune.py')

_NOW = datetime(2026, 9, 13, 12, 0, 0, tzinfo=timezone.utc)


@pytest.fixture(scope='module')
def notarizationPrune():
    spec = importlib.util.spec_from_file_location('notarizationPrune', _SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _sha(seed):
    """A distinct 40-character hexadecimal commit hash."""
    return f"{seed:040x}"


def _asset(name, ageDays, identifier=1):
    stamp = (_NOW - timedelta(days=ageDays)).strftime("%Y-%m-%dT%H:%M:%SZ")
    return {"id": identifier, "name": name, "updated_at": stamp}


def _pair(seed, ageDays):
    """Both bookkeeping assets for one commit, `ageDays` old."""
    sha = _sha(seed)
    return sha, [_asset(f"notarization-requests-{sha}.json", ageDays, seed * 2),
                 _asset(f"notarization-complete-{sha}.json", ageDays, seed * 2 + 1)]


def _groups(module, pairs):
    assets = [asset for _, group in pairs for asset in group]
    return module.group_by_commit(assets)


# --- grouping --------------------------------------------------------------

def test_only_the_notarization_pair_is_grouped(notarizationPrune):
    """Every other asset on the release must be invisible here: the sweep deletes
    whatever it groups, so a loose pattern would delete binaries."""
    assets = [
        _asset(f"notarization-requests-{_sha(1)}.json", 30),
        _asset(f"notarization-complete-{_sha(1)}.json", 30),
        _asset("Galacticus.exe", 0),
        _asset("toolsMacOSM1.tar.zst", 0),
        _asset("SHA256SUMS", 0),
        _asset("notarization-poll-result.json", 0),        # never uploaded, but adjacent
        _asset("notarization-requests-deadbeef.json", 30),  # abbreviated: not a pair
    ]
    groups = notarizationPrune.group_by_commit(assets)
    assert list(groups) == [_sha(1)]
    assert len(groups[_sha(1)]["assets"]) == 2


def test_a_pair_is_as_old_as_its_newest_member(notarizationPrune):
    """The `complete` marker is written when the result is reported, after the
    `requests` file; taking the older member's date would retire the pair early."""
    sha = _sha(1)
    groups = notarizationPrune.group_by_commit([
        _asset(f"notarization-requests-{sha}.json", 30),
        _asset(f"notarization-complete-{sha}.json", 2),
    ])
    assert groups[sha]["updated"] == _NOW - timedelta(days=2)


# --- retention -------------------------------------------------------------

def test_superseded_commits_are_pruned_oldest_first(notarizationPrune):
    pairs = [_pair(seed, ageDays=30 + seed) for seed in range(1, 6)]
    groups = _groups(notarizationPrune, pairs)
    prunable = notarizationPrune.select_prunable(
        groups, head=None, keep=2, min_age_days=7, now=_NOW)
    # Newest two kept; the rest reported oldest first.
    assert prunable == [_sha(5), _sha(4), _sha(3)]


def test_the_head_commit_is_never_pruned(notarizationPrune):
    """`notarizationRetrieve.py` resolves the assets it reads from HEAD, so losing
    this pair turns the next poll into a silent "no-request"."""
    pairs = [_pair(seed, ageDays=100 + seed) for seed in range(1, 4)]
    groups = _groups(notarizationPrune, pairs)
    head = _sha(3)                                  # also the oldest
    prunable = notarizationPrune.select_prunable(
        groups, head=head, keep=0, min_age_days=0, now=_NOW)
    assert head not in prunable
    assert prunable == [_sha(2), _sha(1)]


def test_a_recent_pair_survives_the_retention_count(notarizationPrune):
    """The count is per commit, so a burst of deploys can put many pairs ahead of
    one that is only hours old; the age floor is what keeps it."""
    pairs = [_pair(seed, ageDays=0.1) for seed in range(1, 21)]
    groups = _groups(notarizationPrune, pairs)
    assert notarizationPrune.select_prunable(
        groups, head=None, keep=2, min_age_days=7, now=_NOW) == []


def test_nothing_is_pruned_when_everything_is_within_the_window(notarizationPrune):
    pairs = [_pair(seed, ageDays=seed) for seed in range(1, 4)]
    groups = _groups(notarizationPrune, pairs)
    assert notarizationPrune.select_prunable(
        groups, head=None, keep=10, min_age_days=7, now=_NOW) == []


def test_a_commit_with_only_one_of_the_pair_is_still_pruned(notarizationPrune):
    """A deploy whose notarization never completed leaves a `requests` file with no
    `complete` sibling; it is history like any other once it ages out."""
    sha = _sha(1)
    groups = notarizationPrune.group_by_commit(
        [_asset(f"notarization-requests-{sha}.json", 60)])
    assert notarizationPrune.select_prunable(
        groups, head=None, keep=0, min_age_days=7, now=_NOW) == [sha]


def test_no_assets_is_not_an_error(notarizationPrune):
    assert notarizationPrune.select_prunable(
        {}, head=_sha(1), keep=10, min_age_days=7, now=_NOW) == []


def test_a_negative_keep_does_not_protect_the_whole_list(notarizationPrune):
    """`max(keep, 0)` guards the slice: `ordered[:-1]` would otherwise protect
    everything but the oldest, quietly inverting the policy."""
    pairs = [_pair(seed, ageDays=30 + seed) for seed in range(1, 4)]
    groups = _groups(notarizationPrune, pairs)
    assert notarizationPrune.select_prunable(
        groups, head=None, keep=-1, min_age_days=7, now=_NOW) == [
            _sha(3), _sha(2), _sha(1)]
