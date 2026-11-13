"""Tests for the pull request validation report builder."""

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import buildPrValidationReport as report  # noqa: E402

PROBE = report.BASELINE_PROBE


def recent(low, high, median, scatter):
    return {"count": 30, "low": low, "high": high, "median": median, "scatter": scatter}


def results_file(sha, timestamp):
    body = {"commit": {"id": sha, "timestamp": timestamp}, "results": []}
    return "window.ANALYSES_DATA = " + json.dumps(body)


def git(repository, *arguments, env=None):
    run = subprocess.run(["git", "-C", str(repository), *arguments], capture_output=True, text=True, check=True, env=env)
    return run.stdout.strip()


@pytest.fixture(name="pages")
def fixture_pages(tmp_path):
    """A repository with results published for two master commits, and a master commit between them with no results."""
    repository = tmp_path / "pages"
    (repository / Path(PROBE).parent).mkdir(parents=True)
    git(repository, "init", "-q", "--template=")
    git(repository, "config", "user.email", "test@example.com")
    git(repository, "config", "user.name", "test")
    for sha, timestamp in (("a" * 40, "2026-01-01 00:00:00+00:00"), ("c" * 40, "2026-01-03 00:00:00+00:00")):
        (repository / PROBE).write_text(results_file(sha, timestamp), encoding="utf-8")
        git(repository, "add", ".")
        git(repository, "commit", "-q", "-m", "Update validation results")
    return repository


def commit_dated(repository, date):
    """Create a commit standing in for a master commit made at the given date, and return its hash."""
    environment = {**os.environ, "GIT_COMMITTER_DATE": date, "GIT_AUTHOR_DATE": date}
    git(repository, "commit", "-q", "--allow-empty", "-m", "master", env=environment)
    return git(repository, "rev-parse", "HEAD")


def test_sanitize_replaces_non_finite_numbers():
    assert report.sanitize({"a": [1.0, float("nan")], "b": float("inf")}) == {"a": [1.0, None], "b": None}


def test_roundoff_is_not_flagged_for_a_deterministic_model():
    outside, deviation = report.compare_to_history(-5817.0700001, recent(-5817.07, -5817.07, -5817.07, 0.0))
    assert not outside
    assert abs(deviation) < 1.0


def test_real_change_is_flagged_for_a_deterministic_model():
    outside, deviation = report.compare_to_history(-5800.0, recent(-5817.07, -5817.07, -5817.07, 0.0))
    assert outside
    assert deviation > 100.0


def test_value_within_recent_range_is_not_flagged():
    outside, deviation = report.compare_to_history(-4.0, recent(-6.0, -2.0, -4.5, 1.0))
    assert not outside
    assert deviation == pytest.approx(0.5)


def test_missing_values_are_not_flagged():
    assert report.compare_to_history(None, recent(-6.0, -2.0, -4.0, 1.0)) == (False, None)
    assert report.compare_to_history(-4.0, None) == (False, None)


def test_recent_statistics_ignores_impossible_likelihoods():
    assert report.recent_statistics([1.0, 2.0, 3.0, 1.8e302])["high"] == 3.0
    assert report.recent_statistics([1.8e302]) is None


def test_analysis_rows_convert_history_to_log_likelihood():
    result = {"attributes": {"name": "massFunction", "logLikelihood": "-5.0"}, "data": {}}
    history = {"Model - Likelihood - massFunction": [{"value": value} for value in (4.0, 5.0, 6.0)]}
    row = report.analysis_rows("suffix", "Title", {"results": [result]}, {"results": [result]}, history, {})[0]
    assert (row["recent"]["low"], row["recent"]["high"]) == (-6.0, -4.0)
    assert not row["outside"]


@pytest.mark.parametrize("analysis,recorded", [
    ("massCumulativeDistribution", " (tau) - Likelihood - mass function"),
    ("radialCumulativeDistribution", " (tau) - Likelihood - orbital radius function"),
    ("Symphony environments (4e8 < M < 4e9)", "Halo mass function - Likelihood - Symphony environments"),
])
def test_analysis_rows_find_history_recorded_under_another_name(analysis, recorded):
    result = {"attributes": {"name": analysis, "logLikelihood": "-5.0"}, "data": {}}
    history = {recorded: [{"value": value} for value in (4.0, 5.0, 6.0)]}
    row = report.analysis_rows("suffix", "Title", {"results": [result]}, {"results": [result]}, history, {})[0]
    assert row["recent"]["count"] == 3


def test_analysis_without_history_has_no_recent_range():
    result = {"attributes": {"name": "massFunction", "logLikelihood": "-5.0"}, "data": {}}
    history = {"Model - Likelihood - otherAnalysis": [{"value": 4.0}]}
    row = report.analysis_rows("suffix", "Title", {"results": [result]}, {"results": [result]}, history, {})[0]
    assert row["recent"] is None


def test_roundoff_in_target_is_not_a_changed_target():
    first = {"data": {"xDataset": [1.0, 2.0], "yDatasetTarget": [3.0, 4.0]}}
    second = {"data": {"xDataset": [1.0, 2.0], "yDatasetTarget": [3.0 * (1.0 + 1.0e-8), 4.0]}}
    assert not report.targets_differ(first, second)


@pytest.mark.parametrize("target", [[3.0, 4.1], [3.0], [3.0, None], None])
def test_real_change_in_target_is_detected(target):
    first = {"data": {"xDataset": [1.0, 2.0], "yDatasetTarget": [3.0, 4.0]}}
    second = {"data": {"xDataset": [1.0, 2.0], "yDatasetTarget": target}}
    assert report.targets_differ(first, second)


def comment_for_baseline(baseline, exact):
    model = {"analyses": [], "timings": [], "missing": []}
    meta = {"prCommit": "a" * 40, "baseSha": "b" * 40, "baselineCommit": baseline, "baselineExact": exact}
    return report.render_comment(model, meta)


def test_comment_warns_when_baseline_is_not_exact():
    comment = comment_for_baseline("c" * 40, False)
    assert "**Warning:**" in comment
    assert "other pull requests merged to master" in comment
    assert f"`{'b' * 10}`" in comment and f"`{'c' * 10}`" in comment


def test_comment_has_no_warning_when_baseline_is_exact():
    comment = comment_for_baseline("b" * 40, True)
    assert "**Warning:**" not in comment
    assert "the commit this run was built on" in comment


def test_comment_reports_absent_baseline():
    comment = comment_for_baseline(None, False)
    assert "**Warning:**" not in comment
    assert "No master results were found" in comment


def test_analysis_missing_from_pull_request_is_listed():
    result = {"attributes": {"name": "massFunction", "logLikelihood": "-5.0"}, "data": {}}
    row = report.analysis_rows("suffix", "Title", {"results": []}, {"results": [result]}, {}, {})[0]
    assert row["pr"] is None
    assert row["statusPr"] == "unknown"


def test_baseline_prefers_exact_commit(pages):
    _, metadata, exact = report.find_baseline(pages, "HEAD", "a" * 40)
    assert exact
    assert metadata["id"] == "a" * 40


def test_baseline_falls_back_to_nearest_earlier_results(pages):
    base = commit_dated(pages, "2026-01-02T00:00:00+00:00")
    _, metadata, exact = report.find_baseline(pages, "HEAD", base)
    assert not exact
    assert metadata["id"] == "a" * 40


def test_baseline_absent_when_all_results_are_newer(pages):
    base = commit_dated(pages, "2025-01-01T00:00:00+00:00")
    tree, metadata, exact = report.find_baseline(pages, "HEAD", base)
    assert (tree, metadata, exact) == (None, None, False)


def test_data_cannot_close_the_script_element():
    page = report.render_report("<x>__DATA__</x>", {}, {"note": "</script><script>alert(1)</script>"})
    assert "</script" not in page
    payload = page[len("<x>"):-len("</x>")]
    assert json.loads(payload)["note"] == "</script><script>alert(1)</script>"
