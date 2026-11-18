#!/usr/bin/env python3
"""Build a single, self-contained HTML report comparing a pull request's validation and benchmark results to master.

The report overlays the target data, the results from master, and the results from the pull request on every validation plot, and
tabulates each likelihood and timing metric against the recent history of that metric on master. Everything it needs (data,
Chart.js, styles) is inlined, because it is uploaded as a single unzipped CI artifact which GitHub then displays in the browser.

Inputs are the ``validate-*`` and ``benchmark-*`` artifacts of the pull request's CI run, and a git repository holding the
``gh-pages`` branch (a blob-less partial fetch of that branch is sufficient, as only a few files are ever read). The set of metrics
is taken from ``scripts/build/ghPagesMetricsManifest.yml``, so a metric added there appears in the report with no further work.

Results from master are read from the ``gh-pages`` commit which published the results for the commit that the pull request was
built on, or for the closest earlier commit if no results were published for that one. Using the latest results instead would
credit the pull request with (or blame it for) any change merged to master since its CI run began.

Run with::

    buildPrValidationReport.py --artifacts <dir> --gh-pages-repo <repo> --base-sha <sha> --output <report.html>
"""

import argparse
import hashlib
import json
import math
import statistics
import subprocess
import sys
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import yaml

from buildGhPagesIndex import (
    LOG_LIKELIHOOD_IMPOSSIBLE,
    classify_log_likelihood,
    load_manifest,
    load_thresholds,
    parse_js_data,
)

REPO_URL = "https://github.com/galacticusorg/galacticus"

# The JavaScript libraries inlined into the report. Their digests are checked so that the page cannot silently change.
LIBRARIES = {
    "__CHARTJS__":    ("https://cdn.jsdelivr.net/npm/chart.js@3.9.1/dist/chart.min.js",
                       "fbc45926e6b46845a0f905552a0e0b1331049bff1115ecf94dbe0904d895e710"),
    "__ERRORBARS__":  ("https://cdn.jsdelivr.net/npm/chartjs-chart-error-bars@3.9.1/build/index.umd.min.js",
                       "9490793f2fa25fc975779eb470fcd9a4539c47a00e8f647a1358448ed54543f3"),
}

# The number of most recent master results which define a metric's "recent range".
HISTORY_LENGTH = 30

# The number of gh-pages publications to search when looking for the results matching the commit a pull request was built on.
BASELINE_SEARCH_LIMIT = 40

# A result file which every master run publishes, used to find the commit that each gh-pages publication describes.
BASELINE_PROBE = "dev/valid/darkMatterOnlySubhalos/results.json"

# The relative difference below which two results are treated as numerically identical. Models which are nominally deterministic
# (the halo mass functions, for example) were measured to differ by up to 1e-6 between runs, with recent ranges of up to 9e-6.
ROUND_OFF_TOLERANCE = 1.0e-5

# The deviation reported for a change from a median of exactly zero, for which no scatter can be defined.
LARGE_DEVIATION = 1.0e6

LIKELIHOOD_MARKER = " - Likelihood - "

# Analyses whose benchmark history is recorded under a different name (by testSuite/validate-darkMatterOnlySubhalos-decayingDarkMatter.py).
HISTORY_NAME_ALIASES = {"massCumulativeDistribution": "mass function", "radialCumulativeDistribution": "orbital radius function"}

STATUS_RANK = {"unknown": -1, "ok": 0, "warn": 1, "fail": 2}


def sanitize(value):
    """Replace non-finite numbers with ``None``, so that the data can be embedded as strict JSON."""
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {key: sanitize(item) for key, item in value.items()}
    if isinstance(value, list):
        return [sanitize(item) for item in value]
    return value


def to_float(value):
    """Return ``value`` as a finite float, or ``None`` if it is not one."""
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def parse_timestamp(text):
    """Parse a commit timestamp, returning ``None`` if it cannot be read."""
    try:
        result = datetime.fromisoformat(str(text).replace("Z", "+00:00"))
    except ValueError:
        return None
    return result if result.tzinfo else result.replace(tzinfo=timezone.utc)


class GitTree:
    """Read access to files in one commit of a git repository."""

    def __init__(self, repository, commit):
        self.repository = repository
        self.commit = commit

    def read(self, path):
        """Return the text of ``path`` in this commit, or ``None`` if it does not exist."""
        result = subprocess.run(["git", "-C", str(self.repository), "show", f"{self.commit}:{path}"],
                                capture_output=True, check=False)
        return result.stdout.decode("utf-8") if result.returncode == 0 else None


def find_baseline(repository, ref, base_sha):
    """Find the ``gh-pages`` commit whose results best describe master at ``base_sha``.

    Returns ``(tree, results commit metadata, exact)``, where ``exact`` says whether the results are for ``base_sha`` itself. All
    three are ``None``/``False`` if no publication of results was found."""
    log = subprocess.run(["git", "-C", str(repository), "log", f"-n{BASELINE_SEARCH_LIMIT}", "--format=%H",
                          ref, "--", BASELINE_PROBE], capture_output=True, text=True, check=False)
    show = subprocess.run(["git", "-C", str(repository), "show", "-s", "--format=%cI", base_sha],
                          capture_output=True, text=True, check=False)
    base_time = parse_timestamp(show.stdout.strip()) if show.returncode == 0 else None
    newest = None
    for commit in log.stdout.split():
        tree = GitTree(repository, commit)
        text = tree.read(BASELINE_PROBE)
        data = parse_js_data(text) if text else None
        if not data:
            continue
        metadata = data.get("commit") or {}
        if metadata.get("id") == base_sha:
            return tree, metadata, True
        if newest is None:
            newest = (tree, metadata)
        published = parse_timestamp(metadata.get("timestamp"))
        if base_time is not None and published is not None and published <= base_time:
            return tree, metadata, False
    # Without a usable comparison time, the newest results are the best available.
    if newest is not None and base_time is None:
        return newest[0], newest[1], False
    return None, None, False


def read_bench_history(text):
    """Return ``{name: [{value, commit}, ...]}``, oldest first, from a github-action-benchmark ``data.js`` file."""
    data = parse_js_data(text) if text else None
    history = {}
    for entries in (data or {}).get("entries", {}).values():
        for entry in entries:
            commit = (entry.get("commit") or {}).get("id")
            for bench in entry.get("benches", []):
                value = to_float(bench.get("value"))
                if value is not None:
                    history.setdefault(bench["name"], []).append({"value": value, "commit": commit,
                                                                  "date": entry.get("date", 0),
                                                                  "unit": bench.get("unit", "")})
    for series in history.values():
        series.sort(key=lambda item: item["date"])
    return history


def read_pr_benches(artifacts, suffix):
    """Return ``{name: {value, unit}}`` from the benchmark and validation metric files of a pull request."""
    benches = {}
    for prefix in ("validate", "benchmark"):
        path = artifacts.get(f"{prefix}_{suffix}.json")
        if path is None:
            continue
        try:
            items = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        for item in items:
            value = to_float(item.get("value"))
            if value is not None:
                benches[item["name"]] = {"value": value, "unit": item.get("unit", "")}
    return benches


def recent_statistics(values):
    """Summarize recent master values: their count, range, median, and a robust estimate of their scatter."""
    values = [value for value in values if abs(value) < LOG_LIKELIHOOD_IMPOSSIBLE]
    if not values:
        return None
    median = statistics.median(values)
    return {"count": len(values), "low": min(values), "high": max(values), "median": median,
            "scatter": 1.4826 * statistics.median([abs(value - median) for value in values])}


def compare_to_history(value, recent):
    """Return ``(outside recent range, deviation from the median in units of the scatter)``; the deviation is ``None`` if unknown.

    The scatter is floored at the round-off level, so that deterministic models do not report round-off as a large deviation."""
    if value is None or recent is None:
        return False, None
    tolerance = ROUND_OFF_TOLERANCE * abs(recent["median"])
    outside = value < recent["low"] - tolerance or value > recent["high"] + tolerance
    scatter = max(recent["scatter"], tolerance)
    if scatter <= 0.0:
        return outside, 0.0 if value == recent["median"] else math.copysign(LARGE_DEVIATION, value - recent["median"])
    return outside, (value - recent["median"]) / scatter


def values_differ(first, second):
    """Determine whether two datasets differ by more than round-off. Some targets are computed afresh in each run, and differ by up to 1e-8."""
    if not isinstance(first, list) or not isinstance(second, list) or len(first) != len(second):
        return first != second
    for value_first, value_second in zip(first, second):
        if isinstance(value_first, (int, float)) and isinstance(value_second, (int, float)):
            if abs(value_first - value_second) > ROUND_OFF_TOLERANCE * max(abs(value_first), abs(value_second)):
                return True
        elif value_first != value_second:
            return True
    return False


def targets_differ(pr_analysis, master_analysis):
    """Determine whether the target data of an analysis differ between the pull request and master."""
    first, second = pr_analysis["data"], master_analysis["data"]
    return values_differ(first.get("yDatasetTarget"), second.get("yDatasetTarget")) or values_differ(first.get("xDataset"), second.get("xDataset"))


def history_for_analysis(series, name):
    """Find the benchmark history of an analysis, whose name there may be an alias, or lack the analysis' parenthetical description."""
    for candidate in (name, HISTORY_NAME_ALIASES.get(name), name.split(" (", 1)[0]):
        if candidate in series:
            return series[candidate]
    return []


def analysis_rows(suffix, title, pr, master, history, thresholds):
    """Build one table row for each analysis in a standard validation metric."""
    master_results = {result["attributes"].get("name"): result for result in (master or {}).get("results", [])}
    pr_results = {result["attributes"].get("name"): result for result in (pr or {}).get("results", [])}
    series = {name.split(LIKELIHOOD_MARKER, 1)[1]: items for name, items in history.items() if LIKELIHOOD_MARKER in name}
    rows = []
    for name in list(pr_results) + [name for name in master_results if name not in pr_results]:
        pr_result, master_result = pr_results.get(name), master_results.get(name)
        pr_ll = to_float(pr_result["attributes"].get("logLikelihood")) if pr_result else None
        master_ll = to_float(master_result["attributes"].get("logLikelihood")) if master_result else None
        # History records -log(likelihood), since github-action-benchmark assumes that smaller is better.
        recent = recent_statistics([-item["value"] for item in history_for_analysis(series, name)[-HISTORY_LENGTH:]])
        outside, deviation = compare_to_history(pr_ll, recent)
        threshold = thresholds.get(name)
        rows.append({"metric": suffix, "metricTitle": title, "analysis": name, "pr": pr_ll, "master": master_ll, "recent": recent,
                     "outside": outside, "deviation": deviation,
                     "statusPr": "unknown" if pr_ll is None else classify_log_likelihood(pr_ll, threshold),
                     "statusMaster": "unknown" if master_ll is None else classify_log_likelihood(master_ll, threshold),
                     "targetChanged": bool(pr_result and master_result and targets_differ(pr_result, master_result))})
    return rows


def timing_rows(suffix, title, pr_benches, history, alert):
    """Build one table row for each benchmark metric which is not a likelihood (wall time, memory use, and so on)."""
    rows = []
    for name, bench in pr_benches.items():
        if LIKELIHOOD_MARKER in name:
            continue
        series = history.get(name, [])
        recent = recent_statistics([item["value"] for item in series[-HISTORY_LENGTH:]])
        outside, deviation = compare_to_history(bench["value"], recent)
        rows.append({"metric": suffix, "metricTitle": title, "name": name, "unit": bench["unit"], "pr": bench["value"],
                     "master": series[-1]["value"] if series else None, "recent": recent, "outside": outside,
                     "deviation": deviation, "alert": alert})
    return rows


def read_alert_thresholds(workflow_path):
    """Return ``{suffix: alert ratio}`` for the timing benchmarks configured in the CI workflow."""
    try:
        with open(workflow_path, "r", encoding="utf-8") as stream:
            workflow = yaml.safe_load(stream)
    except (OSError, yaml.YAMLError):
        return {}
    alerts = {}
    for step in workflow.get("jobs", {}).get("Validate", {}).get("steps", []):
        inputs = step.get("with") or {}
        if step.get("uses") == "./.github/actions/benchmark" and inputs.get("type") == "benchmark":
            alerts[inputs["suffix"]] = float(str(inputs["threshold"]).rstrip("%")) / 100.0
    return alerts


def metric_kind(entry, results_file):
    """Return how a metric's results are displayed, or ``None`` if the report does not cover it."""
    if entry.get("aggregator") == "standardLogLikelihood":
        return "standard"
    if results_file == "results_PonosV.json":
        return "ponosV"
    if results_file == "results_baryonicSuppression.json":
        return "baryonicSuppression"
    # The convergence metrics have a bespoke format, and are not produced by CI runs.
    return None


def load_results(text):
    """Decode a results file, returning sanitized data or ``None``."""
    data = parse_js_data(text) if text else None
    return sanitize(data) if data is not None else None


def build_model(arguments, manifest, artifacts, tree):
    """Assemble the data embedded in the report."""
    thresholds = load_thresholds(arguments.thresholds)
    alerts = read_alert_thresholds(arguments.workflow)
    metrics, analyses, timings, missing = {}, [], [], []
    groups = {key: {"key": key, "title": group.get("title", key), "order": group.get("nav_order", 0), "metrics": []}
              for key, group in manifest.get("groups", {}).items()}
    for suffix, entry in manifest.get("metrics", {}).items():
        if entry.get("direct_path"):
            continue
        results_file = entry.get("valid_results_file", "results.json")
        kind = metric_kind(entry, results_file) if entry.get("has_valid") else None
        group = groups.get(entry.get("group"))
        title = f"{group['title']}: {entry.get('display_name', suffix)}" if group else entry.get("display_name", suffix)
        history = {}
        if entry.get("has_bench") and tree is not None:
            history = read_bench_history(tree.read(f"dev/bench/{suffix}/data.js"))
        pr_benches = read_pr_benches(artifacts, suffix)
        timings.extend(timing_rows(suffix, title, pr_benches, history, alerts.get(suffix)))
        if kind is None:
            continue
        source = results_file if results_file != "results.json" else f"results_{suffix}.json"
        path = artifacts.get(source)
        pr = load_results(path.read_text(encoding="utf-8")) if path else None
        master = load_results(tree.read(f"dev/valid/{suffix}/{results_file}")) if tree is not None else None
        if pr is None:
            missing.append(title)
        if pr is None and master is None:
            continue
        rows = []
        if kind == "standard":
            rows = analysis_rows(suffix, title, pr, master, history, thresholds.get(suffix) or {})
            analyses.extend(rows)
        metrics[suffix] = {"title": title, "shortTitle": entry.get("display_name", suffix), "kind": kind, "pr": pr, "master": master,
                           "missing": pr is None,
                           "worstStatus": max((row["statusPr"] for row in rows), key=STATUS_RANK.get, default="unknown"),
                           "flagged": sum(1 for row in rows if row["outside"])}
        if group:
            group["metrics"].append(suffix)
    return {"metrics": metrics, "analyses": analyses, "timings": timings, "missing": missing,
            "groups": sorted((group for group in groups.values() if group["metrics"]), key=lambda group: group["order"]),
            "historyLength": HISTORY_LENGTH}


def fetch_libraries():
    """Download the JavaScript libraries, checking their digests."""
    libraries = {}
    for placeholder, (url, digest) in LIBRARIES.items():
        with urllib.request.urlopen(url, timeout=120) as response:
            content = response.read()
        if hashlib.sha256(content).hexdigest() != digest:
            sys.exit(f"Digest of {url} does not match the expected value.")
        libraries[placeholder] = content.decode("utf-8").replace("</script", "<\\/script")
    return libraries


def render_report(template, libraries, model):
    """Insert the libraries and data into the page template."""
    # Pull request data are not trusted: escaping "<" guarantees that they cannot close the script element which holds them.
    payload = json.dumps(sanitize(model), separators=(",", ":"), allow_nan=False).replace("<", "\\u003c")
    page = template.replace("__DATA__", payload)
    for placeholder, content in libraries.items():
        page = page.replace(placeholder, content)
    return page


def short(sha):
    """Return the abbreviated form of a commit hash."""
    return (sha or "unknown")[:10]


def render_comment(model, meta):
    """Return the Markdown of the pull request comment. The workflow replaces ``__REPORT_URL__`` once the report is uploaded."""
    analyses = model["analyses"]
    outside = [row for row in analyses if row["outside"]]
    worse = [row for row in analyses if STATUS_RANK[row["statusPr"]] > STATUS_RANK[row["statusMaster"]] >= 0]
    better = [row for row in analyses if 0 <= STATUS_RANK[row["statusPr"]] < STATUS_RANK[row["statusMaster"]]]
    slower = [row for row in model["timings"] if row["outside"] and row["recent"] and row["pr"] > row["recent"]["high"]]
    lines = ["<!-- This comment is updated by each CI run of this pull request. -->",
             f"### Validation report for `{short(meta['prCommit'])}`", "",
             "**[Open the report](__REPORT_URL__)**: target data, master, and this pull request on every validation plot, "
             "plus likelihoods and timings against master's recent history. Viewing it requires a GitHub login.", ""]
    if meta["baselineCommit"] and meta["baselineExact"]:
        lines.append(f"Master results are from `{short(meta['baselineCommit'])}`, the commit this run was built on.")
    elif meta["baselineCommit"]:
        lines.append(f"> **Warning:** results for master commit `{short(meta['baseSha'])}`, which this run was built on, are not available yet, "
                     f"so this comparison uses `{short(meta['baselineCommit'])}`, the nearest master commit with results. Differences may come "
                     f"from other pull requests merged to master in between, not from this pull request. Once master's CI/CD run for "
                     f"`{short(meta['baseSha'])}` has finished, re-run this job for an exact comparison.")
    else:
        lines.append("**No master results were found to compare against.**")
    lines += ["", f"- Analyses compared: {len(analyses)}",
              f"- Status worse than master: {len(worse)}, better than master: {len(better)}",
              f"- Likelihoods outside master's last {HISTORY_LENGTH} results: {len(outside)} "
              f"(about {2.0 * len(analyses) / (HISTORY_LENGTH + 1.0):.0f} expected by chance)",
              f"- Timings above master's last {HISTORY_LENGTH} results: {len(slower)}"]
    if model["missing"]:
        lines.append(f"- **Metrics with no result from this run: {len(model['missing'])}** (a job failed or was skipped)")
    if worse:
        lines += ["", "Status changes for the worse:", ""]
        lines += [f"- {row['metricTitle']}, `{row['analysis']}`: {row['statusMaster']} to {row['statusPr']}" for row in worse[:10]]
        if len(worse) > 10:
            lines.append(f"- and {len(worse) - 10} more in the report")
    return "\n".join(lines) + "\n"


def main(argv=None):
    """Build the report and the pull request comment."""
    script_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__.split("\n", maxsplit=1)[0])
    parser.add_argument("--artifacts", required=True, type=Path, help="directory holding the downloaded CI artifacts")
    parser.add_argument("--gh-pages-repo", required=True, type=Path, help="git repository holding the gh-pages branch")
    parser.add_argument("--gh-pages-ref", default="origin/gh-pages", help="git reference of the gh-pages branch")
    parser.add_argument("--base-sha", required=True, help="the master commit which the pull request was built on")
    parser.add_argument("--output", required=True, type=Path, help="file to which to write the report")
    parser.add_argument("--comment", type=Path, help="file to which to write the Markdown of the pull request comment")
    parser.add_argument("--manifest", type=Path, default=script_dir / "ghPagesMetricsManifest.yml")
    parser.add_argument("--thresholds", type=Path, default=script_dir / "ghPagesAnalysisThresholds.yml")
    parser.add_argument("--workflow", type=Path, default=script_dir.parents[1] / ".github/workflows/cicd.yml")
    parser.add_argument("--template", type=Path, default=script_dir / "prValidationReportTemplate.html")
    parser.add_argument("--pr-number", default="", help="pull request number")
    parser.add_argument("--pr-commit", default="", help="head commit of the pull request")
    parser.add_argument("--run-url", default="", help="URL of the CI run")
    arguments = parser.parse_args(argv)

    artifacts = {path.name: path for path in sorted(arguments.artifacts.rglob("*.json"))}
    manifest = load_manifest(arguments.manifest)
    tree, baseline, exact = find_baseline(arguments.gh_pages_repo, arguments.gh_pages_ref, arguments.base_sha)
    model = build_model(arguments, manifest, artifacts, tree)
    meta = {"prNumber": arguments.pr_number, "prCommit": arguments.pr_commit, "baseSha": arguments.base_sha, "runUrl": arguments.run_url,
            "repoUrl": REPO_URL, "baselineCommit": (baseline or {}).get("id"), "baselineMessage": (baseline or {}).get("message", ""),
            "baselineDate": (baseline or {}).get("timestamp"), "baselineExact": exact,
            "generated": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")}
    model["meta"] = meta
    template = arguments.template.read_text(encoding="utf-8")
    arguments.output.write_text(render_report(template, fetch_libraries(), model), encoding="utf-8")
    if arguments.comment:
        arguments.comment.write_text(render_comment(model, meta), encoding="utf-8")
    print(f"Wrote {arguments.output}: {len(model['metrics'])} metrics, {len(model['analyses'])} analyses, "
          f"{len(model['timings'])} timings, {len(model['missing'])} metrics missing from this run; "
          f"baseline {short(meta['baselineCommit'])} (exact={exact}).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
