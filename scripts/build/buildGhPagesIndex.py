#!/usr/bin/env python3
"""Regenerate the gh-pages landing page, dashboard, and sidebar pages.

Reads ``scripts/build/ghPagesMetricsManifest.yml`` and walks a checkout of
the ``gh-pages`` branch to:

  1. Migrate any legacy split-by-suffix index/results files under
     ``dev/valid/idealizedSubhaloSimulations/`` and
     ``dev/valid/darkMatterOnlySubhalosDecayingDarkMatter/`` into
     per-variant directories matching the bench-side convention. Idempotent:
     skips variants already migrated.
  2. Inspect each metric's ``data.js`` (bench) and ``results.json`` /
     ``results_*.json`` (valid) to extract the latest commit and a status
     signal (per the aggregator declared in the manifest).
  3. Write Jekyll site files: ``_config.yml``, ``Gemfile``,
     ``_data/metrics.yml``, ``index.md``, ``dashboard.md``, and one
     ``groups/<key>.md`` page per group.

Run with::

    buildGhPagesIndex.py <gh-pages-checkout-dir>
"""

import argparse
import json
import os
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import yaml

REPO_URL = "https://github.com/galacticusorg/galacticus"

JS_PREFIX_RE = re.compile(
    r"^\s*window\.[A-Z_]+\s*=\s*", flags=re.MULTILINE
)


def parse_js_data(text):
    """Strip ``window.FOO = `` prefix and trailing semicolon from a JS data
    file written by github-action-benchmark or python/validate.py, and
    return the remaining JSON-decoded body. Returns ``None`` on failure."""
    text = text.strip()
    m = JS_PREFIX_RE.match(text)
    if not m:
        return None
    body = text[m.end():].rstrip().rstrip(";").strip()
    try:
        return json.loads(body)
    except json.JSONDecodeError:
        return None


def load_manifest(manifest_path):
    with open(manifest_path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)


def migrate_legacy_shared_dir(gh_root, shared_dir, suffix_extractor):
    """Move ``index_<suffix>.html`` / ``results_<suffix>.json`` out of a
    shared directory into per-variant subdirectories under ``dev/valid``.
    ``suffix_extractor(filename) -> (variant_suffix, full_results_filename)``
    returns the destination directory name and the corresponding results
    filename to look for. Leaves a ``<meta refresh>`` redirect at the old
    path so external links keep working. No-ops if the variant already
    exists at the new location."""
    shared = gh_root / shared_dir
    if not shared.is_dir():
        return []
    moved = []
    for html in sorted(shared.glob("index_*.html")):
        info = suffix_extractor(html.name)
        if info is None:
            continue
        variant, results_name = info
        new_dir = gh_root / "dev" / "valid" / variant
        new_html = new_dir / "index.html"
        new_json = new_dir / "results.json"
        new_dir.mkdir(parents=True, exist_ok=True)
        # Create the per-variant index.html from the old shared-dir template
        # only if it doesn't already exist; rewrite its script src to point
        # at the canonical "results.json" filename.
        if not new_html.exists():
            body = html.read_text(encoding="utf-8")
            body = body.replace(
                f'src="{results_name}"', 'src="results.json"'
            )
            new_html.write_text(body, encoding="utf-8")
            moved.append(variant)
        # Reconcile the old results JSON in the shared dir with any fresh
        # copy already published in the per-variant dir. Never overwrite a
        # fresh per-variant results.json with the stale shared-dir copy.
        old_json = shared / results_name
        if old_json.exists():
            if not new_json.exists():
                old_json.replace(new_json)
            else:
                old_json.unlink()
        # Always replace the old shared-dir index with a redirect stub.
        _write_redirect_stub(html, variant)
    return moved


def _write_redirect_stub(old_path, variant):
    """Overwrite ``old_path`` with a minimal HTML refresh redirect to the
    new per-variant location."""
    target = f"../{variant}/"
    old_path.write_text(
        "<!DOCTYPE html>\n"
        f'<meta http-equiv="refresh" content="0; url={target}">\n'
        f'<title>Moved</title>\n'
        f'<p>This page has moved to <a href="{target}">{target}</a>.</p>\n',
        encoding="utf-8",
    )


def _idealized_extractor(name):
    # index_idealizedSubhaloSimulation_rpra0.005_gamma1.5.html
    m = re.match(r"^index_(idealizedSubhaloSimulation_.+)\.html$", name)
    if not m:
        return None
    variant = m.group(1)
    return variant, f"results_{variant}.json"


def _decaying_extractor(name):
    # index_lifetime10.0_velocityKick20.0.html (note: prefix omitted)
    m = re.match(r"^index_(lifetime[\d.]+_velocityKick[\d.]+)\.html$", name)
    if not m:
        return None
    tail = m.group(1)
    variant = f"darkMatterOnlySubhalos_decayingDarkMatter_{tail}"
    return variant, f"results_{variant}.json"


def inspect_bench(gh_root, suffix):
    """Read ``dev/bench/<suffix>/data.js`` and return a small dict with
    ``last_commit_id``, ``last_commit_short``, ``last_commit_url``,
    ``last_ts``, ``last_ts_iso``, ``benches`` (list of ``{name, value, unit}``).
    Returns ``None`` if the file is missing or unparseable."""
    f = gh_root / "dev" / "bench" / suffix / "data.js"
    if not f.is_file():
        return None
    data = parse_js_data(f.read_text(encoding="utf-8"))
    if not data or "entries" not in data:
        return None
    # entries is a dict of group-name -> list of entries; pick the most
    # recently updated entry across all groups.
    latest = None
    for entries in data["entries"].values():
        for e in entries:
            if latest is None or e.get("date", 0) > latest.get("date", 0):
                latest = e
    if latest is None:
        return None
    commit = latest.get("commit") or {}
    ts_ms = latest.get("date") or 0
    ts = datetime.fromtimestamp(ts_ms / 1000.0, tz=timezone.utc) if ts_ms else None
    return {
        "last_commit_id":    commit.get("id", ""),
        "last_commit_short": (commit.get("id") or "")[:7],
        "last_commit_url":   commit.get("url", ""),
        "last_ts":           ts_ms,
        "last_ts_iso":       ts.isoformat() if ts else "",
        "benches":           latest.get("benches", []),
    }


def inspect_valid(gh_root, suffix, results_filename="results.json"):
    """Read ``dev/valid/<suffix>/<results_filename>`` and return a small
    dict with last commit info plus the parsed payload under ``data``.
    Returns ``None`` if missing or unparseable."""
    f = gh_root / "dev" / "valid" / suffix / results_filename
    if not f.is_file():
        return None
    data = parse_js_data(f.read_text(encoding="utf-8"))
    if data is None:
        return None
    commit = data.get("commit") or {}
    ts_iso = commit.get("timestamp", "")
    ts_ms = 0
    if ts_iso:
        try:
            ts_ms = int(datetime.fromisoformat(ts_iso).timestamp() * 1000)
        except ValueError:
            ts_ms = 0
    return {
        "last_commit_id":    commit.get("id", ""),
        "last_commit_short": (commit.get("id") or "")[:7],
        "last_commit_url":   commit.get("url", ""),
        "last_ts":           ts_ms,
        "last_ts_iso":       ts_iso,
        "data":              data,
    }


# --- Aggregators ---------------------------------------------------------

def _agg_standard_log_likelihood(valid, defaults, metric):
    """Classify each analysis individually against its per-analysis
    ``threshold_warn`` / ``threshold_fail`` (sourced from
    ``ghPagesAnalysisThresholds.yml`` and injected into ``metric["analyses"]``
    by ``build_records``). Aggregate by reporting the worst status across
    analyses, so a single bad analysis can't be hidden by another that
    happens to be much worse on an absolute scale.

    A metric whose analyses all lack thresholds reports as ``unknown``;
    individual analyses without thresholds are noted in the detail string
    but don't influence the overall status as long as at least one
    analysis is thresholded."""
    results = valid["data"].get("results") or []
    analyses_thresh = metric.get("analyses") or {}
    statuses = []  # list of (status, analysis_name, logL)
    for r in results:
        attrs = r.get("attributes") or {}
        name = attrs.get("name", "?")
        v = attrs.get("logLikelihood")
        if v is None:
            continue
        try:
            ll = float(v)
        except (TypeError, ValueError):
            continue
        th = analyses_thresh.get(name) or {}
        warn = th.get("threshold_warn")
        fail = th.get("threshold_fail")
        if warn is None or fail is None:
            statuses.append(("unknown", name, ll))
        elif ll >= warn:
            statuses.append(("ok", name, ll))
        elif ll >= fail:
            statuses.append(("warn", name, ll))
        else:
            statuses.append(("fail", name, ll))
    if not statuses:
        return "unknown", "no logLikelihood found", None
    thresholded = [s for s in statuses if s[0] != "unknown"]
    n_total = len(statuses)
    n_unk = n_total - len(thresholded)
    if not thresholded:
        return "unknown", f"no thresholds set for {n_total} analyses", None
    rank = {"ok": 0, "warn": 1, "fail": 2}
    overall = max((s[0] for s in thresholded), key=lambda s: rank[s])
    if overall in ("warn", "fail"):
        offenders = [n for s, n, _ in statuses if s == overall]
        head = ", ".join(offenders[:2]) + ("…" if len(offenders) > 2 else "")
        detail = f"{len(offenders)}/{n_total} {overall}: {head}"
    else:
        detail = f"all {len(thresholded)} thresholded analyses pass"
    if n_unk:
        detail += f" ({n_unk} without thresholds)"
    return overall, detail, None


def _agg_ponosv(valid, defaults, metric):
    """PonosV: pass criterion is that no less than 5% of realizations lie
    above/below either of the target normalization or slope values."""
    d = valid["data"]
    pcts = []
    for key in ("surfaceDensity", "slope"):
        sub = d.get(key) or {}
        for k in ("percentageBelow", "percentageAbove"):
            if k in sub:
                try:
                    pcts.append(float(sub[k]))
                except (TypeError, ValueError):
                    pass
    if not pcts:
        return "unknown", "no percentile fields", None
    worst = min(pcts)
    th = float(metric.get("threshold_fail", 5.0))
    status = "ok" if worst >= th else "fail"
    return status, f"min realization fraction = {worst:.1f}%", worst


AGGREGATORS = {
    "standardLogLikelihood": _agg_standard_log_likelihood,
    "ponosVPercentile":      _agg_ponosv,
    "none":                  lambda v, d, m: ("unknown", "bespoke schema", None),
}


# --- Threshold sidecar ---------------------------------------------------

def load_thresholds(thresholds_path):
    """Load the per-analysis logLikelihood thresholds emitted by
    ``--bootstrap-thresholds``. Returns a ``{metric_suffix: {analysis_name:
    {threshold_warn, threshold_fail, current}}}`` map, or an empty dict
    when the sidecar is missing."""
    p = Path(thresholds_path) if thresholds_path else None
    if p is None or not p.is_file():
        return {}
    with open(p, "r", encoding="utf-8") as f:
        data = yaml.safe_load(f) or {}
    return data.get("thresholds", {}) or {}


def bootstrap_thresholds(gh_root, manifest, fraction_warn=0.10, fraction_fail=0.20):
    """Read each ``standardLogLikelihood`` metric's currently-published
    logL values and derive per-analysis warn/fail thresholds at
    ``current - fraction * |current|``. Returns a nested dict keyed by
    metric suffix, then analysis name."""
    out = {}
    for suffix, m in (manifest.get("metrics") or {}).items():
        if m.get("aggregator") != "standardLogLikelihood":
            continue
        if not m.get("has_valid"):
            continue
        results_file = m.get("valid_results_file", "results.json")
        valid = inspect_valid(gh_root, suffix, results_file)
        if not valid:
            continue
        analyses = {}
        for r in (valid["data"].get("results") or []):
            attrs = r.get("attributes") or {}
            name = attrs.get("name")
            ll = attrs.get("logLikelihood")
            if name is None or ll is None:
                continue
            try:
                current = float(ll)
            except (TypeError, ValueError):
                continue
            warn = current - fraction_warn * abs(current)
            fail = current - fraction_fail * abs(current)
            analyses[name] = {
                "current":         round(current, 4),
                "threshold_warn":  round(warn,    4),
                "threshold_fail":  round(fail,    4),
            }
        if analyses:
            out[suffix] = analyses
    return out


# --- Record building -----------------------------------------------------

def build_records(gh_root, manifest, thresholds_map=None):
    """Return a list of metric records (dicts) ready to dump as Jekyll
    data, plus the group definitions in declaration order."""
    defaults = (manifest.get("defaults") or {}).get("status") or {}
    groups = manifest.get("groups") or {}
    metrics = manifest.get("metrics") or {}
    thresholds_map = thresholds_map or {}
    records = []
    for suffix, m in metrics.items():
        rec = {
            "suffix":       suffix,
            "group":        m.get("group", "misc"),
            "display_name": m.get("display_name", suffix),
            "sort_key":     m.get("sort_key", suffix),
            "has_bench":    bool(m.get("has_bench")),
            "has_valid":    bool(m.get("has_valid")),
            "direct_path":  m.get("direct_path", ""),
            "bench_path":   "",
            "valid_path":   "",
            "last_commit_id":    "",
            "last_commit_short": "",
            "last_commit_url":   "",
            "last_ts":           0,
            "last_ts_iso":       "",
            "status":            "unknown",
            "status_detail":     "",
            "status_value":      None,
        }
        # Bench inspection.
        bench = inspect_bench(gh_root, suffix) if rec["has_bench"] else None
        if bench:
            rec["bench_path"] = f"dev/bench/{suffix}/"
            rec["last_commit_id"]    = bench["last_commit_id"]
            rec["last_commit_short"] = bench["last_commit_short"]
            rec["last_commit_url"]   = bench["last_commit_url"]
            rec["last_ts"]           = bench["last_ts"]
            rec["last_ts_iso"]       = bench["last_ts_iso"]
        # Valid inspection (overrides commit info if newer).
        if rec["has_valid"]:
            results_file = m.get("valid_results_file", "results.json")
            valid = inspect_valid(gh_root, suffix, results_file)
            if valid:
                rec["valid_path"] = f"dev/valid/{suffix}/"
                if valid["last_ts"] >= rec["last_ts"]:
                    rec["last_commit_id"]    = valid["last_commit_id"]
                    rec["last_commit_short"] = valid["last_commit_short"]
                    rec["last_commit_url"]   = valid["last_commit_url"]
                    rec["last_ts"]           = valid["last_ts"]
                    rec["last_ts_iso"]       = valid["last_ts_iso"]
                agg = AGGREGATORS.get(m.get("aggregator", "none"), AGGREGATORS["none"])
                m_with_thresh = dict(m)
                m_with_thresh["analyses"] = thresholds_map.get(suffix, {})
                status, detail, value = agg(valid, defaults, m_with_thresh)
                rec["status"]        = status
                rec["status_detail"] = detail
                rec["status_value"]  = value
        if rec["direct_path"]:
            rec["status"] = "ok"
            rec["status_detail"] = "static artifact"
        records.append(rec)
    # Sort by group nav_order, then sort_key, then suffix.
    def sort_tuple(r):
        g = groups.get(r["group"], {}) or {}
        return (g.get("nav_order", 1000), r.get("sort_key") or "", r["suffix"])
    records.sort(key=sort_tuple)
    return records, groups


# --- File emitters -------------------------------------------------------

CONFIG_YML = """\
title: Galacticus
description: Benchmark and validation metrics for the Galacticus galaxy formation model.
url: https://galacticusorg.github.io
baseurl: /galacticus
remote_theme: just-the-docs/just-the-docs
plugins:
  - jekyll-remote-theme
search_enabled: true
heading_anchors: true
color_scheme: light
aux_links:
  "Galacticus on GitHub":
    - https://github.com/galacticusorg/galacticus
  "Wiki":
    - https://github.com/galacticusorg/galacticus/wiki
aux_links_new_tab: true
logo: /assets/New_Logo_Galaxy_192_Transparent.png
defaults:
  - scope:
      path: "dev"
    values:
      sitemap: false
exclude:
  - Gemfile
  - Gemfile.lock
  - vendor
include:
  - _data
  - _includes
  - groups
"""

GEMFILE = """\
source "https://rubygems.org"

gem "github-pages", group: :jekyll_plugins
gem "jekyll-remote-theme"
"""


def write_if_changed(path, content):
    """Write ``content`` to ``path`` only if different from the current
    contents. Returns True if a write occurred."""
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if path.read_text(encoding="utf-8") == content:
            return False
    path.write_text(content, encoding="utf-8")
    return True


def emit_metrics_data(gh_root, records, groups):
    payload = {
        "groups":  groups,
        "metrics": records,
        "generated_at": datetime.now(tz=timezone.utc).isoformat(),
    }
    text = yaml.safe_dump(payload, sort_keys=False, allow_unicode=True)
    write_if_changed(gh_root / "_data" / "metrics.yml", text)


def _status_badge_md(status):
    # just-the-docs / kramdown inline attribute list: the IAL follows the
    # element it modifies. ``Status``{: .label .label-COLOR } renders as a
    # colored pill in the just-the-docs theme.
    return {
        "ok":      "`PASS`{: .label .label-green }",
        "warn":    "`WARN`{: .label .label-yellow }",
        "fail":    "`FAIL`{: .label .label-red }",
        "unknown": "`N/A`{: .label .label-grey }",
    }.get(status, "`N/A`")


def emit_index(gh_root, records, groups):
    # Group records for the home page summary cards.
    by_group = {}
    for r in records:
        by_group.setdefault(r["group"], []).append(r)
    parts = [
        "---",
        "layout: default",
        "title: Home",
        "nav_order: 1",
        "---",
        "",
        "# Galacticus",
        "",
        "![logo](assets/New_Logo_Galaxy_192_Transparent.png)",
        "",
        "Benchmark and validation metrics for the "
        "[Galacticus](https://github.com/galacticusorg/galacticus) "
        "galaxy formation model. Documentation lives on the "
        "[wiki](https://github.com/galacticusorg/galacticus/wiki).",
        "",
        "See the [status dashboard](dashboard.html) for a one-page summary "
        "of every metric, or browse by group below.",
        "",
    ]
    for key, g in sorted(groups.items(), key=lambda kv: kv[1].get("nav_order", 1000)):
        items = by_group.get(key, [])
        if not items:
            continue
        parts.append(f"## [{g.get('title', key)}](groups/{key}.html)")
        parts.append("")
        if g.get("description"):
            parts.append(g["description"])
            parts.append("")
        if g.get("citation"):
            url = g.get("citation_url", "")
            cite = f"[{g['citation']}]({url})" if url else g["citation"]
            parts.append(f"_{cite}_")
            parts.append("")
        parts.append(f"{len(items)} metric(s).")
        parts.append("")
    write_if_changed(gh_root / "index.md", "\n".join(parts) + "\n")


def emit_dashboard(gh_root, records, groups):
    parts = [
        "---",
        "layout: default",
        "title: Dashboard",
        "nav_order: 2",
        "---",
        "",
        "# Status dashboard",
        "",
        "Latest pass / warn / fail status for every published validation, "
        "plus the most recent commit that produced each result. Click a "
        "metric name to drill into the detail page.",
        "",
        "Generated " + datetime.now(tz=timezone.utc).strftime("%Y-%m-%d %H:%M UTC") + ".",
        "",
        "| Group | Metric | Status | Detail | Last commit | When |",
        "|-------|--------|--------|--------|-------------|------|",
    ]
    for r in records:
        group_title = (groups.get(r["group"]) or {}).get("title", r["group"])
        if r["valid_path"]:
            link = f"[{r['display_name']}]({r['valid_path']})"
        elif r["bench_path"]:
            link = f"[{r['display_name']}]({r['bench_path']})"
        elif r["direct_path"]:
            link = f"[{r['display_name']}]({r['direct_path']})"
        else:
            link = r["display_name"]
        commit = (
            f"[`{r['last_commit_short']}`]({r['last_commit_url']})"
            if r["last_commit_short"] else "—"
        )
        when = r["last_ts_iso"][:10] if r["last_ts_iso"] else "—"
        parts.append(
            f"| {group_title} | {link} | {_status_badge_md(r['status'])} "
            f"| {r['status_detail'] or '—'} | {commit} | {when} |"
        )
    parts.append("")
    write_if_changed(gh_root / "dashboard.md", "\n".join(parts) + "\n")


def emit_group_pages(gh_root, records, groups):
    by_group = {}
    for r in records:
        by_group.setdefault(r["group"], []).append(r)
    out_dir = gh_root / "groups"
    out_dir.mkdir(parents=True, exist_ok=True)
    for key, g in groups.items():
        items = by_group.get(key, [])
        # nav_order: groups come after Home (1) and Dashboard (2).
        nav = 10 + (g.get("nav_order", 1000) // 10)
        parts = [
            "---",
            "layout: default",
            f"title: \"{g.get('title', key)}\"",
            f"nav_order: {nav}",
            "---",
            "",
            f"# {g.get('title', key)}",
            "",
        ]
        if g.get("description"):
            parts.append(g["description"])
            parts.append("")
        if g.get("citation"):
            url = g.get("citation_url", "")
            cite = f"[{g['citation']}]({url})" if url else g["citation"]
            parts.append(f"_Reference: {cite}_")
            parts.append("")
        if not items:
            parts.append("_No metrics in this group yet._")
        else:
            parts.append("| Metric | Status | Validation | Benchmark | Last commit |")
            parts.append("|--------|--------|------------|-----------|-------------|")
            for r in items:
                valid_link = f"[plots]({{{{ '/{r['valid_path']}' | relative_url }}}})" if r["valid_path"] else "—"
                bench_link = f"[trend]({{{{ '/{r['bench_path']}' | relative_url }}}})" if r["bench_path"] else "—"
                if not r["valid_path"] and not r["bench_path"] and r["direct_path"]:
                    valid_link = f"[open]({{{{ '/{r['direct_path']}' | relative_url }}}})"
                commit = (
                    f"[`{r['last_commit_short']}`]({r['last_commit_url']})"
                    if r["last_commit_short"] else "—"
                )
                parts.append(
                    f"| {r['display_name']} "
                    f"| {_status_badge_md(r['status'])} "
                    f"| {valid_link} | {bench_link} | {commit} |"
                )
        parts.append("")
        write_if_changed(out_dir / f"{key}.md", "\n".join(parts) + "\n")


# --- Entry point ---------------------------------------------------------

def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    p.add_argument("gh_pages_dir", help="Path to a checkout of the gh-pages branch")
    p.add_argument(
        "--manifest",
        default=str(Path(__file__).with_name("ghPagesMetricsManifest.yml")),
        help="Path to the metrics manifest YAML (default: alongside this script).",
    )
    p.add_argument(
        "--thresholds",
        default=str(Path(__file__).with_name("ghPagesAnalysisThresholds.yml")),
        help="Path to the per-analysis logL threshold sidecar YAML.",
    )
    p.add_argument(
        "--bootstrap-thresholds",
        action="store_true",
        help=(
            "Read each metric's currently-published logL values, write a "
            "fresh per-analysis threshold sidecar at --thresholds, and exit "
            "without emitting any site files. Existing entries in the "
            "sidecar are overwritten."
        ),
    )
    p.add_argument(
        "--bootstrap-fraction-warn", type=float, default=0.10,
        help="Bootstrap warn threshold = current - fraction*|current| (default 0.10).",
    )
    p.add_argument(
        "--bootstrap-fraction-fail", type=float, default=0.20,
        help="Bootstrap fail threshold = current - fraction*|current| (default 0.20).",
    )
    args = p.parse_args(argv)
    gh_root = Path(args.gh_pages_dir).resolve()
    if not gh_root.is_dir():
        print(f"error: {gh_root} is not a directory", file=sys.stderr)
        return 2
    manifest = load_manifest(args.manifest)
    if args.bootstrap_thresholds:
        out = bootstrap_thresholds(
            gh_root, manifest,
            fraction_warn=args.bootstrap_fraction_warn,
            fraction_fail=args.bootstrap_fraction_fail,
        )
        target = Path(args.thresholds)
        payload = {
            "generated_at":   datetime.now(tz=timezone.utc).isoformat(),
            "fraction_warn":  args.bootstrap_fraction_warn,
            "fraction_fail":  args.bootstrap_fraction_fail,
            "thresholds":     out,
        }
        target.parent.mkdir(parents=True, exist_ok=True)
        with open(target, "w", encoding="utf-8") as f:
            f.write(
                "# Per-analysis logLikelihood warn/fail thresholds.\n"
                "# Generated by buildGhPagesIndex.py --bootstrap-thresholds; safe to hand-edit\n"
                "# values to taste. ``current`` records the logL at bootstrap time as a hint\n"
                "# for tuning. The standardLogLikelihood aggregator classifies each analysis\n"
                "# individually and reports the worst status across all of a metric's\n"
                "# analyses, so any single failing analysis surfaces as a failing metric.\n\n"
            )
            yaml.safe_dump(
                payload, f, sort_keys=True, allow_unicode=True, default_flow_style=False,
            )
        n = sum(len(v) for v in out.values())
        print(f"wrote {n} analysis thresholds for {len(out)} metrics to {target}")
        return 0
    # Migrate legacy shared-dir layouts to per-variant directories.
    moved_id = migrate_legacy_shared_dir(
        gh_root, "dev/valid/idealizedSubhaloSimulations", _idealized_extractor,
    )
    moved_dd = migrate_legacy_shared_dir(
        gh_root, "dev/valid/darkMatterOnlySubhalosDecayingDarkMatter", _decaying_extractor,
    )
    if moved_id:
        print(f"migrated {len(moved_id)} idealized-subhalo variants")
    if moved_dd:
        print(f"migrated {len(moved_dd)} decaying-DM variants")
    # Build records and emit Jekyll site files.
    thresholds_map = load_thresholds(args.thresholds)
    records, groups = build_records(gh_root, manifest, thresholds_map=thresholds_map)
    write_if_changed(gh_root / "_config.yml", CONFIG_YML)
    write_if_changed(gh_root / "Gemfile", GEMFILE)
    emit_metrics_data(gh_root, records, groups)
    emit_index(gh_root, records, groups)
    emit_dashboard(gh_root, records, groups)
    emit_group_pages(gh_root, records, groups)
    n_ok   = sum(1 for r in records if r["status"] == "ok")
    n_warn = sum(1 for r in records if r["status"] == "warn")
    n_fail = sum(1 for r in records if r["status"] == "fail")
    n_unk  = sum(1 for r in records if r["status"] == "unknown")
    print(
        f"wrote site for {len(records)} metrics "
        f"(ok={n_ok}, warn={n_warn}, fail={n_fail}, unknown={n_unk})"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
