#!/usr/bin/env python3
import json
import os
import re
import subprocess
import sys
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone

# Count source lines of code in Galacticus source files, accounting for embedded
# XML and LaTeX, and report month-over-month deltas to Slack.
# Andrew Benson (04-January-2024)

HISTORY_PATH = ".github/metrics/slocHistory.json"


def count_fortran_xml_latex():
    """Walk source/ and doc/ to count Fortran/embedded-XML/embedded-LaTeX/LaTeX lines."""
    counts = {}
    source_files = []
    for root, _dirs, files in os.walk('source'):
        for f in files:
            if f.endswith('.F90') or f.endswith('.Inc'):
                source_files.append(os.path.join(root, f))
    doc_files = []
    for root, _dirs, files in os.walk('doc'):
        for f in files:
            if f.endswith('.tex'):
                doc_files.append(os.path.join(root, f))
    for file_name in source_files:
        in_xml   = False
        in_latex = False
        with open(file_name, 'r', errors='replace') as f:
            for line in f:
                if   re.match(r'^\s*!!\[\s*$', line):
                    in_xml   = True
                elif re.match(r'^\s*!!\]\s*$', line):
                    in_xml   = False
                elif re.match(r'^\s*!!\{\s*$', line):
                    in_latex = True
                elif re.match(r'^\s*!!\}\s*$', line):
                    in_latex = False
                if re.match(r'^\s*$', line) or re.match(r'^\s*&', line):
                    continue
                if re.match(r'^\s*!', line) and not re.match(r'^\s*!\$', line):
                    continue
                if   in_xml:
                    counts['xml']     = counts.get('xml',     0) + 1
                elif in_latex:
                    counts['latex']   = counts.get('latex',   0) + 1
                else:
                    counts['fortran'] = counts.get('fortran', 0) + 1
    for file_name in doc_files:
        with open(file_name, 'r', errors='replace') as f:
            for line in f:
                if re.match(r'^\s*$', line) or re.match(r'^\s*%', line):
                    continue
                counts['latex'] = counts.get('latex', 0) + 1
    return counts


def count_via_sloccount():
    """Run sloccount and return per-language line counts (excluding f90, handled above)."""
    counts = {}
    result = subprocess.run(
        'sloccount aux constraints parameters parameters.xml perl plots schema scripts testSuite source',
        shell=True, capture_output=True, text=True
    )
    if result.returncode != 0:
        print(f"Error: sloccount failed with exit code {result.returncode}.", file=sys.stderr)
        if result.stderr:
            print(result.stderr, file=sys.stderr)
        raise SystemExit(result.returncode)
    for line in result.stdout.splitlines():
        match = re.match(r'^([a-z0-9]+):\s*(\d+)', line)
        if match:
            language    = match.group(1)
            count_lines = int(match.group(2))
            if language != 'f90':
                counts[language] = counts.get(language, 0) + count_lines
    return counts


def compute_counts():
    counts = count_fortran_xml_latex()
    for language, n in count_via_sloccount().items():
        counts[language] = counts.get(language, 0) + n
    return counts


def load_history():
    if not os.path.exists(HISTORY_PATH):
        return {"schema_version": 1, "entries": []}
    with open(HISTORY_PATH, 'r') as f:
        return json.load(f)


def save_history(history):
    os.makedirs(os.path.dirname(HISTORY_PATH), exist_ok=True)
    with open(HISTORY_PATH, 'w') as f:
        json.dump(history, f, indent=2, sort_keys=True)
        f.write('\n')


def indicator(delta, is_new):
    if is_new:
        return ':new:'
    if delta > 0:
        return ':large_green_circle:'
    if delta < 0:
        return ':red_circle:'
    return ':white_circle:'


def format_report(current, previous_entry, sha, measured):
    previous_counts = previous_entry["counts"] if previous_entry else {}
    previous_date   = previous_entry["timestamp"][:10] if previous_entry else None

    lines = [f"*Galacticus SLOC report — {measured}* (commit `{sha}`)"]
    if previous_date:
        lines.append(f"Previous: {previous_date}")
    else:
        lines.append("Baseline measurement (no prior data).")
    lines.append("")

    sorted_langs = sorted(current.items(), key=lambda kv: kv[1], reverse=True)
    for language, value in sorted_langs:
        if previous_entry is None:
            lines.append(f":new: *{language}*: {value:,}")
            continue
        is_new   = language not in previous_counts
        previous = previous_counts.get(language, 0)
        delta    = value - previous
        if is_new:
            lines.append(f":new: *{language}*: {value:,} (new)")
        else:
            pct = (100.0 * delta / previous) if previous else 0.0
            lines.append(
                f"{indicator(delta, False)} *{language}*: {value:,} "
                f"({delta:+,}, {pct:+.2f}%)"
            )

    total = sum(current.values())
    lines.append("─" * 8)
    if previous_entry is None:
        lines.append(f":new: *total*: {total:,}")
    else:
        previous_total = sum(previous_counts.values())
        delta_total    = total - previous_total
        pct_total      = (100.0 * delta_total / previous_total) if previous_total else 0.0
        lines.append(
            f"{indicator(delta_total, False)} *total*: {total:,} "
            f"({delta_total:+,}, {pct_total:+.2f}%)"
        )
    return "\n".join(lines)


def build_payload(current, previous_entry, sha, measured):
    previous_date = previous_entry["timestamp"][:10] if previous_entry else ""
    return {
        "report":   format_report(current, previous_entry, sha, measured),
        "commit":   sha,
        "measured": measured,
        "previous": previous_date,
    }


def post_to_slack(payload, max_retries=4):
    slack_url = os.environ['SLACK_WEBHOOK_SLOCREPORT_URL']
    data      = json.dumps(payload).encode('utf-8')
    delay     = 1
    for attempt in range(1, max_retries + 1):
        req = urllib.request.Request(
            slack_url,
            data=data,
            headers={'Content-type': 'application/json'},
            method='POST',
        )
        try:
            with urllib.request.urlopen(req, timeout=10) as response:
                status = response.getcode()
                if 200 <= status < 300:
                    return
                transient = status == 429 or 500 <= status < 600
                print(f"Slack webhook returned HTTP {status} (attempt {attempt}).", file=sys.stderr)
                if not transient:
                    raise SystemExit(1)
        except urllib.error.HTTPError as e:
            transient = e.code == 429 or 500 <= e.code < 600
            print(f"Slack webhook HTTPError {e.code} (attempt {attempt}): {e}", file=sys.stderr)
            if not transient:
                raise SystemExit(1)
        except urllib.error.URLError as e:
            print(f"Slack webhook URLError (attempt {attempt}): {e}", file=sys.stderr)
        if attempt == max_retries:
            print("Error: failed to post SLOC report to Slack after retries.", file=sys.stderr)
            raise SystemExit(1)
        time.sleep(delay)
        delay *= 2


def main():
    current   = compute_counts()
    history   = load_history()
    previous  = history["entries"][-1] if history.get("entries") else None
    sha       = os.environ.get("GITHUB_SHA", "unknown")[:8]
    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    measured  = timestamp[:10]
    payload   = build_payload(current, previous, sha, measured)
    post_to_slack(payload)
    history.setdefault("schema_version", 1)
    history.setdefault("entries", [])
    history["entries"].append({"timestamp": timestamp, "sha": sha, "counts": current})
    save_history(history)


if __name__ == "__main__":
    main()
