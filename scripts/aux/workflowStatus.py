#!/usr/bin/env python3
import json
import os
import subprocess
import sys
import urllib.request
import urllib.error

# Report status of GitHub Actions workflows to Slack.
# Andrew Benson

repo = "galacticus"

workflows = [
    {"file": "cicd.yml",      "name": "CI/CD"},
    {"file": "linkCheck.yml", "name": "Link-Check"},
]

# Get the current HEAD SHA once.
head_result = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True)
if head_result.returncode != 0:
    print(f"Error: git rev-parse HEAD failed: {head_result.stderr}", file=sys.stderr)
    raise SystemExit(head_result.returncode)
head_sha = head_result.stdout.strip()

for workflow in workflows:
    result = subprocess.run(
        ["gh", "run", "list", "--workflow", workflow["file"],
         "--branch", "master", "--json", "conclusion,headSha"],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        print(f"Error: gh run list failed: {result.stderr}", file=sys.stderr)
        raise SystemExit(result.returncode)
    data   = json.loads(result.stdout)
    status = ":question:"
    for run in data:
        # Exclude commits that are not ancestors of the current HEAD (e.g. PRs from forks).
        check = subprocess.run(
            ["git", "merge-base", "--is-ancestor", run["headSha"], head_sha],
            capture_output=True
        )
        if check.returncode != 0:
            continue
        conclusion = run.get("conclusion")
        if conclusion in (None, ""):
            status = ":clock2:"
        elif conclusion == "failure":
            status = ":x:"
        elif conclusion == "success":
            status = ":white_check_mark:"
        if status != ":question:":
            break
    payload = json.dumps({
        "repo":     repo,
        "workflow": workflow["name"],
        "status":   status,
        "url":      f"https://github.com/galacticusorg/{repo}/actions/workflows/{workflow['file']}",
    }).encode("utf-8")
    req = urllib.request.Request(
        os.environ["SLACK_WEBHOOK_STATUS_URL"],
        data=payload,
        headers={"Content-type": "application/json"},
        method="POST"
    )
    try:
        with urllib.request.urlopen(req, timeout=10) as response:
            pass
    except urllib.error.URLError as e:
        print(f"Error: failed to post to Slack: {e}", file=sys.stderr)
        raise SystemExit(1)
