#!/usr/bin/env python3
import json
import os
import sys
import urllib.error
import urllib.request

# Report macOS notarization workflow status to Slack.

status_name = os.environ.get("NOTARIZATION_STATUS", "").strip().lower()
if status_name == "success":
    status = ":white_check_mark:"
elif status_name == "failure":
    status = ":x:"
else:
    status = ":clock2:"

repository = os.environ.get("NOTARIZATION_REPOSITORY", "galacticusorg/galacticus")
repo_name = repository.split("/", 1)[-1]
details = os.environ.get("NOTARIZATION_DETAILS", "")
url = os.environ.get("NOTARIZATION_URL", "")

payload = json.dumps(
    {
        "repo": repo_name,
        "workflow": "Notarization (macOS)",
        "status": status,
        "details": details,
        "url": url,
    }
).encode("utf-8")
req = urllib.request.Request(
    os.environ["SLACK_WEBHOOK_STATUS_URL"],
    data=payload,
    headers={"Content-type": "application/json"},
    method="POST",
)
try:
    with urllib.request.urlopen(req, timeout=10):
        pass
except urllib.error.URLError as error:
    print(f"Error: failed to post to Slack: {error}", file=sys.stderr)
    raise SystemExit(1)
