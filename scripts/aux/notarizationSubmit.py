#!/usr/bin/env python3
import datetime
import json
import os
import subprocess
import sys

files = [
    "Galacticus_MacOS.zip",
    "Galacticus_MacOS-M1.zip",
    "toolsMacOS.zip",
    "toolsMacOSM1.zip"
]
for required in ("APPLE_ID", "APPLE_TEAM_ID", "APPLE_APP_SPECIFIC_PASSWORD"):
    if not os.environ.get(required):
        print(f"Error: missing required environment variable '{required}' for notarization submit.", file=sys.stderr)
        raise SystemExit(1)
requests = []
for file_name in files:
    if not os.path.exists(file_name):
        print(f"Error: notarization input '{file_name}' was not found.", file=sys.stderr)
        raise SystemExit(1)
    result = subprocess.run(
        [
            "xcrun", "notarytool", "submit", file_name,
            "--apple-id", os.environ["APPLE_ID"],
            "--password", os.environ["APPLE_APP_SPECIFIC_PASSWORD"],
            "--team-id", os.environ["APPLE_TEAM_ID"],
            "--output-format", "json"
        ],
        capture_output=True,
        text=True
    )
    if result.returncode != 0:
        print(result.stdout)
        print(result.stderr, file=sys.stderr)
        raise SystemExit(result.returncode)
    data = json.loads(result.stdout)
    requests.append(
        {
            "asset": file_name,
            "submissionId": data.get("id"),
            "status": data.get("status"),
        }
    )
output = {
    "repository": os.environ["GITHUB_REPOSITORY"],
    "commitSha": os.environ["GITHUB_SHA"],
    "submittedAtUtc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    "workflowRunUrl": (
        f"{os.environ['GITHUB_SERVER_URL']}/"
        f"{os.environ['GITHUB_REPOSITORY']}/actions/runs/{os.environ['GITHUB_RUN_ID']}"
    ),
    "requests": requests,
}
output_file = f"notarization-requests-{os.environ['GITHUB_SHA']}.json"
with open(output_file, "w", encoding="utf-8") as handle:
    json.dump(output, handle, indent=2)
    handle.write("\n")
print(f"Wrote {output_file}")
