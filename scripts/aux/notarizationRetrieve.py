#!/usr/bin/env python3
import json
import os
import subprocess

master_sha = subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
metadata_file = f"notarization-requests-{master_sha}.json"

download = subprocess.run(
    [
        "gh", "release", "download", "bleeding-edge",
        "--pattern", metadata_file,
        "--dir", "."
    ],
    capture_output=True,
    text=True
)
if download.returncode != 0:
    details = f"No notarization request metadata found for master SHA {master_sha}."
    result = {"status": "pending", "details": details}
    with open("notarization-poll-result.json", "w", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2)
        handle.write("\n")
    with open(os.environ["GITHUB_STEP_SUMMARY"], "a", encoding="utf-8") as summary:
        summary.write("## Notarization (macOS)\n")
        summary.write(f"- Status: pending\n- Details: {details}\n")
    raise SystemExit(0)

with open(metadata_file, "r", encoding="utf-8") as handle:
    metadata = json.load(handle)

requests = metadata.get("requests", [])
statuses = []
failure_details = []
for request in requests:
    submission_id = request.get("submissionId")
    asset = request.get("asset", "unknown-asset")
    info = subprocess.run(
        [
            "xcrun", "notarytool", "info", submission_id,
            "--apple-id", os.environ["APPLE_ID"],
            "--password", os.environ["APPLE_APP_SPECIFIC_PASSWORD"],
            "--team-id", os.environ["APPLE_TEAM_ID"],
            "--output-format", "json"
        ],
        capture_output=True,
        text=True
    )
    if info.returncode != 0:
        failure_details.append(f"{asset}: unable to query status for {submission_id}")
        statuses.append("failure")
        continue
    info_data = json.loads(info.stdout)
    status = str(info_data.get("status", "unknown")).lower()
    if status == "accepted":
        statuses.append("success")
    elif status in {"in progress", "submitted"}:
        statuses.append("pending")
    else:
        statuses.append("failure")
        failure_details.append(f"{asset}: {info_data.get('status', 'unknown')}")
        log_result = subprocess.run(
            [
                "xcrun", "notarytool", "log", submission_id,
                "--apple-id", os.environ["APPLE_ID"],
                "--password", os.environ["APPLE_APP_SPECIFIC_PASSWORD"],
                "--team-id", os.environ["APPLE_TEAM_ID"],
                "--output-format", "json"
            ],
            capture_output=True,
            text=True
        )
        if log_result.returncode == 0:
            log_data = json.loads(log_result.stdout)
            issues = log_data.get("issues", [])
            issue_messages = [issue.get("message") for issue in issues if issue.get("message")]
            for message in issue_messages:
                failure_details.append(f"{asset}: {message}")

if "failure" in statuses:
    overall = "failure"
    details = "; ".join(failure_details) if failure_details else "Notarization rejected."
elif "pending" in statuses:
    overall = "pending"
    details = "At least one notarization request is still in progress."
else:
    overall = "success"
    details = "All notarization requests were accepted."

result = {"status": overall, "details": details}
with open("notarization-poll-result.json", "w", encoding="utf-8") as handle:
    json.dump(result, handle, indent=2)
    handle.write("\n")

with open(os.environ["GITHUB_STEP_SUMMARY"], "a", encoding="utf-8") as summary:
    summary.write("## Notarization (macOS)\n")
    summary.write(f"- Commit: {metadata.get('commitSha', master_sha)}\n")
    summary.write(f"- Status: {overall}\n")
    summary.write(f"- Details: {details}\n")
