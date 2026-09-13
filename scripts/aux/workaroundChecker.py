#!/usr/bin/env python3
import html
import json
import os
import re
import subprocess
import time

# Check for closed GCC PRs in workarounds.
# Andrew Benson (17-September-2021)

WORKAROUNDS = {}
DIRECTORIES = ("./source",)
# Match the opening tag of a `<workaround …>` directive, capturing its
# attributes.  The tag may be preceded by a `#` (workarounds inside
# preprocessor-conditional code).
PATTERN = re.compile(r'^\s*#?\s*<workaround\s([^>]*)')
# Attribute extraction from the captured attribute string.
ATTRIBUTE = r'\b{name}="([^"]*)"'
# Repository against which bare `issue="<number>"` attributes are resolved.
REPOSITORY = os.environ.get("GITHUB_REPOSITORY", "galacticusorg/galacticus")


def attribute(attributes: str, name: str) -> str | None:
    """Return the value of the named attribute, or `None` if it is absent.

    Attribute values in directives are XML-escaped (URLs, for example, have
    their `/` written as `&#x2F;`), so unescape before returning.
    """
    match = re.search(ATTRIBUTE.format(name=name), attributes)
    return html.unescape(match.group(1)) if match else None


def issue_url(issue: str) -> str:
    """Normalize an `issue` attribute (a bare number, or a full URL) to a URL."""
    if re.fullmatch(r'#?\d+', issue):
        return f"https://github.com/{REPOSITORY}/issues/{issue.lstrip('#')}"
    return issue


def file_matcher() -> None:
    for directory in DIRECTORIES:
        for root, _, files in os.walk(directory):
            for file_name in files:
                if not file_name.endswith((".F90", ".Inc")):
                    continue
                full_name = os.path.join(root, file_name)
                print(f"Scanning {full_name}")
                try:
                    with open(full_name, "r", encoding="utf-8", errors="replace") as file_handle:
                        for line in file_handle:
                            match = PATTERN.search(line)
                            if not match:
                                continue
                            attributes = match.group(1)
                            pr = attribute(attributes, "PR")
                            if pr is None:
                                continue
                            if pr not in WORKAROUNDS:
                                WORKAROUNDS[pr] = {"files": set(), "issues": set(), "status": "UNKNOWN"}
                            WORKAROUNDS[pr]["files"].add(full_name)
                            # An `issue` attribute records that removal of this
                            # workaround has already been staged as an issue in
                            # our own repository — see `check_issues()`.
                            issue = attribute(attributes, "issue")
                            if issue:
                                WORKAROUNDS[pr]["issues"].add(issue_url(issue))
                except OSError:
                    continue


def check_links() -> None:
    for pr in WORKAROUNDS:
        WORKAROUNDS[pr]["status"] = "UNKNOWN"
        url = f"https://gcc.gnu.org/bugzilla/show_bug.cgi?id={pr}"
        time.sleep(1)
        result = subprocess.run(
            ["curl", "--silent", "--location", "--fail", url],
            capture_output=True,
            text=True,
            check=False,
        )
        if result.returncode != 0:
            continue
        for line in result.stdout.splitlines():
            match = re.search(r'<span id="static_bug_status">([A-Z]+)', line)
            if match:
                WORKAROUNDS[pr]["status"] = match.group(1)
                break


def issue_state(url: str) -> str:
    """Return the state ("open"/"closed") of one of our issues, or "UNKNOWN".

    A workaround is only "staged" while its issue is *open*: if the issue has
    been closed but the workaround is still present, the workaround has been
    forgotten and should be reported again.  If the state can not be
    determined (no network, API rate limit, a non-GitHub URL, …) we return
    "UNKNOWN" and give the workaround the benefit of the doubt.
    """
    match = re.match(r'https://github\.com/([^/]+/[^/]+)/issues/(\d+)', url)
    if not match:
        return "UNKNOWN"
    api = f"https://api.github.com/repos/{match.group(1)}/issues/{match.group(2)}"
    command = ["curl", "--silent", "--location", "--fail",
               "--header", "Accept: application/vnd.github+json"]
    token = os.environ.get("GH_TOKEN") or os.environ.get("GITHUB_TOKEN")
    if token:
        command += ["--header", f"Authorization: Bearer {token}"]
    result = subprocess.run(command + [api], capture_output=True, text=True, check=False)
    if result.returncode != 0:
        return "UNKNOWN"
    try:
        return json.loads(result.stdout).get("state", "UNKNOWN")
    except json.JSONDecodeError:
        return "UNKNOWN"


def check_issues() -> None:
    """Determine which resolved PRs have their removal staged as an open issue."""
    for pr in WORKAROUNDS:
        states = {url: issue_state(url) for url in sorted(WORKAROUNDS[pr]["issues"])}
        WORKAROUNDS[pr]["issueStates"] = states
        # "Staged" means at least one referenced issue is not known to be
        # closed.  The same PR is often worked around in many files, but a
        # single issue tracks removal of all of them, so we do not require
        # every occurrence to carry the attribute.
        WORKAROUNDS[pr]["staged"] = any(state != "closed" for state in states.values())


def main() -> int:
    file_matcher()
    check_links()
    check_issues()

    resolved = [pr for pr in WORKAROUNDS
                if WORKAROUNDS[pr]["status"] == "RESOLVED" and not WORKAROUNDS[pr]["staged"]]
    staged = [pr for pr in WORKAROUNDS
              if WORKAROUNDS[pr]["status"] == "RESOLVED" and WORKAROUNDS[pr]["staged"]]
    if resolved:
        status = 1
        print("!!! Resolved PRs with workarounds exist !!!\n")
    else:
        status = 0
        print("No resolved PRs with unstaged workarounds exist\n")
    if staged:
        print(f"({len(staged)} resolved PR(s) have their workaround removal staged as an open issue)\n")

    for pr in sorted(WORKAROUNDS, key=int):
        workaround = WORKAROUNDS[pr]
        if workaround["status"] == "RESOLVED":
            print("!!! " if not workaround["staged"] else "--- ", end="")
        print(f"PR{pr} (https://gcc.gnu.org/bugzilla/show_bug.cgi?id={pr}):")
        for file_name in sorted(workaround["files"]):
            print(f" -> {file_name}")
        for url, state in workaround["issueStates"].items():
            label = "staged for removal" if state != "closed" else "issue is CLOSED but workaround remains"
            print(f" => {label}: {url}")
        print("\n")

    return status


if __name__ == "__main__":
    raise SystemExit(main())
