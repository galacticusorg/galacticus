#!/usr/bin/env python3
import os
import re
import subprocess
import time

# Check for closed GCC PRs in workarounds.
# Andrew Benson (17-September-2021) [Perl]; ported to Python 2026

WORKAROUNDS = {}
DIRECTORIES = ("./source", "./perl/Galacticus/Build")
PATTERN = re.compile(r'^\s*#?\s*<workaround\s.*PR="(\d+)"')


def file_matcher() -> None:
    for directory in DIRECTORIES:
        for root, _, files in os.walk(directory):
            for file_name in files:
                if not file_name.endswith((".F90", ".Inc", ".pm")):
                    continue
                full_name = os.path.join(root, file_name)
                print(f"Scanning {full_name}")
                try:
                    with open(full_name, "r", encoding="utf-8", errors="replace") as file_handle:
                        for line in file_handle:
                            match = PATTERN.search(line)
                            if not match:
                                continue
                            pr = match.group(1)
                            if pr not in WORKAROUNDS:
                                WORKAROUNDS[pr] = {"files": set(), "status": "UNKNOWN"}
                            WORKAROUNDS[pr]["files"].add(full_name)
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


def main() -> int:
    file_matcher()
    check_links()

    resolved = sum(1 for pr in WORKAROUNDS if WORKAROUNDS[pr]["status"] == "RESOLVED")
    if resolved:
        status = 1
        print("!!! Resolved PRs with workarounds exist !!!\n")
    else:
        status = 0
        print("No resolved PRs with workarounds exist\n")

    for pr in sorted(WORKAROUNDS, key=int):
        if WORKAROUNDS[pr]["status"] == "RESOLVED":
            print("!!! ", end="")
        print(f"PR{pr} (https://gcc.gnu.org/bugzilla/show_bug.cgi?id={pr}):")
        for file_name in sorted(WORKAROUNDS[pr]["files"]):
            print(f" -> {file_name}")
        print("\n")

    return status


if __name__ == "__main__":
    raise SystemExit(main())
