#!/usr/bin/env python3
import json
import os
import re
import subprocess
import sys
import urllib.request
import urllib.error

# Retrieve and archive copies of all run-time downloaded data files for Galacticus.
# Andrew Benson (19-April-2024)

if len(sys.argv) != 4:
    print("Usage: archive.py <galacticusPath> <archivePath> <slackToken>", file=sys.stderr)
    raise SystemExit(1)

galacticus_path = sys.argv[1]
archive_path    = sys.argv[2]
slack_token     = sys.argv[3]

# Parse the dependencies file.
dependencies = {}
with open(os.path.join(galacticus_path, "aux", "dependencies.yml")) as f:
    for line in f:
        m = re.match(r'^(.*):\s+([\d\.]+)', line)
        if m:
            code    = m.group(1)
            version = m.group(2)
            parts   = version.split('.')
            dependencies[code] = {
                "version":      version,
                "versionMajor": parts[0],
            }

# Scan the source directory for files.
links = []

def link_finder(file_name, full_path):
    # Ignore files not matching expected extensions or Makefile.
    if not (re.search(r'\.(F90|Inc|py)$', file_name) or os.path.basename(file_name) == 'Makefile'):
        return
    try:
        with open(full_path, 'r', errors='replace') as f:
            for line in f:
                # Fortran source.
                if re.search(r'\.(F90|Inc)$', file_name):
                    m = re.match(r'^\s*call\s+download\s*\(\s*(["\'][^,]+)', line)
                    if m:
                        link = re.sub(r'["\']', '', m.group(1))
                        # Replace dependency version placeholders.
                        def replace_major(match):
                            name = match.group(1)
                            return dependencies.get(name, {}).get("versionMajor", match.group(0))
                        link = re.sub(r'//char\(([a-zA-Z]+)VersionMajor\)//', replace_major, link)
                        cloudy_version = dependencies.get('cloudy', {}).get('version', '')
                        link = re.sub(r'//char\(cloudyVersion\)//', 'c' + cloudy_version, link)
                        def replace_version(match):
                            name = match.group(1)
                            return dependencies.get(name, {}).get("version", match.group(0))
                        link = re.sub(r'//char\(([a-zA-Z]+)Version\)//', replace_version, link)
                        # Skip the "backup" ("old") Cloudy path.
                        if not re.search(r'cloudy_releases/c\d+/old/', link):
                            links.append(link)
                # Python scripts.
                if re.search(r'\.py$', file_name):
                    m = re.match(r'^\s*urllib\.request\.urlretrieve\s*\(\s*(["\'][^,]+)', line)
                    if m:
                        link = re.sub(r'["\']', '', m.group(1))
                        links.append(link)
                # Makefiles.
                if os.path.basename(file_name) == 'Makefile':
                    m = re.match(r'^\s*wget\s+(?:--\S+\s+)*(\S+)', line)
                    if m:
                        links.append(m.group(1))
    except (IOError, OSError):
        pass

# Walk source and scripts directories.
for directory in [os.path.join(galacticus_path, "source"), os.path.join(galacticus_path, "scripts")]:
    for root, dirs, files in os.walk(directory):
        for fname in files:
            link_finder(fname, os.path.join(root, fname))

# Process the top-level Makefile.
makefile_path = os.path.join(galacticus_path, "Makefile")
link_finder("Makefile", makefile_path)

# Retrieve links.
report = {"report": ""}
for link in links:
    m = re.match(r'^https?://(.+)/(.+)', link)
    if m:
        path      = m.group(1)
        file_name = m.group(2)
        dest_dir  = os.path.join(archive_path, path)
        os.makedirs(dest_dir, exist_ok=True)
        dest = os.path.join(dest_dir, file_name)
        if not os.path.exists(dest):
            report["report"] += f"RETRIEVING: {link}\n"
            result = subprocess.run(["wget", link, "-O", dest])
            if result.returncode != 0:
                report["report"] += f"\tFAILED: {link}\n"
        else:
            report["report"] += f"SKIPPING: (already archived) {link}\n"

# Report the results.
payload = json.dumps(report).encode("utf-8")
req = urllib.request.Request(
    f"https://hooks.slack.com/triggers/{slack_token}",
    data=payload,
    headers={"Content-type": "application/json"},
    method="POST"
)
try:
    with urllib.request.urlopen(req, timeout=30) as response:
        pass
except urllib.error.URLError as e:
    print(f"Error: failed to post report to Slack: {e}", file=sys.stderr)
    raise SystemExit(1)
