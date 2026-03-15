#!/usr/bin/env python3
import os
import re
import json
import subprocess
import sys
import urllib.request
import urllib.error

# Count source lines of code in Galacticus source files, accounting for embedded XML and LaTeX.
# Andrew Benson (04-January-2024)

count = {}

# Get a list of all source files to process.
source_files = []
for root, dirs, files in os.walk('source'):
    for f in files:
        if f.endswith('.F90') or f.endswith('.Inc'):
            source_files.append(os.path.join(root, f))

doc_files = []
for root, dirs, files in os.walk('doc'):
    for f in files:
        if f.endswith('.tex'):
            doc_files.append(os.path.join(root, f))

# Iterate over source files.
for file_name in source_files:
    in_xml   = False
    in_latex = False
    with open(file_name, 'r', errors='replace') as f:
        for line in f:
            # Detect entry and exit of embedded regions.
            if   re.match(r'^\s*!!\[\s*$', line):
                in_xml   = True
            elif re.match(r'^\s*!!\]\s*$', line):
                in_xml   = False
            elif re.match(r'^\s*!!\{\s*$', line):
                in_latex = True
            elif re.match(r'^\s*!!\}\s*$', line):
                in_latex = False
            # Skip empty lines, continuations, and comments (unless they are an OpenMP directive).
            if re.match(r'^\s*$', line) or re.match(r'^\s*&', line):
                continue
            if re.match(r'^\s*!', line) and not re.match(r'^\s*!\$', line):
                continue
            # Count the line.
            if   in_xml:
                count['xml']     = count.get('xml',     0) + 1
            elif in_latex:
                count['latex']   = count.get('latex',   0) + 1
            else:
                count['fortran'] = count.get('fortran', 0) + 1

# Iterate over documentation files.
for file_name in doc_files:
    with open(file_name, 'r', errors='replace') as f:
        for line in f:
            # Skip empty lines and comments.
            if re.match(r'^\s*$', line) or re.match(r'^\s*%', line):
                continue
            count['latex'] = count.get('latex', 0) + 1

# Count everything else using sloccount.
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
            count[language] = count.get(language, 0) + count_lines

# Report the results.
slack_url = os.environ['SLACK_WEBHOOK_SLOCREPORT_URL']
data = json.dumps(count).encode('utf-8')
req = urllib.request.Request(
    slack_url,
    data=data,
    headers={'Content-type': 'application/json'},
    method='POST'
)
try:
    with urllib.request.urlopen(req, timeout=10) as response:
        status_code = response.getcode()
        if status_code < 200 or status_code >= 300:
            print(f"Error: Slack webhook returned HTTP status {status_code}.", file=sys.stderr)
            raise SystemExit(1)
except urllib.error.URLError as e:
    print(f"Error: failed to post SLOC report to Slack: {e}", file=sys.stderr)
    raise SystemExit(1)
