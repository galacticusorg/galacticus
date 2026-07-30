#!/usr/bin/env python3
import json
import os
import re
import subprocess
import sys
import urllib.request
import urllib.error
import urllib.parse

# Retrieve and archive copies of all run-time downloaded data files for Galacticus.
# Andrew Benson (19-April-2024)

if len(sys.argv) != 4:
    print("Usage: archive.py <galacticusPath> <archivePath> <slackToken>", file=sys.stderr)
    raise SystemExit(1)

galacticus_path = sys.argv[1]
archive_path    = sys.argv[2]
slack_token     = sys.argv[3]

# Parse the dependencies file. A version may be a dotted release number, or - for dependencies whose upstream publishes no
# releases - a git commit SHA. Names are stored folded to lower case, since the corresponding Fortran variable names (e.g.
# `cosmicEmuVersion`) are camel-cased.
dependencies = {}
with open(os.path.join(galacticus_path, "aux", "dependencies.yml")) as f:
    for line in f:
        m = re.match(r'^([^:]+):\s+([0-9A-Za-z._-]+)\s*$', line)
        if m:
            code    = m.group(1).lower()
            version = m.group(2)
            parts   = version.split('.')
            dependencies[code] = {
                "version":      version,
                "versionMajor": parts[0],
            }

# Scan the source directory for files. Links are collected as a list of fallback groups. Each group is a list of alternate
# URLs for the same content, any one of which is sufficient - mirroring `downloadMultiple()` in `source/system/downloader.F90`,
# which tries each URL of an array in turn and stops at the first success.
links = []

def resolve_versions(link):
    """Substitute dependency versions into the placeholders of a link extracted from Fortran source.

    A placeholder is the name of the local variable holding the version, embedded in the URL by string concatenation. It
    appears either wrapped in `char()` (where the URL is assembled from `character` variables) or bare (where it is assembled
    from `varying_string` variables) - normalize to the bare form first so that a single set of patterns handles both.
    """
    link = re.sub(r'char\s*\(\s*([a-zA-Z]+Version(?:Major)?)\s*\)', r'\1', link)
    def replace_major(match):
        name = match.group(1).lower()
        return dependencies.get(name, {}).get("versionMajor", match.group(0))
    link = re.sub(r'//([a-zA-Z]+)VersionMajor//', replace_major, link)
    # Cloudy tarballs are named for the version prefixed by "c".
    cloudy_version = dependencies.get('cloudy', {}).get('version', '')
    link = re.sub(r'//cloudyVersion//', 'c' + cloudy_version, link)
    def replace_version(match):
        name = match.group(1).lower()
        return dependencies.get(name, {}).get("version", match.group(0))
    link = re.sub(r'//([a-zA-Z]+)Version//', replace_version, link)
    return link

# Top-level domains and second-level names reserved by RFC 2606 and RFC 6761 for documentation and testing. These never
# resolve to a real host, so a link using one is a test fixture - e.g. `source/tests/system/download.F90`, which exercises the
# downloader's certificate handling and its defenses against shell metacharacters in URLs - and not a run-time dependency to
# be archived.
reservedTLDs   = ('test', 'example', 'invalid', 'localhost')
reservedHosts  = ('localhost', 'example.com', 'example.net', 'example.org')

def is_reserved(link):
    """Return whether the given link points at a hostname reserved for documentation and testing."""
    try:
        host = urllib.parse.urlsplit(link).hostname or ""
    except ValueError:
        # A malformed URL is not something we can archive in any case.
        return True
    return host in reservedHosts or host.rsplit('.', 1)[-1] in reservedTLDs

def link_finder(file_name, full_path):
    # Ignore files not matching expected extensions or Makefile.
    if not (re.search(r'\.(F90|Inc|py)$', file_name) or os.path.basename(file_name) == 'Makefile'):
        return
    # Fallback groups accumulated within this file, keyed on the name of the array variable holding the alternate URLs.
    groups = {}
    def add_link(link, variable=None):
        # Ignore links to reserved hostnames, which are test fixtures rather than run-time dependencies.
        if is_reserved(link):
            return
        # An indexed assignment (e.g. `urls(1)=...`, `urls(2)=...`) builds an array of alternate URLs, so all elements of a
        # given array in a given file belong to a single fallback group. Anything else stands alone.
        if variable is not None and variable in groups:
            groups[variable].append(link)
        else:
            group = [link]
            links.append(group)
            if variable is not None:
                groups[variable] = group
    try:
        with open(full_path, 'r', errors='replace') as f:
            for line in f:
                # Fortran source. A URL appears either directly as the first argument of a `download()` call, or is built up
                # in a variable which is then passed to `download()` - match both forms.
                if re.search(r'\.(F90|Inc)$', file_name):
                    variable = None
                    m = re.match(r'^\s*call\s+download\s*\(\s*(["\'][^,]+)', line)
                    if m:
                        link = re.sub(r'["\']', '', m.group(1))
                    else:
                        m = re.match(r'^\s*([a-zA-Z_]\w*)\s*(\([^)]*\))?\s*=\s*(?:var_str\s*\(\s*)?(["\']https?://.*)$', line)
                        if m:
                            # Trim any trailing `var_str()` closing parenthesis before removing quotes.
                            link     = re.sub(r'["\']', '', re.sub(r'[\s)]*$', '', m.group(3)))
                            variable = m.group(1) if m.group(2) else None
                    if m:
                        link = resolve_versions(link)
                        # Skip the "backup" ("old") Cloudy path.
                        if not re.search(r'cloudy_releases/c\d+/old/', link):
                            add_link(link, variable)
                # Python scripts.
                if re.search(r'\.py$', file_name):
                    m = re.match(r'^\s*urllib\.request\.urlretrieve\s*\(\s*(["\'][^,]+)', line)
                    if m:
                        add_link(re.sub(r'["\']', '', m.group(1)))
                # Makefiles.
                if os.path.basename(file_name) == 'Makefile':
                    m = re.match(r'^\s*wget\s+(?:--\S+\s+)*(\S+)', line)
                    if m:
                        add_link(m.group(1))
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

def destination(link):
    """Return the path under the archive at which the given link should be stored, or `None` if it is not a retrievable URL."""
    m = re.match(r'^https?://(.+)/(.+)', link)
    if not m:
        return None
    return os.path.join(archive_path, m.group(1), m.group(2))

# Retrieve links. Each group holds alternate URLs for the same content, so it is satisfied as soon as any one of them is
# archived - the remaining alternates are then not retrieved.
report     = {"report": ""}
unresolved = []
failed     = []
for group in links:
    # Guard against unsubstituted placeholders (e.g. a `<name>Version` for which no matching entry was found in the
    # dependencies file), which are left as a bare identifier flanked by Fortran concatenation operators. Attempting to
    # retrieve such a link would silently fail against a non-existent URL.
    for link in group:
        if re.search(r'//[a-zA-Z_]\w*//', link) or re.search(r'char\s*\(', link):
            unresolved.append(link)
            report["report"] += f"UNRESOLVED: {link}\n"
    group = [link for link in group if link not in unresolved and destination(link) is not None]
    if not group:
        continue
    # If any alternate is already archived, the content is held and there is nothing to do.
    archived = [link for link in group if os.path.exists(destination(link))]
    if archived:
        report["report"] += f"SKIPPING: (already archived) {archived[0]}\n"
        continue
    # Try each alternate in turn, stopping at the first success.
    for index, link in enumerate(group):
        dest = destination(link)
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        report["report"] += f"RETRIEVING: {link}\n"
        # `--no-check-certificate` matches the behavior of the run-time downloader (`downloadMultiple()` in
        # `source/system/downloader.F90`), so that the archiver can retrieve exactly what a model run would.
        result = subprocess.run(["wget", "--no-check-certificate", link, "-O", dest])
        if result.returncode == 0:
            break
        # Do not leave a partial file behind, which would otherwise count as "already archived" on the next run.
        if os.path.exists(dest):
            os.remove(dest)
        remaining = len(group) - index - 1
        if remaining > 0:
            report["report"] += f"\tFAILED: {link} (trying fallback URL)\n"
        else:
            report["report"] += f"\tFAILED: {link}\n"
            failed.append(link)

# Report the results. A group counts as failed only if every alternate URL in it failed.
if unresolved or failed:
    report["report"] += f"SUMMARY: {len(unresolved)} unresolved link(s), {len(failed)} unretrievable link(s)\n"
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

# Fail if any link could not be fully resolved.
if unresolved:
    for link in unresolved:
        print(f"Error: unresolved placeholder in link: {link}", file=sys.stderr)
    raise SystemExit(1)
