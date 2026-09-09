#!/usr/bin/env python3
# Check for broken links in Galacticus documentation.
# Andrew Benson (21-September-2020)

import sys
import os
import re
import html
import subprocess
import json
import random
import time
import xml.etree.ElementTree as ET
import requests


def load_failures(filename):
    """Load consecutive-failure records from XML file."""
    failures = {}
    if os.path.exists(filename):
        tree = ET.parse(filename)
        root = tree.getroot()
        for url_elem in root.findall('url'):
            key = url_elem.get('key')
            consec = int(url_elem.get('consecutiveFailures', 0))
            if key:
                failures[key] = {'consecutiveFailures': consec}
    return failures


def save_failures(failures, filename):
    """Save consecutive-failure records to XML file."""
    root = ET.Element('failures')
    for url in sorted(failures):
        url_elem = ET.SubElement(root, 'url')
        url_elem.set('key', url)
        url_elem.set('consecutiveFailures', str(failures[url]['consecutiveFailures']))
    tree = ET.ElementTree(root)
    ET.indent(tree, space='  ')
    with open(filename, 'w', encoding='utf-8') as f:
        f.write('<?xml version="1.0" encoding="utf-8"?>\n')
        tree.write(f, encoding='unicode')
        f.write('\n')


def record_failure(url, failures):
    """Increment consecutive-failure count. Returns True once count reaches 3."""
    if url not in failures:
        failures[url] = {'consecutiveFailures': 0}
    failures[url]['consecutiveFailures'] += 1
    return failures[url]['consecutiveFailures'] >= 3


def record_success(url, failures):
    """Reset consecutive-failure count to zero."""
    if url not in failures:
        failures[url] = {'consecutiveFailures': 0}
    failures[url]['consecutiveFailures'] = 0


# Match a URL up to the first whitespace or delimiter.  This naturally
# terminates the forms the URLs appear in across the docs: RST hyperlinks
# (`` `text <url>`_ ``), XML directive attributes (``referenceURL="url"``,
# ``externalDescription="url"``, ``url="url"``), LaTeX ``\href{url}``, and bare
# URLs in prose.  ``)`` is allowed (Wikipedia-style URLs); ``'`` is allowed
# (apostrophes appear inside real URLs, e.g.
# ``.../Claudia's_Stellar_Population_Model.html``).  Trailing sentence
# punctuation — including a ``'`` used as a closing delimiter around the URL —
# is stripped afterwards.  Backticks are excluded since they delimit RST inline
# literals (`` ``url`` ``) and are never valid inside a URL.
_URL_RE = re.compile(r'https?://[^\s<>"}\]`]+')

# A delimiter that terminates a URL, but which may be hidden behind an entity
# reference in the source -- an RST link target ``<url>`` embedded in an XML
# directive writes its closing ``>`` as ``&gt;``, which the regex above cannot
# see and which unescaping then leaves stuck to the end of the URL.
_URL_TERMINATOR_RE = re.compile(r'[\s<>"`]')

# Fortran string concatenation immediately inside, or immediately following, a
# matched URL: ``'https://…/archive/v'//version//'.tar.gz'``.  The literal is
# only a *fragment* of a URL which is assembled at run time, so there is no
# complete URL here to check -- and the fragment, checked as though it were
# one, is reported as broken.  Such URLs are skipped: missing a link is much
# cheaper than a standing false alarm.
_CONCATENATION_RE       = re.compile(r'[\'"]\s*//')
_CONCATENATION_AFTER_RE = re.compile(r'^[\'"]?\s*//')

# URLs matching any of these patterns are placeholders / illustrative
# examples (not real links) and are skipped during checking.
_EXCLUDED_URL_RES = [
    re.compile(r'^https?://\.+\)?$'),  # placeholder, e.g. "http://......." in code examples
    # Defunct Google Drive link for BPASS data, retained only for the historical record.
    re.compile(r'^https://drive\.google\.com/open\?id=0B7vqPPPgOdtIfjUtb3RsV2JUOTFFX29WV1FZNURPMHAxTEtZQjhJOGtyNXZUTTNVSzFZazQ$'),
    # Defunct Maraston stellar population model link, retained only for the historical record.
    re.compile(r"^http://www\.icg\.port\.ac\.uk/~maraston/Claudia's_Stellar_Population_Model\.html$"),
    # A documentation URL written with a placeholder version, naming the form of
    # a per-release documentation URL rather than any one release's -- see the
    # "Versions and Releases" page of the developer guide.
    re.compile(r'^https://galacticus\.readthedocs\.io/en/vX\.Y\.Z(?:[/?#]|$)'),
    # Defunct link to Inoue's own IGM attenuation code, from which the file that
    # `source/tests/spectra/postprocess/Inoue2014.F90` compares against was
    # extracted.  No archived copy of the tarball could be found, so the URL is
    # retained only to record where that file came from.
    re.compile(r'^http://www\.las\.osaka-sandai\.ac\.jp/~inoue/ANAIGM/ANAIGM\.tar\.gz$'),
    # RFC 2606 reserves the `.invalid` TLD, and the `example.com`/`.net`/`.org`
    # domains, for use in documentation and testing.  The source uses both for
    # URLs which are deliberately not real -- the download system test's
    # shell-injection payloads and its dummy download target, for example -- so
    # such a URL failing is the expected result, not a broken link.
    re.compile(r'^https?://[^/?#]*\.invalid(?:[:/?#]|$)'),
    re.compile(r'^https?://(?:[^/?#]*\.)?example\.(?:com|net|org)(?:[:/?#]|$)'),
]


# Hosts which serve HTTP 403 to the checker -- bot or datacenter-IP blocking by
# a CDN -- while the page is up for ordinary clients.  A browser-like
# user-agent does not help for these (unlike w3schools.com, handled with one
# below), and a genuinely-missing page would answer 404/410, so a 403 from one
# of these hosts is tolerated rather than reported.
_FORBIDDEN_TOLERATED_RES = [
    re.compile(r'^https://www\.openmp\.org/'),
    re.compile(r'^https?://(?:www\.)?math\.stackexchange\.com/'),
    re.compile(r'^https://goldbook\.iupac\.org/'),
    re.compile(r'^https://journals\.aps\.org/'),
]


def tolerates_forbidden(url):
    """Return True if an HTTP 403 from this URL's host is expected, not a fault."""
    return any(pattern.search(url) for pattern in _FORBIDDEN_TOLERATED_RES)


def is_excluded(url):
    """Return True if the URL is a known placeholder that should not be checked."""
    return any(pattern.search(url) for pattern in _EXCLUDED_URL_RES)


def scan_file(file_name, path, urls):
    """Scan a source / RST / TeX file for ``http(s)`` URLs in any markup form."""
    line_number = 0
    try:
        with open(os.path.join(path, file_name), 'r', errors='replace') as f:
            for line in f:
                line_number += 1
                for m in _URL_RE.finditer(line):
                    # Skip a URL which is only one operand of a Fortran string
                    # concatenation — it is a fragment, not a URL.
                    if (_CONCATENATION_RE.search(m.group(0)) or
                            _CONCATENATION_AFTER_RE.match(line[m.end():])):
                        continue
                    # Unescape XML/LaTeX (``&amp;``, ``&#x2F;``, ``\_`` …).
                    url = html.unescape(m.group(0))
                    url = re.sub(r'\\(.)', r'\1', url)
                    # Unescaping can expose a delimiter that terminated the URL
                    # in the source but was entity-encoded there, so truncate
                    # again before trimming punctuation.
                    url = _URL_TERMINATOR_RE.split(url, 1)[0]
                    # Drop trailing sentence punctuation, any closing-quote
                    # delimiter (a ``'`` that wraps the URL, vs. an apostrophe
                    # inside it which is followed by further URL characters),
                    # and a trailing ``)`` that closes an enclosing construct
                    # (shell ``$(curl ...)``, prose parenthetical) rather than
                    # belonging to the URL.  Balanced parens — e.g. Wikipedia's
                    # ``..._(computer_programming)`` — are kept; only an excess
                    # of closing over opening parens is stripped.  Repeat until
                    # stable, since the two forms interleave (``…/meta.)``).
                    previous = None
                    while url != previous:
                        previous = url
                        url = url.rstrip('.,;:\'')
                        while url.endswith(')') and url.count(')') > url.count('('):
                            url = url[:-1]
                    if not re.match(r'^https?://\S', url):
                        continue
                    urls.setdefault(url, []).append(
                        {'file': file_name, 'path': path, 'lineNumber': line_number})
    except OSError as e:
        print(f"Warning: could not read {path}/{file_name}: {e}")


def scan_sources(source_path, urls):
    """Collect URLs from every Fortran source file under `source_path`.

    The tree must be *walked*: `source/` is a directory hierarchy, so all but a
    handful of its files live in subdirectories, and a non-recursive listing
    would see almost none of the URLs embedded in the source.
    """
    for root, _dirs, files in os.walk(source_path):
        for file_name in files:
            if re.search(r'\.(F90|Inc)$', file_name):
                scan_file(file_name, root, urls)


def _find_closing_paren(s, open_pos):
    """Return the index of the ')' that balances the '(' at open_pos."""
    depth = 0
    for i in range(open_pos, len(s)):
        if s[i] == '(':
            depth += 1
        elif s[i] == ')':
            depth -= 1
            if depth == 0:
                return i
    return -1


def scan_wiki(file_name, path, urls):
    """Scan a wiki Markdown file for links."""
    line_number = 0
    # Matches [text]( — the URL portion is then extracted with balanced-paren logic.
    _link_start_re = re.compile(r'\[([^\]]*)\]\(')
    try:
        with open(os.path.join(path, file_name), 'r', errors='replace') as f:
            for line in f:
                line_number += 1
                # Markdown [text](url) patterns — use balanced () to handle
                # URLs that themselves contain parentheses, e.g.
                # https://en.wikipedia.org/wiki/Destructor_(computer_programming)
                while True:
                    m = _link_start_re.search(line)
                    if not m:
                        break
                    open_pos = m.end() - 1   # index of the '(' in line
                    close_pos = _find_closing_paren(line, open_pos)
                    if close_pos == -1:
                        break
                    url = line[open_pos + 1:close_pos]
                    line = line[:m.start()] + line[close_pos + 1:]
                    urls.setdefault(url, []).append(
                        {'file': file_name, 'path': path, 'lineNumber': line_number})
    except OSError as e:
        print(f"Warning: could not read {path}/{file_name}: {e}")


def check_urls(urls, api_token, failures):
    """Check all collected URLs. Returns (exit_status, list_of_bad_urls)."""
    status = 0
    bad_urls = []
    bib_codes = {}

    url_keys = list(urls.keys())
    random.shuffle(url_keys)

    _ads_re = re.compile(r'adsabs\.harvard\.edu/abs/([^/]+)')

    for url_key in url_keys:
        url = url_key.replace('\\#', '#')

        if url.startswith('mailto:'):
            continue
        if url.startswith('#'):
            continue
        if is_excluded(url):
            continue

        # --- NASA ADS links ---
        m = _ads_re.search(url)
        if m:
            bib_code = m.group(1).replace('%26', '&')
            if bib_code not in bib_codes:
                bib_codes[bib_code] = {'sources': [], 'urls': {}}
            bib_codes[bib_code]['sources'].extend(urls[url_key])
            bib_codes[bib_code]['urls'][url_key] = 1
            continue

        # --- General external links ---
        time.sleep(1)
        options = [
            '--max-time', '60', '--insecure', '--location',
            '--output', '/dev/null', '--fail-with-body',
            '--cipher', 'DEFAULT:!DH',
        ]
        no_range = (
            re.search(r'^https://www\.drdobbs\.com/', url) or
            re.search(r'^https://www\.openmp\.org/', url) or
            re.search(r'^https://git-scm\.com/', url)
        )
        if not no_range:
            options += ['--range', '0-0']
        if re.search(r'sharepoint\.com', url):
            options += ['--user-agent', 'Mozilla']
        # w3schools.com rejects curl's default user-agent with HTTP 403 but
        # serves normally to a browser-like user-agent.
        if re.search(r'w3schools\.com', url):
            options += ['--user-agent', 'Mozilla']
        if re.search(r'docker\.com', url):
            options += ['--user-agent', 'Wget/1.21.2']
        if re.search(r'docs\.github\.com', url):
            options += ['--compressed']
        if re.search(r'camb\.info', url):
            options += ['--http1.1']
        if re.search(r'www\.gnu\.org', url):
            options += ['--retry', '5']

        with open('curl.log', 'w') as log_fh:
            result = subprocess.run(
                ['curl'] + options + [url],
                stdout=log_fh, stderr=log_fh)

        error = result.returncode != 0

        if error:
            # Check for known tolerated errors
            with open('curl.log', 'r') as log_fh:
                for line in log_fh:
                    if re.search(r'http://heasarc\.gsfc\.nasa\.gov/xanadu/xspec/', url):
                        if re.search(
                                r'error:0A000126:SSL routines::unexpected eof while reading, errno 0',
                                line):
                            error = False
                            break
                    if re.search(r'www\.gnu\.org', url):
                        if re.search(
                                r'curl: \(28\) Connection timed out after \d+ milliseconds',
                                line):
                            error = False
                            break
                    if tolerates_forbidden(url):
                        if re.search(
                                r'curl: \(22\) The requested URL returned error: 403',
                                line):
                            error = False
                            break

        if error:
            if record_failure(url, failures):
                status = 1
                bad_urls.append(url)
            print(f"Broken link: \"{url}\" "
                  f"(for past {failures[url]['consecutiveFailures']} attempts) in:")
            for source in urls[url_key]:
                print(f"\t{source['path']}/{source['file']} line {source['lineNumber']}")
            print("Log:")
            with open('curl.log', 'r') as log_fh:
                print(log_fh.read(), end='')
        else:
            record_success(url, failures)

    # --- NASA ADS bulk bibcode check ---
    if bib_codes:
        count_records = len(bib_codes)
        ads_url = (
            f'https://api.adsabs.harvard.edu/v1/search/bigquery'
            f'?q=*:*&rows={count_records}'
            f'&fl=bibcode,alternate_bibcode,title,author,year,pub,volume,page'
        )
        headers = {
            'Authorization': f'Bearer {api_token}',
            'Content-Type': 'big-query/csv',
        }
        post_data = 'bibcode\n' + '\n'.join(sorted(bib_codes))

        try:
            response = requests.post(ads_url, headers=headers, data=post_data, timeout=60)
            if response.status_code == 200:
                records = response.json()
                for entry in records.get('response', {}).get('docs', []):
                    found = False
                    if entry['bibcode'] in bib_codes:
                        bib_codes[entry['bibcode']]['found'] = True
                        found = True
                    for alt in entry.get('alternate_bibcode', []):
                        if alt in bib_codes:
                            bib_codes[alt]['found'] = True
                            found = True
                    if not found:
                        status = 1
                        print(f"Received unrequested record for bibcode '{entry['bibcode']}'")
            else:
                status = 1
                print(f"Failed to retrieve record identifiers: "
                      f"{response.status_code} {response.text}")
        except requests.RequestException as e:
            status = 1
            print(f"Failed to retrieve record identifiers: {e}")

        for bib_code, data in bib_codes.items():
            if data.get('found'):
                record_success(bib_code, failures)
            else:
                if record_failure(bib_code, failures):
                    status = 1
                    bad_urls.extend(data['urls'])
                consec = failures.get(bib_code, {}).get('consecutiveFailures', 0)
                print(f"Broken link (for past {consec} attempts): "
                      f"{{bibCode: {bib_code}}} "
                      f"\"{'; '.join(data['urls'])}\" in:")
                for source in data['sources']:
                    print(f"\t{source['path']}/{source['file']} line {source['lineNumber']}")

    return status, bad_urls


def main():
    if len(sys.argv) != 3:
        print("Usage: linkChecker.py <apiToken> <logFile>", file=sys.stderr)
        sys.exit(1)

    api_token    = sys.argv[1]
    log_filename = sys.argv[2]

    # Load existing consecutive-failure records.
    failures_file = 'linkCheckFailures.xml'
    failures = load_failures(failures_file)

    print("Current consecutive failure count (in order of increasing number of failures:")
    for url in sorted(failures, key=lambda u: failures[u]['consecutiveFailures']):
        print(f"{failures[url]['consecutiveFailures']}:\t{url}")

    # Collect URLs from source files.
    urls = {}

    # Embedded docstrings (RST) and constant/workaround directive attributes.
    scan_sources('./source', urls)

    # Committed RST documentation (manuals + landing page) and the glossary.
    for base in ('./docs', './doc'):
        for root, _dirs, files in os.walk(base):
            if '_build' in root.split(os.sep):
                continue
            for file_name in files:
                if file_name.endswith('.rst') or file_name.endswith('.tex'):
                    scan_file(file_name, root, urls)

    wiki_path = './galacticus.wiki'
    subprocess.run(
        ['git', 'clone', 'https://github.com/galacticusorg/galacticus.wiki.git'],
        check=False)
    subprocess.run(['git', 'pull'], cwd=wiki_path, check=False)
    if os.path.isdir(wiki_path):
        for file_name in os.listdir(wiki_path):
            if file_name.endswith('.md'):
                scan_wiki(file_name, wiki_path, urls)

    # Check all collected URLs.
    status, bad_urls = check_urls(urls, api_token, failures)
    
    # Persist updated failure records.
    save_failures(failures, failures_file)

    # Write log of bad URLs if any were found.
    if bad_urls:
        with open(log_filename, 'w') as log_fh:
            log_fh.write('# :warning: Broken links found :warning:\n')
            for bad_url in bad_urls:
                log_fh.write(f'* [`{bad_url}`]({bad_url})\n')

    sys.exit(status)


if __name__ == '__main__':
    main()
