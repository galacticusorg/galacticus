#!/usr/bin/env python3
# Check for broken links in Galacticus documentation.
# Andrew Benson (21-September-2020) [Perl]; ported to Python 2026

import sys
import os
import re
import subprocess
import json
import random
import time
import xml.etree.ElementTree as ET
import requests
from PyPDF2 import PdfReader


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


def scan_file(file_name, path, urls):
    """Scan a TeX or Fortran source file for links."""
    line_number = 0
    try:
        with open(os.path.join(path, file_name), 'r', errors='replace') as f:
            for line in f:
                line_number += 1
                # \href{url} patterns
                while True:
                    m = re.search(r'\\href\{([^\}]+)\}', line)
                    if not m:
                        break
                    url = m.group(1)
                    line = line[:m.start()] + line[m.end():]
                    url = re.sub(r'\\(.)', r'\1', url)  # remove LaTeX escapes
                    urls.setdefault(url, []).append(
                        {'file': file_name, 'path': path, 'lineNumber': line_number})
                # \refPhysics{ref} patterns
                while True:
                    m = re.search(r'\\refPhysics\{([^\}]+)\}', line)
                    if not m:
                        break
                    ref = m.group(1)
                    url = ('https://github.com/galacticusorg/galacticus/releases/'
                           'download/bleeding-edge/Galacticus_Physics.pdf#physics.' + ref)
                    line = line[:m.start()] + line[m.end():]
                    urls.setdefault(url, []).append(
                        {'type': 'refPhysics', 'ref': ref,
                         'file': file_name, 'path': path, 'lineNumber': line_number})
                # \refClass{ref} patterns
                while True:
                    m = re.search(r'\\refClass\{([^\}]+)\}', line)
                    if not m:
                        break
                    ref = m.group(1)
                    url = ('https://github.com/galacticusorg/galacticus/releases/'
                           'download/bleeding-edge/Galacticus_Development.pdf#class.' + ref)
                    line = line[:m.start()] + line[m.end():]
                    urls.setdefault(url, []).append(
                        {'type': 'refClass', 'ref': ref,
                         'file': file_name, 'path': path, 'lineNumber': line_number})
    except OSError as e:
        print(f"Warning: could not read {path}/{file_name}: {e}")


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


def check_urls(urls, pdf_destinations, api_token, failures):
    """Check all collected URLs. Returns (exit_status, list_of_bad_urls)."""
    status = 0
    bad_urls = []
    bib_codes = {}

    url_keys = list(urls.keys())
    random.shuffle(url_keys)

    _pdf_re = re.compile(
        r'^https://github\.com/galacticusorg/galacticus/releases/download/bleeding-edge/'
        r'Galacticus_(Usage|Physics|Development|Source)\.pdf#(.+)')
    _ads_re = re.compile(r'adsabs\.harvard\.edu/abs/([^/]+)')

    for url_key in url_keys:
        url = url_key.replace('\\#', '#')

        if url.startswith('mailto:'):
            continue
        if url.startswith('#'):
            continue

        # --- PDF anchor links ---
        m = _pdf_re.match(url)
        if m:
            suffix = m.group(1)
            anchor = m.group(2)
            if suffix not in pdf_destinations or anchor not in pdf_destinations[suffix]:
                status = 1
                source0 = urls[url_key][0]
                bad_urls.append(source0.get('ref', url))
                print(f"Broken {source0.get('type', 'link')}"
                      f"{{{source0.get('ref', url)}}} link in:")
                for source in urls[url_key]:
                    print(f"\t {source['path']}/{source['file']}"
                          f" line {source['lineNumber']}")
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

    # Download PDFs and extract named destinations.
    pdf_destinations = {}
    for suffix in ('Usage', 'Physics', 'Development', 'Source'):
        pdf_file   = f'Galacticus_{suffix}.pdf'
        dests_file = f'Galacticus_{suffix}.dests'
        subprocess.run([
            'curl', '--silent', '-L', '--insecure', '--output', pdf_file,
            f'https://github.com/galacticusorg/galacticus/releases/download/'
            f'bleeding-edge/{pdf_file}',
        ], check=False)

        with open('Galacticus_'+suffix+'.pdf','rb') as file:
            reader = PdfReader(file)
            pdf_destinations[suffix] = {}
            for dest in reader.named_destinations:
                pdf_destinations[suffix][dest] = 1

    # Load existing consecutive-failure records.
    failures_file = 'linkCheckFailures.xml'
    failures = load_failures(failures_file)

    print("Current consecutive failure count (in order of increasing number of failures:")
    for url in sorted(failures, key=lambda u: failures[u]['consecutiveFailures']):
        print(f"{failures[url]['consecutiveFailures']}:\t{url}")

    # Collect URLs from source files.
    urls = {}

    doc_path = './doc'
    if os.path.isdir(doc_path):
        for file_name in os.listdir(doc_path):
            if re.match(r'^[A-Z].*\.tex$', file_name):
                scan_file(file_name, doc_path, urls)

    source_path = './source'
    if os.path.isdir(source_path):
        for file_name in os.listdir(source_path):
            if re.search(r'\.(F90|Inc)$', file_name):
                scan_file(file_name, source_path, urls)

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
    status, bad_urls = check_urls(urls, pdf_destinations, api_token, failures)
    
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
