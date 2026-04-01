#!/usr/bin/env python3
import json
import os
import re
import shutil
import sys
import urllib.request
import urllib.error

# Update bibliography entries.
# Andrew Benson (18-March-2024)
## Entries that are arXiv records will be pulled from NASA ADS which will return the actual
## journal reference where available.

if len(sys.argv) != 2:
    print("Usage: bibliographyUpdate.py <apiToken>", file=sys.stderr)
    raise SystemExit(1)

api_token = sys.argv[1]
bib_file  = os.path.join(os.environ['GALACTICUS_EXEC_PATH'], 'doc', 'Galacticus.bib')
bib_new   = bib_file + '.new'


def get_bibtex_key(entry_text):
    """Extract the citation key from a BibTeX entry."""
    m = re.match(r'@\w+\{([^,]+),', entry_text.strip())
    return m.group(1).strip() if m else None


def get_bibtex_field(entry_text, field_name):
    """Extract a field value from a BibTeX entry, handling nested braces."""
    pattern = re.compile(
        r'(?:^|[\s,])' + re.escape(field_name) + r'\s*=\s*',
        re.IGNORECASE | re.MULTILINE
    )
    m = pattern.search(entry_text)
    if not m:
        return None
    pos = m.end()
    if pos >= len(entry_text):
        return None
    ch = entry_text[pos]
    if ch == '{':
        depth = 0
        start = pos + 1
        for i in range(pos, len(entry_text)):
            if entry_text[i] == '{':
                depth += 1
            elif entry_text[i] == '}':
                depth -= 1
                if depth == 0:
                    return entry_text[start:i]
        return None
    elif ch == '"':
        end = entry_text.find('"', pos + 1)
        return entry_text[pos + 1:end] if end != -1 else None
    else:
        m2 = re.match(r'(\S+)', entry_text[pos:])
        return m2.group(1).rstrip(',') if m2 else None


# Read all bibliography entries as raw text, splitting on lines that start with '@'.
with open(bib_file, 'r') as f:
    content = f.read()
entries = re.split(r'(?=^@)', content, flags=re.MULTILINE)
entries = [e for e in entries if e.strip()]

count_updated = 0

with open(bib_new, 'w') as out:
    for entry_text in entries:
        # By default, assume we can simply reuse the raw text of the existing entry.
        output_text = entry_text
        # Check for an arXiv entry.
        journal = get_bibtex_field(entry_text, 'journal')
        if journal and re.search(r'arxiv', journal, re.IGNORECASE):
            url = get_bibtex_field(entry_text, 'url') or get_bibtex_field(entry_text, 'adsurl')
            if url is None:
                print("Error: no URL found", file=sys.stderr)
                raise SystemExit(1)
            identifier = None
            m = re.match(r'https?://(ui\.)?adsabs\.harvard\.edu/abs/(.+)', url)
            if m:
                # URL is an ADS arXiv reference. Extract the identifier directly.
                identifier = m.group(2)
            else:
                m = re.match(r'https?://arxiv\.org/abs/([\d\.]+)', url)
                if m:
                    # URL is a direct arXiv reference. Construct the corresponding ADS identifier.
                    identifier = m.group(1)
                    year = get_bibtex_field(entry_text, 'year')
                    if year:
                        identifier = year + 'arXiv' + identifier
                    else:
                        print("no 'year' field exists")
                    author = get_bibtex_field(entry_text, 'author')
                    if author:
                        identifier += author[0]
                    else:
                        print("no 'author' field exists")
                else:
                    print(f"URL '{url}' is not recognized")
            if identifier is None:
                key = get_bibtex_key(entry_text)
                raise SystemExit(f'unable to update entry for arXiv record - key is "{key}"')
            # Query NASA ADS for the BibTeX record.
            payload = json.dumps({"bibcode": [identifier]}).encode('utf-8')
            req = urllib.request.Request(
                'https://api.adsabs.harvard.edu/v1/export/bibtex',
                data=payload,
                headers={
                    'Authorization': f'Bearer {api_token}',
                    'Content-Type':  'application/json',
                },
                method='POST'
            )
            try:
                with urllib.request.urlopen(req, timeout=30) as response:
                    response_body = response.read().decode('utf-8')
            except urllib.error.URLError as e:
                raise SystemExit(f"Failed to retrieve record identifier '{identifier}': {e}")
            bibtex_entry = json.loads(response_body)
            if bibtex_entry.get('export'):
                export = bibtex_entry['export']
                # If the returned record is something other than an arXiv reference, update to this new record.
                is_arxiv = (
                    re.search(
                        r'^\s*adsurl\s*=\s*\{https://ui\.adsabs\.harvard\.edu/abs/\d+arXiv[\d\.]+[A-Z]\}',
                        export, re.MULTILINE
                    ) or
                    re.search(
                        r'^\s*adsurl\s*=\s*\{https://ui\.adsabs\.harvard\.edu/abs/\d+astro\.ph\.\d+[A-Z]\}',
                        export, re.MULTILINE
                    )
                )
                if not is_arxiv:
                    # Replace the BibTeX key with the original key so we don't have to
                    # re-write references in all of our documentation.
                    count_updated += 1
                    key         = get_bibtex_key(entry_text)
                    output_text = re.sub(r'^@ARTICLE\{[^,]+,', f'@article{{{key},', export, count=1)
        out.write(output_text)

# Replace the original bibliography.
shutil.move(bib_new, bib_file)

# Report.
print(f"Updated {count_updated} of {len(entries)} records")
