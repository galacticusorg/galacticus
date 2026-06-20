#!/usr/bin/env python3
import argparse
import json
import os
import re
import shutil
import sys
import urllib.request
import urllib.error

# Update bibliography entries.
# Andrew Benson (18-March-2024)
## Default mode: any entry that looks like an arXiv preprint is queried against
## NASA ADS, and replaced if ADS now has a published (non-preprint) record for
## it. With --refresh-all, every entry for which an ADS bibcode can be derived
## is re-queried, and replaced if the published metadata (volume, pages, DOI,
## year) has changed.

ARXIV_BIBCODE_RE = re.compile(r'^\d{4}(arXiv|astro\.ph)', re.IGNORECASE)


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


def first_author_initial(author_field):
    """Extract the first author's family-name initial for ADS bibcode construction."""
    if not author_field:
        return None
    first = re.split(r'\s+and\s+', author_field, flags=re.IGNORECASE)[0]
    if ',' in first:
        last = first.split(',', 1)[0]
    else:
        parts = first.strip().split()
        last = parts[-1] if parts else ''
    last = re.sub(r'[\\\{\}\'"~^]', '', last).strip()
    return last[:1].upper() if last else None


def clean_bibcode(s):
    """Normalize a bibcode pulled from a free-text source (notes, URLs)."""
    if s is None:
        return None
    s = s.strip().rstrip(',').rstrip('}').strip()
    # ADS bibcodes containing '&' are sometimes BibTeX-escaped as '\&'.
    s = s.replace('\\&', '&')
    return s


def extract_bibcode(entry_text):
    """Try to determine an ADS bibcode for the entry. Returns the bibcode or None."""
    # 1. Direct URL/adsurl pointing at ADS.
    for fld in ('adsurl', 'url'):
        v = get_bibtex_field(entry_text, fld)
        if v:
            m = re.search(r'adsabs\.harvard\.edu/abs/([^/}\s]+)', v)
            if m:
                return clean_bibcode(m.group(1))
    # 2. note = {ADS Bibcode: ...}
    note = get_bibtex_field(entry_text, 'note')
    if note:
        m = re.search(r'ADS Bibcode:\s*(\S+)', note)
        if m:
            return clean_bibcode(m.group(1))
    # 3. arXiv URL of the form arxiv.org/abs/NNNN.NNNNN.
    url = get_bibtex_field(entry_text, 'url') or ''
    m = re.search(r'arxiv\.org/abs/(\d{4}\.\d+)', url, re.IGNORECASE)
    if m:
        arxiv_id = m.group(1)
        year = (get_bibtex_field(entry_text, 'year') or '').strip()
        initial = first_author_initial(get_bibtex_field(entry_text, 'author'))
        if year and initial:
            return f"{year}arXiv{arxiv_id}{initial}"
    # 4. eprint + archivePrefix/eprinttype = arXiv.
    eprint = (get_bibtex_field(entry_text, 'eprint') or '').strip()
    archive = (get_bibtex_field(entry_text, 'archivePrefix')
               or get_bibtex_field(entry_text, 'eprinttype') or '')
    if eprint and re.search(r'arxiv', archive, re.IGNORECASE) \
            and re.match(r'^\d{4}\.\d+$', eprint):
        year = (get_bibtex_field(entry_text, 'year') or '').strip()
        initial = first_author_initial(get_bibtex_field(entry_text, 'author'))
        if year and initial:
            return f"{year}arXiv{eprint}{initial}"
    # 5. note = {arXiv:NNNN.NNNNN ...}.
    if note:
        m = re.search(r'arXiv:\s*(\d{4}\.\d+)', note)
        if m:
            arxiv_id = m.group(1)
            year = (get_bibtex_field(entry_text, 'year') or '').strip()
            initial = first_author_initial(get_bibtex_field(entry_text, 'author'))
            if year and initial:
                return f"{year}arXiv{arxiv_id}{initial}"
    return None


def is_preprint_entry(entry_text, bibcode):
    """Decide whether the entry is currently a preprint that may be upgradable."""
    if bibcode and ARXIV_BIBCODE_RE.match(bibcode):
        return True
    journal = get_bibtex_field(entry_text, 'journal')
    if journal and re.search(r'arxiv', journal, re.IGNORECASE):
        return True
    return False


def query_ads(api_token, bibcode):
    payload = json.dumps({"bibcode": [bibcode]}).encode('utf-8')
    req = urllib.request.Request(
        'https://api.adsabs.harvard.edu/v1/export/bibtex',
        data=payload,
        headers={
            'Authorization': f'Bearer {api_token}',
            'Content-Type':  'application/json',
        },
        method='POST'
    )
    with urllib.request.urlopen(req, timeout=30) as response:
        return json.loads(response.read().decode('utf-8'))


def response_is_preprint(export_text):
    return bool(
        re.search(
            r'^\s*adsurl\s*=\s*\{https://ui\.adsabs\.harvard\.edu/abs/\d+arXiv\S+\}',
            export_text, re.MULTILINE
        ) or re.search(
            r'^\s*adsurl\s*=\s*\{https://ui\.adsabs\.harvard\.edu/abs/\d+astro\.ph\S+\}',
            export_text, re.MULTILINE
        )
    )


def normalize_value(s):
    if s is None:
        return ''
    s = s.strip()
    s = re.sub(r'[{}\\]', '', s)
    s = re.sub(r'[‐-―]|--', '-', s)
    return re.sub(r'\s+', ' ', s).lower()


def field_diff(old_entry, new_entry, fields):
    diffs = []
    for f in fields:
        old = normalize_value(get_bibtex_field(old_entry, f))
        new = normalize_value(get_bibtex_field(new_entry, f))
        if old != new and (old or new):
            diffs.append(f)
    return diffs


def rewrite_with_key(export_text, key):
    return re.sub(r'^@(\w+)\{[^,]+,', lambda m: f'@{m.group(1).lower()}{{{key},',
                  export_text, count=1)


def validate_rewrite(new_path, expected_keys):
    """Structurally validate the rewritten bibliography before swapping it in."""
    with open(new_path, 'r') as f:
        text = f.read()
    if text.count('{') != text.count('}'):
        raise SystemExit("Validation failed: unbalanced braces in rewritten bibliography.")
    parts = re.split(r'(?=^@)', text, flags=re.MULTILINE)
    parts = [p for p in parts if p.strip()]
    found_keys = set()
    for p in parts:
        k = get_bibtex_key(p)
        if k:
            found_keys.add(k)
    missing = expected_keys - found_keys
    if missing:
        sample = sorted(missing)[:5]
        raise SystemExit(
            f"Validation failed: {len(missing)} keys missing after rewrite (e.g. {sample})."
        )
    if len(parts) != len(expected_keys):
        # Duplicate keys are allowed in principle, but don't expect entry count to drop.
        if len(parts) < len(expected_keys):
            raise SystemExit(
                f"Validation failed: entry count dropped from {len(expected_keys)} to {len(parts)}."
            )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('api_token', help='NASA ADS API token.')
    ap.add_argument('--refresh-all', action='store_true',
                    help='Re-query every entry with a discoverable ADS bibcode and '
                         'update if the published metadata has changed.')
    ap.add_argument('--summary', default=None,
                    help='Write a Markdown change-log to this file.')
    args = ap.parse_args()

    bib_file = os.path.join(os.environ['GALACTICUS_EXEC_PATH'], 'docs', 'Galacticus.bib')
    bib_new = bib_file + '.new'

    with open(bib_file, 'r') as f:
        content = f.read()
    entries = re.split(r'(?=^@)', content, flags=re.MULTILINE)
    entries = [e for e in entries if e.strip()]

    expected_keys = set()
    changes = []
    count_queried = 0
    count_updated = 0

    refresh_fields = ['volume', 'pages', 'doi', 'year']

    with open(bib_new, 'w') as out:
        for entry_text in entries:
            key = get_bibtex_key(entry_text)
            if key:
                expected_keys.add(key)
            output_text = entry_text

            bibcode = extract_bibcode(entry_text)
            preprint = is_preprint_entry(entry_text, bibcode)

            should_query = False
            preprint_upgrade = False
            if preprint and bibcode:
                should_query = True
                preprint_upgrade = True
            elif args.refresh_all and bibcode:
                should_query = True

            if should_query:
                count_queried += 1
                try:
                    resp = query_ads(args.api_token, bibcode)
                except urllib.error.URLError as e:
                    print(f"WARN: ADS query failed for '{key}' (bibcode={bibcode}): {e}",
                          file=sys.stderr)
                    out.write(output_text)
                    continue
                except (json.JSONDecodeError, ValueError) as e:
                    print(f"WARN: ADS response unparseable for '{key}' (bibcode={bibcode}): {e}",
                          file=sys.stderr)
                    out.write(output_text)
                    continue
                export = resp.get('export') if isinstance(resp, dict) else None
                if not export:
                    out.write(output_text)
                    continue
                if preprint_upgrade:
                    if not response_is_preprint(export):
                        count_updated += 1
                        output_text = rewrite_with_key(export, key)
                        old = (get_bibtex_field(entry_text, 'journal') or 'arXiv').strip()
                        new = (get_bibtex_field(export, 'journal') or '').strip()
                        old_year = (get_bibtex_field(entry_text, 'year') or '').strip()
                        new_year = (get_bibtex_field(export, 'year') or '').strip()
                        changes.append((key, 'preprint-upgrade',
                                        f"{old} ({old_year})", f"{new} ({new_year})"))
                else:
                    diffs = field_diff(entry_text, export, refresh_fields)
                    if diffs:
                        count_updated += 1
                        output_text = rewrite_with_key(export, key)
                        before = ', '.join(
                            f"{f}={(get_bibtex_field(entry_text, f) or '').strip()}"
                            for f in diffs
                        )
                        after = ', '.join(
                            f"{f}={(get_bibtex_field(export, f) or '').strip()}"
                            for f in diffs
                        )
                        changes.append((key, 'metadata-refresh', before, after))
            out.write(output_text)

    validate_rewrite(bib_new, expected_keys)
    shutil.move(bib_new, bib_file)

    print(f"Queried {count_queried} entries, updated {count_updated} of {len(entries)}.")

    if args.summary:
        with open(args.summary, 'w') as f:
            if changes:
                f.write(f"Updated {len(changes)} bibliography entries "
                        f"(queried {count_queried} of {len(entries)}).\n\n")
                f.write("| Key | Kind | Before | After |\n")
                f.write("|-----|------|--------|-------|\n")
                for key, kind, before, after in changes:
                    safe_before = before.replace('|', '\\|')
                    safe_after = after.replace('|', '\\|')
                    f.write(f"| `{key}` | {kind} | {safe_before} | {safe_after} |\n")
            else:
                f.write(f"No bibliography entries required updates "
                        f"(queried {count_queried} of {len(entries)}).\n")


if __name__ == '__main__':
    main()
