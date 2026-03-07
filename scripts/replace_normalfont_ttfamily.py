#!/usr/bin/env python3
"""Replace {\normalfont \ttfamily ...} with \mono{...} throughout the codebase.

Handles both single-line and split-line occurrences, and Perl files that
emit LaTeX using double-escaped {\normalfont \ttfamily ...}.
"""

import os
import re
import glob

# Matches '{\normalfont \ttfamily' followed by any whitespace (incl. newlines)
OPEN_RE = re.compile(r'\{\\normalfont \\ttfamily\s+')
# Perl files use doubled backslashes in string literals
OPEN_RE_PERL = re.compile(r'\{\\\\normalfont \\\\ttfamily\s+')


def replace_in_text(text, open_re, replacement_cmd='\\mono'):
    """Replace all occurrences of {\\normalfont \\ttfamily CONTENT} with
    \\mono{CONTENT}, correctly handling nested braces and multiline content.
    Content that spans multiple lines is joined with a single space."""
    result = []
    i = 0
    count = 0
    while i < len(text):
        m = open_re.search(text, i)
        if m is None:
            result.append(text[i:])
            break
        result.append(text[i:m.start()])
        # Walk forward from end of match to find the matching closing '}'
        depth = 1
        j = m.end()
        while j < len(text) and depth > 0:
            c = text[j]
            if c == '{':
                depth += 1
            elif c == '}':
                depth -= 1
            if depth > 0:
                j += 1
        if depth != 0:
            # Unmatched brace - leave as-is
            result.append(text[m.start()])
            i = m.start() + 1
            continue
        # Extract content between the opening pattern and the closing '}'
        content = text[m.end():j]
        # Normalise any internal line-break+indent to a single space
        content = re.sub(r'\n[ \t]*', ' ', content).rstrip()
        result.append(replacement_cmd + '{' + content + '}')
        i = j + 1
        count += 1
    return ''.join(result), count


def process_file(filepath, perl=False):
    with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
        original = f.read()
    open_re = OPEN_RE_PERL if perl else OPEN_RE
    if not open_re.search(original):
        return 0
    new_text, count = replace_in_text(original, open_re)
    if new_text != original:
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(new_text)
    return count


def main():
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    total_files = 0
    total_replacements = 0

    # .tex files
    for filepath in sorted(glob.glob(os.path.join(repo_root, 'doc', '**', '*.tex'), recursive=True)):
        count = process_file(filepath)
        if count:
            print(f"  {os.path.relpath(filepath, repo_root)}: {count} replacement(s)")
            total_files += 1
            total_replacements += count

    # .F90 files
    for filepath in sorted(glob.glob(os.path.join(repo_root, '**', '*.F90'), recursive=True)):
        count = process_file(filepath)
        if count:
            print(f"  {os.path.relpath(filepath, repo_root)}: {count} replacement(s)")
            total_files += 1
            total_replacements += count

    # Perl files (.pl and .pm) — backslashes are doubled in string literals
    for pattern in ('**/*.pl', '**/*.pm'):
        for filepath in sorted(glob.glob(os.path.join(repo_root, pattern), recursive=True)):
            count = process_file(filepath, perl=True)
            if count:
                print(f"  {os.path.relpath(filepath, repo_root)}: {count} replacement(s)")
                total_files += 1
                total_replacements += count

    print(f"\nTotal: {total_replacements} replacements across {total_files} files.")


if __name__ == '__main__':
    main()
