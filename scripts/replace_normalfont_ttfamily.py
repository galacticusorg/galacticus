#!/usr/bin/env python3
"""Replace {\normalfont \ttfamily ...} with \source{...} throughout the codebase."""

import os
import sys
import glob

PATTERN = '{\\normalfont \\ttfamily '

def replace_in_text(text):
    """Replace all occurrences of {\\normalfont \\ttfamily CONTENT} with \\source{CONTENT},
    correctly handling nested braces."""
    result = []
    i = 0
    count = 0
    while i < len(text):
        idx = text.find(PATTERN, i)
        if idx == -1:
            result.append(text[i:])
            break
        # Append text before the pattern
        result.append(text[i:idx])
        # Walk from idx+len(PATTERN) to find the matching closing brace
        # At idx we have '{', depth starts at 1
        depth = 1
        j = idx + len(PATTERN)
        while j < len(text) and depth > 0:
            c = text[j]
            if c == '{':
                depth += 1
            elif c == '}':
                depth -= 1
            if depth > 0:
                j += 1
            # if depth==0 we stop (j points to the closing '}')
        if depth != 0:
            # Unmatched brace - leave as-is and move forward one char
            result.append(text[idx])
            i = idx + 1
            continue
        # text[idx+len(PATTERN):j] is the content
        content = text[idx + len(PATTERN):j]
        result.append('\\source{' + content + '}')
        i = j + 1  # skip past the closing '}'
        count += 1
    return ''.join(result), count


def process_file(filepath):
    with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
        original = f.read()
    if PATTERN not in original:
        return 0
    new_text, count = replace_in_text(original)
    if new_text != original:
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(new_text)
    return count


def main():
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    total_files = 0
    total_replacements = 0

    # Process .tex files in doc/
    tex_files = glob.glob(os.path.join(repo_root, 'doc', '**', '*.tex'), recursive=True)
    for filepath in sorted(tex_files):
        count = process_file(filepath)
        if count:
            print(f"  {os.path.relpath(filepath, repo_root)}: {count} replacement(s)")
            total_files += 1
            total_replacements += count

    # Process .F90 files everywhere in the repo
    f90_files = glob.glob(os.path.join(repo_root, '**', '*.F90'), recursive=True)
    for filepath in sorted(f90_files):
        count = process_file(filepath)
        if count:
            print(f"  {os.path.relpath(filepath, repo_root)}: {count} replacement(s)")
            total_files += 1
            total_replacements += count

    print(f"\nTotal: {total_replacements} replacements across {total_files} files.")


if __name__ == '__main__':
    main()
