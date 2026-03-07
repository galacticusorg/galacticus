#!/usr/bin/env python3
"""Scan Fortran90/C++/Perl source code and extract contributor data from '!+' / '#+' lines.

Andrew Benson 25-Mar-2012 (original Perl); Python port 2026.

Usage: Extract_Contributors.py <galacticusDir> <outputFile>
"""

import os
import re
import sys
from collections import defaultdict
from pathlib import Path

# Mapping of accented/special characters to their LaTeX equivalents.
ACCENT_MAP = {
    'Ç': r'\c{C}',
    'ü': r'\"u',
    'é': r"\'e",
    'â': r'\^a',
    'ä': r'\"a',
    'à': r'\`a',
    'ç': r'\c{c}',
    'ê': r'\^e',
    'ë': r'\"e',
    'è': r'\`e',
    'ï': r'\"i',
    'î': r'\^i',
    'ì': r'\`i',
    'Ä': r'\"A',
    'Å': r'\AA',
    'É': r"\'E",
    'ô': r'\^o',
    'ö': r'\"o',
    'ò': r'\`o',
    'û': r'\^u',
    'ù': r'\`u',
    'ÿ': r'\"y',
    'Ö': r'\"O',
    'Ü': r'\"U',
    'á': r"\'a",
    'í': r"\'i",
    'ó': r"\'o",
    'ú': r"\'u",
    'ñ': r'\~n',
    'Ñ': r'\~N',
}

# Special LaTeX characters that must be escaped in arbitrary text (e.g. filenames).
_LATEX_SPECIAL = [
    ('\\', r'\textbackslash{}'),
    ('{',  r'\{'),
    ('}',  r'\}'),
    ('$',  r'\$'),
    ('&',  r'\&'),
    ('%',  r'\%'),
    ('#',  r'\#'),
    ('^',  r'\^{}'),
    ('_',  r'\_'),
    ('~',  r'\~{}'),
]


def latex_encode(text):
    """Escape special LaTeX characters in *text* (equivalent to LaTeX::Encode)."""
    # Backslash must be replaced first to avoid double-escaping.
    for char, replacement in _LATEX_SPECIAL:
        text = text.replace(char, replacement)
    return text


def encode_name_for_latex(name):
    """Replace accented characters in a contributor name with LaTeX accent commands."""
    for char, replacement in ACCENT_MAP.items():
        name = name.replace(char, replacement)
    return name


def extract_contributors(line, marker_re):
    """Return list of contributor names found in *line* using *marker_re*.

    *marker_re* must be a compiled pattern with one capturing group that
    captures everything after the marker (e.g. ``!+`` or ``#+``).
    """
    match = marker_re.match(line)
    if not match:
        return []
    contributors = match.group(1)
    # Strip any leading label such as "Author:" or "Contributors:".
    contributors = re.sub(r'^.*:\s*', '', contributors)
    # Strip trailing period.
    contributors = re.sub(r'\.\s*$', '', contributors)
    return [p.strip() for p in re.split(r'\s*,\s*', contributors) if p.strip()]


def scan_file(file_path, marker_re, contributions, key):
    """Scan *file_path* for contributor markers and record them under *key*."""
    with open(file_path, encoding='utf-8', errors='replace') as fh:
        for line in fh:
            for person in extract_contributors(line, marker_re):
                contributions[person][key] = True


def main():
    if len(sys.argv) != 3:
        print(
            'Usage: Extract_Contributors.py <galacticusDir> <outputFile>',
            file=sys.stderr,
        )
        sys.exit(1)

    galacticus_dir = sys.argv[1]
    output_file    = sys.argv[2]

    # contributions[person][file_key] = True
    contributions = defaultdict(dict)

    fortran_marker = re.compile(r'^\s*!\+\s*(.*)')
    perl_marker    = re.compile(r'^\s*#\+\s*(.*)')

    # Scan Fortran 90, include, and C++ source files.
    source_dir = os.path.join(galacticus_dir, 'source')
    if os.path.isdir(source_dir):
        for file_name in os.listdir(source_dir):
            if (re.search(r'\.(F90|Inc|cpp)$', file_name)
                    and not re.match(r'^\.\#', file_name)):
                scan_file(
                    os.path.join(source_dir, file_name),
                    fortran_marker,
                    contributions,
                    key=file_name,  # store bare filename, matching Perl behaviour
                )

    # Scan Perl modules (full path stored, matching Perl behaviour).
    perl_dir = os.path.join(galacticus_dir, 'perl')
    if os.path.isdir(perl_dir):
        for pm_path in sorted(Path(perl_dir).rglob('*.pm')):
            scan_file(str(pm_path), perl_marker, contributions, key=str(pm_path))

    # Generate LaTeX output.
    with open(output_file, 'w', encoding='utf-8') as out:
        out.write('\\chapter{Contributions}\n\n')
        out.write(
            'Contributions to the \\glc\\ project have been made by the following people:\n\n'
        )
        out.write('\\begin{description}\n')

        for person in sorted(contributions):
            name = encode_name_for_latex(person)
            out.write(f'\\item[{name}] \\hfill\n')
            out.write('\\begin{itemize}\n')
            for file_key in sorted(contributions[person]):
                encoded = latex_encode(file_key)
                if re.search(r'\.(pm|Inc)$', file_key):
                    # Perl modules and include files: no hyperlink.
                    out.write(f'\\item {{\\normalfont \\ttfamily {encoded}}}\n')
                else:
                    out.write(
                        f'\\item \\hyperlink{{{file_key}}}{{\\normalfont \\ttfamily {encoded}}}\n'
                    )
            out.write('\\end{itemize}\n')

        out.write('\\end{description}\n')


if __name__ == '__main__':
    main()
