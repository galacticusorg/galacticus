#!/usr/bin/env python3
"""Extract the Galacticus contributor list from source markers.

Andrew Benson 25-Mar-2012 (original Perl); Python port 2026.

Contributor names are recorded inline in the source as ``!+`` lines (Fortran,
include, and C++ files) and ``#+`` lines (Perl modules).  This module both:

* runs as a script — ``extractContributors.py <galacticusDir> <outputFile>`` —
  writing the contributor list as a reStructuredText paragraph; and
* exposes ``collect_contributor_names()`` / ``render_contributors_rst()`` for
  ``extractDocsRST.py`` to regenerate the Acknowledgments contributor list at
  documentation-build time.

Names are read directly as UTF-8, so accented names need no special handling in
the RST output.
"""

import os
import re
import sys
from pathlib import Path

_FORTRAN_MARKER = re.compile(r'^\s*!\+\s*(.*)')
_PERL_MARKER    = re.compile(r'^\s*#\+\s*(.*)')


def extract_contributors(line, marker_re):
    """Return contributor names found in *line* using *marker_re* (a compiled
    pattern with one group capturing everything after the ``!+``/``#+`` marker).

    Markers may carry a leading label ("Author:") and a trailing period, and may
    list several comma-separated names.  A few markers hold a whole sentence
    rather than a name; those (more than five words) are dropped by callers.
    """
    match = marker_re.match(line)
    if not match:
        return []
    contributors = match.group(1)
    contributors = re.sub(r'^.*:\s*', '', contributors)   # strip "Author:" label
    contributors = re.sub(r'\.\s*$', '', contributors)    # strip trailing period
    return [p.strip() for p in re.split(r'\s*,\s*', contributors) if p.strip()]


def _scan(file_path, marker_re, names):
    with open(file_path, encoding='utf-8', errors='replace') as fh:
        for line in fh:
            for person in extract_contributors(line, marker_re):
                if person and len(person.split()) <= 5:   # real names, not prose
                    names.add(person)


def collect_contributor_names(galacticus_dir='.'):
    """Return the de-duplicated contributor names found in the Galacticus
    source, sorted case-insensitively."""
    names = set()

    source_dir = os.path.join(galacticus_dir, 'source')
    if os.path.isdir(source_dir):
        for file_name in sorted(os.listdir(source_dir)):
            if (re.search(r'\.(F90|Inc|cpp)$', file_name)
                    and not file_name.startswith('.#')):
                _scan(os.path.join(source_dir, file_name), _FORTRAN_MARKER, names)

    perl_dir = os.path.join(galacticus_dir, 'perl')
    if os.path.isdir(perl_dir):
        for pm_path in sorted(Path(perl_dir).rglob('*.pm')):
            _scan(str(pm_path), _PERL_MARKER, names)

    return sorted(names, key=str.casefold)


def render_contributors_rst(names):
    """Render contributor *names* as a single RST paragraph ending in a period."""
    return ', '.join(names) + '.\n'


def main():
    if len(sys.argv) != 3:
        print('Usage: extractContributors.py <galacticusDir> <outputFile>',
              file=sys.stderr)
        sys.exit(1)
    names = collect_contributor_names(sys.argv[1])
    with open(sys.argv[2], 'w', encoding='utf-8') as out:
        out.write(render_contributors_rst(names))


if __name__ == '__main__':
    main()
