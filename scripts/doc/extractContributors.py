#!/usr/bin/env python3
"""Extract the Galacticus contributor list from source markers.

Andrew Benson 25-Mar-2012.

Contributor names are recorded inline in the source as ``!+`` lines (Fortran,
include, and C++ files).  This module both:

* runs as a script — ``extractContributors.py <galacticusDir> <outputFile>`` —
  writing the contributor list as a reStructuredText paragraph; and
* exposes ``collect_contributor_names()`` / ``render_contributors_rst()`` for
  ``extractDocsRST.py`` to regenerate the Acknowledgments contributor list at
  documentation-build time.

Names are read directly as UTF-8, so accented names need no special handling in
the RST output.

A marker records **names only** — the canonical form is::

    !+    Contributions to this file made by: Andrew Benson, Claude.

Everything after the label is parsed as a comma-separated list of contributors,
so a description of *what* was contributed would be rendered as if it were a
person: prose markers previously leaked entries such as "assistance from
Claude" and "per-object OpenMP lock" into the documentation's Acknowledgments.
Where an AI tool contributed, name it as a person would be named (``Claude``,
``Copilot``, …) and describe its contribution in the commit message.

``--check`` enforces this (for CI): it reports every marker which departs from
the canonical form, and exits non-zero if any is found.  See the *Contributor
Attribution* section of ``CONTRIBUTING.md``.
"""

import argparse
import os
import re
import sys

_FORTRAN_MARKER = re.compile(r'^\s*!\+\s*(.*)')

# The canonical label which opens a marker.  A marker may be continued onto
# further `!+` lines (for a long list of names), which carry no label.
MARKER_LABEL = 'Contributions to this file made by:'

# Lowercase particles which may legitimately appear within a name, and so are
# exempt from the "words are capitalized" rule below.
_NAME_PARTICLES = frozenset({
    'al', 'bin', 'da', 'das', 'de', 'del', 'della', 'den', 'der', 'di', 'do',
    'dos', 'du', 'ibn', 'la', 'le', 'van', 'von', 'zu',
})

# Characters which may appear within a word of a name, besides letters.
_NAME_PUNCTUATION = frozenset("'’.-")

_MAX_NAME_WORDS = 5


def extract_contributors(line, marker_re):
    """Return contributor names found in *line* using *marker_re* (a compiled
    pattern with one group capturing everything after the ``!+`` marker).

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


def name_problem(name):
    """Return a description of why *name* is not a personal name, or ``None``
    if it is one.

    This is deliberately shape-based rather than a lookup: it must accept names
    never seen before (including accented and hyphenated ones) while rejecting
    the fragments of prose which result from describing a contribution in the
    marker.
    """
    words = name.split()
    if not words:
        return 'is empty'
    if len(words) > _MAX_NAME_WORDS:
        return f'has more than {_MAX_NAME_WORDS} words, so reads as prose'
    for word in words[:-1]:
        # A period within a name marks an initial or an abbreviation ("S.",
        # "Jr."); a longer word ending in one ends a sentence, so what follows
        # it is a description rather than a name.
        if word.endswith('.') and len(word) > 3:
            return f'continues past the sentence ending "{word}"'
    for word in words:
        if not any(character.isalpha() for character in word):
            return f'contains "{word}", which is not a word of a name'
        for character in word:
            if not character.isalpha() and character not in _NAME_PUNCTUATION:
                return f'contains "{word}", which is not a word of a name'
        if not word[0].isupper() and word.lower() not in _NAME_PARTICLES:
            return f'contains the uncapitalized word "{word}", so reads as prose'
    return None


def iter_source_files(galacticus_dir='.'):
    """Yield the paths of the source files which may carry ``!+`` markers."""
    source_dir = os.path.join(galacticus_dir, 'source')
    if not os.path.isdir(source_dir):
        return
    for root, _dirs, files in sorted(os.walk(source_dir)):
        for file_name in sorted(files):
            if (re.search(r'\.(F90|Inc|cpp)$', file_name)
                    and not file_name.startswith('.#')):
                yield os.path.join(root, file_name)


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

    for file_path in iter_source_files(galacticus_dir):
        _scan(file_path, _FORTRAN_MARKER, names)

    return sorted(names, key=str.casefold)


def check_markers_in_text(text, path='<text>'):
    """Return a list of ``(path, lineNumber, problem)`` for every malformed
    ``!+`` marker in *text*."""
    problems = []
    lines = text.split('\n')
    index = 0
    while index < len(lines):
        if not _FORTRAN_MARKER.match(lines[index]):
            index += 1
            continue
        # Consecutive marker lines form a single marker: only the first carries
        # the label, any others continue its list of names.
        start = index
        while index < len(lines) and _FORTRAN_MARKER.match(lines[index]):
            index += 1
        block = lines[start:index]
        content = _FORTRAN_MARKER.match(block[0]).group(1).strip()
        if not content.startswith(MARKER_LABEL):
            problems.append((path, start + 1,
                             f'does not begin "{MARKER_LABEL}"'))
        for offset, line in enumerate(block):
            names = extract_contributors(line, _FORTRAN_MARKER)
            if not names:
                problems.append((path, start + offset + 1, 'lists no names'))
            for name in names:
                problem = name_problem(name)
                if problem is not None:
                    problems.append((path, start + offset + 1,
                                     f'entry "{name}" {problem}'))
    return problems


def check_markers(galacticus_dir='.'):
    """Return a list of ``(path, lineNumber, problem)`` for every malformed
    ``!+`` marker in the Galacticus source."""
    problems = []
    for file_path in iter_source_files(galacticus_dir):
        with open(file_path, encoding='utf-8', errors='replace') as fh:
            problems.extend(check_markers_in_text(fh.read(), file_path))
    return problems


def render_contributors_rst(names):
    """Render contributor *names* as a single RST paragraph ending in a period."""
    return ', '.join(names) + '.\n'


def _check_main(galacticus_dir):
    problems = check_markers(galacticus_dir)
    for path, line_number, problem in problems:
        print(f'{path}:{line_number}: contributor marker {problem}',
              file=sys.stderr)
    if problems:
        print(f'\n{len(problems)} malformed contributor marker(s) found.\n\n'
              'A `!+` marker records names only, in the form:\n\n'
              f'    !+    {MARKER_LABEL} Andrew Benson, Claude.\n\n'
              'Describe what was contributed - including by any AI tool - in '
              'the commit message,\nnot in the marker: the marker is parsed as '
              'a comma-separated list of names, so any\nother text appears in '
              "the documentation's Acknowledgments as if it were a person.\n"
              'See the "Contributor Attribution" section of CONTRIBUTING.md.',
              file=sys.stderr)
        return 1
    print('All contributor markers are well formed.', file=sys.stderr)
    return 0


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('galacticusDir', nargs='?', default='.',
                        help='the root of the Galacticus source tree '
                             '(default: the current directory)')
    parser.add_argument('outputFile', nargs='?',
                        help='file to which to write the contributor list')
    parser.add_argument('--check', action='store_true',
                        help='do not write; report any contributor marker which '
                             'is not of the canonical, names-only form, exiting '
                             'non-zero if any is found (for CI)')
    arguments = parser.parse_args()

    if arguments.check:
        if arguments.outputFile is not None:
            parser.error('--check takes no output file')
        return _check_main(arguments.galacticusDir)

    if arguments.outputFile is None:
        parser.error('an output file is required (or use --check)')
    names = collect_contributor_names(arguments.galacticusDir)
    with open(arguments.outputFile, 'w', encoding='utf-8') as out:
        out.write(render_contributors_rst(names))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
