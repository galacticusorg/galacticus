#!/usr/bin/env python3
"""Package the built documentation as a single downloadable archive.

Run as ``archiveDocs.py <outputDir>``, where ``<outputDir>`` is the ReadTheDocs
build output directory (``$READTHEDOCS_OUTPUT``).  The HTML site that Sphinx
wrote to ``<outputDir>/html`` is zipped into ``<outputDir>/htmlzip``, from where
ReadTheDocs offers it as the "HTML" download for the version being built — so
every version, and in particular every release, has a self-contained copy of its
documentation that can be kept offline.

This lives in a script rather than inline in ``.readthedocs.yaml`` because the
commands there do not survive nested quoting: the single quotes inside a
``python -c "..."`` command are stripped before the command is run.
"""
import os
import shutil
import sys


def archive(output_directory: str) -> str:
    """Zip the built HTML site, returning the path of the archive written."""
    html    = os.path.join(output_directory, 'html'   )
    htmlzip = os.path.join(output_directory, 'htmlzip')
    if not os.path.isdir(html):
        raise SystemExit(f"archiveDocs.py: no HTML build found in '{html}'")
    os.makedirs(htmlzip, exist_ok=True)
    # ReadTheDocs requires exactly one file in the directory for each format.
    return shutil.make_archive(os.path.join(htmlzip, 'galacticus'), 'zip', html)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise SystemExit('usage: archiveDocs.py <outputDir>')
    print(f'Documentation archived as {archive(sys.argv[1])}')
