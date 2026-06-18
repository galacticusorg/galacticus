#!/usr/bin/env python3
"""One-time migration: convert the hand-written ``doc/*.tex`` manuals to RST.

The four LaTeX manuals are reorganised into three RST guides under
``docs/manuals/`` (committed, then hand-maintained — the ``.tex`` and the PDF
build are retired separately):

* User Guide      — About, Running, Input_Data, Advanced, Python_Interface,
                    Acknowledgments (+ contributor list folded in).
* Developer Guide — Development, Coding, Methods.
* Physics         — Definitions, Components.

Auto-generated includes (autoPhysics/Functions, source_documentation,
contributions, the enumeration/constant/method dumps) are skipped — that content
already lives on the generated Reference pages.  ``contributions`` is folded into
Acknowledgments as a plain contributor list (its per-file links pointed into the
dropped source documentation).

Usage: convertManualsRST.py <docDir> <outDir>   (e.g. doc docs/manuals)
"""
import os
import re
import subprocess
import sys

from latexToRST import glossary_display_map, latex_to_rst, parse_glossary

GUIDES = [
    {'dir': 'user-guide', 'title': 'User Guide',
     'chapters': ['About', 'Running', 'Input_Data', 'Advanced',
                  'Python_Interface', 'Acknowledgments']},
    {'dir': 'developer-guide', 'title': 'Developer Guide',
     'chapters': ['Development', 'Coding', 'Methods']},
    {'dir': 'physics', 'title': 'Physics',
     'chapters': ['Definitions', 'Components']},
]

# Includes that are auto-generated (or non-content) and must not be inlined.
_SKIP_INCLUDES = {
    'autoEnumerationDefinitions', 'constants', 'dataMethods', 'autoPhysics',
    'source_documentation', 'contributions', 'Functions', 'commands',
    'Glossary',
}
_SKIP_RE = re.compile(
    r'\\(?:input|include)\{(' + '|'.join(map(re.escape, _SKIP_INCLUDES))
    + r')\}')


def _page_name(chapter: str) -> str:
    return chapter.lower().replace('_', '-')


def _heading(text: str, char: str) -> str:
    return f'{text}\n{char * max(len(text), 3)}\n'


def _contributors_rst(doc_dir: str, glsmap: dict) -> str:
    """Generate a plain contributor list (names only) from the source markers."""
    gal_dir = os.path.dirname(os.path.abspath(doc_dir))   # repo root (holds source/)
    tmp = os.path.join(gal_dir, '.contributions.tmp.tex')
    here = os.path.dirname(os.path.abspath(__file__))
    try:
        subprocess.run([sys.executable, os.path.join(here, 'extractContributors.py'),
                        gal_dir, tmp], check=True, capture_output=True, cwd=gal_dir)
        raw = open(tmp, encoding='utf-8').read()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return ''
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)
    names = []
    for m in re.findall(r'\\item\[(.*?)\]', raw):
        name = latex_to_rst(m, glsmap).strip()
        # Skip malformed source markers (a few hold a whole sentence rather than
        # a name); real names are a handful of words.
        if not name or len(name.split()) > 5 or name in names:
            continue
        names.append(name)
    if not names:
        return ''
    names.sort(key=str.casefold)
    return ('\n\n' + _heading('Contributors', '-')
            + '\nThe following people have contributed to the Galacticus '
              'project:\n\n'
            + ', '.join(names) + '.\n')


def convert_chapter(doc_dir: str, chapter: str, glsmap: dict) -> str:
    tex = open(os.path.join(doc_dir, f'{chapter}.tex'), encoding='utf-8').read()
    tex = _SKIP_RE.sub('', tex)                 # drop auto-generated includes
    rst = latex_to_rst(tex, glsmap)
    if chapter == 'Acknowledgments':
        rst += _contributors_rst(doc_dir, glsmap)
    return rst.rstrip('\n') + '\n'


def main() -> int:
    if len(sys.argv) != 3:
        print('Usage: convertManualsRST.py <docDir> <outDir>', file=sys.stderr)
        return 1
    doc_dir, out_dir = sys.argv[1], sys.argv[2]
    glsmap = glossary_display_map(parse_glossary(os.path.join(doc_dir, 'Glossary.tex')))

    n_pages = 0
    for guide in GUIDES:
        gdir = os.path.join(out_dir, guide['dir'])
        os.makedirs(gdir, exist_ok=True)
        entries = []
        for chapter in guide['chapters']:
            page = convert_chapter(doc_dir, chapter, glsmap)
            with open(os.path.join(gdir, f'{_page_name(chapter)}.rst'), 'w',
                      encoding='utf-8') as fh:
                fh.write(page)
            entries.append(_page_name(chapter))
            n_pages += 1
        index = (_heading(guide['title'], '=')
                 + '\n.. toctree::\n   :maxdepth: 2\n\n'
                 + '\n'.join(f'   {e}' for e in entries) + '\n')
        with open(os.path.join(gdir, 'index.rst'), 'w', encoding='utf-8') as fh:
            fh.write(index)

    print(f'Wrote {n_pages} manual pages across {len(GUIDES)} guides '
          f'under {out_dir}.', file=sys.stderr)
    return 0


if __name__ == '__main__':
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    main()
