"""Sphinx configuration for the Galacticus documentation on ReadTheDocs.

The physics pages, glossary and references under this directory are generated
at build time by ``scripts/doc/extractDocsRST.py`` (see ``.readthedocs.yaml``)
from the RST docstrings embedded in the Fortran source.
"""
import dataclasses
import datetime

import sphinxcontrib.bibtex.plugin
from sphinxcontrib.bibtex.style.referencing import BracketStyle
from sphinxcontrib.bibtex.style.referencing.author_year import (
    AuthorYearReferenceStyle,
)

project = 'Galacticus'
author = 'Andrew Benson'
copyright = f'2009–{datetime.date.today().year} Andrew Benson'

extensions = [
    'sphinx.ext.mathjax',
    'sphinxcontrib.bibtex',
]


# Author-year citations using round brackets, matching the natbib style of the
# LaTeX manuals: ``\cite``/``\citet`` -> "Author (Year)", ``\citep`` ->
# "(Author Year)" (the default author-year style uses square brackets).
def _round_brackets() -> BracketStyle:
    return BracketStyle(left='(', right=')')


@dataclasses.dataclass
class _RoundAuthorYearStyle(AuthorYearReferenceStyle):
    bracket_parenthetical: BracketStyle = dataclasses.field(
        default_factory=_round_brackets)
    bracket_textual: BracketStyle = dataclasses.field(
        default_factory=_round_brackets)


sphinxcontrib.bibtex.plugin.register_plugin(
    'sphinxcontrib.bibtex.style.referencing',
    'author_year_round', _RoundAuthorYearStyle)

# Bibliography (the same database used by the LaTeX/PDF manuals).
bibtex_bibfiles = ['../doc/Galacticus.bib']
bibtex_default_style = 'plain'
bibtex_reference_style = 'author_year_round'

# Custom math macros mirroring those defined in ``doc/commands.tex`` so the
# embedded equations render identically under MathJax.
mathjax3_config = {
    'tex': {
        'macros': {
            'd':      r'\mathrm{d}',
            'G':      r'\mathrm{G}',
            'clight': r'\mathrm{c}',
            'e':      r'\mathrm{e}',
        },
    },
}

html_theme = 'furo'
html_title = 'Galacticus'

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# The physics pages are generated; do not fail the build on the occasional
# duplicate label that can arise from auto-generated cross-reference anchors.
suppress_warnings = ['epub.unknown_project_files']
