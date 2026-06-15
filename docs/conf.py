"""Sphinx configuration for the Galacticus documentation on ReadTheDocs.

The physics pages, glossary and references under this directory are generated
at build time by ``scripts/doc/extractDocsRST.py`` (see ``.readthedocs.yaml``)
from the RST docstrings embedded in the Fortran source.
"""
import dataclasses
import datetime
import re

from docutils import nodes
from sphinx import addnodes
from sphinx.util.docutils import SphinxRole
from sphinx.util.nodes import make_refnode
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

# Custom math macros (formerly defined in the retired ``doc/commands.tex``) so
# the embedded equations render identically under MathJax.
mathjax3_config = {
    'tex': {
        'macros': {
            'd':      r'\mathrm{d}',
            'G':      r'\mathrm{G}',
            'clight': r'\mathrm{c}',
            'e':      r'\mathrm{e}',
            # \mono is used (rarely) inside equations in the embedded docs.
            'mono':   [r'\texttt{#1}', 1],
        },
    },
}

# Number figures, tables and equations so they can be cross-referenced with
# ``:numref:``/``:eq:`` (the converted docstrings reference them this way).
numfig = True

html_theme = 'furo'
html_title = 'Galacticus'
# The galaxy icon shown at the top of the sidebar on every page (transparent,
# so it works in both light and dark modes).
html_logo = '_static/galacticus-icon.png'
html_favicon = '_static/galacticus-icon.png'

templates_path = ['_templates']
html_static_path = ['_static']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# The physics pages are generated; do not fail the build on the occasional
# duplicate label that can arise from auto-generated cross-reference anchors.
suppress_warnings = ['epub.unknown_project_files']


# --- ``:galacticus-class:`` role -------------------------------------------
# Converted from ``\refClass``/``\refPhysics`` in the docstrings.  Resolves to a
# link to that class's section (anchor ``physics-<name>``, emitted by
# extractDocsRST.py) when one exists, otherwise renders as inline code with no
# warning — so references to abstract/utility classes that have no page degrade
# gracefully.
class _GalacticusClassRole(SphinxRole):
    def run(self):
        name = self.text.strip()
        # Use a private reftype with no domain so Sphinx leaves resolution
        # entirely to the ``missing-reference`` handler below — that keeps the
        # monospace ``<code>`` child (the ``:ref:`` resolver would replace it).
        node = addnodes.pending_xref(
            '', refdomain='', reftype='galacticus-class', reftarget=name,
            refexplicit=True, refwarn=False)
        node += nodes.literal(name, name, classes=['galacticus-class'])
        return [node], []


def _resolve_galacticus_class(app, env, node, contnode):
    if node.get('reftype') != 'galacticus-class':
        return None
    name = node['reftarget']
    # Section labels are stored lower-cased.  Try the class itself, then (for an
    # abstract ``…Class`` base) its family.
    anonlabels = env.get_domain('std').anonlabels
    candidates = [('physics-' + name).lower()]
    if name.endswith('Class'):
        candidates.append(('physics-' + name[:-len('Class')]).lower())
    fromdoc = node.get('refdoc', env.docname)
    for target in candidates:
        if target in anonlabels:
            docname, labelid = anonlabels[target]
            return make_refnode(app.builder, fromdoc, docname, labelid,
                                contnode, name)
    return contnode      # no page for this class: leave the inline code as-is


# --- ``:galacticus-ref:`` role ---------------------------------------------
# Converted from \ref{sec:…} in the manuals.  Resolves to that section (using
# its real title), or — for the few \ref{sec:…} that are dangling in the source
# LaTeX too — degrades to readable text derived from the label, with no warning.
def _decamel(name: str) -> str:
    s = re.sub(r'(?<=[a-z0-9])(?=[A-Z])', ' ', name)
    s = re.sub(r'[\s:_-]+', ' ', s).strip()
    return s[:1].upper() + s[1:] if s else s


class _GalacticusRefRole(SphinxRole):
    def run(self):
        name = self.text.strip()
        slug = re.sub(r'[^A-Za-z0-9]+', '-',
                      'manual-sec-' + name).strip('-').lower()
        display = _decamel(name)
        node = addnodes.pending_xref(
            '', refdomain='', reftype='galacticus-ref', reftarget=slug,
            refexplicit=True, refwarn=False)
        node += nodes.inline(display, display)
        return [node], []


def _resolve_galacticus_ref(app, env, node, contnode):
    if node.get('reftype') != 'galacticus-ref':
        return None
    std = env.get_domain('std')
    target = node['reftarget']
    fromdoc = node.get('refdoc', env.docname)
    if target in std.labels:                      # section target: use its title
        docname, labelid, title = std.labels[target]
        return make_refnode(app.builder, fromdoc, docname, labelid,
                            nodes.inline(title, title), title)
    if target in std.anonlabels:                  # bare target: keep our text
        docname, labelid = std.anonlabels[target]
        return make_refnode(app.builder, fromdoc, docname, labelid, contnode)
    return contnode                               # dangling in source: plain text


def setup(app):
    app.add_role('galacticus-class', _GalacticusClassRole())
    app.add_role('galacticus-ref', _GalacticusRefRole())
    app.connect('missing-reference', _resolve_galacticus_class)
    app.connect('missing-reference', _resolve_galacticus_ref)
    return {'parallel_read_safe': True, 'parallel_write_safe': True}
