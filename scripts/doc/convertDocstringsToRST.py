#!/usr/bin/env python3
r"""One-time converter: rewrite Galacticus embedded docstrings LaTeX -> RST.

Galacticus documents its physics in two embedded forms inside ``source/*.F90``:

* ``!!{ … !!}``  — free LaTeX comment blocks (module/type/function descriptions).
* ``!![ … !!]``  — XML directive blocks whose ``<description>`` (and nested
  ``<inputParameter><description>`` / ``<defaultSource>``) hold LaTeX.

This script converts that embedded LaTeX to reStructuredText *in place* so the
documentation can be served on ReadTheDocs via Sphinx, while the legacy
LaTeX/PDF pipeline learns to skip the converted content.  Two markers signal
"this is already RST" so the conversion is idempotent and the PDF build can opt
out:

* ``!!{`` openers become ``!!{RST``.
* The root element of any directive block carrying a converted description
  gains a ``docformat="rst"`` attribute.

Only text-mode markup is translated; math is preserved verbatim for MathJax
(see :mod:`latexToRST`).

Usage::

    scripts/doc/convertDocstringsToRST.py [options] [PATH ...]

    PATH         Source files or directories (default: ``source``).
    --dry-run    Do not write; print a unified diff to stdout.
    --report     After processing, list files containing constructs (tables,
                 figures, …) that need manual review.
    --glossary   Path to Glossary.tex (default: ``doc/Glossary.tex``).

Andrew Benson / Galacticus — RST documentation migration (2026).
"""
from __future__ import annotations

import argparse
import difflib
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from latexToRST import (                                       # noqa: E402
    parse_glossary, glossary_display_map, latex_to_rst, find_review_constructs,
)


# ``!!{ … !!}`` LaTeX comment block (skip ones already marked ``!!{RST``).
_LATEX_BLOCK_RE = re.compile(
    r'(?P<indent>[ \t]*)!!\{(?!RST)(?P<body>.*?)!!\}', re.DOTALL)

# ``!![ … !!]`` XML directive block.
_XML_BLOCK_RE = re.compile(
    r'(?P<indent>[ \t]*)!!\[(?P<body>.*?)!!\]', re.DOTALL)

# Inner-text capture for description-bearing elements.
_DESC_RE = re.compile(
    r'(?P<indent>[ \t]*)<(?P<tag>description|defaultSource)>'
    r'(?P<body>.*?)</(?P=tag)>', re.DOTALL)

# First real XML tag (root element) in a directive block.
_ROOT_TAG_RE = re.compile(r'<(?P<name>[A-Za-z][\w]*)(?P<attrs>[^>]*?)(?P<slash>/?)>')


def _reindent(rst: str, indent: str) -> str:
    """Indent every (non-empty) line of ``rst`` by ``indent``."""
    out = []
    for line in rst.split('\n'):
        out.append(indent + line if line.strip() else '')
    return '\n'.join(out)


def _convert_latex_block(m: re.Match, glsmap: dict[str, str]) -> str:
    indent = m.group('indent')
    converted = latex_to_rst(m.group('body'), glsmap)
    body = _reindent(converted, indent)
    return f'{indent}!!{{RST\n{body}\n{indent}!!}}'


# A ``&`` that already begins a valid XML entity reference (``&amp;`` ``&lt;``
# ``&gt;`` ``&quot;`` ``&apos;`` or a numeric ``&#NN;`` / ``&#xHH;``).  These came
# from the original (already valid-XML) source — re-escaping them would
# double-encode, e.g. math alignment ``&`` (written ``&amp;``) would become
# ``&amp;amp;`` and render as a literal "amp;" / "Misplaced &" under MathJax.
_EXISTING_ENTITY_RE = re.compile(r'&(?:amp|lt|gt|quot|apos);|&#[0-9]+;|&#[xX][0-9A-Fa-f]+;')


def _xml_escape(text: str) -> str:
    """Escape RST so it is valid XML text inside a ``<description>`` element.

    The converted descriptions can contain ``<``/``>``/``&`` (RST link targets
    ``<url>``, ``\\mono{…<…>…}`` rendered to plain text, …) that must be
    entity-escaped so the build's directive parser (``ET.fromstring``) still
    parses the block; :mod:`extractDocsRST` unescapes them when reading.

    Entity references already present in the source (math alignment ``&amp;``,
    ``&lt;`` inequalities, ``&#x2F;`` in URLs, …) are left intact so the
    escape/unescape round-trip is lossless.
    """
    # Escape only bare ``&`` (those not already starting an entity).
    out = []
    pos = 0
    for m in _EXISTING_ENTITY_RE.finditer(text):
        out.append(text[pos:m.start()].replace('&', '&amp;'))
        out.append(m.group(0))
        pos = m.end()
    out.append(text[pos:].replace('&', '&amp;'))
    text = ''.join(out)
    # Raw ``<`` / ``>`` (entities contain no raw angle brackets, so untouched).
    return text.replace('<', '&lt;').replace('>', '&gt;')


def _xml_escape_attr(text: str) -> str:
    """Like :func:`_xml_escape`, but also escape ``"`` for an attribute value."""
    return _xml_escape(text).replace('"', '&quot;')


# A ``<methods>`` block of type-bound procedures whose ``<method …/>`` carry a
# (LaTeX) ``description="…"`` attribute rather than a ``<description>`` element.
_METHODS_RE = re.compile(r'(<methods\b[^>]*>)(.*?)(</methods>)', re.DOTALL)


def _convert_method_descriptions(block_body: str,
                                 glsmap: dict[str, str]) -> tuple[str, bool]:
    """Convert the LaTeX ``description="…"`` attributes inside ``<methods>``
    blocks to (single-line) RST, marking the block ``docformat="rst"``.

    Returns ``(new_body, changed)``.
    """
    changed = False

    def methods_repl(mb: re.Match) -> str:
        nonlocal changed
        open_tag, inner, close = mb.group(1), mb.group(2), mb.group(3)
        if 'docformat=' in open_tag:
            return mb.group(0)                         # already converted

        def desc_repl(dm: re.Match) -> str:
            nonlocal changed
            rst = _xml_escape_attr(latex_to_rst(dm.group(1), glsmap))
            rst = re.sub(r'\s+', ' ', rst).strip()
            changed = True
            return f'description="{rst}"'

        new_inner = re.sub(r'description="([^"]*)"', desc_repl, inner)
        return f'{open_tag[:-1]} docformat="rst">{new_inner}{close}'

    return _METHODS_RE.sub(methods_repl, block_body), changed


def _convert_descriptions(block_body: str, glsmap: dict[str, str]) -> tuple[str, bool]:
    """Convert every ``<description>`` / ``<defaultSource>`` in a directive body.

    Returns ``(new_body, changed)``.
    """
    changed = False

    def repl(m: re.Match) -> str:
        nonlocal changed
        inner = m.group('body')
        if not inner.strip():
            return m.group(0)
        indent = m.group('indent')
        tag = m.group('tag')
        converted = _xml_escape(latex_to_rst(inner, glsmap))
        changed = True
        body = _reindent(converted, indent)
        return f'{indent}<{tag}>\n{body}\n{indent}</{tag}>'

    return _DESC_RE.sub(repl, block_body), changed


# Any XML tag (open / close / self-closing).  Raw ``>`` inside descriptions is
# always in text (sources write ``$<$`` / ``&lt;`` for less-than), so this never
# matches inside element text.
_ANY_TAG_RE = re.compile(r'<(/)?([A-Za-z][\w]*)([^<>]*?)(/)?>')


def _top_level_element_spans(block_body: str):
    """Yield ``(start, end, opentag_start, opentag_end)`` for each top-level
    element in a directive block.

    A ``!![ … !!]`` block may contain several sibling directives (e.g. many
    ``<inputParameter>`` blocks); the build's directive parser treats each
    depth-zero element as its own directive, so each must be marked.
    """
    depth = 0
    elem_start = None
    otag = None
    for m in _ANY_TAG_RE.finditer(block_body):
        closing = m.group(1) == '/'
        selfclose = m.group(4) == '/'
        if depth == 0 and not closing:
            if selfclose:
                yield (m.start(), m.end(), m.start(), m.end())
                continue
            elem_start, otag, depth = m.start(), m, 1
        elif not closing and not selfclose:
            depth += 1
        elif closing and depth > 0:
            depth -= 1
            if depth == 0 and otag is not None:
                yield (elem_start, m.end(), otag.start(), otag.end())
                otag = None
    # A top-level element opened but never closed within this block — e.g. a
    # `<workaround>` that brackets code and is closed in a later `!![…!!]`.
    # Still mark its opening tag so its description is recognised as RST.
    if depth > 0 and otag is not None:
        yield (elem_start, len(block_body), otag.start(), otag.end())


def _add_docformat(block_body: str) -> str:
    """Add ``docformat="rst"`` to every top-level directive in the block that
    contains a converted ``<description>`` / ``<defaultSource>``."""
    out = []
    cursor = 0
    for start, end, ots, ote in _top_level_element_spans(block_body):
        out.append(block_body[cursor:start])
        elem = block_body[start:end]
        otag = block_body[ots:ote]
        if (('<description>' in elem or '<defaultSource>' in elem)
                and 'docformat=' not in otag):
            new_otag = (otag[:-2] + ' docformat="rst"/>' if otag.endswith('/>')
                        else otag[:-1] + ' docformat="rst">')
            elem = new_otag + block_body[ote:end]
        out.append(elem)
        cursor = end
    out.append(block_body[cursor:])
    return ''.join(out)


def _convert_xml_block(m: re.Match, glsmap: dict[str, str]) -> str:
    indent = m.group('indent')
    body = m.group('body')
    if 'docformat="rst"' in body:
        return m.group(0)                          # already converted (idempotent)
    has_desc = '<description>' in body or '<defaultSource>' in body
    has_methods = '<methods' in body
    if not has_desc and not has_methods:
        return m.group(0)                          # nothing to convert
    new_body, changed = body, False
    if has_desc:
        new_body, changed = _convert_descriptions(new_body, glsmap)
        if changed:
            new_body = _add_docformat(new_body)
    if has_methods:
        new_body, mchanged = _convert_method_descriptions(new_body, glsmap)
        changed = changed or mchanged
    if not changed:
        return m.group(0)
    return f'{indent}!![{new_body}!!]'


def convert_text(text: str, glsmap: dict[str, str]) -> str:
    """Convert all embedded docstrings in one file's text."""
    text = _LATEX_BLOCK_RE.sub(lambda m: _convert_latex_block(m, glsmap), text)
    text = _XML_BLOCK_RE.sub(lambda m: _convert_xml_block(m, glsmap), text)
    return text


def iter_source_files(paths: list[str]):
    for path in paths:
        if os.path.isdir(path):
            for root, _dirs, files in os.walk(path):
                for fn in sorted(files):
                    if re.search(r'\.(F90|Inc)$', fn) and not fn.startswith('.#'):
                        yield os.path.join(root, fn)
        else:
            yield path


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('paths', nargs='*', default=['source'],
                    help='source files or directories (default: source)')
    ap.add_argument('--dry-run', action='store_true',
                    help='print a unified diff instead of writing')
    ap.add_argument('--check', action='store_true',
                    help='do not write; exit non-zero if any file would change '
                         '(for CI: enforces that all docstrings are already RST)')
    ap.add_argument('--report', action='store_true',
                    help='list files with constructs needing manual review')
    ap.add_argument('--glossary', default='doc/Glossary.tex',
                    help='path to Glossary.tex (default: doc/Glossary.tex)')
    args = ap.parse_args()

    glsmap = glossary_display_map(parse_glossary(args.glossary))

    paths = args.paths or ['source']
    write = not (args.dry_run or args.check)
    n_changed = 0
    changed_files: list[str] = []
    review: dict[str, list[str]] = {}
    for path in iter_source_files(paths):
        with open(path, encoding='utf-8', errors='replace') as fh:
            original = fh.read()
        converted = convert_text(original, glsmap)
        constructs = find_review_constructs(original)
        if constructs:
            review[path] = constructs
        if converted == original:
            continue
        n_changed += 1
        changed_files.append(path)
        if args.dry_run:
            diff = difflib.unified_diff(
                original.splitlines(keepends=True),
                converted.splitlines(keepends=True),
                fromfile=path, tofile=path + ' (RST)')
            sys.stdout.writelines(diff)
        elif write:
            with open(path, 'w', encoding='utf-8') as fh:
                fh.write(converted)
            print(f'converted: {path}')

    verb = 'would be ' if (args.dry_run or args.check) else ''
    print(f'\n{n_changed} file(s) {verb}changed.', file=sys.stderr)
    if args.report and review:
        print(f'\n{len(review)} file(s) contain constructs needing manual '
              f'review:', file=sys.stderr)
        for path in sorted(review):
            print(f'  {path}: {", ".join(sorted(review[path]))}', file=sys.stderr)

    if args.check and changed_files:
        print('\nThe following files contain old-style (LaTeX) docstrings; run '
              '`scripts/doc/convertDocstringsToRST.py source` to convert them:',
              file=sys.stderr)
        for path in changed_files:
            print(f'  {path}', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
