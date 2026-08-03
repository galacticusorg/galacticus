#!/usr/bin/env python3
r"""Detect LaTeX left behind inside regions that are supposed to be RST.

:mod:`convertDocstringsToRST` converts embedded LaTeX to RST and stamps each
converted region with a marker (``!!{RST`` / ``docformat="rst"``).  Its
``--check`` mode then asks "would converting change this file?", which is a
guard that a region was *visited* — not that it came out clean: the converter
skips anything already marked, so LaTeX surviving *inside* a marked region is
invisible to it.  This module supplies the complementary check, scanning the
marked regions themselves for constructs that should no longer be there.

Two kinds of region are deliberately **not** scanned, because their content is
still LaTeX by design — :mod:`extractDocsRST` runs it through
:func:`latexToRST.latex_to_rst` at build time:

* ``<entry … description="…"/>``  — enumeration entries (hence the surviving
  ``\gls{…}`` references).
* ``<constant …/>`` attributes    — ``description`` / ``units`` / ``reference``.

``docs/Glossary.tex`` is likewise raw LaTeX by design and is not a source file.

Andrew Benson / Galacticus — RST documentation migration (2026).
Written with assistance from Claude.
"""
from __future__ import annotations

import html
import os
import re
import textwrap


# ===========================================================================
# Region extraction
# ===========================================================================

# A converted free-form docstring: ``!!{RST … !!}``.
_RST_DOCSTRING_RE = re.compile(r'!!\{RST(.*?)!!\}', re.DOTALL)
# An ``!![ … !!]`` XML directive block.
_XML_BLOCK_RE = re.compile(r'^[ \t]*!!\[(.*?)^[ \t]*!!\]', re.DOTALL | re.MULTILINE)
# Description-bearing *elements* (the attribute forms are handled separately so
# that the by-design ``<entry>`` / ``<constant>`` ones can be excluded).
_DESC_ELEMENT_RE = re.compile(
    r'<(description|defaultSource)>(.*?)</\1>', re.DOTALL)
# ``<method … description="…" …/>`` — the only attribute-form description that
# is RST rather than LaTeX.
_METHOD_DESC_RE = re.compile(r'<method\b[^>]*?\bdescription="([^"]*)"')

SOURCE_SUFFIXES = ('.F90', '.Inc')


def iter_source_files(paths):
    """Yield every Fortran source file under ``paths`` (files or directories)."""
    for path in paths:
        if os.path.isdir(path):
            for root, _dirs, files in os.walk(path):
                for fn in sorted(files):
                    if fn.endswith(SOURCE_SUFFIXES) and not fn.startswith('.#'):
                        yield os.path.join(root, fn)
        else:
            yield path


def rst_regions(text: str):
    """Yield ``(offset, body, context)`` for each RST region in ``text``.

    ``offset`` is the index of ``body`` within ``text`` when the body is used
    verbatim.  Bodies taken from XML are entity-unescaped, so their internal
    offsets no longer align with ``text``; those are reported against the start
    of the region, which is accurate to the enclosing element.
    """
    for m in _RST_DOCSTRING_RE.finditer(text):
        yield m.start(1), m.group(1), 'docstring', True

    for bm in _XML_BLOCK_RE.finditer(text):
        block = bm.group(1)
        if 'docformat="rst"' not in block:
            # Unconverted block: ``--check`` already covers it, and flagging it
            # here would double-report the same defect.
            continue
        base = bm.start(1)
        for dm in _DESC_ELEMENT_RE.finditer(block):
            yield base + dm.start(2), html.unescape(dm.group(2)), \
                f'<{dm.group(1)}>', False
        for mm in _METHOD_DESC_RE.finditer(block):
            yield base + mm.start(1), html.unescape(mm.group(1)), \
                '<method description>', False


# ===========================================================================
# Masking
# ===========================================================================

# Regions where a backslash or ``$`` is expected and must not be flagged.
_INLINE_LITERAL_RE = re.compile(r'``.*?``', re.DOTALL)
_ROLE_RE = re.compile(r':[A-Za-z][A-Za-z:\-]*:`[^`]*`')
# OpenMP / OpenACC sentinels quoted in prose (``!$omp atomic read``), whose
# ``$`` is not math.
_SENTINEL_RE = re.compile(r'!\$[A-Za-z]*')
_DIRECTIVE_RE = re.compile(
    r'^(?P<indent>[ \t]*)\.\.[ \t]+(?:math|code-block|code|literalinclude)::')
_LITERAL_BLOCK_RE = re.compile(r'::[ \t]*$')


def _blank(chars, start, end):
    for i in range(start, end):
        if chars[i] != '\n':
            chars[i] = ' '


def mask_verbatim(text: str) -> str:
    """Blank out every span in which LaTeX-looking characters are legitimate.

    Newlines are preserved so that line numbers survive the masking.
    """
    chars = list(text)
    for pattern in (_ROLE_RE, _INLINE_LITERAL_RE, _SENTINEL_RE):
        for m in pattern.finditer(text):
            _blank(chars, *m.span())

    # Explicit markup blocks (``.. math::``, ``.. code-block::``) and literal
    # blocks introduced by a trailing ``::`` own every following line indented
    # more deeply than the introducer.
    lines = text.split('\n')
    offsets, pos = [], 0
    for line in lines:
        offsets.append(pos)
        pos += len(line) + 1

    i = 0
    while i < len(lines):
        line = lines[i]
        if _DIRECTIVE_RE.match(line) or _LITERAL_BLOCK_RE.search(line):
            introducer = len(line) - len(line.lstrip())
            j = i + 1
            while j < len(lines):
                nxt = lines[j]
                if nxt.strip() and (len(nxt) - len(nxt.lstrip())) <= introducer:
                    break
                j += 1
            if j > i + 1:
                _blank(chars, offsets[i + 1],
                       offsets[j - 1] + len(lines[j - 1]))
            i = j
            continue
        i += 1
    return ''.join(chars)


# ===========================================================================
# Detection
# ===========================================================================

# ``{\normalfont \ttfamily …}`` and friends — a LaTeX font group.
_FONT_GROUP_RE = re.compile(
    r'\{\s*\\(?:normalfont|scshape|itshape|bfseries|ttfamily|rmfamily|sffamily|'
    r'boldmath|unboldmath|sc|it|bf|tt|rm|sf|em)\b')
# ``\begin{…}`` / ``\end{…}`` environments.
_ENV_RE = re.compile(r'\\(?:begin|end)\{([A-Za-z*]+)\}')
# Any control word.  RST's own escapes are always ``\`` before a *non*-letter
# (``\ `` for a zero-width separator, ``\*``, ``\|``), so requiring a letter
# here cannot collide with them.
_CONTROL_WORD_RE = re.compile(r'\\([A-Za-z]+)')
# LaTeX character escapes.
_ESCAPE_RE = re.compile(r'\\([_%&#])')
# LaTeX's inter-word space (``e.g.\ foo``).  RST reads ``\ `` as an *escaped*
# space and deletes it, so the words render glued together ("e.g.foo") — the
# one piece of residue that silently corrupts prose rather than showing up as
# visible markup.  It is legitimate, and common, immediately next to inline
# markup, where it is the zero-width separator that lets docutils recognise the
# markup at all (``literal``\ s, cm\ :math:`^2`).
_CONTROL_SPACE_RE = re.compile(r'(?<!\\)\\ ')
_MARKUP_START_RE = re.compile(r'(?:``|:[A-Za-z][A-Za-z:\-]*:`|`|\[[#*][^\]]*\]_|'
                              r'\|[^|]+\||\*\*?\S)')
# Unescaped ``$`` — inline math that was never converted to a ``:math:`` role.
# The paired form is matched first so that the macros *inside* it are reported
# once, as the enclosing ``inline-math``, rather than one finding per macro.
_INLINE_MATH_RE = re.compile(r'(?<!\\)\$.+?(?<!\\)\$', re.DOTALL)
_DOLLAR_RE = re.compile(r'(?<!\\)\$')


class Finding:
    """One residual-LaTeX hit."""

    __slots__ = ('path', 'line', 'context', 'kind', 'snippet')

    def __init__(self, path, line, context, kind, snippet):
        self.path = path
        self.line = line
        self.context = context
        self.kind = kind
        self.snippet = snippet

    @property
    def key(self):
        """Baseline key: deliberately excludes the line number so that the
        baseline does not churn when unrelated edits move code around."""
        return (self.path, self.kind)

    def __repr__(self):                                        # pragma: no cover
        return f'<Finding {self.path}:{self.line} {self.kind}>'


def _snippet(body: str, pos: int, width: int = 60) -> str:
    start = max(0, pos - width // 3)
    return re.sub(r'\s+', ' ', body[start:pos + width]).strip()


def scan_text(text: str, path: str = '<text>'):
    """Return the :class:`Finding` list for one source file's text."""
    findings = []
    for offset, body, context, aligned in rst_regions(text):
        masked = mask_verbatim(body)
        base_line = text.count('\n', 0, offset) + 1

        def add(pos, kind):
            line = base_line + (body.count('\n', 0, pos) if aligned else 0)
            findings.append(Finding(path, line, context, kind,
                                    _snippet(body, pos)))

        # Unconverted inline math is reported as a single finding per ``$…$``
        # span, and the span is then masked so its LaTeX body (``\mathcal``,
        # ``\odot``, …) is not re-reported as loose commands.
        span_chars = list(masked)
        for m in _INLINE_MATH_RE.finditer(masked):
            add(m.start(), 'inline-math')
            _blank(span_chars, *m.span())
        masked = ''.join(span_chars)

        for m in _FONT_GROUP_RE.finditer(masked):
            add(m.start(), 'font-group')
        for m in _ENV_RE.finditer(masked):
            add(m.start(), f'environment:{m.group(1)}')
        for m in _CONTROL_WORD_RE.finditer(masked):
            # ``\begin``/``\end`` are already reported as an environment.
            if m.group(1) in ('begin', 'end'):
                continue
            add(m.start(), f'command:{m.group(1)}')
        for m in _ESCAPE_RE.finditer(masked):
            add(m.start(), f'escape:{m.group(1)}')
        # Adjacency is judged on the *unmasked* body: masking replaces inline
        # markup with spaces, which would hide exactly the neighbours that make
        # a ``\ `` legitimate.
        for m in _CONTROL_SPACE_RE.finditer(masked):
            if body[:m.start()].endswith('`'):
                continue
            if _MARKUP_START_RE.match(body[m.end():]):
                continue
            add(m.start(), 'control-space')
        # Any ``$`` the paired pass did not consume is an unbalanced leftover.
        for m in _DOLLAR_RE.finditer(masked):
            add(m.start(), 'inline-math')
    findings.sort(key=lambda f: (f.line, f.kind))
    return findings


def scan_paths(paths):
    """Scan every source file under ``paths``; returns a flat Finding list."""
    findings = []
    for path in iter_source_files(paths):
        with open(path, encoding='utf-8', errors='replace') as fh:
            findings.extend(scan_text(fh.read(), path))
    return findings


# ===========================================================================
# Baseline
# ===========================================================================

BASELINE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             'residualLatexBaseline.txt')

_BASELINE_HEADER = textwrap.dedent("""\
    # Known residual LaTeX inside regions that should be RST.
    #
    # Generated by:  scripts/doc/convertDocstringsToRST.py --update-baseline
    # Enforced by:   scripts/doc/convertDocstringsToRST.py --check-residual
    #
    # Each line is ``<count>\\t<path>\\t<kind>``.  The check fails when a count
    # rises above its baseline, or when a (path, kind) pair appears that is not
    # listed here — so this file may only ever shrink.  Deleting entries as the
    # migration proceeds is the point; regenerate after fixing a batch.
    #
    # Line numbers are intentionally absent: they would churn on every edit.
    """)


def counts(findings) -> dict:
    """Collapse findings to ``{(path, kind): count}``."""
    out: dict = {}
    for f in findings:
        out[f.key] = out.get(f.key, 0) + 1
    return out


def load_baseline(path: str = BASELINE_PATH) -> dict:
    """Read a baseline file; a missing file is an empty baseline."""
    result: dict = {}
    if not os.path.isfile(path):
        return result
    with open(path, encoding='utf-8') as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line.strip() or line.lstrip().startswith('#'):
                continue
            count, file_path, kind = line.split('\t')
            result[(file_path, kind)] = int(count)
    return result


def write_baseline(findings, path: str = BASELINE_PATH) -> None:
    """Write ``findings`` out as a baseline file."""
    with open(path, 'w', encoding='utf-8') as fh:
        fh.write(_BASELINE_HEADER)
        for (file_path, kind), count in sorted(counts(findings).items()):
            fh.write(f'{count}\t{file_path}\t{kind}\n')


def compare(findings, baseline):
    """Compare findings against a baseline.

    Returns ``(regressions, improvements)`` where ``regressions`` is a list of
    ``(key, allowed, actual)`` that must fail the build, and ``improvements``
    is the same shape for keys that have dropped below (or fallen out of) the
    baseline and could be tightened.
    """
    actual = counts(findings)
    regressions, improvements = [], []
    for key, count in sorted(actual.items()):
        allowed = baseline.get(key, 0)
        if count > allowed:
            regressions.append((key, allowed, count))
        elif count < allowed:
            improvements.append((key, allowed, count))
    for key, allowed in sorted(baseline.items()):
        if key not in actual:
            improvements.append((key, allowed, 0))
    return regressions, improvements
