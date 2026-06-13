#!/usr/bin/env python3
"""Shared LaTeX -> reStructuredText conversion for Galacticus embedded docs.

Two consumers import this module:

* ``convertDocstringsToRST.py`` — the one-time rewriter that converts the
  embedded LaTeX docstrings in ``source/*.F90`` to RST in place.
* ``extractDocsRST.py`` — the build-time extractor, which uses
  :func:`parse_glossary` and :func:`latex_to_rst` to render the glossary page.

The conversion deliberately leaves *math* untouched: inline ``$...$`` becomes an
RST ``:math:`...``` role and display ``\\begin{equation}...`` becomes a
``.. math::`` directive, but the LaTeX *inside* the math is preserved verbatim
for MathJax to render.  Only text-mode markup (``\\mono``, ``\\cite``, ``\\gls``,
lists, …) is translated.

Andrew Benson / Galacticus — RST documentation migration (2026).
"""
from __future__ import annotations

import re
import textwrap
import unicodedata


# ===========================================================================
# Brace / key-value helpers
# ===========================================================================

def extract_braced(text: str, open_pos: int) -> tuple[str, int]:
    """Given ``text`` and the index of an opening ``{``, return ``(inner, end)``
    where ``inner`` is the balanced contents and ``end`` is the index just past
    the matching ``}``.  Honours backslash-escaped braces ``\\{`` / ``\\}``.
    """
    assert text[open_pos] == '{'
    depth = 0
    i = open_pos
    n = len(text)
    while i < n:
        c = text[i]
        if c == '\\' and i + 1 < n:
            i += 2
            continue
        if c == '{':
            depth += 1
        elif c == '}':
            depth -= 1
            if depth == 0:
                return text[open_pos + 1:i], i + 1
        i += 1
    # Unbalanced — return the rest.
    return text[open_pos + 1:], n


def _take_arg(text: str, pos: int) -> tuple[str | None, int]:
    """Read a ``{...}`` argument starting at-or-after ``pos`` (skipping spaces).

    Returns ``(arg, new_pos)``; ``arg`` is ``None`` when no brace group follows.
    """
    j = pos
    while j < len(text) and text[j] in ' \t\n':
        j += 1
    if j < len(text) and text[j] == '{':
        return extract_braced(text, j)
    return None, pos


def _take_optional(text: str, pos: int) -> tuple[str | None, int]:
    """Read an optional ``[...]`` argument at ``pos`` (no space skipping).

    Returns ``(arg, new_pos)`` or ``(None, pos)``.
    """
    if pos < len(text) and text[pos] == '[':
        depth = 0
        i = pos
        while i < len(text):
            if text[i] == '[':
                depth += 1
            elif text[i] == ']':
                depth -= 1
                if depth == 0:
                    return text[pos + 1:i], i + 1
            i += 1
    return None, pos


def _kv_value(body: str, key: str) -> str:
    """Extract the value of ``key=`` from a glossary-entry option ``body``.

    Handles both ``key={braced value}`` and bare ``key=value`` (terminated by a
    top-level comma).  Returns ``''`` when the key is absent.
    """
    m = re.search(r'(?:^|,)\s*' + re.escape(key) + r'\s*=\s*', body)
    if not m:
        return ''
    pos = m.end()
    if pos < len(body) and body[pos] == '{':
        inner, _ = extract_braced(body, pos)
        return inner
    # Bare value: read until the next top-level comma (ignoring braces).
    depth = 0
    i = pos
    while i < len(body):
        c = body[i]
        if c == '{':
            depth += 1
        elif c == '}':
            depth -= 1
        elif c == ',' and depth == 0:
            break
        i += 1
    return body[pos:i].strip()


# ===========================================================================
# Glossary
# ===========================================================================

# LaTeX font-switch wrappers that should collapse to their inner text.
_FONT_GROUP_RE = re.compile(
    r'\{\s*\\(?:normalfont|scshape|itshape|bfseries|ttfamily|rmfamily|sffamily|'
    r'boldmath|unboldmath|sc|it|bf|tt|rm|sf|em)\b\s*([^{}]*)\}'
)


# Bare (brace-less) font-switch declarations, e.g. ``\scshape`` in
# ``{\normalfont \scshape GraphViz}`` or ``\boldmath`` before a math cell.
_FONT_DECL_RE = re.compile(
    r'\\(?:normalfont|scshape|itshape|bfseries|ttfamily|rmfamily|sffamily|'
    r'boldmath|unboldmath|sc|it|bf|tt|rm|sf|em)\b')


def _clean_name(name: str) -> str:
    """Reduce a glossary ``name={...}`` value to plain display text."""
    prev = None
    while prev != name:
        prev = name
        name = _FONT_GROUP_RE.sub(lambda m: m.group(1), name)
    name = _FONT_DECL_RE.sub('', name)
    name = name.replace(r'\&', '&').replace(r'\_', '_').replace(r'\%', '%')
    name = name.replace('{', '').replace('}', '')
    return re.sub(r'\s+', ' ', name).strip()


def parse_glossary(path: str) -> dict[str, dict]:
    """Parse ``doc/Glossary.tex`` into ``{key: entry}``.

    Each ``entry`` is ``{'name', 'description', 'is_acronym'}`` where ``name`` is
    the plain-text display form used for ``\\gls`` references and ``description``
    is the raw LaTeX description (still to be passed through
    :func:`latex_to_rst` by the caller).
    """
    with open(path, encoding='utf-8') as fh:
        text = fh.read()

    entries: dict[str, dict] = {}

    # \newacronym{key}{SHORT}{long form}
    for m in re.finditer(r'\\newacronym\s*\{', text):
        key, p = extract_braced(text, m.end() - 1)
        short, p = _take_arg(text, p)
        long_, p = _take_arg(text, p)
        if short is None:
            continue
        entries[key.strip()] = {
            'name':        _clean_name(short),
            'description': (long_ or '').strip(),
            'is_acronym':  True,
        }

    # \newglossaryentry{key}{ option list }
    for m in re.finditer(r'\\newglossaryentry\s*\{', text):
        key, p = extract_braced(text, m.end() - 1)
        body, p = _take_arg(text, p)
        if body is None:
            continue
        name = _kv_value(body, 'name')
        desc = _kv_value(body, 'description')
        entries[key.strip()] = {
            'name':        _clean_name(name) or key.strip(),
            'description': desc,
            'is_acronym':  r'\acronymtype' in body,
        }

    return entries


def glossary_display_map(glossary: dict[str, dict]) -> dict[str, str]:
    """Map each glossary key to the plain-text term used in ``:term:`` refs."""
    return {key: entry['name'] for key, entry in glossary.items()}


# ===========================================================================
# LaTeX -> RST
# ===========================================================================

# Display-math environments (content preserved verbatim for MathJax).
_DISPLAY_MATH_ENVS = ('equation', 'eqnarray', 'align', 'displaymath', 'multline')
_VERBATIM_ENVS     = ('verbatim', 'lstlisting')

# Simple, argument-less text replacements applied late (outside math/literals).
_SYMBOL_REPLACEMENTS = [
    (r'\textgreater', '>'), (r'\textless', '<'),
    (r'\textbrokenbar', '¦'),
    (r'\textasteriskcentered', '*'),
    (r'\newline', ' '),
    (r'\textbackslash', '\\'),
    (r'\textasciitilde', '~'),
    (r'\LaTeX', 'LaTeX'), (r'\TeX', 'TeX'),
    (r'\AA', 'Å'),
    (r'\ldots', '…'), (r'\dots', '…'),
    (r'\S', 'Section '),
    (r'\copyright', '©'),
    # Standalone letter macros.
    (r'\ss', 'ß'), (r'\ae', 'æ'), (r'\AE', 'Æ'), (r'\oe', 'œ'), (r'\OE', 'Œ'),
    (r'\aa', 'å'), (r'\o', 'ø'), (r'\O', 'Ø'), (r'\l', 'ł'), (r'\L', 'Ł'),
    (r'\i', 'ı'), (r'\j', 'ȷ'),
    (r'\hfill', ''),
    (r'\,', ' '), (r'\;', ' '), (r'\:', ' '), (r'\!', ''),
]

# LaTeX accent commands -> Unicode combining marks (composed via NFC).  Both the
# symbol forms (\"o, \'e, \`a, \^o, \~n, \=o, \.o) and the letter-command forms
# (\c c, \v s, \H o, …), each accepting an optional braced argument.
_ACCENT_COMBINING = {
    '"': '̈', "'": '́', '`': '̀', '^': '̂',
    '~': '̃', '=': '̄', '.': '̇',
    'c': '̧', 'v': '̌', 'u': '̆', 'H': '̋',
    'k': '̨', 'r': '̊', 'd': '̣', 'b': '̱',
}
_ACCENT_SYMBOL_RE = re.compile(r'\\(["\'`^~=.])\s*(?:\{([A-Za-z])\}|([A-Za-z]))')
_ACCENT_LETTER_RE = re.compile(r'\\([cvuHkrdb])(?:\{([A-Za-z])\}|\s+([A-Za-z])\b)')


def _apply_accents(s: str) -> str:
    def repl(m: re.Match) -> str:
        comb = _ACCENT_COMBINING[m.group(1)]
        return unicodedata.normalize('NFC', (m.group(2) or m.group(3)) + comb)
    return _ACCENT_LETTER_RE.sub(repl, _ACCENT_SYMBOL_RE.sub(repl, s))


def _apply_symbols(s: str) -> str:
    """Apply ``_SYMBOL_REPLACEMENTS``.  For control-word macros (alphabetic
    name) also consume the spaces LaTeX swallows after them — so
    ``\\textless fraction`` becomes ``<fraction``, not ``< fraction`` — while
    the trailing ``(?![A-Za-z])`` avoids clipping a longer name (e.g.
    ``\\Sigma``).  A lambda replacement keeps ``dst`` literal (backslashes)."""
    s = _apply_accents(s)
    for src, dst in _SYMBOL_REPLACEMENTS:
        if src[-1].isalpha():
            # ``\cmd{}`` first (the empty group preserves the following space),
            # then ``\cmd`` followed by the spaces LaTeX would swallow.
            s = re.sub(re.escape(src) + r'\{\}', lambda m, d=dst: d, s)
            s = re.sub(re.escape(src) + r'(?![A-Za-z])[ \t]*',
                       lambda m, d=dst: d, s)
        else:
            s = s.replace(src, dst)
    return s

_ESCAPE_REPLACEMENTS = [
    (r'\%', '%'), (r'\_', '_'), (r'\&', '&'), (r'\#', '#'),
    (r'\{', '{'), (r'\}', '}'), (r'\$', '$'),
]


# RST inline markup is only recognised when the start-string is preceded by
# whitespace or one of these characters, and the end-string is followed by
# whitespace or one of these.  When a protected inline fragment lands next to a
# character outside these sets (e.g. ``$x$\mono{y}`` -> two adjacent roles, or
# ``f=$x$``) we splice in an escaped space ``\ `` — a zero-width separator that
# lets docutils recognise the markup without inserting a visible space.
_OK_BEFORE = set(" \t\n-:/'\"<([{")
_OK_AFTER  = set(" \t\n-.,:;!?/'\")]}>")


class _Vault:
    """Holds protected fragments (math, literals, links) behind opaque tokens.

    Inline fragments are spliced back with RST adjacency fix-ups; block
    fragments (display math, code blocks) already carry their own blank-line
    padding and are spliced verbatim.
    """

    def __init__(self) -> None:
        self.items: list[tuple[str, bool]] = []   # (replacement, is_block)

    def stash(self, replacement: str, block: bool = False) -> str:
        token = f'\x00{len(self.items)}\x00'
        self.items.append((replacement, block))
        return token

    def _expand(self) -> list[tuple[str, bool]]:
        """Resolve any tokens nested inside stored fragments."""
        expanded = list(self.items)
        tok_re = re.compile(r'\x00(\d+)\x00')
        for _ in range(5):
            changed = False
            for idx, (val, blk) in enumerate(expanded):
                new = tok_re.sub(lambda m: expanded[int(m.group(1))][0], val)
                if new != val:
                    expanded[idx] = (new, blk)
                    changed = True
            if not changed:
                break
        return expanded

    def restore(self, text: str) -> str:
        expanded = self._expand()
        parts = re.split(r'\x00(\d+)\x00', text)
        seq: list[tuple[str, str]] = []           # (kind, value)
        for i, part in enumerate(parts):
            if i % 2 == 0:
                if part:
                    seq.append(('text', part))
            else:
                val, blk = expanded[int(part)]
                seq.append(('block' if blk else 'inline', val))

        s = ''                                     # output accumulator
        for i, (kind, val) in enumerate(seq):
            if kind == 'text':
                s += val
            elif kind == 'block':
                # Re-indent the block (``.. math::`` / code) to the column where
                # its token sits, so a block inside a list/definition item stays
                # part of that item instead of breaking out at column zero.
                m = re.search(r'(?:^|\n)([ \t]*)$', s)
                indent = m.group(1) if m else ''
                if indent:
                    s = s[:len(s) - len(indent)]
                    val = '\n'.join(
                        (indent + ln if ln.strip() else '')
                        for ln in val.split('\n'))
                s += val
            else:                                  # inline
                prev = s[-1] if s else ''
                if prev and prev not in _OK_BEFORE:
                    s += '\\ '
                s += val
                nxt = seq[i + 1] if i + 1 < len(seq) else ('', '')
                next_char = nxt[1][0] if nxt[1] else ''
                if next_char and next_char not in _OK_AFTER:
                    s += '\\ '
        return s


def _convert_command(text: str, name: str, repl) -> str:
    """Replace every ``\\name{arg}`` with ``repl(arg)`` (brace-balanced)."""
    out = []
    i = 0
    pat = '\\' + name
    while True:
        j = text.find(pat, i)
        if j == -1:
            out.append(text[i:])
            break
        # Ensure the match is the whole command name (not a prefix).
        after = j + len(pat)
        if after < len(text) and (text[after].isalpha()):
            out.append(text[i:after])
            i = after
            continue
        arg, end = _take_arg(text, after)
        if arg is None:
            out.append(text[i:after])
            i = after
            continue
        out.append(text[i:j])
        out.append(repl(arg))
        i = end
    return ''.join(out)


# natbib command -> sphinxcontrib-bibtex role.  With the round-bracket
# author-year style (docs/conf.py) these render as:
#   \cite / \citet      -> "Author (Year)"        (:cite:t:)
#   \citep              -> "(Author Year)"         (:cite:p:)
#   \citealt            -> "Author Year"           (:cite:ts:, no brackets)
#   \citealp            -> "Author, Year"          (:cite:ps:, no brackets)
#   \citeauthor         -> "Author"                (:cite:author:)
#   \citeyear           -> "Year"                  (:cite:year:)
#   \citeyearpar        -> "(Year)"                (:cite:yearpar:)
_CITE_ROLE = {
    'cite': 't', 'citet': 't', 'citep': 'p',
    'citealt': 't', 'citealp': 'p',
    'citeauthor': 'author', 'citeyear': 'year', 'citeyearpar': 'yearpar',
}


def _convert_citations(text: str, vault: '_Vault') -> str:
    r"""Convert natbib ``\cite*`` commands to ``sphinxcontrib-bibtex`` roles.

    The generated roles are stashed in ``vault`` so the adjacency fix-up can
    separate them from neighbouring characters where needed.
    """
    names = sorted(_CITE_ROLE, key=len, reverse=True)
    pat = re.compile(r'\\(' + '|'.join(names) + r')\b')
    out = []
    i = 0
    while True:
        m = pat.search(text, i)
        if not m:
            out.append(text[i:])
            break
        out.append(text[i:m.start()])
        cmd = m.group(1)
        p = m.end()
        # Skip natbib optional arguments [pre][post] / [post].
        _, p = _take_optional(text, p)
        _, p = _take_optional(text, p)
        keys, end = _take_arg(text, p)
        if keys is None:
            out.append(m.group(0))
            i = m.end()
            continue
        key_list = [k.strip() for k in keys.split(',')]
        keys = ','.join(key_list)
        # natbib \citealt is "Author Year" (abbreviated, no brackets).  For a
        # single key, compose :cite:author: + :cite:year: to get exactly that;
        # multi-key \citealt falls back to :cite:t: ("Author (Year)").
        if cmd == 'citealt' and len(key_list) == 1:
            k = key_list[0]
            out.append(vault.stash(f':cite:author:`{k}`') + ' '
                       + vault.stash(f':cite:year:`{k}`'))
        else:
            out.append(vault.stash(f':cite:{_CITE_ROLE[cmd]}:`{keys}`'))
        i = end
    return ''.join(out)


_LIST_ENVS = ('itemize', 'enumerate', 'description')
_ANY_LIST_BEGIN = re.compile(r'\\begin\{(?:itemize|enumerate|description)\}')


def _render_one_list(env: str, body: str) -> str:
    """Render a single (already innermost-flattened) list body as RST.

    ``itemize`` / ``enumerate`` become bullet / numbered lists; ``description``
    becomes an RST definition list (so a ``\\mono`` label can carry inline
    markup, which strong-emphasis ``**…**`` cannot).  The result is padded with
    blank lines so RST recognises the block as separate from its surroundings.
    """
    items = [c.strip() for c in re.split(r'\\item\b', body) if c.strip()]
    lines = []
    for item in items:
        if env == 'description':
            label, p = _take_optional(item, 0)
            text = item[p:].strip() if label is not None else item
            if label is not None:
                # Definition list: term line, then indented definition.  The
                # term may still contain raw markup (e.g. ``\mono{…}``) that
                # later passes convert.
                lines.append(label.strip())
                lines.append(textwrap.indent(text, '   '))
            else:
                lines.append('* ' + textwrap.indent(text, '  ')[2:])
            lines.append('')
            continue
        bullet = '#.' if env == 'enumerate' else '*'
        indented = textwrap.indent(item, '  ')
        indented = indented[2:] if indented.startswith('  ') else indented
        lines.append(f'{bullet} {indented}')
    if env == 'description':
        return '\n\n' + '\n'.join(lines).rstrip('\n') + '\n\n'
    return '\n\n' + '\n'.join(lines) + '\n\n'


def _convert_lists(text: str, glsmap: dict[str, str]) -> str:
    """Convert itemize/enumerate/description environments to RST lists.

    Works innermost-first: an environment whose body contains no further list
    ``\\begin`` is rendered and substituted in place, so that splitting an outer
    body on ``\\item`` never trips over a nested list's items.
    """
    for _ in range(50):                            # bounded; lists never nest deeply
        # Find an innermost list: a \begin whose body has no nested \begin{list}.
        best = None
        for m in re.finditer(r'\\begin\{(itemize|enumerate|description)\}', text):
            env = m.group(1)
            end_m = re.compile(r'\\end\{' + env + r'\}').search(text, m.end())
            if not end_m:
                continue
            inner = text[m.end():end_m.start()]
            if _ANY_LIST_BEGIN.search(inner):
                continue                           # not innermost yet
            best = (m.start(), m.end(), end_m.start(), end_m.end(), env)
            break
        if best is None:
            break
        b0, b1, e0, e1, env = best
        rendered = _render_one_list(env, text[b1:e0])
        # Strip leading spaces the reflow left after ``\end{…}`` (e.g.
        # ``\end{itemize} Here, …``) — otherwise that trailing paragraph is
        # indented one space and docutils renders it as a block quote.
        tail = re.sub(r'^[ \t]+', '', text[e1:])
        text = text[:b0] + rendered + tail
    return text


def _sanitize_label(name: str) -> str:
    """Turn a LaTeX ``\\label`` key (e.g. ``eq:GalliPalla``) into a Sphinx
    cross-reference target (``eq-GalliPalla``).  The ``eq:``/``fig:``/``tb:``/
    ``table:`` prefix is kept so equation/figure/table namespaces stay distinct;
    the original keys are already globally unique."""
    return re.sub(r'[^A-Za-z0-9]+', '-', name.strip()).strip('-')


# Sectioning commands -> a fixed RST underline character per level.  (Heading
# levels in RST are by order of first appearance; a consistent char per level
# keeps the hierarchy stable across a document.)
_SECTION_LEVELS = {
    'chapter': '=', 'section': '-', 'subsection': '~',
    'subsubsection': '^', 'paragraph': '"',
}
_SECTION_RE = re.compile(r'\\(' + '|'.join(_SECTION_LEVELS) + r')\*?\s*\{')


def _convert_sections(text: str, vault: '_Vault', glsmap: dict) -> str:
    """Convert ``\\chapter``/``\\section``/… (with an optional trailing
    ``\\label{…}``) into vaulted RST heading blocks so the prose reflow leaves
    them intact.  The title is converted recursively, then underlined."""
    out, i = [], 0
    while True:
        m = _SECTION_RE.search(text, i)
        if not m:
            out.append(text[i:])
            break
        out.append(text[i:m.start()])
        char = _SECTION_LEVELS[m.group(1)]
        title_raw, p = extract_braced(text, m.end() - 1)
        target = ''
        # The \label may sit after \index/\hyperdef/\hypertarget noise that
        # commonly follows a heading; skip those to find it.
        lm = re.match(
            r'(?:\s*\\(?:index|hyperdef|hypertarget)\b(?:\{[^{}]*\})*)*'
            r'\s*\\label\{([^}]*)\}', text[p:])
        if lm:
            target = '.. _' + _sanitize_label('manual-' + lm.group(1)) + ':\n\n'
            p += lm.end()
        title = re.sub(r'\s+', ' ', latex_to_rst(title_raw, glsmap)).strip()
        block = target + title + '\n' + char * max(len(title), 3)
        out.append(vault.stash('\n\n' + block + '\n\n', block=True))
        i = p
    return ''.join(out)


def _strip_float_envs(text: str, vault: '_Vault', glsmap: dict) -> str:
    """Convert figure/table floats.

    A ``figure`` becomes ``@@FIGURE@@<image path>@@<name>@@`` + caption (finalised
    into a ``.. figure::`` by :func:`_finalise_figures` once the caption is
    converted), or just the caption if it has no ``\\includegraphics``.  A
    ``table`` is unwrapped: the inner list-table (already stashed by
    :func:`_convert_tabular`) gains the float's caption and ``:name:`` so it can
    be ``:numref:``-referenced.
    """
    for env in ('figure', 'figure*'):
        begin_re = re.compile(r'\\begin\{' + re.escape(env) + r'\}(?:\[[^\]]*\])?')
        end_re = re.compile(r'\\end\{' + re.escape(env) + r'\}')
        while True:
            b = begin_re.search(text)
            if not b:
                break
            e = end_re.search(text, b.end())
            if not e:
                break
            inner = text[b.end():e.start()]
            caption = ''
            cm = re.search(r'\\caption\b', inner)
            if cm:
                arg, _ = _take_arg(inner, cm.end())
                caption = (arg or '').strip()
            lm = re.search(r'\\label\{([^}]*)\}', inner)
            name = _sanitize_label(lm.group(1)) if lm else ''
            image = ''
            im = re.search(r'\\includegraphics(?:\[[^\]]*\])?\s*\{', inner)
            if im:
                image, _ = extract_braced(inner, im.end() - 1)
            vm = re.search(r'\\begin\{verbatim\}(.*?)\\end\{verbatim\}',
                           inner, re.DOTALL)
            if image.strip():
                replacement = (f'\n\n@@FIGURE@@{image.strip()}@@{name}@@'
                               f'\n\n{caption}\n\n')
            elif vm and name:
                # An imageless figure wrapping a verbatim block (e.g. an ASCII
                # diagram): render as a captioned, named code-block so it is
                # numbered and ``:numref:``-referenceable.
                cap = _conv_cell(caption, glsmap) if caption else ''
                block = ('.. code-block:: none\n   :name: ' + name
                         + (('\n   :caption: ' + cap) if cap else '')
                         + '\n\n' + textwrap.indent(vm.group(1).strip('\n'), '   '))
                replacement = vault.stash('\n\n' + block + '\n\n', block=True)
            else:
                replacement = f'\n\n{caption}\n\n' if caption else '\n\n'
            text = text[:b.start()] + replacement + text[e.end():]

    # Table floats: fold the caption + ``:name:`` into the stashed list-table.
    def _table_repl(m: re.Match) -> str:
        inner = m.group(1)
        cap = ''
        cm = re.search(r'\\caption\b', inner)
        if cm:
            arg, _ = _take_arg(inner, cm.end())
            cap = _conv_cell(arg or '', glsmap)
        lm = re.search(r'\\label\{([^}]*)\}', inner)
        name = _sanitize_label(lm.group(1)) if lm else ''
        tm = re.search('\x00(\\d+)\x00', inner)
        if tm:
            idx = int(tm.group(1))
            block, is_block = vault.items[idx]
            head = '.. list-table::'
            new_head = head + ((' ' + cap) if cap else '')
            if name:
                new_head += '\n   :name: ' + name
            if new_head != head:
                vault.items[idx] = (block.replace(head, new_head, 1), is_block)
            return '\n\n\x00' + str(idx) + '\x00\n\n'
        return ('\n\n' + cap + '\n\n') if cap else '\n\n'

    text = re.sub(r'\\begin\{table\*?\}(?:\[[^\]]*\])?(.*?)\\end\{table\*?\}',
                  _table_repl, text, flags=re.DOTALL)
    text = _convert_command(text, 'caption', lambda a: '\n\n' + a.strip() + '\n\n')
    return text


def _finalise_figures(text: str) -> str:
    """Turn ``@@FIGURE@@path@@`` + caption placeholders (caption already
    converted) into ``.. figure::`` directives."""
    def _fig(path: str, name: str) -> str:
        body = '.. figure:: ' + path.strip()
        if name.strip():
            body += '\n   :name: ' + name.strip()
        return body

    def with_caption(m: re.Match) -> str:
        body = _fig(m.group(1), m.group(2))
        cap = m.group(3).strip()
        if cap:
            body += '\n\n' + textwrap.indent(cap, '   ')
        return body

    text = re.sub(r'@@FIGURE@@(.+?)@@(.*?)@@\n\n([^\n]+)', with_caption, text)
    text = re.sub(r'@@FIGURE@@(.+?)@@(.*?)@@',
                  lambda m: _fig(m.group(1), m.group(2)), text)
    return text


def _conv_cell(text: str, glsmap: dict) -> str:
    """Convert one tabular cell to inline RST by running the full pipeline on
    the fragment (handles ``\\mono``, ``\\gls``, ``$math$``, font/boldmath …)."""
    return re.sub(r'\s+', ' ', latex_to_rst(text, glsmap)).strip()


def _parse_colspec(spec: str):
    """Parse a tabular column spec into ``(n_logical, phys_to_logical,
    phys_sep)``.  Physical columns joined by ``@{sep}`` (e.g. the ``r@{.}l`` of
    a decimal-aligned number) collapse to one logical column."""
    seps, pending, i = [], '', 0
    while i < len(spec):
        ch = spec[i]
        if ch in 'lcrLCR':
            seps.append(pending); pending = ''; i += 1
        elif ch in 'pmb' and i + 1 < len(spec) and spec[i + 1] == '{':
            _, i = extract_braced(spec, i + 1)
            seps.append(pending); pending = ''
        elif ch == '@' and i + 1 < len(spec) and spec[i + 1] == '{':
            pending, i = extract_braced(spec, i + 1)
        else:
            i += 1
    phys_to_logical, phys_sep, logical = [], [], -1
    for sep in seps:
        if sep == '':
            logical += 1
        phys_to_logical.append(logical)
        phys_sep.append(sep)
    return logical + 1, phys_to_logical, phys_sep


_MULTICOL_RE = re.compile(r'^\s*\\multicolumn\s*\{(\d+)\}\s*\{[^{}]*\}\s*\{')


def _row_cells(row, n_logical, phys_to_logical, phys_sep, glsmap):
    logical = [''] * n_logical
    phys = 0
    for raw in row.split('&amp;'):
        mc = _MULTICOL_RE.match(raw)
        if mc:
            content, _ = extract_braced(raw, mc.end() - 1)
            span = int(mc.group(1))
        else:
            content, span = raw, 1
        lc = phys_to_logical[phys] if phys < len(phys_to_logical) else n_logical - 1
        conv = _conv_cell(content, glsmap)
        if logical[lc]:
            sep = phys_sep[phys] if phys < len(phys_sep) and phys_sep[phys] else ' '
            logical[lc] += sep + conv
        elif conv:
            logical[lc] = conv
        phys += span
    return logical


def _convert_tabular(text: str, vault: '_Vault', glsmap: dict) -> str:
    """Convert ``\\begin{tabular}{spec} … \\end{tabular}`` to an RST list-table.

    Handles ``@{.}`` decimal-aligned columns, ``\\multicolumn`` spans, multi-row
    headers (the first ``\\hline``-delimited group), and per-cell markup.
    Vaulted as a block so the reflow/list passes leave it intact.
    """
    out, i = [], 0
    while True:
        b = text.find(r'\begin{tabular}', i)
        if b == -1:
            out.append(text[i:])
            break
        j = b + len(r'\begin{tabular}')
        while j < len(text) and text[j] in ' \t\n':
            j += 1
        spec, body_start = (extract_braced(text, j) if j < len(text)
                            and text[j] == '{' else ('', j))
        e = text.find(r'\end{tabular}', body_start)
        if e == -1:
            out.append(text[i:])
            break
        out.append(text[i:b])
        out.append(_build_listtable(spec, text[body_start:e], vault, glsmap))
        i = e + len(r'\end{tabular}')
    return ''.join(out)


def _build_listtable(spec, body, vault, glsmap):
    n_logical, p2l, psep = _parse_colspec(spec)
    if n_logical == 0:
        n_logical, p2l, psep = 1, [0], ['']
    segments = re.split(r'\\hline', body)
    seg_rows = [[r.strip() for r in re.split(r'\\\\', s) if r.strip()]
                for s in segments]
    nonempty = [rs for rs in seg_rows if rs]
    if not nonempty:
        return ''
    rows = [_row_cells(r, n_logical, p2l, psep, glsmap)
            for rs in nonempty for r in rs]
    # Drop trailing columns empty in every row (e.g. an unused decimal column).
    while n_logical > 1 and all(not r[n_logical - 1] for r in rows):
        n_logical -= 1
        rows = [r[:n_logical] for r in rows]
    header = len(nonempty[0]) if len(nonempty) >= 2 else 0
    lines = ['.. list-table::']
    if header:
        lines.append(f'   :header-rows: {header}')
    lines.append('')
    for row in rows:
        for k in range(n_logical):
            cell = row[k] if k < len(row) else ''
            lines.append(('   * - ' if k == 0 else '     - ')
                         + (cell or '​'))      # zero-width space for empties
    return vault.stash('\n\n' + '\n'.join(lines) + '\n\n', block=True)


# One step of a ``\noindent\hspace{N mm} $\rightarrow$ \parbox[..]{W}{TEXT}``
# decision tree (the ``{TEXT}`` opening brace is captured by extract_braced).
_ARROW_RE = re.compile(
    r'\\noindent\s*\\hspace\*?\s*\{\s*(\d+)\s*mm\s*\}\s*'
    r'\$\\rightarrow\$\s*\\parbox\s*(?:\[[^\]]*\])?\s*\{[^{}]*\}\s*\{')


def _arrows_to_itemize(text: str) -> str:
    """Rewrite a run of indented ``→`` / ``\\parbox`` steps (indentation set by
    ``\\hspace{N mm}``) into nested ``\\begin{itemize}`` so the normal list
    machinery renders the nesting."""
    items = []
    i = 0
    while True:
        m = _ARROW_RE.search(text, i)
        if not m:
            break
        body, end = extract_braced(text, m.end() - 1)
        items.append((int(m.group(1)), body.strip(), m.start(), end))
        i = end
    if not items:
        return text
    lines, stack = [], []
    for depth, body, _s, _e in items:
        while stack and stack[-1] > depth:
            lines.append('\\end{itemize}')
            stack.pop()
        if not stack or stack[-1] < depth:
            lines.append('\\begin{itemize}')
            stack.append(depth)
        lines.append('\\item ' + body)
    while stack:
        lines.append('\\end{itemize}')
        stack.pop()
    start, finish = items[0][2], items[-1][3]
    tail = re.sub(r'^\s*\\\\', '', text[finish:])          # drop trailing row break
    return text[:start] + '\n'.join(lines) + '\n' + tail


def _strip_parbox(text: str) -> str:
    """Replace ``\\parbox[opt]{width}{text}`` with just ``text``."""
    out, i = [], 0
    while True:
        j = text.find('\\parbox', i)
        if j == -1:
            out.append(text[i:])
            break
        p = j + len('\\parbox')
        _, p = _take_optional(text, p)
        _, p = _take_arg(text, p)                          # width
        body, p2 = _take_arg(text, p)                      # text
        out.append(text[i:j])
        if body is not None:
            out.append(body)
            i = p2
        else:
            out.append('\\parbox')
            i = j + len('\\parbox')
    return ''.join(out)


def latex_to_rst(text: str, glsmap: dict[str, str] | None = None) -> str:
    """Convert one LaTeX description string to reStructuredText.

    ``glsmap`` maps glossary keys to display terms (see
    :func:`glossary_display_map`).  Math is preserved verbatim; only text-mode
    markup is translated.  The result is dedented to column zero so callers can
    re-indent it as needed.
    """
    glsmap = glsmap or {}
    text = textwrap.dedent(text.replace('\r\n', '\n')).strip('\n')
    vault = _Vault()

    # --- Sectioning -> RST heading blocks (vaulted) ---------------------
    text = _convert_sections(text, vault, glsmap)

    # --- Footnotes: capture bodies now (so they are not processed inline);
    # the definitions are converted and appended at the very end. -----------
    footnotes: list[str] = []

    def _footnote_repl(text_in: str) -> str:
        out, i = [], 0
        while True:
            fm = re.compile(r'\\footnote\s*\{').search(text_in, i)
            if not fm:
                out.append(text_in[i:])
                break
            body, p = extract_braced(text_in, fm.end() - 1)
            out.append(text_in[i:fm.start()])
            footnotes.append(body)
            out.append(vault.stash('[#]_'))
            i = p
        return ''.join(out)

    text = _footnote_repl(text)

    # --- tabular -> list-table (vaulted block) ---------------------------
    # Before float handling, so a tabular inside a ``table`` float survives.
    text = _convert_tabular(text, vault, glsmap)

    # --- Float environments ----------------------------------------------
    # figure -> ".. figure::" (image + caption); table -> unwrapped, keeping the
    # list-table above plus its caption.
    text = _strip_float_envs(text, vault, glsmap)

    # --- \noindent\hspace{N} $\rightarrow$ \parbox{…} decision trees -----
    # Rewrite into nested \begin{itemize} so the normal list machinery renders
    # the nesting; then drop any leftover layout-only commands.
    text = _arrows_to_itemize(text)
    text = re.sub(r'\\hyperdef(?:\{[^{}]*\}){3}', '', text)
    text = re.sub(r'\\noindent\b', '', text)
    text = re.sub(r'\\hspace\*?\s*\{[^{}]*\}', '', text)
    text = _strip_parbox(text)                             # \parbox[..]{W}{T} -> T

    # --- Protect display math --------------------------------------------
    for env in _DISPLAY_MATH_ENVS:
        env_re = re.compile(
            r'\\begin\{' + env + r'\}(.*?)\\end\{' + env + r'\}', re.DOTALL)
        is_eqnarray = env.startswith('eqnarray')

        def repl_disp(m: re.Match, is_eqnarray=is_eqnarray) -> str:
            # Sphinx wraps a multi-row body (containing ``\\``) in ``split`` and
            # a single-row body in plain ``\[ … \]``.  ``eqnarray`` is ``rcl``
            # (two ``&`` per row: ``lhs &op& rhs``) which, left as-is, right-
            # aligns the RHS — and a single-row ``&`` has no alignment env at
            # all ("Misplaced &").  Convert each row to a single ``&`` so
            # ``split`` left-aligns the RHS next to the relation, and drop the
            # ``&`` entirely from a single-row equation.
            inner = textwrap.dedent(m.group(1).strip('\n')).strip('\n')
            # Pull out \label{…} -> directive ``:label:`` so the equation can be
            # referenced with ``:eq:``.  Keep the first label for the block.
            eq_label = None
            eq_labels = re.findall(r'\\label\{([^}]*)\}', inner)
            if eq_labels:
                eq_label = _sanitize_label(eq_labels[0])
                inner = re.sub(r'\\label\{[^}]*\}', '', inner)
            if is_eqnarray:
                # Alignment ampersands are XML-escaped as ``&amp;`` in the
                # source description (``&lt;``/``&gt;`` inequalities are left
                # untouched).
                if re.search(r'\\\\', inner):                  # multi-row
                    rows = []
                    for row in re.split(r'\\\\', inner):
                        parts = row.split('&amp;')
                        if len(parts) > 1:
                            row = (parts[0].rstrip() + ' &amp; '
                                   + ' '.join(p.strip() for p in parts[1:]))
                        rows.append(row.strip())
                    inner = ' \\\\\n'.join(rows)
                else:                                          # single row
                    inner = inner.replace('&amp;', ' ')
            directive = '.. math::'
            if eq_label:
                directive += '\n   :label: ' + eq_label
            block = directive + '\n\n' + textwrap.indent(inner, '   ')
            return vault.stash('\n\n' + block + '\n\n', block=True)

        text = env_re.sub(repl_disp, text)

    # --- Protect verbatim / listings -------------------------------------
    for env in _VERBATIM_ENVS:
        env_re = re.compile(
            r'\\begin\{' + env + r'\}(?:\[[^\]]*\])?(.*?)\\end\{' + env + r'\}',
            re.DOTALL)

        def repl_verb(m: re.Match) -> str:
            body = m.group(1).strip('\n')
            block = '.. code-block:: none\n\n' + textwrap.indent(body, '   ')
            return vault.stash('\n\n' + block + '\n\n', block=True)

        text = env_re.sub(repl_verb, text)

    # \lstset{…} is listings configuration with no rendered output.
    text = _convert_command(text, 'lstset', lambda a: '')
    # Drop \index{…} now (brace-balanced, so a nested \mono in an index key like
    # ``\index{foo@\mono{foo}}`` goes too) — before \mono is turned into a
    # literal that the later, single-level drop could not reach into.
    text = _convert_command(text, 'index', lambda a: '')
    # \verb<delim>code<delim> -> inline literal (delimiter is any char).
    text = re.sub(
        r'\\verb\*?(\S)(.*?)\1',
        lambda m: vault.stash('``' + m.group(2).replace('`', "'") + '``'),
        text)

    # --- Protect inline math ---------------------------------------------
    def repl_inline(m: re.Match) -> str:
        # Collapse internal whitespace so a ``$…$`` that wrapped across source
        # lines becomes a single-line ``:math:`…``` role (newlines inside an
        # inline role confuse docutils).
        content = re.sub(r'\s+', ' ', m.group(1)).strip()
        return vault.stash(':math:`' + content + '`')

    text = re.sub(r'(?<!\\)\$(.+?)(?<!\\)\$', repl_inline, text, flags=re.DOTALL)

    # --- Reflow prose -----------------------------------------------------
    # First isolate block-level protected fragments (display math, code) onto
    # their own paragraph, stripping the surrounding spaces — otherwise the
    # reflow below would leave a stray leading space on the line after a math
    # block, which docutils would treat as a continuation of the ``.. math::``
    # directive (any indent > 0 after a column-0 directive is its content).
    def _isolate_blocks(m: re.Match) -> str:
        idx = int(m.group(1))
        if vault.items[idx][1]:                    # is_block
            return f'\n\n\x00{idx}\x00\n\n'
        return m.group(0)

    text = re.sub(r'[ \t]*\x00(\d+)\x00[ \t]*', _isolate_blocks, text)

    # Collapse intra-paragraph newlines + source indentation to single spaces,
    # preserving blank-line paragraph breaks.  Remaining (inline) protected
    # fragments are opaque tokens, so they are untouched; this gives list/quote
    # conversion clean, indentation-free input to work with.
    paragraphs = re.split(r'\n[ \t]*\n', text)
    text = '\n\n'.join(re.sub(r'\s*\n\s*', ' ', p).strip() for p in paragraphs)

    # --- Hyperlinks (before literals: a label may be ``\mono{…}``) --------
    def _clean_link_label(lbl: str) -> str:
        prev = None
        while prev != lbl:
            prev = lbl
            lbl = _FONT_GROUP_RE.sub(lambda m: m.group(1), lbl)
        lbl = _FONT_DECL_RE.sub('', lbl)
        lbl = re.sub(r'\\[a-zA-Z]+\*?\s*\{([^{}]*)\}', r'\1', lbl)
        for src, dst in _ESCAPE_REPLACEMENTS:
            lbl = lbl.replace(src, dst)
        return re.sub(r'\s+', ' ', lbl).strip()

    def href_repl(text_in: str) -> str:
        out = []
        i = 0
        while True:
            j = text_in.find('\\href', i)
            if j == -1:
                out.append(text_in[i:])
                break
            url, p = _take_arg(text_in, j + len('\\href'))
            label, q = _take_arg(text_in, p)
            # Tolerate the non-standard ``\href{url}\mono{label}`` form, where
            # the link text is a macro argument rather than a brace group.
            if label is None and url is not None:
                mm = re.match(r'\s*\\[a-zA-Z]+\*?\s*\{', text_in[p:])
                if mm:
                    label, q = extract_braced(text_in, p + mm.end() - 1)
            out.append(text_in[i:j])
            if url is None or label is None:
                out.append(text_in[j:j + len('\\href')])
                i = j + len('\\href')
                continue
            out.append(vault.stash(
                f'`{_clean_link_label(label)} <{url.strip()}>`_'))
            i = q
        return ''.join(out)

    text = href_repl(text)
    text = _convert_command(
        text, 'url', lambda a: vault.stash(f'`{a.strip()} <{a.strip()}>`_'))

    # --- Lists ------------------------------------------------------------
    text = _convert_lists(text, glsmap)

    # --- Drop layout-only environments, keep their content ---------------
    for env in ('center', 'flushleft', 'flushright', 'quotation', 'quote'):
        text = re.sub(r'\\begin\{' + env + r'\}', '', text)
        text = re.sub(r'\\end\{' + env + r'\}', '', text)

    # --- Inline literal / code-like commands -> ``...`` ------------------
    # Stash the generated literals so later quote/tilde fix-ups cannot disturb
    # their double-backtick delimiters or contents.  Inner markup (a nested
    # ``\gls``/``\refClass``, inline ``$…$`` math, a ``\_`` escape, …) is
    # flattened to plain text because RST inline literals cannot contain
    # further markup.
    def _math_to_text(mth: str) -> str:
        mth = mth.replace('&lt;', '<').replace('&gt;', '>').replace('&amp;', '&')
        mth = mth.replace(r'\langle', '<').replace(r'\rangle', '>')
        mth = re.sub(r'\\mathrm\{([^{}]*)\}', r'\1', mth)
        mth = re.sub(r'\\[a-zA-Z]+', '', mth)
        return mth.replace('{', '').replace('}', '').strip()

    def lit(a: str) -> str:
        # Resolve protected inline-math tokens to plain text (a literal cannot
        # contain a ``:math:`` role).
        def _resolve(mt: re.Match) -> str:
            val = vault.items[int(mt.group(1))][0]
            if val.startswith(':math:`') and val.endswith('`'):
                return _math_to_text(val[len(':math:`'):-1])
            return val
        a = re.sub(r'\x00(\d+)\x00', _resolve, a)
        a = re.sub(r'\\[a-zA-Z]+\*?\s*\{([^{}]*)\}', r'\1', a)
        a = _FONT_DECL_RE.sub('', a)
        # Resolve bare symbol macros (e.g. \textless inside \mono{<fraction>}),
        # which the brace-flatten above leaves untouched.
        a = _apply_symbols(a)
        for src, dst in _ESCAPE_REPLACEMENTS:
            a = a.replace(src, dst)
        return vault.stash('``' + a.replace('`', "'").strip() + '``')

    for cmd in ('mono', 'texttt', 'refMethod', 'path'):
        text = _convert_command(text, cmd, lit)

    # --- Emphasis ---------------------------------------------------------
    for cmd in ('emph', 'textit', 'textsl'):
        text = _convert_command(text, cmd, lambda a: '*' + a.strip() + '*')
    for cmd in ('textbf', 'textsc'):
        text = _convert_command(text, cmd, lambda a: '**' + a.strip() + '**')
    for cmd in ('textrm', 'textnormal', 'mbox', 'text'):
        text = _convert_command(text, cmd, lambda a: a)

    # --- Glossary ---------------------------------------------------------
    # Stash the generated roles so the adjacency fix-up inserts ``\ `` when a
    # term abuts a word character (e.g. ``\glspl{imf}`` -> ``:term:`IMF`\ s``).
    def gls_repl(key: str) -> str:
        term = glsmap.get(key.strip())
        if term is None:
            return key.strip()
        return vault.stash(f':term:`{term}`')

    for cmd in ('glsplural', 'glspl', 'Glspl', 'gls', 'Gls', 'acrshort',
                'acrlong', 'acrfull'):
        text = _convert_command(text, cmd, gls_repl)
    # \glslink{target}{shown text} -> keep shown text.
    def glslink_repl(text_in: str) -> str:
        out = []
        i = 0
        while True:
            j = text_in.find('\\glslink', i)
            if j == -1:
                out.append(text_in[i:])
                break
            _, p = _take_arg(text_in, j + len('\\glslink'))
            shown, q = _take_arg(text_in, p)
            out.append(text_in[i:j])
            out.append(shown if shown is not None else '')
            i = q if shown is not None else j + len('\\glslink')
        return ''.join(out)

    text = glslink_repl(text)
    # \hyperlink{target}{shown text} -> keep the shown text (the PDF-internal
    # targets, e.g. into the dropped source documentation, do not exist on RTD).
    def hyperlink_repl(text_in: str) -> str:
        out, i = [], 0
        while True:
            j = text_in.find('\\hyperlink', i)
            if j == -1:
                out.append(text_in[i:])
                break
            _, p = _take_arg(text_in, j + len('\\hyperlink'))
            shown, q = _take_arg(text_in, p)
            out.append(text_in[i:j])
            out.append(shown if shown is not None else '')
            i = q if shown is not None else j + len('\\hyperlink')
        return ''.join(out)

    text = hyperlink_repl(text)

    # --- Cross references ------------------------------------------------
    # \refClass/\refPhysics -> a custom role resolved at build time to a link to
    # that class's page (or rendered as inline code when no page exists).
    def refclass_repl(a: str) -> str:
        name = re.sub(r'\\[a-zA-Z]+\*?\s*\{([^{}]*)\}', r'\1', a)
        name = name.replace('{', '').replace('}', '').strip()
        return vault.stash(':galacticus-class:`' + name + '`')

    for cmd in ('refClass', 'refPhysics'):
        text = _convert_command(text, cmd, refclass_repl)

    # \ref{eq:…} -> ``:eq:`` (equation number); \ref{fig:…}/\ref{tb:…}/
    # \ref{table:…} -> ``:numref:`` ("Fig. N"/"Table N"); \ref{sec:…} -> ``:ref:``
    # to the manual section anchor; others are dropped.
    def ref_repl(a: str) -> str:
        # The source prose already supplies the leading word ("Table~\ref{…}",
        # "Figure~\ref{…}"), so ``:numref:`` emits just the number ("{number}")
        # to mirror LaTeX ``\ref``; equations use ``:eq:`` ("(N)") after "eq.".
        key = a.strip()
        if key.startswith('eq:'):
            return vault.stash(':eq:`' + _sanitize_label(key) + '`')
        if re.match(r'(?:fig|tb|table):', key):
            return vault.stash(':numref:`{number} <' + _sanitize_label(key) + '>`')
        if key.startswith('sec:'):
            # Custom role: resolves to the manual section (anchor manual-sec-…),
            # else degrades to readable plain text — some source \ref{sec:…} are
            # dangling in the original LaTeX too.
            return vault.stash(':galacticus-ref:`' + key[4:] + '`')
        return ''

    # A leftover \label (one not consumed by a heading or a math/figure env);
    # keep ``sec:`` ones as targets so cross-references still resolve.
    def label_repl(a: str) -> str:
        key = a.strip()
        if key.startswith('sec:'):
            return vault.stash(
                '\n\n.. _' + _sanitize_label('manual-' + key) + ':\n\n',
                block=True)
        return ''

    text = _convert_command(text, 'label', label_repl)
    text = _convert_command(text, 'ref', ref_repl)
    text = _convert_command(text, 'pageref', lambda a: '')
    text = _convert_command(text, 'index', lambda a: '')
    # Drop \includegraphics — figures are handled separately/by review.
    text = _convert_command(text, 'includegraphics', lambda a: '')

    # --- Citations (stashed for adjacency fix-up) ------------------------
    text = re.sub(r'\\protect\b\s*', '', text)     # \protect\cite{…} -> \cite{…}
    text = _convert_citations(text, vault)

    # --- Product name and font groups ------------------------------------
    text = re.sub(r'\\glc\\?(?:\s|\{\})', 'Galacticus ', text)
    text = text.replace(r'\glc', 'Galacticus')
    prev = None
    while prev != text:
        prev = text
        text = _FONT_GROUP_RE.sub(lambda m: m.group(1), text)
    text = _FONT_DECL_RE.sub('', text)

    # --- Symbols ---------------------------------------------------------
    # Before brace-stripping so a macro's empty group (``\LaTeX{} at`` — used to
    # protect the following space) is handled, not silently dropped.
    text = _apply_symbols(text)

    # --- Strip leftover LaTeX grouping braces ----------------------------
    # Any ``{`` / ``}`` still present are LaTeX grouping (e.g. ``{\glspl{x}}``);
    # they have no RST meaning and would break adjacent roles.  Escaped braces
    # (``\{`` / ``\}``) are preserved here and turned into literals just below.
    text = re.sub(r'(?<!\\)[{}]', '', text)

    # --- Escapes ---------------------------------------------------------
    # LaTeX quotes -> straight quotes.
    text = text.replace('``', '"').replace("''", '"')
    # Non-breaking space and explicit spaces.
    text = text.replace('~', ' ')
    for src, dst in _ESCAPE_REPLACEMENTS:
        text = text.replace(src, dst)

    # --- Restore protected fragments -------------------------------------
    text = vault.restore(text)

    # --- Finalise figures (caption now converted) ------------------------
    text = _finalise_figures(text)

    # --- Footnote definitions (auto-numbered, matched to ``[#]_`` in order) -
    if footnotes:
        defs = '\n'.join(
            '.. [#] ' + re.sub(r'\s*\n\s*', ' ',
                               latex_to_rst(b, glsmap)).strip()
            for b in footnotes)
        text = text.rstrip('\n') + '\n\n' + defs + '\n'

    # --- Tidy whitespace --------------------------------------------------
    text = re.sub(r'[ \t]+\n', '\n', text)
    text = re.sub(r'\n{3,}', '\n\n', text)
    return text.strip('\n')


# Constructs this converter does not translate; callers may warn on these.
# (``array`` / ``cases`` are math-internal and rendered by MathJax, so they are
# intentionally not flagged here.)
REVIEW_PATTERNS = {
    'tabular':  re.compile(r'\\begin\{tabular\}'),
    'table':    re.compile(r'\\begin\{table\}'),
    'figure':   re.compile(r'\\begin\{figure\}'),
    'tabbing':  re.compile(r'\\begin\{tabbing\}'),
}


def find_review_constructs(text: str) -> list[str]:
    """Return the names of constructs in ``text`` that need manual review."""
    return [name for name, pat in REVIEW_PATTERNS.items() if pat.search(text)]
