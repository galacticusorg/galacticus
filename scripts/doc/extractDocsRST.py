#!/usr/bin/env python3
"""Build-time extractor: embedded RST docstrings -> Sphinx pages.

Scans ``source/*.F90`` for ``!![ … !!]`` directive blocks that carry
``docformat="rst"`` (written by :mod:`convertDocstringsToRST`) and emits a
Sphinx documentation tree:

* ``<out>/physics/<family>.rst``  — one page per ``functionClass`` family, with
  a section per implementation class, its description, and the input parameters
  declared in the same source file.
* ``<out>/physics/index.rst``      — a toctree of all family pages.
* ``<out>/glossary.rst``           — the glossary, from ``doc/Glossary.tex``.
* ``<out>/references.rst``         — a single project-wide bibliography.

Descriptions are read **straight from the raw source** (not via the directive
parser, which left-strips each line and would destroy the indentation- and
blank-line-sensitive RST).  No Fortran compilation is required.

Usage::

    scripts/doc/extractDocsRST.py <sourceDir> <outDir>

Andrew Benson / Galacticus — RST documentation migration (2026).
"""
from __future__ import annotations

import html
import os
import re
import subprocess
import sys
import textwrap

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from latexToRST import (                                       # noqa: E402
    parse_glossary, glossary_display_map, latex_to_rst,
)

# ``!![ … !!]`` directive block.
_BLOCK_RE = re.compile(r'^[ \t]*!!\[(.*?)^[ \t]*!!\]', re.DOTALL | re.MULTILINE)
# Root (first) element opening tag of a block.
_ROOT_RE = re.compile(r'<([A-Za-z][\w]*)((?:\s+[\w:]+="[^"]*")*)\s*/?>')
# First <description>…</description> inside a block (the class/family/parameter
# description; nested method/entry descriptions are handled separately).
_DESC_RE = re.compile(r'<description>(.*?)</description>', re.DOTALL)
# A <method …>…</method> inside a functionClass block.
_METHOD_RE = re.compile(r'<method\b([^>]*?)>(.*?)</method>', re.DOTALL)
# An <entry label="…" [description="…"]/> inside an enumeration block.
_ENTRY_RE = re.compile(r'<entry\b([^>]*?)/?>')
_ARG_RE = re.compile(r'<argument>(.*?)</argument>', re.DOTALL)


def _attr(attrs: str, name: str) -> str | None:
    m = re.search(rf'\b{name}="([^"]*)"', attrs)
    return m.group(1) if m else None


def _child(block: str, tag: str) -> str | None:
    m = re.search(rf'<{tag}>(.*?)</{tag}>', block, re.DOTALL)
    return m.group(1).strip() if m else None


def _extract_methods(block: str) -> list[dict]:
    """Extract the ``<method>`` entries (name, description, type, arguments)
    that a ``functionClass`` directive declares for its class interface."""
    methods = []
    for m in _METHOD_RE.finditer(block):
        name = _attr(m.group(1), 'name')
        if not name:
            continue
        inner = m.group(2)
        dm = _DESC_RE.search(inner)
        methods.append({
            'name':        name,
            'description': dm.group(1) if dm else None,
            'type':        _child(inner, 'type'),
            'arguments':   [re.sub(r'\s+', ' ', a).strip()
                            for a in _ARG_RE.findall(inner)],
        })
    return methods


# A ``<methods>`` block of type-bound procedures (one per implementation), each
# ``<method method="…" description="…"/>``.  Unlike the functionClass
# ``<method name=…>`` interface, the description is a still-LaTeX attribute.
_TYPE_METHODS_RE = re.compile(r'<methods\b[^>]*>(.*?)</methods>', re.DOTALL)
_TYPE_METHOD_RE = re.compile(r'<method\b([^>]*?)/?>')


def _extract_type_methods(text: str) -> list[dict]:
    """Extract ``<method method="…" description="…"/>`` entries from every
    ``<methods>`` block in a file — the methods a concrete implementation's
    derived type declares, beyond the base-class interface."""
    out = []
    for mb in _TYPE_METHODS_RE.finditer(text):
        for m in _TYPE_METHOD_RE.finditer(mb.group(1)):
            name = _attr(m.group(1), 'method')
            if not name:
                continue
            out.append({'name': name,
                        'description': _attr(m.group(1), 'description')})
    return out


def _extract_entries(block: str) -> list[dict]:
    """Extract ``<entry label="…" description="…"/>`` items of an enumeration."""
    entries = []
    for m in _ENTRY_RE.finditer(block):
        label = _attr(m.group(1), 'label')
        if label is None:
            continue
        entries.append({'label': label,
                        'description': _attr(m.group(1), 'description')})
    return entries


def _desc_to_rst(raw: str | None, glsmap: dict) -> str:
    """Dedent a raw ``<description>`` body and pass it through the converter.

    The body is already RST (converted in the source).  It is XML-entity
    escaped there (so the build's directive parser accepts it), so we unescape
    it back to plain RST and dedent.  It is *not* re-converted from LaTeX — that
    already happened in the source.
    """
    if not raw:
        return ''
    return _process_figures(
        textwrap.dedent(html.unescape(raw).strip('\n')).strip('\n'))


# Build-time figure handling: each ``.. figure:: <path>.pdf`` (the figure's
# \includegraphics path, relative to doc/) is rasterised to a PNG under
# <out>/_figures and rewritten to a web path.  Set by main().
_FIGURES_DIR: str | None = None
_FIGURE_CACHE: dict[str, str | None] = {}


def _process_figures(text: str) -> str:
    if _FIGURES_DIR is None:
        return text

    def repl(m: re.Match) -> str:
        path = m.group(1).strip()
        if not path.lower().endswith('.pdf'):
            return m.group(0)
        pdf = os.path.normpath(os.path.join('doc', path))   # LaTeX runs in doc/
        if pdf not in _FIGURE_CACHE:
            _FIGURE_CACHE[pdf] = _rasterise(pdf)
        web = _FIGURE_CACHE[pdf]
        return f'.. figure:: {web}' if web else m.group(0)

    return re.sub(r'\.\. figure:: (\S+)', repl, text)


def _rasterise(pdf: str) -> str | None:
    """PDF -> PNG (150 dpi) under _FIGURES_DIR; returns the web path or None."""
    if not os.path.isfile(pdf):
        print(f'warning: figure not found: {pdf}', file=sys.stderr)
        return None
    name = re.sub(r'[^A-Za-z0-9_.-]', '_',
                  os.path.splitext(os.path.basename(pdf))[0])
    out = os.path.join(_FIGURES_DIR, name)                   # pdftocairo adds .png
    try:
        subprocess.run(['pdftocairo', '-png', '-r', '150', '-singlefile',
                        pdf, out], check=True, capture_output=True)
    except (OSError, subprocess.CalledProcessError) as exc:
        print(f'warning: pdftocairo failed for {pdf}: {exc}', file=sys.stderr)
        return None
    return f'/_figures/{name}.png'


def _heading(text: str, char: str) -> str:
    return f'{text}\n{char * max(3, len(text))}\n'


def _lcfirst(s: str) -> str:
    return s[:1].lower() + s[1:] if s else s


# ---------------------------------------------------------------------------
# Raw source scan
# ---------------------------------------------------------------------------

def scan_source(source_dir: str):
    """Return ``(families, implementations, params_by_file)`` from raw source.

    * ``families``: ``{name: {descriptiveName, description, default}}`` for every
      ``functionClass``.
    * ``implementations``: ``{family: [{name, description, file}, …]}`` keyed by
      the root element name.
    * ``params_by_file``: ``{file_base: [{name, default, description}, …]}``.
    * ``methods_by_file``: ``{file_base: [{name, description}, …]}`` — the
      type-bound methods declared in that file's ``<methods>`` block(s).
    * ``modules``: ``[{name, file, description, classRef}, …]`` — one per
      documented module (``classRef`` links class modules to their page).
    * ``workarounds``: ``[{type, pr, url, description, seeAlso, file}, …]``.
    """
    families: dict[str, dict] = {}
    by_root: dict[str, list[dict]] = {}
    params_by_file: dict[str, list[dict]] = {}
    methods_by_file: dict[str, list[dict]] = {}
    enumerations: list[dict] = []
    modules: list[dict] = []
    workarounds: list[dict] = []
    class_by_file: dict[str, str] = {}

    for root, _dirs, files in os.walk(source_dir):
        for fn in sorted(files):
            if not fn.endswith('.F90') or fn.startswith('.#'):
                continue
            base = os.path.splitext(fn)[0]
            with open(os.path.join(root, fn), encoding='utf-8',
                      errors='replace') as fh:
                text = fh.read()
            mmod = re.search(r'(?im)^[ \t]*module[ \t]+([a-z0-9_]+)[ \t]*$', text)
            module_name = mmod.group(1) if mmod else base
            # Module-level description: the !!{RST …!!} block right after the
            # ``module`` line (already RST; not XML-escaped like <description>).
            if mmod:
                md = re.match(r'\s*!!\{RST\b(.*?)!!\}',
                              text[mmod.end():mmod.end() + 800], re.DOTALL)
                if md and md.group(1).strip():
                    modules.append({'name': module_name, 'file': base,
                                    'description': md.group(1).strip()})
            tmeths = _extract_type_methods(text)
            if tmeths:
                methods_by_file[base] = tmeths
            for bm in _BLOCK_RE.finditer(text):
                block = bm.group(1)
                if 'docformat="rst"' not in block:
                    continue
                rm = _ROOT_RE.search(block)
                if not rm:
                    continue
                rtype, attrs = rm.group(1), rm.group(2)
                desc_m = _DESC_RE.search(block)
                desc_raw = desc_m.group(1) if desc_m else None

                if rtype == 'functionClass':
                    name = _child(block, 'name')
                    if not name:
                        continue
                    families[name] = {
                        'descriptiveName': _child(block, 'descriptiveName') or name,
                        'description':     desc_raw,
                        'default':         _child(block, 'default'),
                        'methods':         _extract_methods(block),
                    }
                    class_by_file[base] = name
                elif rtype == 'workaround':
                    if desc_raw is None:
                        continue
                    workarounds.append({
                        'type':        _attr(attrs, 'type'),
                        'pr':          _attr(attrs, 'PR'),
                        'url':         _attr(attrs, 'url'),
                        'description': desc_raw,
                        'seeAlso':     [{'type': _attr(s, 'type'),
                                         'pr':   _attr(s, 'PR'),
                                         'url':  _attr(s, 'url')}
                                        for s in re.findall(
                                            r'<seeAlso\b([^>]*?)/?>', block)],
                        'file':        base,
                    })
                elif rtype == 'enumeration':
                    name = _child(block, 'name')
                    if not name:
                        continue
                    enumerations.append({
                        'name':        name,
                        'description': desc_raw,
                        'entries':     _extract_entries(block),
                        'module':      module_name,
                    })
                elif rtype == 'inputParameter':
                    name = _child(block, 'name')
                    if not name:
                        continue
                    params_by_file.setdefault(base, []).append({
                        'name':        name,
                        'default':     _child(block, 'defaultValue'),
                        'description': desc_raw,
                    })
                else:
                    name = _attr(attrs, 'name')
                    if not name or desc_raw is None:
                        continue
                    by_root.setdefault(rtype, []).append({
                        'name':        name,
                        'description': desc_raw,
                        'file':        base,
                    })

    implementations = {fam: by_root[fam] for fam in by_root if fam in families}
    enumerations.sort(key=lambda e: e['name'].lower())
    # Tag each module that is a (page-bearing) functionClass file so the modules
    # index can link to its physics page.
    paged = set(implementations)
    for m in modules:
        cls = class_by_file.get(m['file'])
        m['classRef'] = cls if cls in paged else None
    modules.sort(key=lambda m: m['name'].lower())
    # The same compiler workaround is applied at many call sites with identical
    # text; collapse to one entry per (type, PR, description), unioning seeAlso.
    deduped: dict[tuple, dict] = {}
    for w in workarounds:
        key = (w.get('type'), w.get('pr'),
               ' '.join((w.get('description') or '').split()))
        entry = deduped.setdefault(key, {**w, 'seeAlso': []})
        have = {(s.get('type'), s.get('pr'), s.get('url'))
                for s in entry['seeAlso']}
        for s in w.get('seeAlso') or []:
            if (s.get('type'), s.get('pr'), s.get('url')) not in have:
                entry['seeAlso'].append(s)
    workarounds = sorted(deduped.values(),
                         key=lambda w: (str(w.get('type')), str(w.get('pr'))))
    return (families, implementations, params_by_file, methods_by_file,
            enumerations, modules, workarounds)


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------

def render_parameter(p: dict, glsmap: dict) -> str:
    name = p.get('name', '')
    desc = _desc_to_rst(p.get('description'), glsmap)
    # Parameter descriptions are short; collapse to a single logical line.
    desc = re.sub(r'\s*\n\s*', ' ', desc).strip()
    default = p.get('default')
    head = f'``[{name}]``'
    if default not in (None, ''):
        head += f' (default ``{str(default).strip()}``)'
    return f'* {head} — {desc}' if desc else f'* {head}'


def render_family(fam: str, families: dict, implementations: dict,
                  params_by_file: dict, methods_by_file: dict,
                  glsmap: dict) -> str:
    directive = families[fam]
    descriptive = directive.get('descriptiveName', fam)
    out = [f'.. _physics-{fam}:\n', _heading(descriptive, '=')]

    desc = _desc_to_rst(directive.get('description'), glsmap)
    if desc:
        out.append(desc + '\n')

    default = directive.get('default')
    if default:
        target = f'{fam}{default[:1].upper()}{default[1:]}'
        out.append(f'**Default implementation:** ``{target}``\n')

    methods = directive.get('methods') or []
    if methods:
        out.append(_heading('Methods', '-'))
        for meth in methods:
            sig = f'``{meth["name"]}``'
            rtype = (meth.get('type') or '').strip()
            if rtype:
                sig += f' → ``{rtype}``'
            mdesc = re.sub(r'\s*\n\s*', ' ',
                           _desc_to_rst(meth.get('description'), glsmap)).strip()
            body = mdesc or '—'
            for arg in meth.get('arguments') or []:
                body += f'\n\n* ``{arg}``'
            out.append(sig)
            out.append(textwrap.indent(body, '   '))
            out.append('')

    impls = sorted(implementations.get(fam, []),
                   key=lambda d: str(d.get('name', '')).lower())
    # (No ``.. contents::`` — the Furo theme provides a right-hand on-page nav
    # and errors on an explicit local table of contents.)

    for impl in impls:
        name = impl['name']
        out.append(f'.. _physics-{name}:\n')
        out.append(_heading(f'``{name}``', '-'))
        idesc = _desc_to_rst(impl.get('description'), glsmap)
        if idesc:
            out.append(idesc + '\n')
        if default and f'{fam}{default[:1].upper()}{default[1:]}' == name:
            out.append('**(Default implementation)**\n')
        imethods = methods_by_file.get(impl.get('file', ''), [])
        if imethods:
            out.append('**Methods**\n')
            for meth in imethods:
                # The description attribute is already RST (converted in source);
                # just unescape the XML entities, as for <description> elements.
                mdesc = re.sub(
                    r'\s*\n\s*', ' ',
                    html.unescape(meth.get('description') or '')).strip()
                out.append(f'* ``{meth["name"]}`` — {mdesc}' if mdesc
                           else f'* ``{meth["name"]}``')
            out.append('')
        params = params_by_file.get(impl.get('file', ''), [])
        if params:
            out.append('**Parameters**\n')
            for p in params:
                out.append(render_parameter(p, glsmap))
            out.append('')
    return '\n'.join(out) + '\n'


def render_glossary(glossary: dict, glsmap: dict) -> str:
    # Deduplicate by display name; keep the longest description on collision.
    by_name: dict[str, str] = {}
    for entry in glossary.values():
        name = entry['name']
        desc = entry.get('description', '')
        if name not in by_name or len(desc) > len(by_name[name]):
            by_name[name] = desc

    out = [_heading('Glossary', '='), '.. glossary::\n']
    for name in sorted(by_name, key=str.lower):
        desc_rst = latex_to_rst(by_name[name], glsmap).strip()
        # Glossary definitions are single-paragraph prose; keep them tidy.
        desc_rst = re.sub(r'\s*\n\s*', ' ', desc_rst) or '—'
        out.append(textwrap.indent(name, '   '))
        out.append(textwrap.indent(desc_rst, '      '))
        out.append('')
    return '\n'.join(out) + '\n'


def render_enumerations(enumerations: list[dict], glsmap: dict) -> str:
    out = [_heading('Enumerations', '='),
           'Enumerated types defined in the Galacticus source code.\n']
    # Enumeration names are not globally unique (several modules define their
    # own e.g. ``propertyType``); disambiguate any repeated name by appending
    # its defining module so the headings (and Sphinx anchors) stay distinct.
    name_counts: dict[str, int] = {}
    for en in enumerations:
        name_counts[en['name']] = name_counts.get(en['name'], 0) + 1
    for en in enumerations:
        title = f'``{en["name"]}``'
        if name_counts[en['name']] > 1:
            title += f' ({en.get("module") or "?"})'
        out.append(_heading(title, '-'))
        desc = re.sub(r'\s*\n\s*', ' ',
                      _desc_to_rst(en.get('description'), glsmap)).strip()
        if desc:
            out.append(desc + '\n')
        entries = en.get('entries') or []
        if entries:
            out.append('**Values**\n')
            for e in entries:
                ed = e.get('description')
                if ed:
                    ed = re.sub(r'\s*\n\s*', ' ', latex_to_rst(ed, glsmap)).strip()
                    out.append(f'* ``{e["label"]}`` — {ed}')
                else:
                    out.append(f'* ``{e["label"]}``')
            out.append('')
    return '\n'.join(out) + '\n'


def render_modules(modules: list[dict]) -> str:
    """A reference index of every documented source module + its summary."""
    out = [_heading('Modules', '='),
           'Every documented Fortran module in the Galacticus source, with the '
           'summary from its embedded documentation.  Modules implementing a '
           'pluggable physics class link to that class.\n']
    for m in modules:
        desc = re.sub(r'\s*\n\s*', ' ',
                      textwrap.dedent(m['description']).strip()).strip()
        body = desc or '—'
        if m.get('classRef'):
            body += f'\n\nSee :ref:`physics-{m["classRef"]}`.'
        out.append(f'``{m["name"]}``')
        out.append(textwrap.indent(body, '   '))
        out.append('')
    return '\n'.join(out) + '\n'


def render_workarounds(workarounds: list[dict], glsmap: dict) -> str:
    """A developer reference of the compiler workarounds documented in source."""
    def link(typ, pr, url):
        label = ' '.join(x for x in (typ, f'PR {pr}' if pr else None) if x)
        label = label or 'reference'
        url = html.unescape(url) if url else None
        return f'`{label} <{url}>`_' if url else label

    out = [_heading('Source-Code Workarounds', '='),
           'Workarounds for compiler bugs and limitations, documented inline in '
           'the Galacticus source.\n']
    for w in workarounds:
        body = textwrap.dedent(_desc_to_rst(w.get('description'), glsmap)).strip()
        body = body or '—'
        for s in w.get('seeAlso') or []:
            body += ('\n\nSee also: '
                     + link(s.get('type'), s.get('pr'), s.get('url')) + '.')
        out.append(link(w.get('type'), w.get('pr'), w.get('url')))
        out.append(textwrap.indent(body, '   '))
        out.append('')
    return '\n'.join(out) + '\n'


# ``<constant …/>`` directives in the source define the physical/mathematical/…
# constants.  (constants.py reads the build's *.constants.xml; the same data is
# in the source directives, so we read those directly — no compiled build.)
_CONSTANT_RE = re.compile(r'<constant\b([^>]*?)/?>')
_CONSTANT_GROUPS = {
    'astrophysical': 'Astrophysical constants',
    'atomic':        'Atomic physics constants',
    'physical':      'Physical constants',
    'math':          'Mathematical constants',
    'units':         'Units',
    'prefixes':      'SI prefixes',
    'GSL':           'GNU Scientific Library constants',
    'Kernel':        'Kernel constants',
    'misc':          'Miscellaneous constants',
}


def scan_constants(source_dir: str) -> list[dict]:
    constants = []
    for root, _dirs, files in os.walk(source_dir):
        for fn in sorted(files):
            if not fn.endswith('.F90') or fn.startswith('.#'):
                continue
            with open(os.path.join(root, fn), encoding='utf-8',
                      errors='replace') as fh:
                text = fh.read()
            for m in _CONSTANT_RE.finditer(text):
                attrs = {k: html.unescape(v)
                         for k, v in re.findall(r'(\w+)="([^"]*)"', m.group(1))}
                if 'variable' in attrs:
                    constants.append(attrs)
    return constants


def _format_constant_value(value: str) -> str:
    value = re.sub(r'd([+\-0-9]+)$', r'e\1', value)   # Fortran d-exponent -> e
    value = re.sub(r'_[_a-zA-Z0-9]+$', '', value)     # strip Fortran kind suffix
    return value


def render_constants(constants: list[dict], glsmap: dict) -> str:
    """A reference of the constants defined via ``<constant>`` source directives."""
    groups: dict[str, list[dict]] = {}
    for c in constants:
        for g in (c.get('group') or 'misc').split(':'):
            groups.setdefault(g if g in _CONSTANT_GROUPS else 'misc', []).append(c)

    out = [_heading('Constants', '='),
           'Physical, mathematical and astrophysical constants available in the '
           'Galacticus source (defined via ``<constant>`` directives). Constants '
           'whose value is provided by the GNU Scientific Library note the GSL '
           'symbol rather than a literal value.\n']
    for g in sorted(groups, key=lambda k: _CONSTANT_GROUPS[k].lower()):
        out.append(_heading(_CONSTANT_GROUPS[g], '-'))
        for c in sorted(groups[g], key=lambda c: c.get('variable', '').lower()):
            term = f'``{c["variable"]}``'
            if c.get('value'):
                term += f' = ``{_format_constant_value(c["value"])}``'
            if c.get('symbol'):
                term += f' (:math:`{c["symbol"]}`)'
            units = c.get('units', '')
            if units and units != 'dimensionless':
                term += f' [{units}]'
            desc = re.sub(r'\s*\n\s*', ' ',
                          latex_to_rst(c.get('description', ''), glsmap)).strip()
            extras = []
            if c.get('gslSymbol'):
                extras.append(f'value from GSL ``{c["gslSymbol"]}``')
            if c.get('reference'):
                extras.append(f'`{c["reference"]} <{c["referenceURL"]}>`_'
                              if c.get('referenceURL') else c['reference'])
            if c.get('externalDescription'):
                extras.append(f'`more <{c["externalDescription"]}>`_')
            body = (desc or '—') + ('\n\n' + '; '.join(extras) if extras else '')
            out.append(term)
            out.append(textwrap.indent(body, '   '))
            out.append('')
    return '\n'.join(out) + '\n'


def render_physics_index(families: dict, implementations: dict) -> str:
    fams = sorted((f for f in families if f in implementations),
                  key=lambda f: families[f].get('descriptiveName', f).lower())
    out = [_heading('Physics Classes', '='),
           'Documentation for each pluggable physics class, generated from the '
           'source code.\n',
           '.. toctree::\n   :maxdepth: 1\n']
    for fam in fams:
        out.append(f'   {fam}')
    return '\n'.join(out) + '\n'


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> int:
    if len(sys.argv) != 3:
        print('Usage: extractDocsRST.py <sourceDir> <outDir>', file=sys.stderr)
        return 1
    source_dir, out_dir = sys.argv[1], sys.argv[2]

    global _FIGURES_DIR
    _FIGURES_DIR = os.path.join(out_dir, '_figures')
    os.makedirs(_FIGURES_DIR, exist_ok=True)

    glossary = parse_glossary(os.path.join('doc', 'Glossary.tex'))
    glsmap = glossary_display_map(glossary)

    (families, implementations, params_by_file, methods_by_file, enumerations,
     modules, workarounds) = scan_source(source_dir)

    physics_dir = os.path.join(out_dir, 'physics')
    os.makedirs(physics_dir, exist_ok=True)

    for fam in implementations:
        page = render_family(fam, families, implementations,
                             params_by_file, methods_by_file, glsmap)
        with open(os.path.join(physics_dir, f'{fam}.rst'), 'w',
                  encoding='utf-8') as fh:
            fh.write(page)

    with open(os.path.join(physics_dir, 'index.rst'), 'w',
              encoding='utf-8') as fh:
        fh.write(render_physics_index(families, implementations))

    with open(os.path.join(out_dir, 'glossary.rst'), 'w',
              encoding='utf-8') as fh:
        fh.write(render_glossary(glossary, glsmap))

    with open(os.path.join(out_dir, 'enumerations.rst'), 'w',
              encoding='utf-8') as fh:
        fh.write(render_enumerations(enumerations, glsmap))

    with open(os.path.join(out_dir, 'references.rst'), 'w',
              encoding='utf-8') as fh:
        # Only the works actually cited in the documentation — not the whole
        # (400+ entry) Galacticus.bib database.
        fh.write(_heading('References', '=') +
                 '\n.. bibliography::\n   :filter: cited\n')

    with open(os.path.join(out_dir, 'modules.rst'), 'w',
              encoding='utf-8') as fh:
        fh.write(render_modules(modules))

    with open(os.path.join(out_dir, 'workarounds.rst'), 'w',
              encoding='utf-8') as fh:
        fh.write(render_workarounds(workarounds, glsmap))

    constants = scan_constants(source_dir)
    with open(os.path.join(out_dir, 'constants.rst'), 'w',
              encoding='utf-8') as fh:
        fh.write(render_constants(constants, glsmap))

    n_impl = sum(len(v) for v in implementations.values())
    n_meth = sum(len(f.get('methods') or []) for f in families.values())
    impl_files = {i.get('file') for v in implementations.values() for i in v}
    n_tmeth = sum(len(m) for f, m in methods_by_file.items() if f in impl_files)
    print(f'Wrote {len(implementations)} physics pages ({n_impl} '
          f'implementations, {n_meth} interface + {n_tmeth} type methods), '
          f'{len(enumerations)} enumerations, {len(modules)} modules, '
          f'{len(workarounds)} workarounds, {len(constants)} constants, glossary '
          f'({len(glossary)} entries), references.', file=sys.stderr)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
