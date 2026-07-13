"""Scan a Fortran source file for `!![…!!]` XML directive blocks on disk.

Andrew Benson (ported to Python 2026)

Mirrors perl/Galacticus/Build/Directives.pm — specifically Extract_Directive
(incremental single-directive iterator) and Extract_Directives (return-all
wrapper).  This is distinct from SourceTree.Parse.Directives.parse_directives,
which walks an already-parsed AST; the functions here open a file on disk and
extract its directives without building a SourceTree.
"""
from __future__ import annotations

import os
import re
import xml.etree.ElementTree as ET

from typing import Any

from XML.Utils import xml_to_dict

__all__ = ['extract_directives', 'extract_directive']


_END_MARKER_RE    = re.compile(r'^\s*!!\]')
_START_MARKER_RE  = re.compile(r'^\s*!!\[')
_INSTRUMENT_RE    = re.compile(r'^!-->')
_STRIP_BANG_LT_RE = re.compile(r'^(!<)?\s*')
_OPEN_TAG_RE      = re.compile(r'<([a-zA-Z0-9]+)[^/>]*>')
_CLOSE_TAG_RE     = re.compile(r'</([a-zA-Z0-9]+)>')


def _matches_conditions(directive: dict, conditions: dict | None) -> bool:
    """Return True when every key/value in `conditions` is present and equal in
    `directive`.  Mirrors the inner filter loop at Directives.pm:55-65.
    """
    if not conditions:
        return True
    for k, v in conditions.items():
        if directive.get(k) != v:
            return False
    return True


def extract_directives(file_name: str, directive_name,
                       conditions: dict | None = None,
                       set_root_element_type: bool = False) -> list[dict]:
    """Return every directive in `file_name` whose root element matches
    `directive_name` — a single name, a set/tuple/list of names, or `'*'`
    to match any root element.

    Passing several names in one call reads and parses the file once
    (each call performs a full line-by-line read, so N single-name calls
    cost N reads); use `set_root_element_type=True` to recover which name
    each returned directive matched.

    Mirrors Perl Extract_Directives (Directives.pm:80-98) plus the underlying
    state machine from Extract_Directive (Directives.pm:11-78):

    - Lines matching `^\\s*!![ ` start an XML section.
    - Lines matching `^\\s*!!]` end the current XML section.
    - Lines starting with `!-->` (SourceIntrospection instrumentation) are skipped.
    - Inside an XML section, `!<` prefixes and trailing whitespace are stripped
      and `&nbsp;` entities are converted to spaces.
    - XML parsing kicks in when tag-depth returns to zero and we have accumulated
      non-empty text — handling both inline `<foo .../>` forms and multi-line
      `<foo> ... </foo>` blocks.

    Parameters
    ----------
    file_name : str
        Path to a Fortran source file.  If it does not exist, an empty list
        is returned (matches Perl's `return unless -e $fileName`).
    directive_name : str
        Either a specific root element name, or `'*'` to match any.
    conditions : dict, optional
        Attribute/value pairs a directive must match to be included.
    set_root_element_type : bool, optional
        When True, each returned dict gains a `rootElementType` key carrying
        the actual root element name (useful when `directive_name == '*'`).

    Returns
    -------
    list[dict]
        Each directive's XML, converted via `xml_to_dict`, without the root
        element wrapper (matching Perl's `$directive->{$rootName}` unwrap).
    """
    if not os.path.exists(file_name):
        return []

    results   = []
    in_xml    = False
    depth     = 0
    xml_text  = ''

    with open(file_name, 'r', errors='replace') as fh:
        for line in fh:
            if _END_MARKER_RE.match(line):
                in_xml = False
                # Fall through to accumulate tail of line if currently parsing.
            elif _INSTRUMENT_RE.match(line):
                continue

            if in_xml or depth > 0:
                body = _STRIP_BANG_LT_RE.sub('', line) if in_xml else line
                body = body.replace('&nbsp;', ' ')
                xml_text += body
                depth += len(_OPEN_TAG_RE.findall(line))
                depth -= len(_CLOSE_TAG_RE.findall(line))
                if xml_text.strip() and depth == 0:
                    directive = _parse_xml_block(
                        xml_text, file_name, directive_name,
                        conditions, set_root_element_type)
                    if directive is not None:
                        results.append(directive)
                    xml_text = ''

            if _START_MARKER_RE.match(line):
                in_xml = True

    return results


def extract_directive(file_name: str, directive_name: str, **kwargs: Any) -> dict | None:
    """Return the first directive matching `directive_name` (or None).

    Thin convenience wrapper around `extract_directives`.
    """
    directives = extract_directives(file_name, directive_name, **kwargs)
    return directives[0] if directives else None


def _parse_xml_block(xml_text: str, file_name: str, directive_name,
                     conditions: dict | None,
                     set_root_element_type: bool) -> dict | None:
    """Parse one accumulated XML block and return a dict if the root matches
    `directive_name`, else None.

    Matches the Perl error-handling idiom: any XML parse failure is fatal
    (raised as RuntimeError) since the caller has no way to recover; Perl
    `die`s in the same spot.
    """
    try:
        elem = ET.fromstring(xml_text)
    except ET.ParseError as exc:
        raise RuntimeError(
            f"extract_directives: failed parsing XML while extracting directive "
            f"'{directive_name}' from '{file_name}': {exc}\nXML content was:\n{xml_text}"
        )

    if isinstance(directive_name, str):
        if directive_name != '*' and elem.tag != directive_name:
            return None
    elif elem.tag not in directive_name:
        return None

    directive = xml_to_dict(elem)
    if not isinstance(directive, dict):
        # A directive whose root is text-only (`<foo>hello</foo>`) — wrap it
        # so callers can still attach rootElementType / conditions uniformly.
        directive = {'content': directive}

    if set_root_element_type:
        directive['rootElementType'] = elem.tag

    if not _matches_conditions(directive, conditions):
        return None

    return directive
