# Processes XML directives embedded in Fortran/C source files.
# Python port of perl/Galacticus/Build/Directives.pm.
# Andrew Benson (ported to Python 2026)

import os
import re
import xml.etree.ElementTree as ET


_XML_START     = re.compile(r'^\s*!!\[')
_XML_END       = re.compile(r'^\s*!!\]')
_INSTRUMENT    = re.compile(r'^!-->')
_IN_XML_PREFIX = re.compile(r'^(!<)?\s*')
_OPEN_TAG      = re.compile(r'<([a-zA-Z0-9]+)[^/>]*>')
_CLOSE_TAG     = re.compile(r'</([a-zA-Z0-9]+)>')

# Attribute names used by Perl's XML::Simple to auto-fold a group of
# same-tagged siblings into a single dict keyed by that attribute.
_KEY_ATTRS = ('name', 'key', 'id')


def xml_to_dict_simple(element):
    """Convert an XML element to a nested dict, mirroring Perl XML::Simple's
    default behaviour (KeepRoot not applied; the caller handles that if it
    wants the root element as a dict key).

    Rules:
    - Element attributes are merged into the result dict.
    - Multiple same-tag children are folded into a dict keyed by 'name', 'key'
      or 'id' (in that order) if every sibling in the group carries the same
      attribute; otherwise they become a list.
    - A single same-tag child becomes the child's value directly.
    - Text-only elements (no attributes, no children) return a plain string.
    - For elements with child elements and additional text ("mixed content"),
      every text fragment (element.text and each child's .tail) that is not
      whitespace-only contributes verbatim to a 'content' field.
    - Self-closing elements with no attributes return an empty dict.
    """
    result = dict(element.attrib)
    children_by_tag = {}
    for child in element:
        children_by_tag.setdefault(child.tag, []).append(child)
    for tag, children in children_by_tag.items():
        converted = [xml_to_dict_simple(child) for child in children]
        if len(children) > 1:
            key_attr = None
            if all(isinstance(c, dict) for c in converted):
                for candidate in _KEY_ATTRS:
                    if all(candidate in c for c in converted):
                        key_attr = candidate
                        break
            if key_attr is not None:
                folded = {}
                for value in converted:
                    key_val = value[key_attr]
                    folded[key_val] = {k: v for k, v in value.items() if k != key_attr}
                result[tag] = folded
            else:
                result[tag] = converted
        else:
            result[tag] = converted[0]
    has_children = len(element) > 0
    if has_children:
        pieces = []
        if element.text and element.text.strip():
            pieces.append(element.text)
        for child in element:
            if child.tail and child.tail.strip():
                pieces.append(child.tail)
        if pieces:
            result['content'] = ''.join(pieces)
        return result
    if result:
        if element.text and element.text.strip():
            result['content'] = element.text
        return result
    # No attributes and no children: text-only or empty element.
    if element.text is None or element.text == '':
        return {}
    return element.text


def extract_directive(line_iter, directive_name, state,
                      set_root_element_type=False, conditions=None,
                      file_name=None):
    """Extract a single named directive from an iterator of source-file lines.

    Mirrors Perl Galacticus::Build::Directives::Extract_Directive().

    line_iter is any iterator yielding lines (typically an open file handle).
    state is a dict used to carry the in-XML flag between calls (so repeated
    invocations against the same iterator can parse consecutive directives).
    directive_name is the root-element name to match, or '*' to match any.
    If set_root_element_type is True, the returned directive dict will contain
    a 'rootElementType' field holding the matched root element name.
    conditions, if given, is a dict of key/value pairs that must all be
    present and equal in the directive for it to be returned.
    file_name is used in error messages if XML parsing fails.

    Returns the parsed directive as a dict, or None if the iterator is
    exhausted without finding a matching directive.
    """
    state.setdefault('inXML', False)
    state.setdefault('lineNumber', 0)
    xml_text = ''
    depth = 0
    for line in line_iter:
        state['lineNumber'] += 1
        if _XML_END.match(line):
            state['inXML'] = False
        if _INSTRUMENT.match(line):
            continue
        if state['inXML'] or depth > 0:
            processed = line
            if state['inXML']:
                processed = _IN_XML_PREFIX.sub('', processed, count=1)
            processed = processed.replace('&nbsp;', ' ')
            xml_text += processed
            depth += len(_OPEN_TAG .findall(line))
            depth -= len(_CLOSE_TAG.findall(line))
            if xml_text.strip() and depth == 0:
                try:
                    root = ET.fromstring(xml_text)
                except ET.ParseError as err:
                    where = file_name if file_name is not None else '<unknown>'
                    raise RuntimeError(
                        "extract_directive: while extracting directive '"
                        + directive_name + "' from line "
                        + str(state['lineNumber']) + " of file '" + where
                        + "' failed parsing with message:\n" + str(err)
                        + "\n XML content was:\n" + xml_text
                    )
                matched = root.tag
                if directive_name == '*' or matched == directive_name:
                    directive = xml_to_dict_simple(root)
                    if not isinstance(directive, dict):
                        directive = {'content': directive}
                    if set_root_element_type:
                        directive['rootElementType'] = matched
                    if conditions:
                        mismatch = False
                        for key, value in conditions.items():
                            if directive.get(key) != value:
                                mismatch = True
                                break
                        if mismatch:
                            xml_text = ''
                            continue
                    return directive
                xml_text = ''
        if _XML_START.match(line):
            state['inXML'] = True
    return None


def extract_directives(file_name, directive_name,
                       set_root_element_type=False, conditions=None):
    """Extract all matching directives from the given file.

    Mirrors Perl Galacticus::Build::Directives::Extract_Directives().

    Returns a list of directive dicts (empty list if the file does not exist
    or contains no matching directives).
    """
    if not os.path.exists(file_name):
        return []
    directives = []
    state = {}
    with open(file_name, 'r', encoding='utf-8', errors='replace') as handle:
        while True:
            directive = extract_directive(
                handle, directive_name, state,
                set_root_element_type=set_root_element_type,
                conditions=conditions,
                file_name=file_name,
            )
            if directive is None:
                break
            directives.append(directive)
    return directives
