# Provides XML element-to-dict conversion utilities.
# Mirrors Perl XML::Simple default behaviour, with optional KeyAttr and ForceArray support.
# Andrew Benson (2026)

import xml.etree.ElementTree as ET


def xml_to_dict(element, keyed_tags=None, force_array=None):
    """Recursively convert an XML element to nested Python dicts/lists/strings.

    Mirrors Perl XML::Simple default behaviour:
    - Element attributes are merged into the result dict.
    - Multiple same-tag children -> list; single same-tag child -> dict.
    - force_array: set of tag names always returned as a list (even if only one element).
    - keyed_tags:  set of tag names whose children are keyed by their 'name' attribute
                   (mirrors XML::Simple KeyAttr behaviour).
    - Text-only elements (no attributes, no children) -> plain string.
    - Elements with text plus attributes/children -> text stored under 'content' key.

    Works with both lxml.etree and stdlib xml.etree.ElementTree elements.
    """
    result = dict(element.attrib)
    children_by_tag = {}
    for child in element:
        children_by_tag.setdefault(child.tag, []).append(child)
    for tag, children in children_by_tag.items():
        if keyed_tags and tag in keyed_tags:
            result[tag] = {
                child.get('name'): xml_to_dict(child, keyed_tags, force_array)
                for child in children
            }
        else:
            converted = [xml_to_dict(child, keyed_tags, force_array) for child in children]
            is_forced = force_array is not None and tag in force_array
            result[tag] = converted if (is_forced or len(converted) > 1) else converted[0]
    text = (element.text or '').strip()
    if text:
        if result:
            result['content'] = text
        else:
            return text
    return result


def dict_to_xml_string(root_name, data):
    """Serialise `data` to a pretty-printed XML string, XML::Simple-style.

    Mirrors Perl `XML::Simple::XMLout($data, RootName => root_name, NoAttr => 1)`:
    every dict key becomes a child element (no attributes); list values are
    rendered as repeated sibling elements; scalars become the element text.
    Dict keys are emitted in sorted order for deterministic output.
    """
    root = ET.Element(root_name)
    _fill_element(root, data)
    _indent(root)
    return ET.tostring(root, encoding='unicode') + "\n"


def _fill_element(elem, data):
    if isinstance(data, dict):
        for key in sorted(data.keys()):
            value = data[key]
            if isinstance(value, list):
                for item in value:
                    child = ET.SubElement(elem, key)
                    _fill_element(child, item)
            elif isinstance(value, dict):
                child = ET.SubElement(elem, key)
                _fill_element(child, value)
            else:
                child = ET.SubElement(elem, key)
                if value is not None:
                    child.text = str(value)
    elif isinstance(data, list):
        for item in data:
            child = ET.SubElement(elem, 'item')
            _fill_element(child, item)
    else:
        if data is not None:
            elem.text = str(data)


def _indent(elem, level=0):
    """In-place two-space pretty-printer."""
    pad = "\n" + "  " * level
    if len(elem):
        if not (elem.text and elem.text.strip()):
            elem.text = pad + "  "
        for i, child in enumerate(elem):
            _indent(child, level + 1)
            child.tail = (pad + "  ") if (i + 1 < len(elem)) else pad
    elif level and (not elem.tail or not elem.tail.strip()):
        elem.tail = pad
