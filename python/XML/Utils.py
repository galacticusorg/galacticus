# Provides XML element-to-dict conversion utilities.
# Mirrors Perl XML::Simple default behaviour, with optional KeyAttr and ForceArray support.
# Andrew Benson (2026)


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
