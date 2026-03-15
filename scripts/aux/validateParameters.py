#!/usr/bin/env python3
import sys
import os
import re
import argparse
import xml.etree.ElementTree as ET
from collections import defaultdict

# Validate a Galacticus XML parameter file.
# Andrew Benson (1-December-2013)

parser = argparse.ArgumentParser(prog='validateParameters.py', description='Validate a Galacticus XML parameter file.')
parser.add_argument('file', help='Parameter file to validate')
args = parser.parse_args()

# Parse the file.
try:
    tree = ET.parse(args.file)
except ET.ParseError as e:
    print(f"Failed to parse XML: {e}", file=sys.stderr)
    sys.exit(1)

root = tree.getroot()

# Determine the format.
format_version = 2  # Best guess if no other information.
format_elem = root.find('formatVersion')
if format_elem is not None:
    format_text = format_elem.text
    if format_text is None or not format_text.strip():
        print("Invalid or missing 'formatVersion' value in XML: expected a non-empty integer.", file=sys.stderr)
        sys.exit(1)
    try:
        format_version = int(format_text.strip())
    except ValueError:
        print(f"Invalid 'formatVersion' value '{format_text}': expected an integer.", file=sys.stderr)
        sys.exit(1)
elif root.find('parameter') is not None:
    format_version = 1

valid = 0  # 0 = valid, 1 = invalid (matches Perl convention).

if format_version == 1:
    # Handle format version 1.
    # Validate using XML schema if lxml is available.
    try:
        from lxml import etree as lxml_etree
        exec_path = os.environ.get('GALACTICUS_EXEC_PATH', '')
        schema_file = os.path.join(exec_path, 'schema', 'parameters.xsd')
        with open(schema_file, 'rb') as f:
            schema_doc = lxml_etree.parse(f)
        schema = lxml_etree.XMLSchema(schema_doc)
        doc = lxml_etree.parse(args.file)
        if not schema.validate(doc):
            print(f"Parameter file fails XML schema validation\n{schema.error_log}", file=sys.stderr)
            sys.exit(1)
    except ImportError:
        # lxml not available; skip schema validation.
        pass
    except (OSError, lxml_etree.XMLSchemaParseError, lxml_etree.XMLSyntaxError) as e:
        # Fail with a clear, actionable error if the schema file cannot be opened
        # or if the schema/document cannot be parsed.
        print(
            f"Failed to validate parameters using schema '{schema_file}' "
            f"for parameter file '{args.file}': {e}",
            file=sys.stderr,
        )
        sys.exit(1)

    # Check for duplicated entries.
    names = defaultdict(int)
    for param in root.findall('parameter'):
        name_elem = param.find('name')
        if name_elem is not None and name_elem.text:
            names[name_elem.text.strip()] += 1
    for name, count in names.items():
        if name == 'xi:include':
            continue
        if count > 1:
            valid = 1
            print(f"Parameter '{name}' appears {count} times - should appear only once")

elif format_version == 2:
    # Handle format version 2.
    # Validate all elements using a stack-based traversal, mirroring the Perl implementation.
    # Stack entries are dicts: {'name': tag_path, 'elements': [list of ET.Element]}
    # At the top level, each distinct tag maps to a list of matching children.
    # - If more than one child shares a tag, and it's not 'xi:include' or a nested path, it's a duplicate.
    # - Nested paths (containing '->') are sub-parameters where repeats are allowed.

    METADATA_TAGS = {'formatVersion', 'version', 'lastModified'}

    def _split_tag(tag):
        """
        Split an ElementTree tag in Clark notation into (namespace_uri, localname).
        Returns (None, tag) if the tag is not namespaced.
        """
        if tag and tag.startswith('{'):
            uri, local = tag[1:].split('}', 1)
            return uri, local
        return None, tag

    # Group top-level children by tag.
    top_level = defaultdict(list)
    for child in root:
        top_level[child.tag].append(child)

    # Build initial stack from top-level groups.
    stack = [{'name': tag, 'elements': elems} for tag, elems in top_level.items()]

    while stack:
        entry    = stack.pop()
        tag_name = entry['name']
        elements = entry['elements']
        is_nested = '->' in tag_name

        # Skip metadata elements.
        if tag_name in METADATA_TAGS:
            continue

        # Check for duplicate top-level parameters.
        # Treat XInclude <include> elements specially: they may appear multiple times.
        is_xinclude = False
        if elements:
            ns_uri, localname = _split_tag(elements[0].tag)
            if ns_uri == 'http://www.w3.org/2001/XInclude' and localname == 'include':
                is_xinclude = True
        if len(elements) > 1 and not is_nested and not is_xinclude:
            valid = 1
            print(f"Parameter '{tag_name}' appears {len(elements)} times - should appear only once")
            # Still validate each occurrence.

        for elem in elements:
            has_value_attr     = 'value' in elem.attrib
            has_children       = len(elem) > 0
            value_child_elems  = elem.findall('value')

            if not has_value_attr and not has_children:
                # Empty element - missing value.
                valid = 1
                print(f"Parameter '{tag_name}' has no value")
            elif len(value_child_elems) > 1:
                # Multiple <value> child elements.
                valid = 1
                print(f"Parameter '{tag_name}' has multiple values")

            # Push sub-elements onto the stack (grouped by their tag).
            sub_tags = defaultdict(list)
            for child in elem:
                if child.tag == 'value':
                    continue  # Skip <value> text elements.
                sub_tags[child.tag].append(child)
            for sub_tag, sub_elems in sub_tags.items():
                stack.append({'name': tag_name + '->' + sub_tag, 'elements': sub_elems})

sys.exit(valid)
