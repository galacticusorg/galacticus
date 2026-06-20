#!/usr/bin/env python3
import os
import re
import sys
import xml.etree.ElementTree as ET

# Validate embedded XML directives: schema-check each `!![ … !!]` directive block.
# Andrew Benson (28-February-2023) [Python port].
#
# The LaTeX checks (fragment compilation + spell-check) were retired with the
# migration of the embedded documentation to reStructuredText; prose validity is
# now covered by the ReadTheDocs build and the RST spell-checker.

# Optional lxml for XSD validation; fall back to skipping validation if absent.
try:
    from lxml import etree as lxml_etree
    HAS_LXML = True
except ImportError:
    HAS_LXML = False

if len(sys.argv) != 3:
    print("Usage: embeddedAnalyzer.py <fileName> <warningFile>", file=sys.stderr)
    raise SystemExit(1)

file_name         = sys.argv[1]
warning_file_name = sys.argv[2]   # retained for CLI compatibility (unused)

galacticus_path = os.environ.get('GALACTICUS_EXEC_PATH', '.')
schema_dir      = os.path.join(galacticus_path, 'schema')

status = 0


# ── XML directive processing ─────────────────────────────────────────────────

def process_directive(stripped_directive, line_number):
    """Parse and (if a schema exists) XSD-validate one XML directive block."""
    global status

    try:
        root_el = ET.fromstring(stripped_directive)
    except ET.ParseError as exc:
        print(f"Parsing XML fragment failed ({file_name}:{line_number})")
        print(str(exc))
        print(stripped_directive)
        status = 1
        return

    xsd_path = os.path.join(schema_dir, root_el.tag + '.xsd')
    if HAS_LXML and os.path.exists(xsd_path):
        try:
            schema   = lxml_etree.XMLSchema(file=xsd_path)
            document = lxml_etree.fromstring(stripped_directive.encode())
            schema.assertValid(document)
        except (lxml_etree.XMLSchemaError,
                lxml_etree.DocumentInvalid,
                lxml_etree.XMLSyntaxError) as exc:
            print(f"XML fragment validation failed ({file_name}:{line_number}):")
            print(str(exc))
            status = 1


# ── Main file parsing ────────────────────────────────────────────────────────

in_xml        = False
in_directive  = False
line_number   = 0
directive_root     = None
raw_directive      = ''
stripped_directive = ''

with open(file_name, 'r', errors='replace') as code:
    for line in code:

        # Exit an XML directive block.
        if re.match(r'^\s*!!\]', line):
            in_xml = False

        # ── XML directive accumulation / processing ──────────────────────────
        is_directive  = False
        end_directive = False

        # Strip the optional '!< ' prefix used on some XML-comment directive lines.
        stripped_line = re.sub(r'^\s*\!<\s*', '', line)
        stripped_line = stripped_line.replace('&nbsp;', ' ')

        if in_xml:
            m_open = re.match(r'^\s*<([^\s>/]+)', stripped_line)
            if m_open or in_directive:
                is_directive = True
            if is_directive and not in_directive and m_open:
                directive_root = m_open.group(1)
            if is_directive and directive_root:
                # Detect closing of the directive's root element.
                if re.search(r'</' + re.escape(directive_root) + r'>', stripped_line):
                    end_directive = True
                # Self-closing root element.
                if not in_directive and (
                    re.search(r'<' + re.escape(directive_root) + r'\s[^>]*/>', stripped_line)
                    or re.search(r'<' + re.escape(directive_root) + r'/>', stripped_line)
                ):
                    end_directive = True
            if is_directive:
                in_directive = True

        # Accumulate directive text.
        if in_directive:
            raw_directive      += line
            stripped_directive += stripped_line

        # Process a complete directive block.
        if raw_directive and (not in_directive or end_directive):
            process_directive(stripped_directive.strip(), line_number)
            in_directive       = False
            raw_directive      = ''
            stripped_directive = ''
            directive_root     = None

        # Enter an XML directive block.
        if re.match(r'^\s*!!\[', line):
            in_xml = True

        line_number += 1

raise SystemExit(status)
