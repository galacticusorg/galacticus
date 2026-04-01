#!/usr/bin/env python3
import os
import re
import sys
import xml.etree.ElementTree as ET

# Monitor for updates to certain files, outputting warnings if so.
# Andrew Benson (09-April-2025) [Python port]

galacticus_path = os.environ.get('GALACTICUS_EXEC_PATH', '.')
files_changed   = sys.argv[1:]


# ── Directive extraction (mirrors Galacticus::Build::Directives) ─────────────

def _extract_directives(file_name):
    """Yield every XML directive found inside !![ ... !!] blocks in file_name.
    Each yielded value is a dict with keys 'rootElementType' and any XML
    attributes of the root element (notably 'name')."""
    if not os.path.exists(file_name):
        return
    in_xml    = False
    xml_lines = []
    depth     = 0
    try:
        with open(file_name, 'r', errors='replace') as fh:
            for line in fh:
                # End of XML section — flip state before accumulating.
                if re.match(r'^\s*!!\]', line):
                    in_xml = False

                if in_xml or depth > 0:
                    # Strip leading instrumentation prefix and normalise.
                    processed = re.sub(r'^(\!<)?\s*', '', line) if in_xml else line
                    processed = processed.replace('&nbsp;', ' ')
                    xml_lines.append(processed)

                    # Track element depth (opening minus closing tags; self-closing = 0).
                    depth += len(re.findall(r'<[a-zA-Z0-9]+[^/>]*>', processed))
                    depth -= len(re.findall(r'</[a-zA-Z0-9]+>',       processed))

                    if xml_lines and depth == 0:
                        xml_text = ''.join(xml_lines).strip()
                        xml_lines = []
                        if xml_text:
                            # Wrap in a synthetic root so ET can parse multiple
                            # sibling elements from one block.
                            try:
                                wrapper = ET.fromstring(
                                    f'<_root_>{xml_text}</_root_>'
                                )
                                for child in wrapper:
                                    directive = dict(child.attrib)
                                    directive['rootElementType'] = child.tag
                                    yield directive
                            except ET.ParseError:
                                pass

                # Start of XML section.
                if re.match(r'^\s*!!\[', line):
                    in_xml    = True
                    xml_lines = []
                    depth     = 0
    except (IOError, OSError):
        pass


# ── Load watches ─────────────────────────────────────────────────────────────

watches_path = os.path.join(galacticus_path, 'scripts', 'aux', 'watches.xml')
watches      = []
try:
    tree = ET.parse(watches_path)
    for watch_el in tree.getroot().findall('watch'):
        watches.append(dict(watch_el.attrib))
except (ET.ParseError, OSError):
    pass


# ── Check each changed file ───────────────────────────────────────────────────

warnings = ''

for file_name in files_changed:
    # File-level watches.
    for watch in watches:
        if 'file' in watch and watch['file'] == file_name:
            warnings += (
                f":warning: File `{file_name}` has changed. "
                f"{watch['message']}\n"
            )

    # Directive-level watches.
    for directive in _extract_directives(file_name):
        for watch in watches:
            if 'type' not in watch:
                continue
            if (directive.get('rootElementType') == watch['type']
                    and directive.get('name') == watch.get('name')):
                warnings += (
                    f":warning: File `{file_name}` has changed. "
                    f"{watch['message']}\n"
                )

# ── Write output ─────────────────────────────────────────────────────────────

if warnings:
    out_path = os.path.join(galacticus_path, 'warnings.md')
    with open(out_path, 'w') as fh:
        fh.write('# Watched files have changed\n')
        fh.write(warnings)
