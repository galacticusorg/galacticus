#!/usr/bin/env python3
import os
import re
import sys
import xml.etree.ElementTree as ET

# Format constant definitions for incorporation into the documentation.
# Andrew Benson (29-January-2025)

if len(sys.argv) != 3:
    print("Usage: constants.py <buildPath> <outputFile>", file=sys.stderr)
    raise SystemExit(1)

build_path  = sys.argv[1]
output_file = sys.argv[2]

# Initialize a list to contain constants.
constants = []

# Process all files.
for root, dirs, files in os.walk(build_path):
    for fname in files:
        if not fname.endswith('.constants.xml'):
            continue
        full_path = os.path.join(root, fname)
        try:
            tree   = ET.parse(full_path)
        except ET.ParseError:
            continue
        root_el = tree.getroot()
        for constant_el in root_el.findall('constant'):
            # Build a dict from attributes (XML::Simple-style flattening).
            constant = dict(constant_el.attrib)
            # Also include any child element text.
            for child in constant_el:
                if child.text:
                    constant[child.tag] = child.text
            constants.append(constant)

# Define known groups.
known_groups = {
    "astrophysical": "Astrophysical constants",
    "atomic":        "Atomic physics constants",
    "physical":      "Physical constants",
    "math":          "Mathematical constants",
    "units":         "Units",
    "prefixes":      "SI Prefixes",
    "GSL":           "Gnu Scientific Library constants",
    "Kernel":        "Kernel constants",
    "misc":          "Miscellaneous constants",
}

# Find all groups.
groups = {}
for constant in constants:
    if 'value' not in constant:
        continue
    if 'group' in constant:
        for group_name in constant['group'].split(':'):
            if group_name not in known_groups:
                raise SystemExit(f"Unknown group name '{group_name}'")
            groups.setdefault(group_name, []).append(constant)
    else:
        groups.setdefault('misc', []).append(constant)

# Iterate over groups and constants, formatting into LaTeX.
with open(output_file, 'w') as out:
    for group_name in sorted(groups.keys()):
        out.write(f"\\subsection{{{known_groups[group_name]}}}\n")
        sorted_constants = sorted(groups[group_name], key=lambda c: c.get('variable', ''))
        for constant in sorted_constants:
            if 'value' not in constant:
                continue
            ref_text = constant.get('reference', '')
            if 'referenceURL' in constant and ref_text:
                reference = f"\\href{{{constant['referenceURL']}}}{{{ref_text}}}"
            else:
                reference = ref_text
            reference = reference.replace('&', '\\&')
            variable  = constant.get('variable', '').replace('_', '\\_')
            external  = (
                f" (See \\href{{{constant['externalDescription']}}}{{here}}.)"
                if 'externalDescription' in constant else ""
            )
            symbol    = f"${constant['symbol']}$: " if 'symbol' in constant else ""
            units_val = constant.get('units', '')
            units     = f" [{units_val}]" if units_val and units_val != "dimensionless" else ""
            value     = constant.get('value', '')
            value     = re.sub(r'd([\+\-0-9]+)$', r'e\1', value)
            value     = re.sub(r'_[_a-zA-Z0-9]+$', '', value)
            module    = constant.get('module', '').replace('_', '\\_')
            file_name = constant.get('fileName', '').replace('.', '_')
            module_url = (
                "https://github.com/galacticusorg/galacticus/releases/download/masterRelease/"
                f"Galacticus_Source.pdf\\#source.{file_name}:{constant.get('module', '').lower()}"
            )
            out.write(f"\\noindent \\mono{{{variable} = {value}}}{units}\\\\\n")
            out.write(f"\\indent {symbol}{constant.get('description', '')}{external}\\\\\n")
            out.write(f"\\indent Module: \\href{{{module_url}}}\\mono{{{module}}}\\\\\n")
            out.write(f"\\indent Reference: {reference}\n\n\\medskip\n")
