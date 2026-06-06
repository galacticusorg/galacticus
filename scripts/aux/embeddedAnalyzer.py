#!/usr/bin/env python3
import os
import re
import subprocess
import sys
import xml.etree.ElementTree as ET
import latex_spellcheck

# Perform checks on embedded XML and LaTeX.
# Andrew Benson (28-February-2023) [Python port]

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
warning_file_name = sys.argv[2]

galacticus_path = os.environ.get('GALACTICUS_EXEC_PATH', '.')
schema_dir      = os.path.join(galacticus_path, 'schema')

status   = 0
warnings = ''


# ── LaTeX fragment compilation test ─────────────────────────────────────────

def test_latex(raw_latex):
    """Compile a LaTeX fragment.  Returns the log text on failure, None on success."""
    raw_latex = raw_latex.replace('&amp;', '&').replace('&lt;', '<')
    frag_tex = 'doc/frag.tex'
    with open(frag_tex, 'w') as fh:
        fh.write('\\documentclass[letterpaper,10pt,headsepline]{scrbook}\n')
        fh.write('\\usepackage{natbib}\n')
        fh.write('\\usepackage{epsfig}\n')
        fh.write('\\usepackage[acronym]{glossaries}\n')
        fh.write('\\usepackage[backref,colorlinks]{hyperref}\n')
        fh.write('\\usepackage{amssymb}\n')
        fh.write('\\usepackage{amsmath}\n')
        fh.write('\\usepackage{color}\n')
        fh.write('\\usepackage{listings}\n')
        fh.write('\\usepackage{tensor}\n')
        fh.write(f'\\input{{{galacticus_path}/doc/commands}}\n')
        fh.write(f'\\input{{{galacticus_path}/doc/Glossary}}\n')
        fh.write('\\newcommand{\\docname}{tmp}\n')
        fh.write('\\begin{document}\n')
        fh.write(raw_latex)
        fh.write('\\end{document}\n')

    result = subprocess.run(
        'cd doc; pdflatex -halt-on-error frag > frag.tmp',
        shell=True
    )
    log = None
    if result.returncode != 0:
        try:
            with open('doc/frag.log', 'r', errors='replace') as fh:
                log = fh.read()
        except OSError:
            log = '(log unavailable)'

    for ext in ('tex', 'pdf', 'log', 'aux', 'tmp', 'glo'):
        path = f'doc/frag.{ext}'
        try:
            os.unlink(path)
        except OSError:
            pass
    return log


# ── XML directive processing ─────────────────────────────────────────────────

def process_directive(stripped_directive, line_number):
    """Parse, validate, and spell-check one XML directive block."""
    global status, warnings

    # Parse with ElementTree.
    try:
        root_el = ET.fromstring(stripped_directive)
    except ET.ParseError as exc:
        print(f"Parsing XML fragment failed ({file_name}:{line_number})")
        print(str(exc))
        print(stripped_directive)
        status = 1
        return

    directive_name = root_el.tag

    # Validate against XSD schema if one exists.
    xsd_path = os.path.join(schema_dir, directive_name + '.xsd')
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

    # Check the 'description' attribute/child for LaTeX validity.
    description = root_el.get('description')
    if description is None:
        desc_el = root_el.find('description')
        if desc_el is not None:
            description = (desc_el.text or '')
    if description:
        log = test_latex(description)
        if log is not None:
            print(f"XML LaTeX description compilation failed ({file_name}:{line_number}):")
            print(log)
            status = 1
        warnings += latex_spellcheck.spell_check(description, 'latex', file_name)


# ── Main file parsing ────────────────────────────────────────────────────────

in_directive  = False
in_xml        = False
in_latex      = False
in_comment    = False
line_number   = 0
directive_root        = None
raw_directive         = ''
stripped_directive    = ''
raw_latex             = None
raw_comment           = None

with open(file_name, 'r', errors='replace') as code:
    for line in code:

        # ── State exits (checked before accumulation) ────────────────────────
        if re.match(r'^\s*!!\}', line):
            in_latex = False
        if re.match(r'^\s*!!\]', line):
            in_xml = False
        if not re.match(r'^\s*!\s', line):
            in_comment = False

        # ── LaTeX block accumulation / processing ────────────────────────────
        if in_latex:
            raw_latex = (raw_latex or '') + line
        if raw_latex is not None and not in_latex:
            log = test_latex(raw_latex)
            if log is not None:
                print(f"LaTeX fragment compilation failed ({file_name}:{line_number}):")
                print(log)
                status = 1
            warnings += latex_spellcheck.spell_check(raw_latex, 'latex', file_name)
            raw_latex = None

        # ── Comment block accumulation / processing ──────────────────────────
        if in_comment:
            raw_comment = (raw_comment or '') + line
        if raw_comment is not None and not in_comment:
            warnings += latex_spellcheck.spell_check(raw_comment, 'text', file_name)
            raw_comment = None

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

        # Accumulate text.
        if in_directive:
            raw_directive      += line
            stripped_directive += stripped_line
        elif not re.match(r'^\s*!![\[\]]', line):
            pass  # raw_code accumulation (not needed for our checks)

        # Process a complete directive block.
        if raw_directive and (not in_directive or end_directive):
            process_directive(stripped_directive.strip(), line_number)
            in_directive       = False
            raw_directive      = ''
            stripped_directive = ''
            directive_root     = None

        # ── State entries (checked after processing) ─────────────────────────
        if re.match(r'^\s*!!\[', line):
            in_xml = True
        if re.match(r'^\s*!!\{', line):
            in_latex = True
        if re.match(r'^\s*!\s', line):
            in_comment  = True
            raw_comment = (raw_comment or '') + line

        line_number += 1

# ── Write accumulated warnings ───────────────────────────────────────────────
if warnings:
    with open(warning_file_name, 'a') as fh:
        fh.write(warnings)

raise SystemExit(status)
