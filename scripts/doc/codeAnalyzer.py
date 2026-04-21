#!/usr/bin/env python3
"""Scan Fortran90 source code and output a LaTeX file describing it.

Outputs lists of modules used/used-by with hyperlinks, descriptions, and code
line counts.  Collects much more information on subroutine calls, function
calls, module procedures, etc. — could be made to output that data if required.

Andrew Benson (01-May-2010 original Perl); Python port 2026.

Usage: Code_Analyzer.py <sourceDir> <outputFile>
"""

import os
import re
import sys

_exec_path = os.environ.get(
    'GALACTICUS_EXEC_PATH',
    os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')),
)
sys.path.insert(0, os.path.join(_exec_path, 'python'))
from latex_utils import latex_encode                                   # noqa: E402
from build.fortran_utils import (get_fortran_line, extract_bracketed,  # noqa: E402
                                  extract_variables)

# ---------------------------------------------------------------------------
# Global data structures (populated by process_file, consumed by
# build_modules_hash and output_data).
# ---------------------------------------------------------------------------

# unitID → {unitType, unitName, belongsTo?, codeLines?, contains?,
#            modulesUsed?, usedBy?, subroutinesCalled?, functionCalls?,
#            methods?, comments?, <intrinsic-type lists>…}
units = {}

# moduleName → unitID  (filled by build_modules_hash)
modules = {}

# Directory searched for Fortran include files (set by main before os.walk).
include_file_dir = ''

# ---------------------------------------------------------------------------
# Regex constants
# ---------------------------------------------------------------------------

# use statements:  group(3) is the module name.
USE_RE = re.compile(
    r'^\s*use\s*(,\s*intrinsic\s*)*(::)??\s*([a-z0-9_]+)',
    re.IGNORECASE,
)

# Direct subroutine calls:  group(2) is the subroutine name.
CALL_RE = re.compile(
    r'(^|\W)call\s+([a-z0-9_]+)\s*(\(|\s*$)',
    re.IGNORECASE,
)

# Type-bound subroutine calls:  group(2) variable, group(3) method name.
CALL_TYPE_BOUND_RE = re.compile(
    r'(^|\W)call\s+([a-z0-9_]+)\s*%\s*([a-z0-9_]+)\s*(\(|\s*$)',
    re.IGNORECASE,
)

# Module procedure declarations:  group(1) is the procedure name.
MODULE_PROCEDURE_RE = re.compile(
    r'^\s*module\s+procedure\s+([a-z0-9_]+)\s*$',
    re.IGNORECASE,
)

# Derived-type variable declarations:  group(1) type name, group(3) variable list.
DERIVED_TYPE_RE = re.compile(
    r'^\s*type\s*\(\s*([a-z0-9_]+)\s*\)([\sa-z0-9_,:\+\-\*\/\(\)]*::)*\s*([a-z0-9_,:\+\-\*\/\(\)]+)\s*$',
    re.IGNORECASE,
)

# Type-bound procedure declarations:  group(2) method name, group(3) procedure list.
TYPE_BOUND_RE = re.compile(
    r'^\s*(procedure|generic)\s*::\s*([a-z0-9_]+)\s*=>\s*([a-z0-9_,\s]+)$',
    re.IGNORECASE,
)

# Intrinsic variable declarations.
# Each entry: intrinsic label, compiled regex, and 1-based group indices for
# the kind spec (type_grp), attribute list (attrs_grp), and variable list (vars_grp).
#
# character is special: its outer kind/len spec is group 1, but there are two
# nested capturing groups inside it (groups 2 and 3), which pushes attributes
# to group 4 and variables to group 5.
INTRINSIC_DECLARATIONS = {
    'integer': {
        'intrinsic': 'integer',
        'type_grp':  1,
        'attrs_grp': 2,
        'vars_grp':  3,
        'regex': re.compile(
            r'^\s*integer\s*(\(\s*kind\s*=\s*[a-zA-Z0-9_]+\s*\))*([\sa-zA-Z0-9_,:\+\-\*\/\(\)]*)??::\s*([\sa-zA-Z0-9_,:=>\+\-\*\/\(\)\[\]]+)\s*$',
            re.IGNORECASE,
        ),
    },
    'real': {
        'intrinsic': 'real',
        'type_grp':  1,
        'attrs_grp': 2,
        'vars_grp':  3,
        'regex': re.compile(
            r'^\s*real\s*(\(\s*kind\s*=\s*[a-zA-Z0-9_]+\s*\))*([\sa-zA-Z0-9_,:\+\-\*\/\(\)]*)??::\s*([\sa-zA-Z0-9\._,:=>\+\-\*\/\(\)\[\]]+)\s*$',
            re.IGNORECASE,
        ),
    },
    'double': {
        'intrinsic': 'double precision',
        'type_grp':  1,
        'attrs_grp': 2,
        'vars_grp':  3,
        'regex': re.compile(
            r'^\s*double\s+precision\s*(\(\s*kind\s*=\s*[a-zA-Z0-9_]+\s*\))*([\sa-zA-Z0-9_,:=\+\-\*\/\(\)]*)??::\s*([\sa-zA-Z0-9\._,:=>\+\-\*\/\(\)\[\]]+)\s*$',
            re.IGNORECASE,
        ),
    },
    'logical': {
        'intrinsic': 'logical',
        'type_grp':  1,
        'attrs_grp': 2,
        'vars_grp':  3,
        'regex': re.compile(
            r'^\s*logical\s*(\(\s*kind\s*=\s*[a-zA-Z0-9_]+\s*\))*([\sa-zA-Z0-9_,:\+\-\*\/\(\)]*)??::\s*([\sa-zA-Z0-9_,:=>\+\-\*\/\(\)\[\]]+)\s*$',
            re.IGNORECASE,
        ),
    },
    'character': {
        'intrinsic': 'character',
        'type_grp':  1,
        'attrs_grp': 4,  # two nested groups inside the kind/len spec shift attrs to 4
        'vars_grp':  5,
        'regex': re.compile(
            r'^\s*character\s*(\((\s*(len|kind)\s*=\s*[a-zA-Z0-9_,\+\-\*\(\)]+\s*)+\))*([\sa-zA-Z0-9_,:\+\-\*\/\(\)]*)??::\s*([\sa-zA-Z0-9_,:=>\+\-\*\/\(\)\[\]]+)\s*$',
            re.IGNORECASE,
        ),
    },
    'procedure': {
        'intrinsic': 'procedure',
        'type_grp':  1,
        'attrs_grp': 2,
        'vars_grp':  3,
        'regex': re.compile(
            r'^\s*procedure\s*(\(\s*[a-zA-Z0-9_]*\s*\))*([\sa-zA-Z0-9_,:\+\-\*\/\(\)]*)??::\s*([\sa-zA-Z0-9_,:=>\+\-\*\/\(\)]+)\s*$',
            re.IGNORECASE,
        ),
    },
}

# Unit-opening patterns.  unit_name_grp is the 1-based capturing group that
# holds the unit name.
UNIT_OPENERS = {
    'module': {
        'unit_name_grp': 1,
        'regex': re.compile(
            r'^\s*module\s+(?!procedure\s)([a-z0-9_]+)$',
            re.IGNORECASE,
        ),
    },
    'program': {
        'unit_name_grp': 1,
        'regex': re.compile(
            r'^\s*program\s+([a-z0-9_]+)',
            re.IGNORECASE,
        ),
    },
    'subroutine': {
        'unit_name_grp': 2,
        'regex': re.compile(
            r'^\s*(impure\s+|pure\s+|elemental\s+|recursive\s+|module\s+)*\s*subroutine\s+([a-z0-9_]+)',
            re.IGNORECASE,
        ),
    },
    # Groups: 1=modifiers, 2=return-type keyword, 3=kind-spec, 4=kind/len=, 5=kind|len, 6=name
    'function': {
        'unit_name_grp': 6,
        'regex': re.compile(
            r'^\s*(impure\s+|pure\s+|elemental\s+|recursive\s+)*\s*'
            r'(real\s*|integer\s*|double\s+precision\s*|double\s+complex\s*|character\s*|logical\s*|module\s*)*\s*'
            r'(\(((kind|len)=)??\w*\))*\s*function\s+([a-z0-9_]+)',
            re.IGNORECASE,
        ),
    },
    # Groups: 1=abstract?, 2=interface name (may be empty)
    'interface': {
        'unit_name_grp': 2,
        'regex': re.compile(
            r'^\s*(abstract\s+)??interface\s+([a-z0-9_()/\+\-\*\.=]*)',
            re.IGNORECASE,
        ),
    },
    # Groups: 1=attributes?, 2=::?, 3=type name
    'type': {
        'unit_name_grp': 3,
        'regex': re.compile(
            r'^\s*type\s*(,\s*abstract\s*|,\s*public\s*|,\s*private\s*|,\s*extends\s*\([a-zA-Z0-9_]+\)\s*|,\s*bind\(c\)\s*)*(::)??\s*([a-z0-9_{}\u00a6]+)\s*$',
            re.IGNORECASE,
        ),
    },
}

# Unit-closing patterns.  unit_name_grp is the 1-based capturing group for the
# unit name; interface has no capturing group so unit_name_grp=0 signals that.
UNIT_CLOSERS = {
    'module': {
        'unit_name_grp': 1,
        'regex': re.compile(r'^\s*end\s+module\s+([a-z0-9_]+)', re.IGNORECASE),
    },
    'program': {
        'unit_name_grp': 1,
        'regex': re.compile(r'^\s*end\s+program\s+([a-z0-9_]+)', re.IGNORECASE),
    },
    'subroutine': {
        'unit_name_grp': 1,
        'regex': re.compile(r'^\s*end\s+subroutine\s+([a-z0-9_]+)', re.IGNORECASE),
    },
    'function': {
        'unit_name_grp': 1,
        'regex': re.compile(r'^\s*end\s+function\s+([a-z0-9_]+)', re.IGNORECASE),
    },
    'interface': {
        'unit_name_grp': 0,  # no capturing group; name check skipped via unitType match
        'regex': re.compile(r'^\s*end\s+interface', re.IGNORECASE),
    },
    'type': {
        'unit_name_grp': 1,
        'regex': re.compile(r'^\s*end\s+type\s+([a-z0-9_{}\u00a6]+)', re.IGNORECASE),
    },
}

# ---------------------------------------------------------------------------
# Core functions
# ---------------------------------------------------------------------------

def process_file(file_path):
    """Scan one source file and populate the global `units` dict."""
    global units

    # File::Find delivers the bare filename; replicate that as the unit ID so
    # that nested unit IDs (parent:child) stay short and LaTeX-friendly.
    file_name = os.path.basename(file_path)

    # Only process Fortran (.F90 / .f90) and C++ (.cpp) files; skip editor
    # lock files (.\#...).
    is_fortran = bool(re.search(r'\.[fF]90$', file_name))
    is_cpp     = file_name.endswith('.cpp')
    if not (is_fortran or is_cpp) or file_name.startswith('.#'):
        return

    # Root unit entry for this file.
    unit_id_list = [file_name]
    units[file_name] = {'unitType': 'file', 'unitName': file_name}

    # C++ files get an entry but no further parsing.
    if not is_fortran:
        return

    # File stack: each frame is {name, position, in_xml, in_latex}.
    # position=None means "open from the beginning".
    file_stack = [{'name': file_path, 'position': None,
                   'in_xml': False, 'in_latex': False}]

    while file_stack:
        frame = file_stack[0]

        with open(frame['name'], 'r', errors='replace') as fh:
            if frame['position'] is not None:
                fh.seek(frame['position'])

            while True:
                raw, processed, comments = get_fortran_line(fh)
                if not raw:
                    # EOF — remove this frame and resume the parent (if any).
                    file_stack.pop(0)
                    break

                # Detect leaving LaTeX / XML documentation blocks.
                if re.match(r'^\s*!!\]', raw):
                    frame['in_xml']   = False
                if re.match(r'^\s*!!\}', raw):
                    frame['in_latex'] = False

                # Every physical line counts toward code-line totals for all
                # currently open units.
                for uid in unit_id_list:
                    units[uid]['codeLines'] = units[uid].get('codeLines', 0) + 1

                # Detect Fortran include directives and recurse into the file.
                m = re.match(r'''^\s*include\s*['"]([^'"]+)['"]\s*$''', processed)
                if m:
                    inc_path = os.path.join(include_file_dir, m.group(1))
                    if os.path.exists(inc_path):
                        # Save position so we resume after the include line
                        # when we return to this file.
                        frame['position'] = fh.tell()
                        file_stack.insert(0, {'name': inc_path, 'position': None,
                                              'in_xml': False, 'in_latex': False})
                        break  # break inner loop; outer while opens the include

                line_processed = False

                # ---- Unit opening ----------------------------------------
                if not line_processed and not frame['in_xml'] and not frame['in_latex']:
                    for unit_type, opener in UNIT_OPENERS.items():
                        m = opener['regex'].match(processed)
                        if m:
                            grp = opener['unit_name_grp']
                            unit_name = (m.group(grp).lower()
                                         if grp > 0 and m.lastindex and m.lastindex >= grp
                                         else '')
                            parent_id = unit_id_list[-1]
                            unit_id   = f'{parent_id}:{unit_name}'
                            # Strip placeholder chars used in generated names.
                            unit_id = re.sub(r'[¦{}]', '', unit_id)
                            units.setdefault(parent_id, {}).setdefault('contains', {})[unit_id] = True
                            unit_id_list.append(unit_id)
                            units[unit_id] = {
                                'unitType':  unit_type,
                                'unitName':  unit_name,
                                'belongsTo': parent_id,
                            }
                            line_processed = True
                            break

                # ---- Unit closing ----------------------------------------
                if not line_processed and not frame['in_xml'] and not frame['in_latex']:
                    for unit_type, closer in UNIT_CLOSERS.items():
                        m = closer['regex'].match(processed)
                        if m:
                            grp = closer['unit_name_grp']
                            unit_name   = (m.group(grp).lower()
                                           if grp > 0 and m.lastindex and m.lastindex >= grp
                                           else '')
                            opener_id   = unit_id_list[-1]
                            opener_type = units.get(opener_id, {}).get('unitType', '')
                            opener_name = units.get(opener_id, {}).get('unitName', '')
                            if not (unit_type == opener_type
                                    and (unit_name == opener_name
                                         or 'interface' in unit_type.lower())):
                                print('Unit close does not match unit open:')
                                print(f' Closing with: {unit_type} {unit_name}')
                                print(f'  Opened with: {opener_type} {opener_name}')
                                print(f'      In file: {file_name}')
                                sys.exit(1)
                            unit_id_list.pop()
                            line_processed = True
                            break

                # ---- Description comment block  !!{ ... !!} ---------------
                # buffered_comments starts with '!{' when the current line is
                # '!!{...': the comment content (everything after the first '!')
                # begins with another '!' followed by '{'.
                if not line_processed and not frame['in_xml']:
                    if re.match(r'^!\{', comments):
                        uid = unit_id_list[-1]
                        units.setdefault(uid, {})
                        while True:
                            comment_line = fh.readline()
                            if not comment_line:
                                break
                            if re.match(r'^\s*!!\}', comment_line):
                                break
                            comment_line = re.sub(r'^\s*!', '', comment_line)
                            units[uid]['comments'] = (
                                units[uid].get('comments', '') + comment_line
                            )
                        frame['in_latex'] = False
                        line_processed = True

                # ---- use statement ---------------------------------------
                if not line_processed and not frame['in_xml'] and not frame['in_latex']:
                    m = USE_RE.match(processed)
                    if m:
                        module_name = m.group(3).lower()
                        uid = unit_id_list[-1]
                        units.setdefault(uid, {}).setdefault('modulesUsed', {})[module_name] = True
                        line_processed = True

                # ---- Direct subroutine call ------------------------------
                # CALL_RE is not anchored so use search(); must come before
                # type-bound check (the two are mutually exclusive by regex
                # design — CALL_RE requires name followed by '(' or EOL).
                if not line_processed and not frame['in_xml'] and not frame['in_latex']:
                    m = CALL_RE.search(processed)
                    if m:
                        sub_name = m.group(2).lower()
                        uid = unit_id_list[-1]
                        units.setdefault(uid, {}).setdefault('subroutinesCalled', {})[sub_name] = -1
                        line_processed = True

                # ---- Type-bound subroutine call --------------------------
                if not line_processed and not frame['in_xml'] and not frame['in_latex']:
                    m = CALL_TYPE_BOUND_RE.search(processed)
                    if m:
                        sub_name  = m.group(3).lower()
                        var_name  = m.group(2).lower()
                        uid = unit_id_list[-1]
                        units.setdefault(uid, {}).setdefault('subroutinesCalled', {})[sub_name] = var_name
                        line_processed = True

                # ---- Type-bound procedure declaration --------------------
                if not line_processed and not frame['in_xml'] and not frame['in_latex']:
                    if units.get(unit_id_list[-1], {}).get('unitType') == 'type':
                        m = TYPE_BOUND_RE.match(processed)
                        if m:
                            method_name    = m.group(2).lower()
                            procedure_list = re.sub(r'\s', '', m.group(3).lower())
                            procedures     = procedure_list.split(',')
                            uid = unit_id_list[-1]
                            units.setdefault(uid, {}).setdefault('methods', {})[method_name] = procedures
                            line_processed = True

                # ---- Module procedure ------------------------------------
                if not line_processed and not frame['in_xml'] and not frame['in_latex']:
                    m = MODULE_PROCEDURE_RE.match(processed)
                    if m:
                        proc_name = m.group(1).lower()
                        uid = unit_id_list[-1]
                        units.setdefault(uid, {}).setdefault('moduleProcedures', {})[proc_name] = True
                        line_processed = True

                # ---- Intrinsic variable declarations --------------------
                if not line_processed and not frame['in_xml'] and not frame['in_latex']:
                    for intr_type, intr in INTRINSIC_DECLARATIONS.items():
                        m = intr['regex'].match(processed)
                        if m:
                            vars_grp = intr['vars_grp']
                            vars_str = m.group(vars_grp) if m.lastindex and m.lastindex >= vars_grp else None
                            if vars_str:
                                uid = unit_id_list[-1]
                                variables = extract_variables(vars_str.lower())
                                units.setdefault(uid, {}).setdefault(intr_type, []).extend(variables)
                            line_processed = True
                            break

                # ---- Derived-type variable declarations ------------------
                if not line_processed and not frame['in_xml'] and not frame['in_latex']:
                    m = DERIVED_TYPE_RE.match(processed)
                    if m:
                        derived_type = m.group(1).lower()
                        vars_str     = m.group(3)
                        if vars_str:
                            uid = unit_id_list[-1]
                            variables = extract_variables(vars_str.lower())
                            (units.setdefault(uid, {})
                                  .setdefault('derivedTypesUsed', {})
                                  .setdefault(derived_type, [])
                                  .extend(variables))
                        line_processed = True

                # ---- Function call extraction ----------------------------
                # Not mutually exclusive with the checks above: does NOT set
                # line_processed (matching the Perl behaviour).
                if not line_processed and not frame['in_xml'] and not frame['in_latex']:
                    func_seek = processed.lower()
                    # Strip string literals to avoid false matches inside them.
                    func_seek = re.sub(r"''", '', func_seek)
                    func_seek = re.sub(r'""', '', func_seek)
                    func_seek = re.sub(r"'[^']+'", '', func_seek)
                    func_seek = re.sub(r'"[^"]+"', '', func_seek)
                    count_iters = 0
                    while '(' in func_seek:
                        count_iters += 1
                        if count_iters > 1000:
                            print('Code_Analyzer.py: exceeded 1000 iterations '
                                  'extracting function calls — possible parse failure',
                                  file=sys.stderr)
                            sys.exit(1)
                        extracted, remainder, prefix = extract_bracketed(func_seek, "()")
                        if prefix is None:
                            print('failed to find function name:')
                            print(f'\traw line: {raw}')
                            print(f'\tprocessed line: {processed}')
                            print(f'\tbuffered comments: {comments}')
                            print(f'\tremaining line: {func_seek}')
                            sys.exit(1)
                        m = re.search(r'(([a-z0-9_]+)\s*%\s*)([a-z0-9_]+)$', prefix)
                        if m:
                            func_name    = m.group(3)
                            derived_type = re.sub(r'%\s*$', '', m.group(1))
                            if not derived_type:
                                derived_type = -1
                            uid = unit_id_list[-1]
                            units.setdefault(uid, {}).setdefault('functionCalls', {})[func_name] = derived_type
                        func_seek = remainder

                # ---- Detect entering XML / LaTeX documentation blocks -----
                # (must run after all other checks so line_processed is final)
                if re.match(r'^\s*!!\[', raw):
                    frame['in_xml'] = True
                if not line_processed and re.match(r'^\s*!!\{', raw):
                    frame['in_latex'] = True


def build_modules_hash():
    """Build `modules` dict and record usedBy back-references in `units`."""
    global modules
    # Map each module name to its unit ID.
    for unit_id, unit in units.items():
        if unit.get('unitType') == 'module':
            modules[unit['unitName']] = unit_id
    # Record which units each module is used by.
    for unit_id, unit in units.items():
        for mod_name in unit.get('modulesUsed', {}):
            if mod_name in modules:
                mod_unit_id = modules[mod_name]
                units[mod_unit_id].setdefault('usedBy', {})[unit_id] = True


def print_two_column(fh, items):
    """Write items in a two-column supertabular layout."""
    col = 0
    for item in items:
        fh.write(f' & \\RaggedRight {item}')
        if col == 1:
            fh.write('\\\\\n')
        col = 1 - col
    if col == 1:
        fh.write(' & \\\\\n')


def output_data(output_file):
    """Write the LaTeX output file."""
    table_open_str = '\\begin{supertabular}{lp{70mm}p{70mm}}\n'

    with open(output_file, 'w', encoding='utf-8') as fh:
        fh.write('\\section{Program units}\n')

        for unit_id in sorted(units):
            unit      = units[unit_id]
            unit_type = unit.get('unitType', '')

            # Encode the unit name for LaTeX; add soft hyphen before each
            # escaped underscore so long identifiers can line-wrap.
            unit_name = latex_encode(unit.get('unitName', ''))
            unit_name = unit_name.replace('\\_', '\\-\\_')
            unit_name = unit_name.replace('¦', '\\textbrokenbar{}')

            parent_id = unit.get('belongsTo')

            # Skip abstract interfaces (unnamed interface units) and their
            # direct children.
            is_abstract = (unit_type == 'interface' and unit_name == '')
            parent_is_abstract = False
            if parent_id:
                p = units.get(parent_id, {})
                parent_is_abstract = (p.get('unitType') == 'interface'
                                      and latex_encode(p.get('unitName', '')) == '')
            if is_abstract or parent_is_abstract:
                continue

            # Header line.
            hyperdef_target = unit_id.replace('.', '_')
            fh.write(
                f'\\noindent{{\\normalfont \\bfseries {unit_type}:}}'
                f' \\hypertarget{{{unit_id}}}{{\\mono{{{unit_name}}}}}'
                f'\\hyperdef{{source}}{{{hyperdef_target}}}{{}}'
                f'\\index[code]{{{unit_name}@\\mono{{{unit_name}}} ({unit_type})}}\n\n'
            )

            table_is_open = False

            def open_table():
                nonlocal table_is_open
                if not table_is_open:
                    fh.write(table_open_str)
                    table_is_open = True

            # Description (from !!{ ... !!} comment blocks).
            if 'comments' in unit:
                comments = unit['comments']
                if '{verbatim}' not in comments:
                    open_table()
                    fh.write(
                        '\\emph{Description:} & \\multicolumn{2}{l}{\n'
                        '\\begin{minipage}[t]{140mm}\n'
                        f'{comments}\\end{{minipage}}\n}}\\\\\n'
                    )
                else:
                    fh.write(f'\\emph{{Description:}} {comments}')
                    if not re.search(r'\}\s*$', comments):
                        fh.write('\\\\')
                    fh.write('\n')

            # Code line count.
            open_table()
            if 'codeLines' in unit:
                fh.write(f'\\emph{{Code lines:}} & \\multicolumn{{2}}{{l}}{{{unit["codeLines"]}}} \\\\\n')

            # Parent unit.
            if parent_id:
                open_table()
                p = units.get(parent_id, {})
                p_type = p.get('unitType', '')
                p_name = latex_encode(p.get('unitName', ''))
                fh.write(
                    f'\\emph{{Contained by:}} & \\multicolumn{{2}}{{l}}'
                    f'{{{p_type} \\hyperlink{{{parent_id}}}{{\\mono{{{p_name}}}}}}}'
                    f' \\\\ \n'
                )

            # Modules used.
            if unit.get('modulesUsed'):
                open_table()
                fh.write('\\emph{Modules used:} ')
                mod_items = []
                for mod_name in sorted(unit['modulesUsed']):
                    enc = latex_encode(mod_name)
                    if mod_name in modules and modules[mod_name]:
                        mod_items.append(f'\\hyperlink{{{modules[mod_name]}}}{{\\mono{{{enc}}}}}')
                    else:
                        mod_items.append(f'\\mono{{{enc}}}')
                print_two_column(fh, mod_items)

            # Units that use this one.
            if unit.get('usedBy'):
                open_table()
                fh.write('\\emph{Used by:} ')
                used_items = []
                for uid in sorted(unit['usedBy']):
                    u2 = units.get(uid, {})
                    u2_type = u2.get('unitType', '')
                    u2_name = latex_encode(u2.get('unitName', ''))
                    if uid:
                        used_items.append(
                            f'{u2_type} \\hyperlink{{{uid}}}{{\\mono{{{u2_name}}}}}'
                        )
                    else:
                        used_items.append(f'{u2_type} \\mono{{{u2_name}}}')
                print_two_column(fh, used_items)

            # Close table or emit blank line.
            if table_is_open:
                fh.write(' & & \\\\\n\\end{supertabular}\n\\\\\n')
            else:
                fh.write('\n')


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    if len(sys.argv) != 3:
        print('Usage: Code_Analyzer.py <sourceDir> <outputFile>', file=sys.stderr)
        sys.exit(1)

    source_dir  = sys.argv[1]
    output_file = sys.argv[2]

    build_path = os.environ.get('BUILDPATH')
    if not build_path:
        print('Error: BUILDPATH environment variable is not set.', file=sys.stderr)
        sys.exit(1)

    # Replicates Perl: $currentDirectory/$sourceDir/../$BUILDPATH
    global include_file_dir
    include_file_dir = os.path.normpath(
        os.path.join(os.getcwd(), source_dir, '..', build_path)
    )

    # Walk the source tree and process every eligible file.
    for root, _dirs, files in os.walk(source_dir):
        for fname in files:
            process_file(os.path.join(root, fname))

    build_modules_hash()
    output_data(output_file)


if __name__ == '__main__':
    main()
