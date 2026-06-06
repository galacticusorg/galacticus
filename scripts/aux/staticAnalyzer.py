#!/usr/bin/env python3
import re
import sys

# Perform static analysis of Fortran files.
# Andrew Benson (28-February-2023) [Python port]
#
# Checks performed (mirrors staticAnalyzer.pl):
#  1. class/type pointer members in derived types that are not null-initialized.
#  2. Duplicated variable assignments in <constructorAssign> directives.
#  3. Empty constructor functions (no-argument constructors for known types).
#  4. Empty finalizer subroutines.

if len(sys.argv) != 2:
    print("Usage: staticAnalyzer.py <fileName>", file=sys.stderr)
    raise SystemExit(1)

file_name = sys.argv[1]
status    = 0

# ── Fortran regex patterns ────────────────────────────────────────────────────

# Variable declaration (from Fortran::Utils::variableDeclarationRegEx).
_VAR_DECL = re.compile(
    r'^\s*(?:!\$\s*)?'
    r'(?:integer|real|double\s+precision|logical|character|type|class|'
    r'complex|procedure|generic)'
    r'(?:\s*\(\s*[a-zA-Z0-9_=\*]+\s*\))*'
    r'(?:[\sa-zA-Z0-9_,%:\+\-\*\/\(\)]*)?'
    r'::\s*[\sa-zA-Z0-9\._,:=>\+\-\*\/\(\)\[\]]+\s*$',
    re.IGNORECASE,
)

# Derived-type definition opener: "type [,attrs] [::] name"  (not "type(...)").
_TYPE_DEF = re.compile(
    r'^\s*type'
    r'(?:\s*,\s*(?:abstract|public|private|extends\s*\([^)]+\))\s*)*'
    r'(?:\s*::)?\s*([a-zA-Z0-9_]+)\s*$',
    re.IGNORECASE,
)
_END_TYPE      = re.compile(r'^\s*end\s+type\b',      re.IGNORECASE)
_INTERFACE     = re.compile(r'^\s*interface\s+([a-zA-Z0-9_]+)\s*$', re.IGNORECASE)
_END_INTERFACE = re.compile(r'^\s*end\s+interface\b', re.IGNORECASE)
_MODULE_PROC   = re.compile(r'^\s*module\s+procedure\s+(.*)', re.IGNORECASE)

# Function opener — captures name and argument list.
_FUNCTION = re.compile(
    r'^\s*(?:(?:pure|elemental|recursive|module)\s+)*'
    r'(?:(?:real|integer|double\s+precision|double\s+complex|character|logical)'
    r'\s*(?:\([^)]*\))?\s*)?'
    r'function\s+([a-zA-Z0-9_]+)\s*(\([^)]*\))',
    re.IGNORECASE,
)
_MODULE_FUNCTION = re.compile(r'^\s*module\s+function\b', re.IGNORECASE)
_END_FUNCTION    = re.compile(r'^\s*end\s+function\b',    re.IGNORECASE)

_SUBROUTINE    = re.compile(
    r'^\s*(?:(?:pure|elemental|recursive|module)\s+)*subroutine\s+([a-zA-Z0-9_]+)\s*\(',
    re.IGNORECASE,
)
_END_SUBROUTINE = re.compile(r'^\s*end\s+subroutine\b', re.IGNORECASE)

# Pointer member in a derived-type body.
_POINTER_MEMBER = re.compile(
    r'^\s*(?:type|class)\s*\([^)]+\)'      # type(foo) or class(foo)
    r'(?:\s*,\s*[^:,]+)*'                  # optional attributes before pointer
    r'\s*,\s*pointer'                       # , pointer
    r'(?:\s*,\s*[^:,]+)*'                  # optional attributes after pointer
    r'\s*::\s*(.+)',                        # :: variable_list
    re.IGNORECASE,
)
_NULL_INIT    = re.compile(r'=>\s*null\s*\(\s*\)', re.IGNORECASE)
_FINAL_DECL   = re.compile(r'^\s*final\s*::\s*(.+)', re.IGNORECASE)

# Lines that make a function/subroutine body non-trivially empty.
_BLANK        = re.compile(r'^\s*$')
_COMMENT      = re.compile(r'^\s*![^\$]')          # comment, not OpenMP
_USE_STMT     = re.compile(r'^\s*use(\s*,\s*intrinsic)??\s*::', re.IGNORECASE)
_IMPLICIT     = re.compile(r'^\s*implicit\s+none\s*$', re.IGNORECASE)
_RETURN       = re.compile(r'^\s*return\s*$',          re.IGNORECASE)


# ── Helpers ───────────────────────────────────────────────────────────────────

def _join_continuations(lines):
    """Join Fortran continuation lines (trailing &) into single logical lines."""
    result  = []
    current = None
    for raw in lines:
        line = raw.rstrip('\n')
        if current is not None:
            # The next line may optionally start with & (ignored).
            stripped = line.lstrip()
            if stripped.startswith('&'):
                stripped = stripped[1:]
            current = current + stripped
        else:
            current = line
        if current.rstrip().endswith('&'):
            current = current.rstrip()[:-1]    # drop trailing &
        else:
            result.append(current)
            current = None
    if current is not None:
        result.append(current)
    return result


def _is_trivial_body_line(line):
    """Return True if a line in a function/subroutine body is ignorable for the
    purposes of the 'empty function' check."""
    return (
        _BLANK.match(line)
        or _COMMENT.match(line)
        or _USE_STMT.match(line)
        or _IMPLICIT.match(line)
        or _RETURN.match(line)
        or _VAR_DECL.match(line)
        or _END_SUBROUTINE.match(line)
        or _END_FUNCTION.match(line)
    )


def _function_is_empty(body_lines):
    """Return True if the body (a list of logical lines) contains no meaningful
    code — i.e. every line is trivial AND there are no XML directive blocks."""
    in_docstring = False
    for line in body_lines:
        if re.match(r'^\s*!!\[', line):
            # Any XML directive block makes the function non-empty.
            return False
        if re.match(r'^\s*!!\}', line):
            in_docstring = False
            continue
        if re.match(r'^\s*!!\{', line):
            in_docstring = True
            continue
        if in_docstring:
            continue
        if not _is_trivial_body_line(line):
            return False
    return True


# ── Parse the file ────────────────────────────────────────────────────────────

with open(file_name, 'r', errors='replace') as fh:
    raw_lines = fh.readlines()

lines = _join_continuations(raw_lines)

# Collected facts.
known_types       = []      # lowercase type names defined in this file
constructors      = []      # lowercase names of constructor functions
destructors       = []      # lowercase names of finalizer subroutines

# Parsing state.
in_type_block      = False
current_type_name  = None
in_interface_block = False
current_iface_name = None   # lowercase
in_function        = False
current_func_name  = None   # lowercase
func_is_module     = False
func_args_empty    = False
func_body          = []
in_subroutine      = False
current_sub_name   = None   # lowercase
sub_body           = []
in_xml_block       = False

for line in lines:
    # ── XML block state ──────────────────────────────────────────────────────
    if re.match(r'^\s*!!\[', line):
        in_xml_block = True
    if re.match(r'^\s*!!\]', line):
        in_xml_block = False

    # ── Check 2: constructorAssign duplicate variables ───────────────────────
    # These appear inside !![ ... !!] blocks.
    if in_xml_block:
        m = re.match(r'^\s*<constructorAssign\s+variables\s*=\s*"([^"]*)"', line,
                     re.IGNORECASE)
        if m:
            raw_vars = m.group(1)
            # Strip pointer sigils (*) before the variable name.
            var_names = [v.strip().lstrip('*') for v in raw_vars.split(',')]
            counts = {}
            for v in var_names:
                counts[v] = counts.get(v, 0) + 1
            for v, cnt in counts.items():
                if cnt > 1:
                    print(f"Duplicated assignment of '{v}' in `constructorAssign`"
                          f" directive in file '{file_name}'")
                    status = 1

    # ── Collect body lines for open function/subroutine ──────────────────────
    if in_function:
        func_body.append(line)
    if in_subroutine:
        sub_body.append(line)

    # ── End-of-block detection ───────────────────────────────────────────────

    # end type
    if _END_TYPE.match(line):
        in_type_block     = False
        current_type_name = None

    # end interface
    if _END_INTERFACE.match(line):
        in_interface_block = False
        current_iface_name = None

    # end function
    if _END_FUNCTION.match(line) and in_function:
        # Check 3: empty constructor.
        if (current_func_name in constructors
                and func_args_empty
                and not func_is_module
                and _function_is_empty(func_body)):
            print(f"Empty constructor function '{current_func_name}'"
                  f" in file '{file_name}'")
            status = 1
        in_function       = False
        current_func_name = None
        func_body         = []

    # end subroutine
    if _END_SUBROUTINE.match(line) and in_subroutine:
        # Check 4: empty finalizer.
        if (current_sub_name in destructors
                and _function_is_empty(sub_body)):
            print(f"Empty finalizer function '{current_sub_name}'"
                  f" in file '{file_name}'")
            status = 1
        in_subroutine    = False
        current_sub_name = None
        sub_body         = []

    # ── Block-opener detection ───────────────────────────────────────────────

    # type definition
    m = _TYPE_DEF.match(line)
    if m and not in_function and not in_subroutine:
        type_name = m.group(1).lower()
        known_types.append(type_name)
        in_type_block     = True
        current_type_name = type_name

    # interface block
    m = _INTERFACE.match(line)
    if m:
        iface_name = m.group(1).lower()
        if iface_name in known_types:
            in_interface_block = True
            current_iface_name = iface_name

    # function
    m = _FUNCTION.match(line)
    if m and not in_function and not in_subroutine:
        func_name      = m.group(1).lower()
        arg_list       = m.group(2)           # e.g. "()" or "(a, b)"
        in_function      = True
        current_func_name = func_name
        func_is_module   = bool(_MODULE_FUNCTION.match(line))
        func_args_empty  = re.match(r'^\(\s*\)$', arg_list.strip()) is not None
        func_body        = []

    # subroutine
    m = _SUBROUTINE.match(line)
    if m and not in_function and not in_subroutine:
        sub_name         = m.group(1).lower()
        in_subroutine     = True
        current_sub_name  = sub_name
        sub_body          = []

    # ── In-block content ─────────────────────────────────────────────────────

    # Check 1: pointer members not null-initialized (only inside type blocks).
    if in_type_block and not in_function and not in_subroutine:
        m = _POINTER_MEMBER.match(line)
        if m:
            var_list = m.group(1)
            # Split on commas that are not inside parentheses.
            depth     = 0
            variables = []
            current   = ''
            for ch in var_list:
                if ch == '(':
                    depth += 1
                elif ch == ')':
                    depth -= 1
                if ch == ',' and depth == 0:
                    variables.append(current.strip())
                    current = ''
                else:
                    current += ch
            if current.strip():
                variables.append(current.strip())

            for var in variables:
                if not _NULL_INIT.search(var):
                    # Extract the clean variable name (before any =>).
                    var_name = var.split('=>')[0].strip()
                    # Strip array dimensions if present.
                    var_name = re.sub(r'\(.*\)', '', var_name).strip()
                    # Infer type name from the type opener.
                    type_name_display = current_type_name or '(unknown)'
                    m2 = re.match(r'^\s*type\s*(?:[^:]*::)?\s*([a-zA-Z0-9_]+)',
                                  line, re.IGNORECASE)
                    if m2:
                        type_name_display = m2.group(1)
                    print(f"Pointer variable '{var_name}' in type"
                          f" '{type_name_display}' in file '{file_name}'"
                          f" is not null initialized")
                    status = 1

        # final :: subroutine_name
        m = _FINAL_DECL.match(line)
        if m:
            names = [n.strip().lower() for n in m.group(1).split(',')]
            destructors.extend(names)

    # module procedure inside an interface block for a known type.
    if in_interface_block:
        m = _MODULE_PROC.match(line)
        if m:
            names = [n.strip().lower() for n in m.group(1).split(',')]
            constructors.extend(names)

raise SystemExit(status)
