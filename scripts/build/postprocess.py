#!/usr/bin/env python3
import bisect
import re
import sys

# Remap line numbers in preprocessed files to their original sources.
# Reads compiler output from stdin and writes remapped output to stdout.
# Andrew Benson (ported to Python 2026)

if len(sys.argv) != 2:
    print("Usage: postprocess.py <preprocessedSource>", file=sys.stderr)
    sys.exit(1)

preprocessed_source = sys.argv[1]

# Fast path: a clean compile produces no diagnostics at all. Read the first
# input line before anything else — if there is none, skip the (relatively
# expensive) pre-scan of the preprocessed source and the .lmap parse
# entirely. Every compile pays this script's cost, and most compiles are
# clean.
_first_line = sys.stdin.readline()
if not _first_line:
    sys.exit(0)

# Detect whether stdout is a terminal (for colour output).
_have_color = sys.stdout.isatty()
if _have_color:
    # ANSI colour codes.
    _BRIGHT_MAGENTA_BOLD = '\033[1;35m'
    _BRIGHT_GREEN_BOLD   = '\033[1;32m'
    _BOLD                = '\033[1m'
    _RESET               = '\033[0m'

# --- Parse the line-number map (.lmap file) ---
_map = [
    {'source': preprocessed_source, 'line': 1, 'lineOriginal': 1}
]
lmap_path = preprocessed_source + '.lmap'
try:
    with open(lmap_path, 'r') as fh:
        for line in fh:
            m = re.match(
                r'^!-->\s+(\d+)\s+(\d+)\s+"([a-zA-Z0-9_\-.\/()\:]+)"\s*$',
                line
            )
            if m:
                _map.append({
                    'source'      : m.group(3),
                    'line'        : int(m.group(1)),
                    'lineOriginal': int(m.group(2)),
                })
except OSError:
    pass  # If there's no .lmap file, the map stays minimal.

_map_keys = [entry['lineOriginal'] for entry in _map]

# --- Pre-scan the preprocessed source for compiler directives ---
unused_functions       = {}   # function name → True
initialized_variables  = {}   # variable name → True
ignore_outlives        = {}   # pointer name  → True
ignore_unused          = {}   # module var    → True
unused_variables       = {}   # {unit: {var: True}}
interoperable_variables= {}   # {unit: {var: True}}

unit_name = None
try:
    with open(preprocessed_source, 'r', errors='replace') as fh:
        # Regexes from Fortran::Utils for unit openers.
        # Unit-opener regexes — used only to extract the unit NAME so we can
        # later associate `!$GLC attributes unused :: …` directives with the
        # correct enclosing routine.  We deliberately do NOT match the
        # argument list: the original Perl regex used a nested
        # `(?:\(\s*(?:[a-zA-Z0-9_{}¦,\s]*)\))*` for the argument list, and
        # the equivalent Python pattern has Python's NFA engine backtracking
        # exponentially when a long generated line of identifier+whitespace
        # characters never produces a `function`/`subroutine` match.  The
        # argument-list capture isn't read anywhere downstream, so dropping
        # it is purely a simplification, not a behaviour change.  Each
        # alternation has its trailing `\s+` folded inside the group so
        # `(?:opt\s+)*` cannot fight with a separate `\s*` for the same
        # whitespace, and the type-spec is `?` (one or zero) rather than
        # `*` (a Fortran function has at most one return-type qualifier).
        _label = r'[a-zA-Z0-9_{}¦]+'
        _unit_openers = {
            'subroutine': re.compile(
                r'^\s*(?:(?:impure|pure|elemental|recursive|module)\s+)*'
                r'subroutine\s+(' + _label + r')',
                re.IGNORECASE),
            'function': re.compile(
                r'^\s*(?:(?:impure|pure|elemental|recursive|module)\s+)*'
                r'(?:(?:real|integer|double\s+precision|double\s+complex|'
                r'character|logical)\s*(?:\([^)]*\))?\s*)?'
                r'function\s+(' + _label + r')',
                re.IGNORECASE),
            'moduleProcedure': re.compile(
                r'^\s*module\s+procedure\s+(' + _label + r')',
                re.IGNORECASE),
        }
        for line in fh:
            for regex in _unit_openers.values():
                matches = regex.match(line)
                if matches:
                    unit_name = matches.group(1).lower()
                    break   # don't run the remaining regexes — they can't match too

            m = re.match(
                r'^\s*!\$GLC\s+function\s+attributes\s+unused\s*::\s*([a-zA-Z0-9_,\s]+)\s*$',
                line)
            if m:
                for fn in re.split(r'\s*,\s*', m.group(1).rstrip()):
                    unused_functions[fn.lower()] = True

            m = re.match(
                r'^\s*!\$GLC\s+attributes\s+initialized\s*::\s*([a-zA-Z0-9_,\s]+)\s*$',
                line)
            if m:
                for v in re.split(r'\s*,\s*', m.group(1).rstrip()):
                    initialized_variables[v.lower()] = True

            m = re.match(
                r'^\s*!\$GLC\s+ignore\s+outlive\s*::\s*([a-zA-Z0-9_,\s]+)\s*$',
                line)
            if m:
                for v in re.split(r'\s*,\s*', m.group(1).rstrip()):
                    ignore_outlives[v.lower()] = True

            m = re.match(
                r'^\s*!\$GLC\s+ignore\s+unused\s*::\s*([a-zA-Z0-9_,\s]+)\s*$',
                line)
            if m:
                for v in re.split(r'\s*,\s*', m.group(1).rstrip()):
                    ignore_unused[v.lower()] = True

            m = re.match(
                r'^\s*!\$GLC\s+attributes\s+unused\s*::\s*([a-zA-Z0-9_,\s]+)\s*$',
                line)
            if m:
                u = unused_variables.setdefault(unit_name, {})
                for v in re.split(r'\s*,\s*', m.group(1).rstrip()):
                    u[v.lower()] = True

            m = re.match(
                r'^\s*!\$GLC\s+attributes\s+interoperable\s*::\s*([a-zA-Z0-9_,\s]+)\s*$',
                line)
            if m:
                u = interoperable_variables.setdefault(
                    unit_name.lower() if unit_name else '', {}
                )
                for v in re.split(r'\s*,\s*', m.group(1).rstrip()):
                    u[v.lower()] = True
except OSError:
    pass

# --- Process stdin ---
status       = 0
buffer       = None
last_dropped = False
function_name = None
procedure_name = None
pointer_name   = None
bogus_uninit   = {}

# Unit-opener regexes for matching inside compiler output lines.  Same
# simplification as the source-file pre-scan above — see the comment there
# for why nested unbounded quantifiers were dropped.
_label_p = r'[a-zA-Z0-9_{}¦]+'
_opener_re = {
    'subroutine': re.compile(
        r'^\s*(?:(?:impure|pure|elemental|recursive|module)\s+)*'
        r'subroutine\s+(' + _label_p + r')',
        re.IGNORECASE),
    'function': re.compile(
        r'^\s*(?:(?:impure|pure|elemental|recursive|module)\s+)*'
        r'(?:(?:real|integer|double\s+precision|double\s+complex|'
        r'character|logical)\s*(?:\([^)]*\))?\s*)?'
        r'function\s+(' + _label_p + r')',
        re.IGNORECASE),
}

import itertools

for line in itertools.chain([_first_line], sys.stdin):
    # --- Remap preprocessed file:line references ---
    m = re.match(r'^([a-zA-Z0-9_./]+\.p\.F90):(\d+):([\d\-]+):\s*$', line)
    if m:
        line_original = int(m.group(2))
        flag          = m.group(3)
        # The .lmap entries are in ascending lineOriginal order; find the
        # last entry at or before this line.
        idx = bisect.bisect_right(_map_keys, line_original) - 1
        source_entry = _map[idx]
        src = source_entry['source']
        if src.endswith('()'):
            line_descriptor = "auto-generated code (no line number)"
        else:
            ln = line_original - source_entry['lineOriginal'] + source_entry['line']
            line_descriptor = f"line {ln}"
        line = f"{src}; {line_descriptor} [preprocessed line {line_original}]; code {flag}\n"
        if buffer is not None:
            sys.stdout.write(buffer)
        buffer = None

    drop_buffer = False

    # gfortran PR#58175
    if re.search(r'Only array FINAL procedures declared for derived type', line):
        drop_buffer = True

    # gfortran PR#86117 — struct MEM[] uninitialized
    m = re.search(
        r"Warning:\s+['\u2018\u2019]MEM\[\(struct\s+([a-zA-Z0-9_]+)[a-z0-9_\s*]+\)[a-z0-9_&]+\s+\+\s+\d+B\]['\u2018\u2019] (?:is|may be) used uninitialized",
        line)
    if m:
        drop_buffer = True
        bogus_uninit[m.group(1)] = True

    m = re.search(
        r"Warning:\s+['\u2018\u2019](MEM\s<real[a-zA-Z0-9_.<>()\[\]=:\s*]+>\s[a-zA-Z0-9_.<>()\[\]=:\s*]+)['\u2018\u2019] (?:is|may be) used uninitialized",
        line)
    if m:
        drop_buffer = True
        bogus_uninit[m.group(1)] = True

    m = re.search(
        r"Warning:\s+['\u2018\u2019]([a-zA-Z0-9_.<>()\[\]=:\s*]+)['\u2018\u2019]\s+may be used uninitialized",
        line)
    if m:
        elements = m.group(1).split('.')
        for sym in bogus_uninit:
            if sym in elements:
                drop_buffer = True
                break

    # Versions of gfortran ≥ 12 (which include a complete implementation of finalization) emit spurious warnings related to
    # internal finalizers - ignore them.
    m = re.search(
        r"^_F\.DA\d+",
        line)
    if m:
        drop_buffer = True

    # Unused function attributes.
    m = re.match(r'^\s*\d+\s*\|\s*subroutine\s+([a-z0-9_]+)', line, re.IGNORECASE)
    if m:
        function_name = m.group(1).lower()
    if re.search(r'\[-Wunused-function\]', line) and function_name is not None:
        if function_name in unused_functions:
            drop_buffer = True
        function_name = None

    # Unused variable attributes.
    m = re.match(r'^\s*\d+\s*\|\s*(.+)', line)
    if m:
        opener_text = m.group(1)
        for regex in _opener_re.values():
            mm = regex.match(opener_text)
            if mm:
                procedure_name = mm.group(1).lower()
                break

    if re.match(r'^<stdin>:\d+:\d+:', line):
        procedure_name = 'stdin'

    m = re.match(
        r"^Warning: Unused dummy argument ['\u2018\u2019]([a-zA-Z0-9_]+)['\u2018\u2019] at \(\d+\) \[-Wunused-dummy-argument\]",
        line)
    if m and procedure_name is not None:
        var = m.group(1)
        if procedure_name == 'stdin' or \
           var in unused_variables.get(procedure_name, {}):
            drop_buffer = True

    m = re.match(
        r"^Warning: Dummy argument ['\u2018\u2019]([a-zA-Z0-9_]+)['\u2018\u2019] at \(\d+\) was declared INTENT\(OUT\) but was not set \[-Wunused-dummy-argument\]",
        line)
    if m and procedure_name is not None:
        var = m.group(1)
        if procedure_name == 'stdin' or \
           var in unused_variables.get(procedure_name, {}):
            drop_buffer = True

    m = re.match(
        r"^Warning: Derived-type dummy argument ['\u2018\u2019]([a-zA-Z0-9_]+)['\u2018\u2019] at \(\d+\) was declared INTENT\(OUT\) but was not set and does not have a default initializer \[-Wunused-dummy-argument\]",
        line)
    if m and procedure_name is not None:
        var = m.group(1)
        if procedure_name == 'stdin' or \
           var in unused_variables.get(procedure_name, {}):
            drop_buffer = True

    # C-interoperable dummy argument warnings.
    m = re.match(
        r"^Warning: Variable ['\u2018\u2019]([a-zA-Z0-9_]+)['\u2018\u2019] at \(\d+\) is a dummy argument of the BIND\(C\) procedure ['\u2018\u2019]([a-z_]+)['\u2018\u2019] but may not be C interoperable \[-Wc-binding-type\]",
        line)
    if m and procedure_name is not None:
        var_name  = m.group(1)
        proc_name = m.group(2)
        if proc_name.lower() in interoperable_variables and \
           var_name.lower() in interoperable_variables[proc_name.lower()]:
            drop_buffer = True

    # Uninitialized variable attributes.
    m = re.search(
        r"Warning: ['\u2018\u2019](?:\(\*\))?([a-zA-Z0-9_]+)[a-zA-Z0-9_.[\]\s{}):]*['\u2018\u2019] may be used uninitialized(?:\s+in this function)? \[-Wmaybe-uninitialized\]",
        line)
    if m and m.group(1).lower() in initialized_variables:
        drop_buffer = True

    m = re.search(
        r"Warning: ['\u2018\u2019]([a-zA-Z0-9_]+)[a-zA-Z0-9_.[\]]*['\u2018\u2019] is used uninitialized(?:\s+in this function)? \[-Wuninitialized\]",
        line)
    if m and m.group(1).lower() in initialized_variables:
        drop_buffer = True

    m = re.search(
        r"note: ['\u2018\u2019]([a-zA-Z0-9_]+)[a-zA-Z0-9_.[\]]*['\u2018\u2019](?: was)? declared here",
        line)
    if m and m.group(1).lower() in initialized_variables:
        drop_buffer = True

    # "Unused PRIVATE module variable" warnings.
    m = re.match(
        r"^\s*Warning: Unused PRIVATE module variable ['\u2018\u2019]([a-zA-Z0-9_]+)['\u2018\u2019] declared at \(1\) \[-Wunused-value\]",
        line)
    if m and m.group(1).lower() in ignore_unused:
        drop_buffer = True

    # "pointer may outlive target" warnings.
    m = re.match(r'^\s*\d+\s*\|\s*([a-z0-9_]+)\s*=>\s*[a-z0-9_]+', line, re.IGNORECASE)
    if m:
        pointer_name = m.group(1).lower()
    if re.search(r'\[-Wtarget-lifetime\]', line):
        if pointer_name is not None and pointer_name in ignore_outlives:
            drop_buffer = True
        pointer_name = None

    # "note:" following a dropped warning.
    if re.match(r'^note:', line) and last_dropped:
        drop_buffer = True

    # Apply colour highlighting if interactive.
    if _have_color:
        if re.match(r'^Warning:\s', line):
            line = re.sub(r'^Warning:\s',
                          _BRIGHT_MAGENTA_BOLD + 'Warning: ' + _RESET, line)
            line = re.sub(
                r"(['\u2018\u2019])([^'\u2018\u2019]+)(['\u2018\u2019])",
                lambda mo: mo.group(1) + _BOLD + mo.group(2) + _RESET + mo.group(3) + "'",
                line)
        if re.match(r'^\s*\^\s*$', line):
            line = re.sub(r'\^', _BRIGHT_GREEN_BOLD + '^' + _RESET, line)
        if re.match(r'^\s*1\s*$', line):
            line = re.sub(r'1', _BRIGHT_GREEN_BOLD + '1' + _RESET, line)

    # Accumulate into buffer.
    if buffer is None:
        buffer = line
    else:
        buffer += line

    if drop_buffer:
        buffer      = None
        last_dropped = True
    elif re.match(r'^(Error|Warning):', line):
        if re.match(r'^Error:', line):
            status = 1
        sys.stdout.write(buffer)
        buffer      = None
        last_dropped = False

if buffer is not None:
    sys.stdout.write(buffer)

sys.exit(status)
