"""Provides utility functions and patterns for parsing Fortran source code.

Andrew Benson (ported to Python 2026)

Mirrors the subset of perl/Fortran/Utils.pm used by the Galacticus build system.
The get_fortran_line() function lives in Galacticus.Build.FortranUtils; import
it from there when needed.
"""

import re

__all__ = [
    'extract_variables', 'unformat_variables',
    'LABEL', 'ARGUMENT_LIST',
    'UNIT_OPENERS', 'UNIT_CLOSERS', 'INTRINSIC_DECLARATIONS',
]

# ---------------------------------------------------------------------------
# Basic regex fragments
# ---------------------------------------------------------------------------

LABEL         = r'[a-zA-Z0-9_{}¦]+'
ARGUMENT_LIST = r'[a-zA-Z0-9_{}¦,\s]*'

# ---------------------------------------------------------------------------
# Unit opener / closer patterns
# Each entry: {'unit_name': <capture-group index (0-based)>, 'regex': compiled re}
# Optional extra keys: 'intrinsic', 'kind', 'arguments'
# ---------------------------------------------------------------------------

UNIT_OPENERS = {
    # Module — avoid matching "module procedure / function / subroutine"
    'module': {
        'unit_name': 0,
        'regex': re.compile(
            r'^\s*module\s+(' + LABEL + r')\s*$',
            re.IGNORECASE,
        ),
    },
    # Submodule
    'submodule': {
        'unit_name': 0,
        'regex': re.compile(
            r'^\s*submodule\s*\(\s*[a-zA-Z0-9_:{}¦]+\s*\)\s+(' + LABEL + r')\s*$',
            re.IGNORECASE,
        ),
    },
    # Program
    'program': {
        'unit_name': 0,
        'regex': re.compile(
            r'^\s*program\s+(' + LABEL + r')',
            re.IGNORECASE,
        ),
    },
    # Subroutine (pure/elemental/recursive/impure/module prefix allowed).
    # Groups: 0=name, 1=arg-parens (optional), 2=arg-list (optional).
    'subroutine': {
        'unit_name': 0,
        'arguments': 2,
        'regex': re.compile(
            r'^\s*(?:impure\s+|pure\s+|elemental\s+|recursive\s+|module\s+)*'
            r'\s*subroutine\s+(' + LABEL + r')'
            r'\s*(\(\s*(' + ARGUMENT_LIST + r')\))?',
            re.IGNORECASE,
        ),
    },
    # Function (with optional type prefix)
    # Capture groups (0-based): 0=intrinsic, 1=kind, 2=function_name, 3=args_with_parens, 4=args_content
    'function': {
        'unit_name': 2,
        'intrinsic': 0,
        'kind': 1,
        'arguments': 4,
        'regex': re.compile(
            r'^\s*(?:impure\s+|pure\s+|elemental\s+|recursive\s+|module\s+)*'
            r'\s*(real|integer|double\s+precision|double\s+complex|character|logical)?'
            r'\s*(\((?:(?:kind|len)=)?[\w\d]*\))?'
            r'\s*function\s+(' + LABEL + r')'
            r'\s*(\(\s*(' + ARGUMENT_LIST + r')\))?',
            re.IGNORECASE,
        ),
    },
    # Module procedure (inside interface blocks)
    'moduleProcedure': {
        'unit_name': 0,
        'regex': re.compile(
            r'^\s*module\s+procedure\s+(' + LABEL + r')',
            re.IGNORECASE,
        ),
    },
    # Interface block
    'interface': {
        'unit_name': 0,
        'regex': re.compile(
            r'^\s*(?:abstract\s+)?interface\s*([a-zA-Z0-9_()/+\-*.=]*)',
            re.IGNORECASE,
        ),
    },
    # Derived type
    'type': {
        'unit_name': 0,
        'regex': re.compile(
            r'^\s*type\s*'
            r'(?:,\s*(?:abstract|public|private|extends\s*\(' + LABEL + r'\))\s*)*'
            r'(?:::)?\s*(' + LABEL + r')\s*$',
            re.IGNORECASE,
        ),
    },
    # `contains` marker — self-closing.  Perl's parser makes `contains` a
    # container whose children are the post-contains subprograms; here it
    # stays a sibling marker which keeps the parser logic simple, yet lets
    # `insert_pre_contains` / `insert_post_contains` anchor on it.
    'contains': {
        'unit_name': -1,
        'regex': re.compile(r'^\s*contains\s*$', re.IGNORECASE),
    },
}

# Note: The function regex group numbering needs verification against actual patterns.

UNIT_CLOSERS = {
    'module':          re.compile(r'^\s*end\s+module\s+('           + LABEL + r')' , re.IGNORECASE),
    'submodule':       re.compile(r'^\s*end\s+submodule\s+('        + LABEL + r')' , re.IGNORECASE),
    'program':         re.compile(r'^\s*end\s+program\s+('          + LABEL + r')' , re.IGNORECASE),
    'subroutine':      re.compile(r'^\s*end\s+subroutine\s+('       + LABEL + r')' , re.IGNORECASE),
    'function':        re.compile(r'^\s*end\s+function\s+('         + LABEL + r')' , re.IGNORECASE),
    'moduleProcedure': re.compile(r'^\s*end\s+procedure\s+('        + LABEL + r')' , re.IGNORECASE),
    'interface':       re.compile(r'^\s*end\s+interface\s*([a-zA-Z0-9_()/+\-*.=]*)', re.IGNORECASE),
    'type':            re.compile(r'^\s*end\s+type\s+('             + LABEL + r')' , re.IGNORECASE),
}

# ---------------------------------------------------------------------------
# Intrinsic declaration regex patterns
# Each entry: {'intrinsic': str, 'openmp': group-index, 'type': group-index,
#              'attributes': group-index, 'variables': group-index, 'regex': compiled re}
# ---------------------------------------------------------------------------

_KV  = r'[a-zA-Z0-9_=]+'          # kind/len value
_KVX = r'[a-zA-Z0-9_=,+\-*()\s]+' # extended kind value (character)
_ATR = r'[\sa-zA-Z0-9_,%:+\-*/()]*'  # attributes text
_VAR = r'[\sa-zA-Z0-9._,:=>+\-*/()[\]]*'  # variable list
_VAR2 = r'[\sa-zA-Z0-9_.:,=>+\-*/()[\]]*'  # variable list variant 2

INTRINSIC_DECLARATIONS = {
    'integer': {
        'intrinsic':   'integer',
        'openmp':      0,
        'type':        1,
        'attributes':  2,
        'variables':   3,
        'regex': re.compile(
            r'^\s*(!\$)?\s*integer'
            r'(?:\s*\(\s*(?:' + _KV + r'|kind\s*\(\s*[a-zA-Z0-9_]+\s*\))\s*\))?'
            r'([\sa-zA-Z0-9_,%:+\-*/()]*)?'
            r'::\s*([\sa-zA-Z0-9._,:=>+\-*/()[\]]+)\s*$',
            re.IGNORECASE,
        ),
    },
    'real': {
        'intrinsic':   'real',
        'openmp':      0,
        'type':        1,
        'attributes':  2,
        'variables':   3,
        'regex': re.compile(
            r'^\s*(!\$)?\s*real'
            r'(?:\s*\(\s*' + _KV + r'\s*\))?'
            r'([\sa-zA-Z0-9_,%:+\-*/()]*)?'
            r'::\s*([\sa-zA-Z0-9._,:=>+\-*/()[\]]+)\s*$',
            re.IGNORECASE,
        ),
    },
    'double': {
        'intrinsic':   'double precision',
        'openmp':      0,
        'type':        1,
        'attributes':  2,
        'variables':   3,
        'regex': re.compile(
            r'^\s*(!\$)?\s*double\s+precision'
            r'(?:\s*\(\s*' + _KV + r'\s*\))?'
            r'([\sa-zA-Z0-9_,%:=+\-*/()]*)?'
            r'::\s*([\sa-zA-Z0-9._,:=>+\-*/()[\]]+)\s*$',
            re.IGNORECASE,
        ),
    },
    'complex': {
        'intrinsic':   'complex',
        'openmp':      0,
        'type':        1,
        'attributes':  2,
        'variables':   3,
        'regex': re.compile(
            r'^\s*(!\$)?\s*complex'
            r'(?:\s*\(\s*' + _KV + r'\s*\))?'
            r'([\sa-zA-Z0-9_,%:+\-*/()]*)?'
            r'::\s*([\sa-zA-Z0-9._,:=>+\-*/()[\]]+)\s*$',
            re.IGNORECASE,
        ),
    },
    'doubleComplex': {
        'intrinsic':   'double complex',
        'openmp':      0,
        'type':        1,
        'attributes':  2,
        'variables':   3,
        'regex': re.compile(
            r'^\s*(!\$)?\s*double\s+complex'
            r'(?:\s*\(\s*' + _KV + r'\s*\))?'
            r'([\sa-zA-Z0-9_,%:=+\-*/()]*)?'
            r'::\s*([\sa-zA-Z0-9._,:=>+\-*/()[\]]+)\s*$',
            re.IGNORECASE,
        ),
    },
    'logical': {
        'intrinsic':   'logical',
        'openmp':      0,
        'type':        1,
        'attributes':  2,
        'variables':   3,
        'regex': re.compile(
            r'^\s*(!\$)?\s*logical'
            r'(?:\s*\(\s*' + _KV + r'\s*\))?'
            r'([\sa-zA-Z0-9_,%:+\-*/()]*)?'
            r'::\s*([\sa-zA-Z0-9_.:,=>+\-*/()[\]]+)\s*$',
            re.IGNORECASE,
        ),
    },
    'character': {
        'intrinsic':   'character',
        'openmp':      0,
        'type':        1,
        'attributes':  2,
        'variables':   3,
        'regex': re.compile(
            r'^\s*(!\$)?\s*character'
            r'(?:\s*\(\s*' + _KVX + r'\s*\))?'
            r'([\sa-zA-Z0-9_,%:+\-*/()]*)?'
            r'::\s*([\sa-zA-Z0-9._,:=>+\-*/()[\]]+)\s*$',
            re.IGNORECASE,
        ),
    },
    'type': {
        'intrinsic':   'type',
        'openmp':      0,
        'type':        1,
        'attributes':  2,
        'variables':   3,
        'regex': re.compile(
            r'^\s*(!\$)?\s*type'
            r'(?:\s*\(\s*' + LABEL + r'\s*\))?'
            r'([\sa-zA-Z0-9_,%:+\-*/()]*)?'
            r'::\s*([\sa-zA-Z0-9._,:=>+\-*/()[\]]+)\s*$',
            re.IGNORECASE,
        ),
    },
    'class': {
        'intrinsic':   'class',
        'openmp':      0,
        'type':        1,
        'attributes':  2,
        'variables':   3,
        'regex': re.compile(
            r'^\s*(!\$)?\s*class'
            r'(?:\s*\(\s*[a-zA-Z0-9_*]+\s*\))?'
            r'([\sa-zA-Z0-9_,%:+\-*/()]*)?'
            r'::\s*([\sa-zA-Z0-9._,:=>+\-*/()[\]]+)\s*$',
            re.IGNORECASE,
        ),
    },
    'procedure': {
        'intrinsic':   'procedure',
        'openmp':      0,
        'type':        1,
        'attributes':  2,
        'variables':   3,
        'regex': re.compile(
            r'^\s*(!\$)?\s*procedure'
            r'(?:\s*\([a-zA-Z0-9_\s]*\))?'
            r'([\sa-zA-Z0-9_,%:+\-*/()]*)?'
            r'::\s*([a-zA-Z0-9_,:=>+\-*/()<>.{}¦\s]+)\s*$',
            re.IGNORECASE,
        ),
    },
    'generic': {
        'intrinsic':   'generic',
        'openmp':      0,
        'type':        1,
        'attributes':  2,
        'variables':   3,
        'regex': re.compile(
            r'^\s*(!\$)?\s*generic'
            r'(?:\s*\([a-zA-Z0-9_\s]*\))?'
            r'([\sa-zA-Z0-9_,%:+\-*/()]*)?'
            r'::\s*([a-zA-Z0-9_,:=+\-*/()<>.{}¦\s]+)\s*$',
            re.IGNORECASE,
        ),
    },
    'final': {
        'intrinsic':   'final',
        'openmp':      0,
        'type':        1,
        'attributes':  2,
        'variables':   3,
        'regex': re.compile(
            r'^\s*(!\$)?\s*final'
            r'(?:\s*\([a-zA-Z0-9_\s]*\))?'
            r'([\sa-zA-Z0-9_,%:+\-*/()]*)?'
            r'::\s*([a-zA-Z0-9_,\s]+)\s*$',
            re.IGNORECASE,
        ),
    },
}

# ---------------------------------------------------------------------------
# extract_variables
# ---------------------------------------------------------------------------

def extract_variables(variable_list, keep_qualifiers=False, lower_case=True, remove_spaces=True):
    """Parse the post-'::' section of a Fortran declaration line into variable names.

    Mirrors Perl Fortran::Utils::Extract_Variables().

    Parameters
    ----------
    variable_list : str or None
        The comma-separated variable text (after '::').
    keep_qualifiers : bool
        If True, dimension/kind qualifiers like '(3)' are preserved.
    lower_case : bool
        If True (default), variable names are lowercased.
    remove_spaces : bool
        If True (default), whitespace is stripped.

    Returns
    -------
    list of str
        The individual variable names (possibly with qualifiers if keep_qualifiers).
    """
    if variable_list is None:
        return []

    if '::' in variable_list:
        raise ValueError(
            f"extract_variables: variable list '{variable_list}' contains '::' — "
            "most likely regex matching failed"
        )

    # Stash Fortran string literals before lowercasing or whitespace removal.
    # Strings inside `=...` initializers (filenames, descriptions, …) MUST
    # keep their original case — Fortran is case-insensitive for identifiers
    # but the bytes inside a string literal are not, and a lowercased
    # filename like `'cooling/cooling_function_atomic_cie_cloudy.hdf5'`
    # won't open on a case-sensitive filesystem when the original was
    # `'cooling/cooling_function_Atomic_CIE_Cloudy.hdf5'`.  Strings can
    # also contain commas, parens, and asterisks that the downstream
    # bracket-walk and comma-split would otherwise mishandle — stashing
    # them as opaque tokens fixes that too.  We support both single- and
    # double-quoted forms and Fortran's doubled-quote escape (`''` inside
    # `'…'`).
    string_literals = []

    def _stash(m):
        idx = len(string_literals)
        string_literals.append(m.group(0))
        # Lowercase placeholder so the subsequent `.lower()` is a no-op on
        # it — otherwise the placeholder mutates and the restore regex
        # below can't find it.
        return f'%%strlit{idx}%%'

    variable_list = re.sub(
        r"'(?:[^']|'')*'|\"(?:[^\"]|\"\")*\"", _stash, variable_list)

    if lower_case:
        variable_list = variable_list.lower()

    if remove_spaces:
        variable_list = re.sub(r'\s', '', variable_list)
    else:
        variable_list = variable_list.rstrip()

    # Handle '*' (e.g. character*80)
    if remove_spaces:
        if not keep_qualifiers:
            variable_list = variable_list.replace('*', '')
        else:
            variable_list = variable_list.replace('*', '%%ASTERISK%%')

    # Remove bracketed content — round and square brackets.
    # We iterate because brackets can be nested.
    MAX_ITER = 10000
    for iteration in range(MAX_ITER):
        m = re.search(r'[\(\[]', variable_list)
        if not m:
            break
        # Find the matching closing bracket.
        start = m.start()
        open_char  = variable_list[start]
        close_char = ')' if open_char == '(' else ']'
        depth = 0
        end   = start
        for i in range(start, len(variable_list)):
            if variable_list[i] == open_char:
                depth += 1
            elif variable_list[i] == close_char:
                depth -= 1
                if depth == 0:
                    end = i
                    break
        extracted = variable_list[start:end + 1]
        prefix    = variable_list[:start]
        remainder = variable_list[end + 1:]
        if not keep_qualifiers:
            variable_list = prefix + remainder
        else:
            # Escape bracket/comma content so splitting still works.
            escaped = extracted
            escaped = escaped.replace('(', '%%OPEN%%')
            escaped = escaped.replace(')', '%%CLOSE%%')
            escaped = escaped.replace('[', '%%OPENSQ%%')
            escaped = escaped.replace(']', '%%CLOSESQ%%')
            escaped = escaped.replace(',', '%%COMMA%%')
            variable_list = prefix + escaped + remainder
    else:
        import sys
        print(f"extract_variables: maximum iterations exceeded for input: '{variable_list}'",
              file=sys.stderr)
        return []

    # Strip initialisation/association suffixes (e.g. "x=1.0" → "x") unless keeping qualifiers.
    if not keep_qualifiers:
        variable_list = re.sub(r'=[^,]*(,|$)', r'\1', variable_list)

    # Split on commas.
    variables = [v.strip() for v in variable_list.split(',') if v.strip()]

    # Unescape if qualifiers were kept.
    if keep_qualifiers:
        unescaped = []
        for v in variables:
            v = v.replace('%%OPEN%%',    '(')
            v = v.replace('%%CLOSE%%',   ')')
            v = v.replace('%%OPENSQ%%',  '[')
            v = v.replace('%%CLOSESQ%%', ']')
            v = v.replace('%%COMMA%%',   ',')
            v = v.replace('%%ASTERISK%%','*')
            unescaped.append(v)
        variables = unescaped

    # Restore stashed string literals (preserving original case).
    if string_literals:
        def _restore(token):
            return re.sub(
                r'%%strlit(\d+)%%',
                lambda m: string_literals[int(m.group(1))],
                token)
        variables = [_restore(v) for v in variables]

    return variables


def unformat_variables(variable_string):
    """Parse a Fortran variable declaration string into a structured dict.

    Mirrors Perl `Fortran::Utils::Unformat_Variables`.  Returns None if
    the input doesn't parse as a valid declaration.

    The returned dict carries `intrinsic`, `variables` (with qualifiers
    preserved), `variableNames` (without qualifiers), and optionally
    `type` and `attributes`.

    The Perl original walks the `INTRINSIC_DECLARATIONS` regex table
    and indexes capture groups by position, but the Python ports of
    those regexes use slightly different capture-group structures
    (non-capturing type brackets) than the Perl versions.  Rather than
    keep two versions of the regex table in sync we just split the
    string at the `::` separator, find the matching intrinsic prefix,
    extract any `(...)` type bracket, and treat whatever remains as
    attributes.  Equivalent semantics, simpler implementation.
    """
    if variable_string is None or '::' not in variable_string:
        return None
    before_decl, variables_part = variable_string.split('::', 1)

    # Detect optional !$ OpenMP marker.
    m = re.match(r'^\s*(!\$)?\s*(.*)$', before_decl, re.DOTALL)
    if not m:
        return None
    head = m.group(2).rstrip()

    # Match against each intrinsic in turn — the first that matches
    # wins.  Use the same intrinsic names as the Perl table; multi-word
    # intrinsics (`double precision`, `double complex`) need a flexible
    # whitespace pattern.
    intrinsic_patterns = [
        ('double precision', r'double\s+precision'),
        ('double complex',   r'double\s+complex'),
        ('integer',          r'integer'),
        ('real',             r'real'),
        ('complex',          r'complex'),
        ('logical',          r'logical'),
        ('character',        r'character'),
        ('type',             r'type'),
        ('class',            r'class'),
        ('procedure',        r'procedure'),
        ('generic',          r'generic'),
        ('final',            r'final'),
    ]

    matched_intrinsic = None
    rest = None
    for name, pattern in intrinsic_patterns:
        m = re.match(r'^(' + pattern + r')\s*(.*)$', head, re.IGNORECASE)
        if m:
            matched_intrinsic = name
            rest = m.group(2).strip()
            break
    if matched_intrinsic is None:
        return None

    # Pull off the optional `(...)` type bracket if present.  Uses
    # `extract_bracketed` to walk balanced parens (so `kind(c_int)`-style
    # nesting doesn't trip us up).  Local import avoids a top-level
    # circular dep with `build.fortran_utils`.
    type_text = None
    if rest.startswith('('):
        from Galacticus.Build.FortranUtils import extract_bracketed
        bracket, after, _ = extract_bracketed(rest, "()")
        if bracket is not None:
            type_text = bracket
            rest = after.strip()

    # Whatever remains before the leading comma (if any) is leftover
    # after-type whitespace; whatever is after a comma is attributes.
    attributes_string = re.sub(r'^\s*,\s*', '', rest) if rest else ''

    if type_text is not None:
        # Strip the outer parens then any internal whitespace.
        type_text = re.sub(r'^\((.*)\)$', r'\1', type_text)
        type_text = re.sub(r'\s', '', type_text)

    variables_string  = variables_part
    variables         = extract_variables(variables_string,  keep_qualifiers=True,  remove_spaces=True)
    variable_names    = extract_variables(variables_string,  keep_qualifiers=False, remove_spaces=True)
    attributes        = extract_variables(attributes_string, keep_qualifiers=True,  remove_spaces=True) \
        if attributes_string else []

    result = {
        'intrinsic':     matched_intrinsic,
        'variables':     variables,
        'variableNames': variable_names,
    }
    if type_text:
        result['type'] = type_text
    if attributes:
        result['attributes'] = attributes
    return result
