"""Type inference for Galacticus input parameters.

Andrew Benson (2026)

Only ~1.5% of `<inputParameter>` directives carry an explicit `<type>`, so the
type of nearly every parameter must be inferred.  Two mechanical sources are
used, in order of decreasing confidence:

  1. the explicit ``<type>`` element, when present;
  2. the Fortran literal in ``<defaultValue>`` (covers ~2/3 of parameters);
  3. the declared type of the variable the parameter is read into, looked up in
     the enclosing constructor's declarations.

Each parameter records the canonical type plus the *provenance* of that type
(``explicit`` / ``default`` / ``declaration`` / ``unknown``) so consumers can
weight the confidence and so we can measure coverage.

Canonical types mirror the `<inputParameter><type>` schema enumeration --
``boolean``, ``integer``, ``real``, ``string`` -- with an extra ``object`` type
for the (rare) case of a parameter whose variable is a derived type, and
``unknown`` when no source resolves it.  A finer-grained Fortran ``kind`` (e.g.
``double``, ``long``) is recorded separately when detectable, without changing
the canonical type.
"""

import re

# Canonical parameter types.
TYPE_BOOLEAN = 'boolean'
TYPE_INTEGER = 'integer'
TYPE_REAL    = 'real'
TYPE_STRING  = 'string'
TYPE_OBJECT  = 'object'
TYPE_UNKNOWN = 'unknown'

# Provenance tags.
PROV_EXPLICIT    = 'explicit'
PROV_DEFAULT     = 'default'
PROV_DECLARATION = 'declaration'
PROV_UNKNOWN     = 'unknown'


def _scalar(value):
    """Reduce an `xml_to_dict` value to a stripped string, or None.

    `xml_to_dict` may yield a str, a dict (element with attributes/children),
    or a list (repeated elements).  For the simple text-bearing elements we
    care about (`type`, `defaultValue`, `cardinality`, ...) we want the text.
    """
    if value is None:
        return None
    if isinstance(value, dict):
        # An element rendered as a dict -- look for its text content.
        text = value.get('content') if isinstance(value.get('content'), str) else None
        return text.strip() if text and text.strip() else None
    if isinstance(value, list):
        return None
    text = str(value).strip()
    return text or None


def normalize_explicit_type(type_text):
    """Map an explicit `<type>` value to a canonical type, or None."""
    if not type_text:
        return None
    t = type_text.strip().lower()
    return {
        'boolean': TYPE_BOOLEAN,
        'logical': TYPE_BOOLEAN,
        'integer': TYPE_INTEGER,
        'real':    TYPE_REAL,
        'double':  TYPE_REAL,
        'string':  TYPE_STRING,
    }.get(t)


# Fortran literal patterns.
_RE_LOGICAL    = re.compile(r'^\.(true|false)\.$', re.IGNORECASE)
_RE_VARSTR     = re.compile(r'^var_str\s*\(', re.IGNORECASE)
_RE_QUOTED     = re.compile(r"""^['"]""")
_RE_INTEGER    = re.compile(r'^[+-]?\d+(_[a-zA-Z0-9_]+)?$')
# A real/double literal: requires a decimal point or an exponent marker.
_RE_REAL       = re.compile(
    r'^[+-]?(\d+\.\d*|\.\d+|\d+)([deDE][+-]?\d+)?(_[a-zA-Z0-9_]+)?$')
_RE_HAS_POINT_OR_EXP = re.compile(r'[.dDeE]')
# Long-integer kinds.  In a literal the kind follows an underscore (`0_c_long`),
# anchored at the end; in a declaration type-spec it appears bare (`integer(c_long)`
# -> type-spec `c_long`, or `kind=c_long`).
_RE_LONG_KIND     = re.compile(r'_(c_long|kind_int8|int64|c_size_t)$', re.IGNORECASE)
_RE_LONG_TYPESPEC = re.compile(r'(c_long|kind_int8|int64|c_size_t)', re.IGNORECASE)
_RE_DOUBLE_EXP = re.compile(r'[dD][+-]?\d+$')


def classify_literal(literal):
    """Classify a Fortran literal default value.

    Returns ``(canonical_type, kind)`` where ``kind`` is a finer Fortran kind
    descriptor (``'double'``, ``'long'``) or None, or ``(None, None)`` when the
    literal is not recognised (e.g. an expression referencing other entities).
    """
    if literal is None:
        return None, None
    text = literal.strip()
    if not text:
        return None, None

    if _RE_LOGICAL.match(text):
        return TYPE_BOOLEAN, None
    if _RE_VARSTR.match(text) or _RE_QUOTED.match(text):
        return TYPE_STRING, None
    if _RE_INTEGER.match(text):
        kind = 'long' if _RE_LONG_KIND.search(text) else None
        return TYPE_INTEGER, kind
    if _RE_REAL.match(text):
        # Distinguish double (has a `d` exponent) from single-precision real.
        kind = 'double' if _RE_DOUBLE_EXP.search(text) else None
        # A bare integer would have matched _RE_INTEGER above, so anything here
        # with no point/exponent is unusual -- treat as real defensively.
        return TYPE_REAL, kind
    return None, None


def canonical_from_declaration(declaration):
    """Map a parsed declaration dict (from Parse/Declarations) to a canonical
    type and Fortran kind.

    Returns ``(canonical_type, kind)`` or ``(None, None)`` when the declaration
    does not correspond to a scalar input-parameter type.
    """
    if not declaration:
        return None, None
    intrinsic = (declaration.get('intrinsic') or '').lower()
    type_spec = (declaration.get('type') or '').lower()

    if intrinsic == 'logical':
        return TYPE_BOOLEAN, None
    if intrinsic == 'integer':
        kind = 'long' if _RE_LONG_TYPESPEC.search(type_spec) else None
        return TYPE_INTEGER, kind
    if intrinsic == 'double precision':
        return TYPE_REAL, 'double'
    if intrinsic == 'real':
        # `real(kind=...)` with a double kind is still effectively double.
        kind = 'double' if ('c_double' in type_spec or 'kind_double' in type_spec) else None
        return TYPE_REAL, kind
    if intrinsic == 'character':
        return TYPE_STRING, None
    if intrinsic == 'type':
        # A varying_string is our string type; anything else is a derived type.
        if 'varying_string' in type_spec:
            return TYPE_STRING, None
        return TYPE_OBJECT, type_spec or None
    if intrinsic == 'class':
        return TYPE_OBJECT, type_spec or None
    return None, None


def infer_parameter_type(directive, declaration_lookup=None):
    """Infer the canonical type of a single `<inputParameter>` directive.

    Parameters
    ----------
    directive : dict
        The parsed directive dict (``node['directive']``) with possible keys
        ``name``, ``variable``, ``type``, ``defaultValue``, ``cardinality``.
    declaration_lookup : callable, optional
        ``f(variable_name) -> declaration_dict`` returning the parsed Fortran
        declaration for a variable in the enclosing constructor, or None.

    Returns
    -------
    dict
        ``{'type', 'kind', 'provenance'}``.
    """
    # 1. Explicit <type>.
    explicit = normalize_explicit_type(_scalar(directive.get('type')))
    if explicit:
        return {'type': explicit, 'kind': None, 'provenance': PROV_EXPLICIT}

    # 2. Default-value literal.
    default = _scalar(directive.get('defaultValue'))
    if default is not None:
        ctype, kind = classify_literal(default)
        if ctype:
            return {'type': ctype, 'kind': kind, 'provenance': PROV_DEFAULT}

    # 3. Declaration of the variable the parameter is read into.
    variable = _scalar(directive.get('variable')) or _scalar(directive.get('name'))
    if variable and '(' not in variable and declaration_lookup is not None:
        declaration = declaration_lookup(variable)
        ctype, kind = canonical_from_declaration(declaration)
        if ctype:
            return {'type': ctype, 'kind': kind, 'provenance': PROV_DECLARATION}

    return {'type': TYPE_UNKNOWN, 'kind': None, 'provenance': PROV_UNKNOWN}
