"""Renders Fortran declarations in the reStructuredText notation used to
document the arguments and return types of class methods.

Andrew Benson (with assistance from Claude 2026)

An argument is written as an inline literal holding its type spec and name,
followed by a ``(optional)`` marker for optional arguments and an ``[in]``/
``[out]``/``[inout]`` direction marker::

    ``type(treeNode) node`` [inout]
    ``double precision(3) position`` (optional) [in]

Keeping the name inside the literal makes it verbatim, so neither a trailing
underscore (which reStructuredText would read as a reference) nor the
documentation's spell-checker has to be worked around.

The direction markers replace the ``\\argin``/``\\argout`` LaTeX macros of the
retired manual; optionality is spelled out rather than bracketed so that it
does not read as a third direction marker.

Two producers share this notation, which is why it lives here rather than in
either of them:

* :mod:`Galacticus.Build.SourceTree.Process.FunctionClass`, which embeds the
  notation in the ``<methods>`` metadata it generates for each
  ``functionClass``; and
* ``scripts/doc/extractDocsRST.py``, which renders the same notation onto the
  class pages — both for the declarations carried by hand-written
  ``<argument>`` elements and for the type-bound methods whose signatures it
  recovers from the source.
"""

import re

from Galacticus.Build.SourceTree.Parse.Declarations import parse_declaration


__all__ = [
    'arguments_to_rst', 'declaration_to_rst', 'format_argument',
    'is_declaration', 'normalize_argument', 'procedure_signature',
    'return_type_to_rst', 'split_arguments', 'type_spec',
]

# Direction markers, keyed by the `intent(…)` attribute that implies them.
_DIRECTIONS = {
    'intent(in)':    '[in]',
    'intent(out)':   '[out]',
    'intent(inout)': '[inout]',
}

# Return types that carry no information and so are rendered as no return type
# at all: a subroutine has nothing to return, and the component generator spells
# that `void`.
_VOID_TYPES = {'', 'void', '``void``'}


def type_spec(intrinsic, type_=None, rank=None):
    """Return the plain-text type spec for a declaration's parts.

    `type_` is the kind/type spec with its parentheses already stripped (as
    `parse_declaration` returns it), so they are restored here — without this
    a `type(treeNode)` argument renders as the meaningless `typetreeNode`.
    `rank` is the array shape, parentheses included (e.g. `(3)`, `(:,:)`).
    """
    spec = intrinsic
    if type_ is not None and type_ != '':
        spec += f'({type_})'
    if rank:
        spec += rank
    return spec


def normalize_argument(argument):
    """Return an argument already written in the notation, with its name inside
    the inline literal.

    Metadata written by hand keeps the two apart (```` ``double precision``
    timeStep [in] ````), which is easier to author; rendering wants them
    together, so that the name is verbatim like the type spec.  Arguments whose
    name is already inside the literal — or that have no name at all — are
    returned unchanged, so this may be applied freely.
    """
    return re.sub(r'^``([^`]+)`` ([A-Za-z_][A-Za-z0-9_]*)', r'``\1 \2``',
                  argument.strip())


def _rank_of(variable, attributes):
    """Return the array shape of a declared variable, parentheses included.

    The shape may be attached to the variable itself (``:: position(3)``) or
    declared once for the whole statement (``, dimension(:) ::``); either way
    it belongs on the type spec, not on the argument name.
    """
    m = re.search(r'(\(.*\))\s*$', variable)
    if m:
        return m.group(1)
    for attribute in attributes:
        m = re.match(r'dimension\s*(\(.*\))\s*$', attribute, re.IGNORECASE)
        if m:
            return m.group(1)
    return None


def format_argument(declaration, variable, variable_name):
    """Render one variable of a parsed `declaration` in the argument notation.

    `variable` is the declared form (which may carry an array shape) and
    `variable_name` the name alone, matching the `variables` and
    `variableNames` entries of `parse_declaration`.
    """
    attributes = declaration.get('attributes') or []
    spec       = type_spec(declaration['intrinsic'], declaration.get('type'),
                           _rank_of(variable, attributes))
    # Type spec and name share one inline literal, which is verbatim: nothing in
    # either needs escaping, and the documentation's spell-checker skips the
    # whole of it rather than reporting every argument name as a misspelling.
    rendered   = f'``{spec} {variable_name}``'
    if any(attribute.lower() == 'optional' for attribute in attributes):
        rendered += ' (optional)'
    for attribute in attributes:
        direction = _DIRECTIONS.get(attribute.lower())
        if direction:
            rendered += ' ' + direction
    return rendered


def declaration_to_rst(line):
    """Render every variable declared by a single declaration `line`.

    Returns a list — one entry per declared variable, since a declaration may
    name several (``double precision, intent(in) :: a, b``) — or an empty list
    if the line is not a declaration we recognize.
    """
    declaration = parse_declaration(line)
    if declaration is None:
        return []
    variables      = declaration.get('variables')      or []
    variable_names = declaration.get('variableNames')  or []
    return [format_argument(declaration, variable, variable_name)
            for variable, variable_name in zip(variables, variable_names)]


def arguments_to_rst(lines):
    """Render a method's argument declarations, in declaration order."""
    rendered = []
    for line in lines:
        rendered.extend(declaration_to_rst(line))
    return rendered


def procedure_signature(node, skip_passed_object=True):
    """Recover a parsed procedure's signature in the documentation's notation.

    `node` is a `function` or `subroutine` node of a source tree.  Returns
    `(return_type, arguments)`, the first an inline literal (or None for a
    subroutine, or where the type cannot be determined) and the second a list of
    rendered arguments in declaration order.

    `skip_passed_object` drops the leading dummy argument, which for a type-bound
    procedure is the object itself — documenting a method's own object as one of
    its arguments would be noise.  Pass false for a `nopass` binding.

    Arguments whose declarations cannot be found are skipped rather than guessed
    at: a partial signature is still useful, an invented one is not.
    """
    opener = node.get('opener') or ''
    name   = node.get('name') or ''
    m = re.search(r'\b(?:function|subroutine)\s+' + re.escape(name) +
                  r'\s*\(([^)]*)\)', opener, re.IGNORECASE)
    dummies = [d.strip() for d in m.group(1).split(',')] if m else []
    dummies = [d for d in dummies if d]
    if skip_passed_object and dummies:
        dummies = dummies[1:]

    # Index every declaration in the procedure by variable name.
    declared = {}
    child = node.get('firstChild')
    while child:
        for declaration in (child.get('declarations') or []):
            variables      = declaration.get('variables')     or []
            variable_names = declaration.get('variableNames') or []
            for variable, variable_name in zip(variables, variable_names):
                declared[variable_name.lower()] = (declaration, variable,
                                                   variable_name)
        child = child.get('sibling')

    arguments = []
    for dummy in dummies:
        entry = declared.get(dummy.lower())
        if entry is None:
            continue
        arguments.append(format_argument(*entry))

    return_type = None
    if node.get('type') == 'function':
        # A function's type is either a prefix of the `function` statement or a
        # declaration of its result variable.
        m = re.match(r'\s*(?:(?:pure|impure|elemental|recursive|non_recursive'
                     r'|module)\s+)*(.*?)\s*\bfunction\b', opener,
                     re.IGNORECASE)
        prefix = m.group(1).strip() if m else ''
        if prefix:
            return_type = return_type_to_rst(prefix)
        else:
            m = re.search(r'\bresult\s*\(\s*([A-Za-z0-9_]+)\s*\)', opener,
                          re.IGNORECASE)
            entry = declared.get((m.group(1) if m else name).lower())
            if entry is not None:
                declaration, variable, _ = entry
                return_type = return_type_to_rst(type_spec(
                    declaration['intrinsic'], declaration.get('type'),
                    _rank_of(variable, declaration.get('attributes') or [])))
    return return_type, arguments


def split_arguments(text):
    """Split a comma-separated argument list written in the notation.

    The component generator stores a method's arguments as one string; split it
    at the commas that separate arguments rather than at every comma, since an
    array shape may contain them (``double precision(:,:)``).  Every argument
    begins with its inline-literal type spec, which is what marks the boundary.
    """
    if not text:
        return []
    text = text.strip()
    if not text:
        return []
    return [normalize_argument(argument)
            for argument in re.split(r',\s*(?=``)', text) if argument.strip()]


def is_declaration(text):
    """Return true if `text` is a Fortran declaration rather than text already
    written in the argument notation.

    Both forms reach the documentation: hand-written ``<argument>`` elements
    carry raw declarations, while the generators emit the notation directly.
    The `::` separator is what distinguishes them — the notation never contains
    one.
    """
    return '::' in text and parse_declaration(text) is not None


def return_type_to_rst(text):
    """Render a return type as an inline literal, or None if there is none.

    A `void` return (how the component generator spells a subroutine) yields
    None so that the caller renders no return type at all.
    """
    if text is None:
        return None
    text = text.strip()
    if text.lower() in _VOID_TYPES:
        return None
    # Already an inline literal (the component generator's metadata is written
    # that way); do not double-quote it.  A leading `*` inside the literal is
    # how that metadata marks a pointer result; spell it out instead, so that it
    # reads as Fortran and parallels the `(optional)` marker on arguments.
    if text.startswith('``') and text.endswith('``'):
        if text.startswith('``*'):
            return '``' + text[3:] + ' (pointer)'
        return text
    # A rank declared as an attribute (`double precision, dimension(3)`) belongs
    # on the type spec, exactly as for an argument.
    m = re.match(r'(.*?)\s*,\s*dimension\s*(\(.*\))\s*$', text, re.IGNORECASE)
    if m:
        text = m.group(1) + m.group(2)
    return f'``{text}``'
