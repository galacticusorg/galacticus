"""Helper functions used during Fortran code generation for components.

Andrew Benson (ported to Python 2026)
"""

import re

# Recognises a `kind=foo` or just plain `foo` inside a type spec.
_KIND_RE = re.compile(r'\s*(?:kind\s*=\s*)?([a-zA-Z0-9_]+)')


def function_arguments(data):
    """Return the names of every argument-bearing variable in `data`.

    A descriptor contributes its `variables` to the result iff it has at
    least one `intent(...)` attribute.
    """
    arguments = []
    for datum in data or []:
        attributes = datum.get('attributes') or []
        if any('intent' in a for a in attributes):
            arguments.extend(datum.get('variables') or [])
    return arguments


def importables(data):
    """Return the names of every type that must be `import`ed inside an
    abstract interface using `data` as its variable list.

    Rules:
    * `class(<T>)` (T not "*") and `type(<T>)` always import T.
    * Any other intrinsic that carries an explicit kind (e.g.
      `integer(kind=c_int)`) imports the kind name, *unless* the kind text
      contains `len=` (in which case it's a character length and there's
      nothing module-scoped to import).
    """
    out = []
    for datum in data or []:
        intrinsic = datum.get('intrinsic')
        type_text = datum.get('type')
        if intrinsic == 'class' and type_text and type_text != '*':
            out.append(type_text)
        elif intrinsic == 'type' and type_text:
            out.append(type_text)
        elif (intrinsic != 'class'
              and type_text
              and 'len=' not in type_text):
            m = _KIND_RE.search(type_text)
            if m:
                out.append(m.group(1))
    return out
