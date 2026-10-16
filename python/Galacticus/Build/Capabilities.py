"""Parse the capability markers shared by the functionClass build and the
parameter catalog.

Some implementations of a functionClass can only answer a method when given one
of that method's optional arguments, and some consumers never supply such an
argument. Selecting one for the other then fails part-way through a run, when the
method is first called. These markers record both sides, so that the combination
can instead be rejected before a run (by `scripts/build/parameterValidate.py`)
and when the object is built (by the generated `requires` method, which the
`objectBuilder` directive calls).

An implementation declares what it requires:

    <haloMassFunction name="haloMassFunctionPressSchechter">
     <description>...</description>
     <requires method="differential" argument="node"/>
    </haloMassFunction>

An implementation which passes the arguments it is given on to another object of
the same class declares that it inherits that object's requirements:

    <haloMassFunction name="haloMassFunctionMultiplier">
     <description>...</description>
     <forwards object="massFunction_"/>
    </haloMassFunction>

A consumer declares, on the directive which builds the object, what it never
supplies:

    <objectBuilder class="haloMassFunction" name="haloMassFunction_" source="parameters" withholds="differential:node"/>

Andrew Benson (2026).
"""

import re

from List.ExtraUtils import as_array


_PAIR_RE = re.compile(r'^([A-Za-z][A-Za-z0-9_]*):([A-Za-z][A-Za-z0-9_]*)$')


def parse_withholds(text):
    """Return the ``(method, argument)`` pairs listed in a ``withholds``
    attribute: ``method:argument`` tokens separated by commas or whitespace."""
    pairs = []
    for token in re.split(r'[\s,]+', (text or '').strip()):
        if not token:
            continue
        match = _PAIR_RE.match(token)
        if match is None:
            raise ValueError(
                f"malformed withholds entry '{token}' - expected 'method:argument'")
        pairs.append((match.group(1), match.group(2)))
    return pairs


def parse_requires(directive):
    """Return ``(requires, forwards)`` for an implementation's directive: the
    ``(method, argument)`` pairs it requires, and the names of the objects whose
    requirements it inherits."""
    requires = []
    for entry in as_array(directive.get('requires')):
        if not isinstance(entry, dict) or not entry.get('method') or not entry.get('argument'):
            raise ValueError(
                "a <requires> element needs both 'method' and 'argument' attributes")
        requires.append((entry['method'], entry['argument']))
    forwards = []
    for entry in as_array(directive.get('forwards')):
        if not isinstance(entry, dict) or not entry.get('object'):
            raise ValueError("a <forwards> element needs an 'object' attribute")
        forwards.append(entry['object'])
    return requires, forwards
