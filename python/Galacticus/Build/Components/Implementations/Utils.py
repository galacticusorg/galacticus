"""Per-implementation utility hooks.

Andrew Benson (ported to Python 2026)

Mirrors perl/Galacticus/Build/Components/Implementations/Utils.pm.
Two hooks on `functions`:
  - Implementation_Is_Active      — emits one `<member>IsActive`
                                    boolean accessor per (class, member).
  - Implementation_Function_Iterator — drives the
                                    `implementationIteratedFunctions`
                                    sub-iteration phase.
"""



import logging

from Galacticus.Build.Components.Utils import (
    register,
    component_utils,
    _component_properties,
)

logger = logging.getLogger(__name__)


def Implementation_Is_Active(build):
    """Generate one `<class><Member>IsActive` per (class, member).

    Mirrors `Implementation_Is_Active`.  Each function reads the
    matching `IsActiveValue` module variable.
    """
    for class_dict in (build.get('componentClasses') or {}).values():
        cap_class = _ucfirst(class_dict['name'])
        for member in class_dict.get('members') or []:
            cap_member = _ucfirst(member['name'])
            impl_type  = 'nodeComponent' + cap_class + cap_member
            function = {
                'type':        'logical',
                'name':        impl_type + 'IsActive',
                'description': (
                    f"Return true if the {member['name']} implementation of "
                    f"the {class_dict['name']} component is the active choice."
                ),
                'content':     f"{impl_type}IsActive={impl_type}IsActiveValue\n",
            }
            build.setdefault('types', {}) \
                 .setdefault('nodeComponent' + cap_class, {}) \
                 .setdefault('boundFunctions', []) \
                 .append({
                'type':       'procedure',
                'descriptor': function,
                'pass':       'nopass',
                'name':       member['name'] + 'IsActive',
            })


def Implementation_Function_Iterator(build):
    """Drive the `implementationIteratedFunctions` phase.

    Mirrors `Implementation_Function_Iterator`.  Walks the
    `component_utils` registry; for every owner that registered any
    `implementationIteratedFunctions`, calls each function once per
    `(class, member)` pair.
    """
    for owner_name in sorted(component_utils.keys()):
        owner = component_utils[owner_name]
        functions = owner.get('implementationIteratedFunctions')
        if not functions:
            continue
        if not isinstance(functions, list):
            functions = [functions]
        for fn in functions:
            marker = (
                f" {{{getattr(fn, '__name__', '<fn>')}}}"
                if len(functions) > 1 else ''
            )
            logger.info(f"         --> {owner_name}{marker}")
            for class_dict in (build.get('componentClasses') or {}).values():
                for member in class_dict.get('members') or []:
                    fn(build, class_dict, member)


# ---------------------------------------------------------------------------
# Free helpers shared with sister modules — exact ports of the Perl
# `hasRealEvolvers` / `hasRealNonTrivialEvolvers` / `listRealEvolvers`.
# ---------------------------------------------------------------------------

def has_real_evolvers(member):
    """Return True if `member` has at least one non-virtual, evolvable
    property.  Mirrors `hasRealEvolvers`.
    """
    return any(
        not (p.get('attributes') or {}).get('isVirtual')
        and (p.get('data') or {}).get('isEvolvable')
        for p in _component_properties(member)
    )


def has_real_non_trivial_evolvers(member):
    """Return True if `member` has at least one non-virtual, evolvable
    property whose data is rank > 0 or non-double.  Mirrors
    `hasRealNonTrivialEvolvers`.
    """
    for p in _component_properties(member):
        attrs = p.get('attributes') or {}
        data  = p.get('data')       or {}
        if attrs.get('isVirtual'):
            continue
        if not data.get('isEvolvable'):
            continue
        if int(data.get('rank') or 0) > 0:
            return True
        if data.get('type') != 'double':
            return True
    return False


def list_real_evolvers(member):
    """Return the list of non-virtual, evolvable properties on `member`.

    Mirrors `listRealEvolvers`.
    """
    return [
        p for p in _component_properties(member)
        if (
            not (p.get('attributes') or {}).get('isVirtual')
            and (p.get('data')      or {}).get('isEvolvable')
        )
    ]


def _ucfirst(text):
    return text[:1].upper() + text[1:] if text else text


# ---------------------------------------------------------------------------
# Hook registration
# ---------------------------------------------------------------------------

register('implementationUtils', 'functions', Implementation_Is_Active)
register('implementationUtils', 'functions', Implementation_Function_Iterator)
