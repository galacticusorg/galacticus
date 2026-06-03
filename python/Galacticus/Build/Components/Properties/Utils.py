"""Per-property utilities + the `propertyIteratedFunctions` phase iterator.

Andrew Benson (ported to Python 2026)

Mirrors perl/Galacticus/Build/Components/Properties/Utils.pm.
"""



import logging

from Galacticus.Build.Components.Utils import (
    register,
    component_utils,
    _component_properties,
)

logger = logging.getLogger(__name__)


# Maps an action verb (`get` / `set` / `rate`) to the matching property
# attribute key.  Mirrors `%attributeAdjective` at Properties/Utils.pm:30.
attribute_adjective = {
    'get':  'isGettable',
    'set':  'isSettable',
    'rate': 'isEvolvable',
}


def Property_Function_Iterator(build):
    """Drive the `propertyIteratedFunctions` phase.

    Mirrors `Property_Function_Iterator`.  Walks the `component_utils`
    registry; for every owner that registered any
    `propertyIteratedFunctions`, calls each function once per
    `(class, member, property)` triple.
    """
    for owner_name in sorted(component_utils.keys()):
        owner = component_utils[owner_name]
        functions = owner.get('propertyIteratedFunctions')
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
                    for prop in _component_properties(member):
                        fn(build, class_dict, member, prop)


register('propertyUtils', 'functions', Property_Function_Iterator)
