"""Descriptor-parameter utilities for the functionClass pipeline.

Andrew Benson (ported to Python 2026)
"""

import re


from List.ExtraUtils                                           import as_array
from Galacticus.Build.StateStorables                           import (
    function_class_names    as _shared_function_class_names,
    function_class_instances as _shared_function_class_instances,
)
from Galacticus.Build.SourceTree.Process.FunctionClass.Utils   import (
    lctrim, trimlc, strip_variable_name,
)


def potential_descriptor_parameters(declarations, non_abstract_class,
                                    class_record, state_storables,
                                    potential_names):
    """Classify every declaration of the class as a potential descriptor
    parameter, mutating `potential_names` in place.

    Categories (keys of `potential_names`):
      - `objects`           — pointer members whose type is a functionClass
                              (instance or class name) listed in
                              `stateStorables.xml`.
      - `statefulTypes`     — `type(stateful{Integer,Double,Logical})` members.
      - `enumerations`      — `type(enumerationXType)` members.
      - `parameters`        — intrinsic-type members (integer / logical /
                              double precision / character / varying_string).
      - `linkedListObjects` / `linkedLists` — populated from the class's
                              `<linkedList>` metadata block.

    `non_abstract_class['hasCustomDescriptor']` is set True when a
    `procedure :: descriptor => ...` bind is present.
    """
    fc_names = _function_class_names(state_storables)
    fc_instances = _function_class_instances(state_storables)

    for declaration in as_array(declarations):
        if not isinstance(declaration, dict):
            continue
        intrinsic  = declaration.get('intrinsic') or ''
        type_text  = declaration.get('type') or ''
        attributes = declaration.get('attributes') or []
        variables  = declaration.get('variables') or []

        # Object pointers — functionClass members with `pointer` attribute.
        if (intrinsic == 'class'
                and (lctrim(type_text) in fc_names
                     or lctrim(type_text) in fc_instances)
                and any(a == 'pointer' for a in attributes)):
            potential_names.setdefault('objects', []).extend(
                strip_variable_name(v) for v in variables
            )

        # Stateful types.
        if (intrinsic == 'type'
                and re.match(r'^stateful(Integer|Double|Logical)\s*$',
                             type_text, re.IGNORECASE)):
            potential_names.setdefault('statefulTypes', []).append(declaration)

        # Enumerations.
        if (intrinsic == 'type'
                and re.match(r'^enumeration[a-z0-9_]+type\s*$',
                             type_text, re.IGNORECASE)):
            potential_names.setdefault('enumerations', []).append(declaration)

        # Regular intrinsic-type parameters.
        if (intrinsic in ('integer', 'logical', 'double precision', 'character')
                or (intrinsic == 'type'
                    and trimlc(type_text) == 'varying_string')):
            potential_names.setdefault('parameters', []).append(declaration)

        # Custom-descriptor type-bound procedure.
        if (intrinsic == 'procedure'
                and variables
                and re.match(r'^descriptor=>', variables[0])):
            non_abstract_class['hasCustomDescriptor'] = True

    # Linked-list objects declared via the class's `<linkedList>` block.
    if class_record is not None and 'linkedList' in class_record:
        linked = class_record['linkedList']
        for obj in (linked.get('object') or '').split():
            ll_objects = potential_names.setdefault('linkedListObjects', [])
            if obj not in ll_objects:
                ll_objects.append(obj)
            potential_names.setdefault('linkedLists', {})[obj] = linked


def _function_class_names(state_storables):
    """Return the set of lowercased functionClass names."""
    return {n.lower() for n in _shared_function_class_names(state_storables)}


def _function_class_instances(state_storables):
    """Return the set of lowercased functionClass instance names."""
    return {n.lower() for n in _shared_function_class_instances(state_storables)}
