"""Consumer-facing queries over the parameter catalog.

Helpers that resolve an implementation's full parameter set (its own parameters
plus those inherited from ancestor implementations) and attach enumeration
allowed-values, for use by documentation, schema, and configuration generators.
The catalog stores each implementation's *own* parameters and a `parent` link;
inheritance is resolved here (mirroring the validator).

Andrew Benson (2026)
"""

__all__ = ['inheritance_chain', 'enumeration_values', 'resolved_parameters']


def inheritance_chain(catalog, impl_type):
    """Return `impl_type` followed by each ancestor implementation it extends,
    stopping when the parent is not a known implementation (the base class)."""
    implementations = catalog.get('implementations', {})
    chain, seen, current = [], set(), impl_type
    while current in implementations and current not in seen:
        seen.add(current)
        chain.append(current)
        current = implementations[current].get('parent')
    return chain


def enumeration_values(catalog, enumeration):
    """Return the allowed labels for an enumeration, or None."""
    if not enumeration:
        return None
    return catalog.get('enumerations', {}).get(enumeration)


def resolved_parameters(catalog, impl_type):
    """Return the ordered, de-duplicated parameter list for an implementation:
    its own parameters first, then inherited ones.  Each entry is a copy of the
    catalog parameter with two added keys: ``inheritedFrom`` (the ancestor
    implementation type, or None if declared directly) and ``allowedValues``
    (the enumeration's labels, or None)."""
    implementations = catalog.get('implementations', {})
    resolved, seen = [], set()
    for type_name in inheritance_chain(catalog, impl_type):
        for parameter in implementations[type_name].get('parameters', []):
            name = parameter.get('name')
            if not name or name in seen:
                continue
            seen.add(name)
            entry = dict(parameter)
            entry['inheritedFrom'] = None if type_name == impl_type else type_name
            entry['allowedValues'] = enumeration_values(
                catalog, parameter.get('enumeration'))
            resolved.append(entry)
    return resolved
