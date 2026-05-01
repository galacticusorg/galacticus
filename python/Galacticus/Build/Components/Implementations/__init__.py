# Components/Implementations — per-implementation generators.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/Components/Implementations.pm.



from Sort.Topo                                import sort as topo_sort
from Galacticus.Build.Components.Utils        import (
    register,
    apply_defaults,
    verbosity_level,
    boolean_label,
)
from Galacticus.Build.Components.DataTypes    import data_object_definition


def Default_Full_Name(build):
    """Set `fullyQualifiedName = <Class><Name>` on every component.

    Mirrors `Default_Full_Name`.
    """
    for impl in (build.get('components') or {}).values():
        impl['fullyQualifiedName'] = (
            _ucfirst(impl['class']) + _ucfirst(impl['name'])
        )


def Implementation_ID_List(build):
    """Populate `build['componentIdList']` with a sorted list of component
    implementation IDs (e.g. `["BasicStandard", "BlackHoleSimple", …]`).

    Mirrors `Implementation_ID_List`.
    """
    components = build.get('components') or {}
    build['componentIdList'] = sorted(components.keys())


# Default-attribute table for each implementation.  Mirrors the
# `%defaults` hash inside `Implementation_Defaults` (Implementations.pm:64+).
_IMPLEMENTATION_DEFAULTS = {
    'isDefault':      'booleanFalse',
    'bindings':       {
        'binding':    {
            'isDeferred': 'booleanFalse',
            'interface':  {
                'self':   {
                    'pass': 'booleanFalse',
                },
            },
        },
    },
    'createFunction': {
        'isDeferred': 'booleanFalse',
    },
    'output':         {
        'instances': 'all',
    },
}


def Implementation_Defaults(build):
    """Apply per-implementation default attributes and record the default
    implementation name on each class.

    Mirrors the relevant top portion of `Implementation_Defaults`.
    """
    for impl in (build.get('components') or {}).values():
        for key, default in _IMPLEMENTATION_DEFAULTS.items():
            apply_defaults(impl, key, default)

    # Record `defaultImplementation` on each component class.
    component_classes = build.setdefault('componentClasses', {})
    for impl in (build.get('components') or {}).values():
        if impl.get('isDefault'):
            component_classes.setdefault(impl['class'], {})['defaultImplementation'] = impl['name']


def Null_Implementations(build):
    """Synthesise a `null` implementation for any class that lacks one.

    Mirrors `Null_Implementations`.  When a class has no `null` member,
    we add one; if no other implementation was marked `isDefault`, the
    new null member becomes the default.
    """
    components        = build.setdefault('components', {})
    component_classes = build.setdefault('componentClasses', {})

    classes = {}
    for impl in components.values():
        c = classes.setdefault(impl['class'], {'hasNull': False, 'hasDefault': False})
        if impl.get('name') == 'null':
            c['hasNull'] = True
        if impl.get('isDefault'):
            c['hasDefault'] = True

    for class_name in sorted(classes.keys()):
        info = classes[class_name]
        if info['hasNull']:
            continue
        impl_name = _ucfirst(class_name) + 'Null'
        is_default = not info['hasDefault']
        components[impl_name] = {
            'class':     class_name,
            'name':      'null',
            'isDefault': is_default,
        }
        if is_default:
            component_classes.setdefault(class_name, {})['defaultImplementation'] = 'null'
        # Append to the ID list (Implementation_ID_List ran in
        # `preValidate`, so we keep the trailing-append pattern).
        build.setdefault('componentIdList', []).append(impl_name)
        if verbosity_level >= 1:
            prefix = "         --> Adding null implementation "
            qualifier = "as default " if not info['hasDefault'] else ""
            print(f"{prefix}{qualifier}for {class_name} class")


def Implementation_Dependencies(build):
    """Order each class's members in parent->child order, then rebuild
    `componentIdList` to match.

    Mirrors `Implementation_Dependencies`.
    """
    if verbosity_level >= 1:
        print("         --> Sorting implentations into parent->child order:")

    component_classes = build.get('componentClasses') or {}
    components        = build.get('components')        or {}
    component_class_list = build.get('componentClassList') or []

    for class_name in component_class_list:
        class_dict   = component_classes.get(class_name) or {}
        member_names = list(class_dict.get('memberNames') or [])

        dependencies = {}
        for impl_name in member_names:
            impl_id = _ucfirst(class_name) + _ucfirst(impl_name)
            impl    = components.get(impl_id) or {}
            ext     = impl.get('extends')
            if isinstance(ext, dict) and 'name' in ext:
                dependencies.setdefault(impl_name, []).append(ext['name'])

        sorted_member_names = topo_sort(member_names, dependencies)
        class_dict['memberNames'] = sorted_member_names

        members_by_name = {m['name']: m for m in class_dict.get('members') or []}
        class_dict['members'] = [
            members_by_name[n] for n in sorted_member_names if n in members_by_name
        ]

        if verbosity_level >= 1:
            print(f"            --> {class_name}:")
            for n in sorted_member_names:
                print(f"               --> {n}")

    # Rebuild componentIdList in dependency order: outer iterates classes,
    # inner iterates the freshly-sorted memberNames.  Mirrors the Perl
    # nestedmap at Implementations.pm:178-183.
    new_id_list = []
    for class_name in component_class_list:
        for impl_name in component_classes.get(class_name, {}).get('memberNames') or []:
            new_id_list.append(_ucfirst(class_name) + _ucfirst(impl_name))
    build['componentIdList'] = new_id_list


def Implementation_Parents(build):
    """Resolve each implementation's `extends` link to the actual parent
    implementation dict (key `extends.implementation`).

    Mirrors `Implementation_Parents`.
    """
    components = build.get('components') or {}
    for impl in components.values():
        ext = impl.get('extends')
        if not isinstance(ext, dict):
            continue
        parent_id = _ucfirst(ext.get('class', '')) + _ucfirst(ext.get('name', ''))
        if parent_id in components:
            ext['implementation'] = components[parent_id]


def Build_Component_Implementations(build):
    """Generate one `nodeComponent<FullyQualified>` Fortran type per
    component implementation.

    Mirrors `Build_Component_Implementations`.  Linked-data variables come
    from `implementation['content']['data']` (built by
    Properties.Construct_Data).  Each binding becomes a type-bound
    procedure (deferred bindings are wrapped later by Deferred.pm).
    Each property emits two boolean stubs (`<prop>IsGettable` /
    `<prop>IsSettable`) for the type's contract.
    """
    components = build.get('components') or {}
    for impl_id in sorted(components.keys()):
        impl = components[impl_id]

        # Determine extends: either parent implementation's nodeComponent
        # type, or just the class type (when no explicit extends).
        ext = impl.get('extends')
        if isinstance(ext, dict) and isinstance(ext.get('implementation'), dict):
            extension_of = (
                'nodeComponent' + ext['implementation']['fullyQualifiedName']
            )
        else:
            extension_of = 'nodeComponent' + _ucfirst(impl['class'])

        # Linked-data fields: one Fortran variable per `linkedData` entry.
        data_content = []
        content_data = (impl.get('content') or {}).get('data') or {}
        for linked_name in sorted(content_data.keys()):
            type_definition, _label = data_object_definition(
                content_data[linked_name]
            )
            type_definition['variables'] = [linked_name]
            data_content.append(type_definition)

        type_bound_functions = []

        # Bindings: emit non-deferred ones directly; deferred bindings get
        # a wrapper attached later by Implementations/Deferred.
        bindings = (impl.get('bindings') or {}).get('binding') or []
        if not isinstance(bindings, list):
            bindings = [bindings]
        for binding in bindings:
            if not isinstance(binding, dict):
                continue
            if binding.get('isDeferred'):
                continue
            entry = {
                'type': 'procedure',
                'name': binding['method'],
            }
            entry['function'] = binding['function']
            for attr in ('description', 'returnType', 'arguments'):
                if attr in binding:
                    entry[attr] = binding[attr]
            type_bound_functions.append(entry)

        # IsGettable / IsSettable stubs per property — bound to a single
        # `Boolean_True` / `Boolean_False` free function.
        properties = _component_properties(impl)
        for prop in sorted(properties, key=lambda p: p.get('name', '')):
            prop_name = prop['name']
            attrs     = prop.get('attributes') or {}
            for verb in ('isGettable', 'isSettable'):
                value = bool(attrs.get(verb))
                # boolean_label is ('false', 'true') — index 0/1.
                label = boolean_label[1 if value else 0]
                # Perl: `$propertyName."IsGettable"` (matches the
                # `<prop>IsGettable` accessor declared by Defaults).
                method_suffix = 'I' + verb[1:]
                type_bound_functions.append({
                    'type':     'procedure',
                    'pass':     'nopass',
                    'name':     prop_name + method_suffix,
                    'function': 'Boolean_' + _ucfirst(label),
                })

        full = _ucfirst(impl['fullyQualifiedName'])
        type_name = 'nodeComponent' + full
        build.setdefault('types', {})[type_name] = {
            'name':           type_name,
            'comment':        (
                f"Class for the {impl['name']} implementation of the "
                f"{impl['class']} component."
            ),
            'isPublic':       True,
            'extends':        extension_of,
            'boundFunctions': type_bound_functions,
            'dataContent':    data_content,
        }


def _component_properties(impl):
    """Yield every property dict declared on `impl`."""
    props = (impl.get('properties') or {}).get('property')
    if props is None:
        return []
    if isinstance(props, list):
        return props
    if isinstance(props, dict):
        if all(isinstance(v, dict) for v in props.values()):
            return list(props.values())
        return [props]
    return []


# ---------------------------------------------------------------------------
# Hook registration
# ---------------------------------------------------------------------------

# `Implementation_ID_List` registers under `preValidate` so it runs before
# `default`-phase hooks that walk `components` and may mutate the dict.
register('implementations', 'preValidate', Implementation_ID_List)
register('implementations', 'default',     Implementation_Defaults)
register('implementations', 'default',     Null_Implementations)
register('implementations', 'default',     Default_Full_Name)
register('implementations', 'gather',      Implementation_Dependencies)
register('implementations', 'gather',      Implementation_Parents)
register('implementations', 'types',       Build_Component_Implementations)


def _ucfirst(text):
    return text[:1].upper() + text[1:] if text else text
