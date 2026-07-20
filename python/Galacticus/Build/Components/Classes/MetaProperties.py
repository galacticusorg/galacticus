"""Class-level meta-property accessors (get/set) plus the shared
`meta_property_types` data table.

Andrew Benson (ported to Python 2026)
"""



from Galacticus.Build.Components.Utils import register, offset_name
from Galacticus.Build.Components.Properties.Evolve import _is_debugging

# The order of entries fixes the generated identifier suffixes (`Float`,
# `LongInteger`, `Integer` × ranks 0/1) that the rest of the pipeline
# expects — do not reorder.
meta_property_types = [
    {'label': 'float',       'intrinsic': 'double precision',                          'rank': 0},
    {'label': 'float',       'intrinsic': 'double precision',                          'rank': 1},
    {'label': 'longInteger', 'intrinsic': 'integer',          'type': 'kind_int8',     'rank': 0},
    {'label': 'longInteger', 'intrinsic': 'integer',          'type': 'kind_int8',     'rank': 1},
    {'label': 'integer',     'intrinsic': 'integer',                                   'rank': 0},
    {'label': 'integer',     'intrinsic': 'integer',                                   'rank': 1},
]


def Class_Meta_Property_Get(build, class_dict):
    """Generate one `<class><Label>Rank<N>MetaPropertyGet` per
    meta-property type, returning the stored value.
    """
    name      = class_dict['name']
    cap       = _ucfirst(name)
    type_name = 'nodeComponent' + cap
    active    = name in (build.get('componentClassListActive') or [])

    for mpt in meta_property_types:
        rank = mpt['rank']
        cap_label = _ucfirst(mpt['label'])
        rank_attr = (
            ', allocatable, dimension(' + ','.join([':'] * rank) + ')'
            if rank > 0 else ''
        )

        modules = ['ISO_Varying_String']
        if active:
            kind_suffix = (
                f"%values" if rank > 0 else ""
            )
            content = (
                f"if (.not.{name}{cap_label}Rank{rank}MetaPropertyCreator(metaPropertyID))"
                f" call metaPropertyNoCreator('{name}',char({name}{cap_label}"
                f"Rank{rank}MetaPropertyLabels(metaPropertyID)),'{mpt['label']}',"
                f"{rank})\n"
            )
            # In debugging builds, guard reads of inactive float rank-0 meta-properties (the only meta-property kind
            # that supports inactivity) while active-property derivatives are being evaluated - their stored values
            # are stale during that phase, so any such read is a bug (see issue #128).
            if _is_debugging() and mpt['label'] == 'float' and rank == 0:
                offset_all = offset_name('all', name, 'floatRank0MetaProperties')
                content += (
                    "if (evaluationActiveRHS) then\n"
                    f" if (nodeInactives({offset_all}(metaPropertyID))) call Error_Report("
                    f"'value of inactive float rank-0 meta-property \"'//char({name}"
                    f"FloatRank0MetaPropertyLabels(metaPropertyID))//'\" of the {name} "
                    "component was read during active property evolution'"
                    "//{introspection:location})\n"
                    "end if\n"
                )
                modules.append('Error')
            if rank > 0:
                size_args = ','.join(
                    f"size(self%{mpt['label']}Rank{rank}"
                    f"MetaProperties(metaPropertyID)%values,dim={d})"
                    for d in range(1, rank + 1)
                )
                content += f"allocate(value_({size_args}))\n"
            content += (
                f"value_=self%{mpt['label']}Rank{rank}"
                f"MetaProperties(metaPropertyID){kind_suffix}\n"
            )
        else:
            content = ''

        return_type_text = (
            f"{mpt['intrinsic']}"
            + (f"({mpt['type']})" if 'type' in mpt else "")
            + f"{rank_attr} => value_"
        )

        function = {
            'type':        return_type_text,
            'name':        f"{type_name}{cap_label}Rank{rank}MetaPropertyGet",
            'description': (
                f"Return the value of a rank-{rank} {mpt['label']} "
                f"meta-property of a {type_name} component given its ID."
            ),
            'modules':     modules,
            'variables':   [
                {
                    'intrinsic':  'class',
                    'type':       type_name,
                    'attributes': ['intent(in   )'],
                    'variables':  ['self'],
                },
                {
                    'intrinsic':  'integer',
                    'attributes': ['intent(in   )'],
                    'variables':  ['metaPropertyID'],
                },
            ],
            'content':     content,
        }

        build.setdefault('types', {}).setdefault(type_name, {}) \
                                      .setdefault('boundFunctions', []) \
                                      .append({
            'type':       'procedure',
            'descriptor': function,
            'name':       f"{mpt['label']}Rank{rank}MetaPropertyGet",
        })


def Class_Meta_Property_Get_Reference(build, class_dict):
    """Generate one `<class><Label>Rank<N>MetaPropertyGetReference` per
    rank>=1 meta-property type, returning a POINTER to the stored array
    (no allocate, no copy).

    Purely additive sibling of `Class_Meta_Property_Get`.  A pointer to a
    scalar would be pointless, so rank-0 meta-properties are skipped.  The
    returned pointer is read-only by contract and is valid only while the
    component and its allocation persist (see the correctness contract in
    the position-interpolation operator).
    """
    name      = class_dict['name']
    cap       = _ucfirst(name)
    type_name = 'nodeComponent' + cap
    active    = name in (build.get('componentClassListActive') or [])

    for mpt in meta_property_types:
        rank = mpt['rank']
        if rank == 0:                       # a pointer to a scalar is pointless here
            continue
        cap_label = _ucfirst(mpt['label'])
        rank_attr = ', pointer, dimension(' + ','.join([':'] * rank) + ')'

        if active:
            content = (
                f"if (.not.{name}{cap_label}Rank{rank}MetaPropertyCreator(metaPropertyID))"
                f" call metaPropertyNoCreator('{name}',char({name}{cap_label}"
                f"Rank{rank}MetaPropertyLabels(metaPropertyID)),'{mpt['label']}',"
                f"{rank})\n"
                f"value_=>self%{mpt['label']}Rank{rank}"
                f"MetaProperties(metaPropertyID)%values\n"
            )
        else:
            content = ''

        return_type_text = (
            f"{mpt['intrinsic']}"
            + (f"({mpt['type']})" if 'type' in mpt else "")
            + f"{rank_attr} => value_"
        )

        function = {
            'type':        return_type_text,
            # The procedure name is abbreviated to `...MetaPropertyRef` (rather
            # than the full `...MetaPropertyGetReference` used for the public,
            # caller-facing type-bound name below) so that it stays within the
            # Fortran 63-character identifier limit: the longest class+label
            # combination (`darkMatterProfile`+`longInteger`) already reaches
            # 61 characters with the existing `...MetaPropertyGet` suffix.
            'name':        f"{type_name}{cap_label}Rank{rank}MetaPropertyRef",
            'description': (
                f"Return a pointer to the stored rank-{rank} {mpt['label']} "
                f"meta-property of a {type_name} component given its ID (no "
                f"copy; read-only - valid only while the component and its "
                f"allocation persist)."
            ),
            'modules':     ['ISO_Varying_String'],
            'variables':   [
                {
                    'intrinsic':  'class',
                    'type':       type_name,
                    'attributes': ['intent(in   )', 'target'],
                    'variables':  ['self'],
                },
                {
                    'intrinsic':  'integer',
                    'attributes': ['intent(in   )'],
                    'variables':  ['metaPropertyID'],
                },
            ],
            'content':     content,
        }

        build.setdefault('types', {}).setdefault(type_name, {}) \
                                      .setdefault('boundFunctions', []) \
                                      .append({
            'type':       'procedure',
            'descriptor': function,
            'name':       f"{mpt['label']}Rank{rank}MetaPropertyGetReference",
        })


def Class_Meta_Property_Set(build, class_dict):
    """Generate one `<class><Label>Rank<N>MetaPropertySet` per
    meta-property type, storing the provided value.
    """
    name      = class_dict['name']
    cap       = _ucfirst(name)
    type_name = 'nodeComponent' + cap
    active    = name in (build.get('componentClassListActive') or [])

    for mpt in meta_property_types:
        rank = mpt['rank']
        cap_label = _ucfirst(mpt['label'])

        if active:
            kind_suffix = "%values" if rank > 0 else ""
            content = (
                f"if (.not.{name}{cap_label}Rank{rank}MetaPropertyCreator(metaPropertyID))"
                f" call metaPropertyNoCreator('{name}',char({name}{cap_label}"
                f"Rank{rank}MetaPropertyLabels(metaPropertyID)),'{mpt['label']}',"
                f"{rank})\n"
                f"self%{mpt['label']}Rank{rank}MetaProperties(metaPropertyID)"
                f"{kind_suffix}=metaPropertyValue\n"
            )
        else:
            content = ''

        attributes = ['intent(in   )']
        if rank > 0:
            attributes.append('dimension(' + ','.join([':'] * rank) + ')')

        meta_property_variable = {
            'intrinsic':  mpt['intrinsic'],
            'attributes': attributes,
            'variables':  ['metaPropertyValue'],
        }
        if 'type' in mpt:
            meta_property_variable['type'] = mpt['type']

        function = {
            'type':        'void',
            'name':        f"{type_name}{cap_label}Rank{rank}MetaPropertySet",
            'description': (
                f"Set the value of a rank-{rank} {mpt['label']} "
                f"meta-property of a {type_name} component given its ID."
            ),
            'modules':     ['ISO_Varying_String'],
            'variables':   [
                {
                    'intrinsic':  'class',
                    'type':       type_name,
                    'attributes': ['intent(inout)'],
                    'variables':  ['self'],
                },
                {
                    'intrinsic':  'integer',
                    'attributes': ['intent(in   )'],
                    'variables':  ['metaPropertyID'],
                },
                meta_property_variable,
            ],
            'content':     content,
        }

        build.setdefault('types', {}).setdefault(type_name, {}) \
                                      .setdefault('boundFunctions', []) \
                                      .append({
            'type':       'procedure',
            'descriptor': function,
            'name':       f"{mpt['label']}Rank{rank}MetaPropertySet",
        })


def _ucfirst(text):
    return text[:1].upper() + text[1:] if text else text


register('classMetaProperties', 'classIteratedFunctions', Class_Meta_Property_Get          )
register('classMetaProperties', 'classIteratedFunctions', Class_Meta_Property_Get_Reference)
register('classMetaProperties', 'classIteratedFunctions', Class_Meta_Property_Set          )
