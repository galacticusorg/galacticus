"""Class-level meta-property accessors (get/set) plus the shared
`meta_property_types` data table.

Andrew Benson (ported to Python 2026)

Mirrors perl/Galacticus/Build/Components/Classes/MetaProperties.pm.
"""



from Galacticus.Build.Components.Utils import register

# The order of entries is preserved verbatim from the Perl module so that
# generated identifier suffixes (`Float`, `LongInteger`, `Integer` × ranks
# 0/1) keep matching what the rest of the pipeline expects.
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

    Mirrors `Class_Meta_Property_Get`.
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
            'modules':     ['ISO_Varying_String'],
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


def Class_Meta_Property_Set(build, class_dict):
    """Generate one `<class><Label>Rank<N>MetaPropertySet` per
    meta-property type, storing the provided value.

    Mirrors `Class_Meta_Property_Set`.
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


register('classMetaProperties', 'classIteratedFunctions', Class_Meta_Property_Get)
register('classMetaProperties', 'classIteratedFunctions', Class_Meta_Property_Set)
