# Class-default property accessors.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/Components/Classes/Defaults.pm.  Three
# `classIteratedFunctions` hooks emit per-class default-value getters:
# `<prop>IsGettable`, `<prop>` (the default-value lookup), and
# `<prop>RateGet` (a zero-rate fallback).

import os
import re
import sys


from Galacticus.Build.Components.Utils      import (
    register,
    is_intrinsic,
    intrinsic_nulls,
    _component_properties,
)
from Galacticus.Build.Components.DataTypes  import data_object_definition


_SELF_REFERENCE_RE = re.compile(r'self([a-zA-Z]+)\s*%')


def _class_members_by_property(class_dict):
    """Group every property name to a triple: `{name, property, members}`
    where `members` is the subset of class members that supply a
    `classDefault` for that property.

    Mirrors `Class_Members_By_Property`.
    """
    properties = {}
    for member in class_dict.get('members') or []:
        for prop in _component_properties(member):
            entry = properties.setdefault(prop['name'], {
                'name':     prop['name'],
                'property': prop,
                'members':  [],
            })
            if 'classDefault' in prop:
                entry['members'].append(member)
    return properties


def Class_Property_Is_Gettable(build, class_dict):
    """Generate `<class><Prop>IsGettable` per class property —
    returns true when at least one member that supplies that
    property is active at run-time.

    Mirrors `Class_Property_Is_Gettable`.
    """
    class_name = class_dict['name']
    cap_class  = _ucfirst(class_name)
    type_name  = 'nodeComponent' + cap_class

    for prop_entry in _class_members_by_property(class_dict).values():
        prop      = prop_entry['property']
        prop_name = prop_entry['name']
        cap_prop  = _ucfirst(prop_name)

        function = {
            'type':        'logical',
            'name':        f"{class_name}{cap_prop}IsGettable",
            'description': (
                f"Returns true if the \\mono{{{prop_name}}} property is "
                f"gettable for the \\mono{{{class_name}}} component class."
            ),
        }
        content = f"{class_name}{cap_prop}IsGettable=.false.\n"
        for member in prop_entry['members']:
            cap_member = _ucfirst(member['name'])
            content += (
                f"if (nodeComponent{cap_class}{cap_member}IsActiveValue) "
                f"{class_name}{cap_prop}IsGettable=.true.\n"
            )
        function['content'] = content

        build.setdefault('types', {}).setdefault(type_name, {}) \
                                      .setdefault('boundFunctions', []) \
                                      .append({
            'type':       'procedure',
            'descriptor': function,
            'name':       prop_name + 'IsGettable',
            'pass':       'nopass',
        })


def Class_Property_Default(build, class_dict):
    """Generate `<class><Prop>` per class property — returns the
    default value for the property.  Walks active class members that
    declare a `classDefault` and falls back to type-zero when none
    matches.

    Mirrors `Class_Property_Default`.
    """
    class_name = class_dict['name']
    cap_class  = _ucfirst(class_name)
    type_name  = 'nodeComponent' + cap_class

    for prop_entry in _class_members_by_property(class_dict).values():
        prop      = prop_entry['property']
        prop_name = prop_entry['name']
        cap_prop  = _ucfirst(prop_name)
        data      = prop.get('data') or {}

        type_descriptor, type_label = data_object_definition(data)
        return_type = _format_return_type(type_descriptor)

        function = {
            'type':        return_type,
            'name':        f"{class_name}{cap_prop}",
            'description': (
                f"Returns the default value for the \\mono{{{prop_name}}} "
                f"property for the \\mono{{{class_name}}} component class."
            ),
            'variables':   [
                {
                    'intrinsic':  'class',
                    'type':       type_name,
                    'attributes': ['intent(inout)'],
                    'variables':  ['self'],
                },
            ],
        }

        # Discover any cross-component dependencies referenced via
        # `selfXxx%…` inside each member's `classDefault.code`.
        required_components = {'all': set()}
        required_modules    = set()
        for member in prop_entry['members']:
            impl_property = (
                (member.get('properties') or {})
                .get('property') or {}
            )
            if isinstance(impl_property, list):
                impl_property = next(
                    (p for p in impl_property if p.get('name') == prop_name),
                    None,
                )
            elif isinstance(impl_property, dict) and prop_name in impl_property \
                 and isinstance(impl_property[prop_name], dict):
                impl_property = impl_property[prop_name]
            class_default = (impl_property or {}).get('classDefault') or {}
            if 'code' not in class_default:
                continue
            for module in class_default.get('modules') or []:
                required_modules.add(module)
            code = class_default['code']
            for match in _SELF_REFERENCE_RE.findall(code):
                required_components['all'].add(match)
                required_components.setdefault(member['name'], set()).add(match)

        if required_modules:
            function['modules'] = sorted(required_modules)
        if required_components['all']:
            function['variables'].append({
                'intrinsic':  'type',
                'type':       'treeNode',
                'attributes': ['pointer'],
                'variables':  ['node'],
            })
            for component in sorted(required_components['all']):
                function['variables'].append({
                    'intrinsic':  'class',
                    'type':       'nodeComponent' + _ucfirst(component),
                    'attributes': ['pointer'],
                    'variables':  [f'self{_ucfirst(component)}'],
                })

        content = "!$GLC attributes unused :: self\n"
        for member in prop_entry['members']:
            impl_property = _find_member_property(member, prop_name)
            class_default = (impl_property or {}).get('classDefault') or {}
            if 'code' not in class_default:
                continue
            cap_member = _ucfirst(member['name'])
            content += (
                f"if (nodeComponent{cap_class}{cap_member}IsActiveValue) then\n"
            )
            if required_components.get(member['name']):
                content += "   node => self%host()\n"
            for component in sorted(required_components.get(member['name']) or set()):
                content += (
                    f"   self{component} => node%{component}()\n"
                )
            if 'count' in class_default:
                content += f"   allocate(classDefault({class_default['count']}))\n"
            content += (
                f"   classDefault={class_default['code']}\n"
                f"   return\n"
                f"end if\n"
            )

        rank = int(data.get('rank') or 0)
        if rank == 0:
            if is_intrinsic(data.get('type')):
                content += f"classDefault={intrinsic_nulls[data['type']]}\n"
            else:
                content += "call classDefault%reset()\n"
        else:
            content += f"    classDefault=null{type_label}{rank}d\n"
        function['content'] = content

        build.setdefault('types', {}).setdefault(type_name, {}) \
                                      .setdefault('boundFunctions', []) \
                                      .append({
            'type':       'procedure',
            'descriptor': function,
            'name':       prop_name,
        })


def Class_Property_Rate_Default(build, class_dict):
    """Generate `<class><Prop>RateGet` returning a zero rate for the
    property.  Mirrors `Class_Property_Rate_Default`.
    """
    class_name = class_dict['name']
    cap_class  = _ucfirst(class_name)
    type_name  = 'nodeComponent' + cap_class

    for prop_entry in _class_members_by_property(class_dict).values():
        prop      = prop_entry['property']
        prop_name = prop_entry['name']
        cap_prop  = _ucfirst(prop_name)
        data      = prop.get('data') or {}
        rank      = int(data.get('rank') or 0)

        type_descriptor, type_label = data_object_definition(data)
        return_type = _format_return_type(type_descriptor)

        function = {
            'type':        return_type,
            'name':        f"{class_name}{cap_prop}RateGet",
            'description': (
                f"Returns a zero rate for the \\mono{{{prop_name}}} property "
                f"for the \\mono{{{class_name}}} component class."
            ),
            'variables':   [
                {
                    'intrinsic':  'class',
                    'type':       type_name,
                    'attributes': ['intent(inout)'],
                    'variables':  ['self'],
                },
            ],
        }

        content = "!$GLC attributes unused :: self\n"
        if rank == 0:
            if is_intrinsic(data.get('type')):
                content += f"classDefault={intrinsic_nulls[data['type']]}\n"
            else:
                content += "call classDefault%reset()\n"
        else:
            content += f"    classDefault=null{type_label}{rank}d\n"
        function['content'] = content

        build.setdefault('types', {}).setdefault(type_name, {}) \
                                      .setdefault('boundFunctions', []) \
                                      .append({
            'type':       'procedure',
            'descriptor': function,
            'name':       prop_name + 'RateGet',
        })


def _format_return_type(descriptor):
    """Build the `<type>(<inner>), <attrs> => classDefault` return-type
    string from a `data_object_definition` descriptor.  Mirrors the
    Perl ternary that walks `intrinsic` / `type` / `attributes`.
    """
    parts = [descriptor['intrinsic']]
    if 'type' in descriptor:
        parts.append(f"({descriptor['type']})")
    rendered = ''.join(parts)
    if descriptor.get('attributes'):
        rendered += ', ' + ', '.join(descriptor['attributes'])
    rendered += ' => classDefault'
    return rendered


def _find_member_property(member, name):
    """Return the property dict on `member` that has `prop['name'] == name`,
    or None.  Tolerates the three shapes properties take across the
    pipeline (list, dict-by-name, single dict).
    """
    props = (member.get('properties') or {}).get('property')
    if props is None:
        return None
    if isinstance(props, list):
        for p in props:
            if isinstance(p, dict) and p.get('name') == name:
                return p
        return None
    if isinstance(props, dict):
        if name in props and isinstance(props[name], dict):
            return props[name]
        if all(isinstance(v, dict) for v in props.values()):
            for p in props.values():
                if p.get('name') == name:
                    return p
        if props.get('name') == name:
            return props
    return None


def _ucfirst(text):
    return text[:1].upper() + text[1:] if text else text


register('classesDefaults', 'classIteratedFunctions', Class_Property_Is_Gettable)
register('classesDefaults', 'classIteratedFunctions', Class_Property_Default)
register('classesDefaults', 'classIteratedFunctions', Class_Property_Rate_Default)
