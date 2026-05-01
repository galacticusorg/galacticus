# Components/Properties — per-property generators.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/Components/Properties.pm.  Six hooks:
#
#   preValidate → Class_Defaults_Validate, Data_Validate
#   default     → Property_Defaults
#   validate    → Property_Output_Validate
#   gather      → Class_Defaults_Gather
#   scatter     → Class_Defaults_Scatter
#   content     → Construct_Data

import logging
import re

from Galacticus.Build.Components       import Utils as _Utils
from Galacticus.Build.Components.Utils import (
    register,
    apply_defaults,
    _component_properties,
)

logger = logging.getLogger(__name__)


# Default-attribute table for each property.  Mirrors the `%defaults`
# hash inside Property_Defaults (Properties.pm:52-69).
_PROPERTY_DEFAULTS = {
    'properties': {
        'property': {
            'ALL': {
                'attributes': {
                    'isVirtual':      'booleanFalse',
                    'createIfNeeded': 'booleanFalse',
                    'isDeferred':     'false',
                    'isNonNegative':  'booleanFalse',
                },
            },
        },
    },
}


def Property_Defaults(build):
    """Apply per-property default attributes.

    Mirrors `Property_Defaults`.  Walks the nested defaults table and
    coerces missing attributes (`isVirtual` / `createIfNeeded` /
    `isDeferred` / `isNonNegative`) to the matching defaults.
    """
    for impl in (build.get('components') or {}).values():
        for key, default in _PROPERTY_DEFAULTS.items():
            apply_defaults(impl, key, default)


def Data_Validate(build):
    """Verify each property has `type` and `rank`, and forbid
    `createIfNeeded` on rank > 0 properties.

    Mirrors `Data_Validate`.
    """
    for component in (build.get('components') or {}).values():
        for prop in _component_properties(component):
            for required in ('type', 'rank'):
                if required not in prop:
                    raise ValueError(
                        f"Data_Validate: no {required} was specified for "
                        f"the '{prop.get('name', '?')}' property of the '"
                        f"{_lcfirst(component.get('name', '?'))}' component "
                        f"of the '{_lcfirst(component.get('class', '?'))}' "
                        "class"
                    )
            attrs = prop.get('attributes') or {}
            if int(prop.get('rank') or 0) > 0 and attrs.get('createIfNeeded'):
                raise ValueError(
                    "Data_Validate: auto-creation of rank > 0 properties "
                    "is not supported"
                )


_LABELS_BRACKET_RE = re.compile(r'^\[(.*)\]$')


def Property_Output_Validate(build):
    """Validate the `output` block on every property.

    Mirrors `Property_Output_Validate`.  Strips whitespace from
    `labels` and `modules` after validation.
    """
    for component in (build.get('components') or {}).values():
        for prop in _component_properties(component):
            output = prop.get('output')
            if not isinstance(output, dict):
                continue
            attrs = prop.get('attributes') or {}
            if attrs.get('isVirtual') and not attrs.get('isGettable'):
                raise ValueError(
                    "Property_Output_Validate: non-gettable, virtual "
                    "properties can not be output"
                )
            data = prop.get('data') or {}
            rank = int(data.get('rank') or 0)
            if rank > 1:
                raise ValueError(
                    "Property_Output_Validate: output of rank>1 arrays "
                    "is not supported"
                )
            if rank > 0 and 'labels' not in output:
                raise ValueError(
                    "Property_Output_Validate: output of rank>0 objects "
                    "requires a labels attribute"
                )
            labels = output.get('labels')
            if (
                rank > 0
                and labels is not None
                and not _LABELS_BRACKET_RE.match(labels)
                and 'count' not in output
            ):
                raise ValueError(
                    "Property_Output_Validate: output of rank>0 objects "
                    "requires parseable labels or explicit count"
                )
            if (
                labels is not None
                and not _LABELS_BRACKET_RE.match(labels)
                and 'count' not in output
            ):
                raise ValueError(
                    "Property_Output_Validate: no method to determine "
                    "output size for rank-1 property"
                )
            # Strip whitespace from labels and modules.
            for key in ('labels', 'modules'):
                if key in output and isinstance(output[key], str):
                    output[key] = re.sub(r'\s', '', output[key])


def Class_Defaults_Validate(build):
    """Verify that all implementations of a class declare consistent
    `classDefault` blocks for any property they share.

    Mirrors `Class_Defaults_Validate`.
    """
    class_defaults = {}
    for component in (build.get('components') or {}).values():
        for prop in _component_properties(component):
            cd = prop.get('classDefault')
            if cd is None:
                continue
            if isinstance(cd, dict) and 'content' in cd:
                code = cd['content']
            else:
                code = cd

            class_entry = class_defaults.setdefault(component['class'], {}) \
                                         .setdefault(prop['name'], {})
            if 'code' in class_entry:
                if code != class_entry['code']:
                    raise ValueError(
                        f"Class_Defaults_Validate: inconsistent class "
                        f"defaults for property '{prop['name']}' of "
                        f"class '{component['class']}'"
                    )
            else:
                class_entry['code'] = code

            default_count = None
            if isinstance(cd, dict) and 'count' in cd:
                default_count = cd['count']
            elif (
                isinstance(cd, str)
                and cd.startswith('[')
                and cd.endswith(']')
            ):
                default_count = len(cd.split(','))

            if default_count is not None:
                if 'count' in class_entry:
                    if str(default_count) != str(class_entry['count']):
                        raise ValueError(
                            f"Class_Defaults: inconsistent class defaults "
                            f"count for property '{prop['name']}' of class "
                            f"'{component['class']}'"
                        )
                else:
                    class_entry['count'] = default_count


def Class_Defaults_Gather(build):
    """Collect per-class default settings into `build['classDefaults']`.

    Mirrors `Class_Defaults_Gather`.
    """
    class_defaults = build.setdefault('classDefaults', {})
    for component in (build.get('components') or {}).values():
        for prop in _component_properties(component):
            cd = prop.get('classDefault')
            if cd is None:
                continue
            entry = class_defaults.setdefault(component['class'], {}) \
                                   .setdefault(prop['name'], {})
            # Code.
            if isinstance(cd, dict) and 'content' in cd:
                entry['code'] = cd['content']
            else:
                entry['code'] = cd
            # Count.
            if isinstance(cd, dict) and 'count' in cd:
                entry['count'] = cd['count']
            elif (
                isinstance(cd, str)
                and cd.startswith('[')
                and cd.endswith(']')
            ):
                entry['count'] = len(cd.split(','))
            # Modules.
            if isinstance(cd, dict) and 'modules' in cd:
                existing = entry.get('modules') or []
                fresh = re.split(r'\s*,\s*', cd['modules'])
                merged = sorted(set(existing) | set(fresh))
                entry['modules'] = merged


def Class_Defaults_Scatter(build):
    """Push collected class defaults back onto every component property
    that doesn't already declare one.  Mirrors `Class_Defaults_Scatter`.
    """
    class_defaults = build.get('classDefaults') or {}
    for component in (build.get('components') or {}).values():
        for prop in _component_properties(component):
            class_entry = class_defaults.get(component['class'], {}) \
                                         .get(prop['name'])
            if class_entry is None:
                continue
            if 'classDefault' not in prop:
                logger.warning(
                    "         --> property '"
                    f"{prop['name']}' of component '{component['name']}' "
                    f"of class '{component['class']}' is being assigned "
                    "a class default even though it does not have one "
                    "explicitly declared"
                )
            prop['classDefault'] = class_entry


def Construct_Data(build):
    """Seed each property's `data` sub-dict, validate it against any
    parent implementation, and create the per-component `content.data`
    linked-data registry.

    Mirrors `Construct_Data` (Properties.pm:260-318).
    """
    for component in (build.get('components') or {}).values():
        logger.info(
            "         --> Creating linked data objects for "
            f"implementation '{_lcfirst(component.get('name', ''))}'"
            f" of '{_lcfirst(component.get('class', ''))}' class"
        )
        component_content = component.setdefault('content', {})
        component_content_data = component_content.setdefault('data', {})

        for prop in _component_properties(component):
            attributes = prop.get('attributes') or {}
            prop['data'] = {
                'type':        prop.get('type'),
                'rank':        prop.get('rank'),
                'isEvolvable': attributes.get('isEvolvable'),
            }
            prop.setdefault('definedInParent', False)

            # Walk parent implementations and validate matching
            # type / rank / attribute presence.
            parent_implementation = (component.get('extends') or {}).get('implementation')
            while isinstance(parent_implementation, dict):
                parent_props = (parent_implementation.get('properties') or {}).get('property') or {}
                if isinstance(parent_props, list):
                    parent_prop = next(
                        (p for p in parent_props if p.get('name') == prop['name']),
                        None,
                    )
                elif isinstance(parent_props, dict):
                    parent_prop = parent_props.get(prop['name']) \
                        if isinstance(parent_props.get(prop['name']), dict) else None
                else:
                    parent_prop = None
                if parent_prop is not None:
                    prop['definedInParent'] = True
                    for key in ('type', 'rank'):
                        if key not in prop:
                            raise ValueError(
                                "Galacticus::Build::Components::Properties::"
                                f"Construct_Data: property '{prop['name']}' "
                                f"in component '{component['name']}' of "
                                f"class '{component['class']}' lacks "
                                f"'{key}' data characteristic which is "
                                "present in parent implementation"
                            )
                        if str(prop[key]) != str(parent_prop.get(key)):
                            raise ValueError(
                                "Galacticus::Build::Components::Properties::"
                                f"Construct_Data: property '{prop['name']}' "
                                f"in component '{component['name']}' of "
                                f"class '{component['class']}' does not "
                                f"match '{key}' data characteristic in "
                                "parent implementation"
                            )
                    for attr_name in (parent_prop.get('attributes') or {}).keys():
                        if attr_name not in (prop.get('attributes') or {}):
                            raise ValueError(
                                "Galacticus::Build::Components::Properties::"
                                f"Construct_Data: property '{prop['name']}' "
                                f"in component '{component['name']}' of "
                                f"class '{component['class']}' lacks "
                                f"'{attr_name}' attribute which is present "
                                "in parent implementation"
                            )
                parent_implementation = (parent_implementation.get('extends') or {}).get('implementation')

            # Skip linked-data construction for virtual / inherited
            # properties.
            if attributes.get('isVirtual') or prop['definedInParent']:
                continue

            logger.info(f"            --> '{prop['name']}'")

            linked_data_name = prop['name'] + 'Data'
            prop['linkedData'] = linked_data_name
            component_content_data[linked_data_name] = prop['data']
            if len(linked_data_name) > _Utils.linked_data_name_length_max:
                _Utils.linked_data_name_length_max = len(linked_data_name)


def _lcfirst(text):
    return text[:1].lower() + text[1:] if text else text


# ---------------------------------------------------------------------------
# Hook registration.  Order matches Perl Properties.pm:21-44.
# ---------------------------------------------------------------------------

register('properties', 'preValidate', Class_Defaults_Validate)
register('properties', 'preValidate', Data_Validate)
register('properties', 'default',     Property_Defaults)
register('properties', 'validate',    Property_Output_Validate)
register('properties', 'gather',      Class_Defaults_Gather)
register('properties', 'scatter',     Class_Defaults_Scatter)
register('properties', 'content',     Construct_Data)
