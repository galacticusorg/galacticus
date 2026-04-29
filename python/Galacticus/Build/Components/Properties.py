# Components/Properties — per-property generators.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/Components/Properties.pm.  This file
# currently ports only the small `Construct_Data` hook needed to seed
# the `property['data']` sub-dict; that sub-dict is later read by
# Hierarchy/Utils, Properties/Evolve, and several other downstream
# modules.  The remaining hooks (Class_Defaults_*, Property_Defaults,
# Property_Output_Validate, Data_Validate, Class_Defaults_Validate,
# linked-data construction in Construct_Data) come with the full
# Properties stage.

import os
import sys

sys.path.insert(0, os.path.join(os.environ['GALACTICUS_EXEC_PATH'], 'python'))

from Galacticus.Build.Components       import Utils as _Utils
from Galacticus.Build.Components.Utils import register, _component_properties


def Construct_Data(build):
    """Seed each property's `data` sub-dict and create the per-component
    `content.data` linked-data registry.

    Mirrors `Construct_Data` (Properties.pm:260-318).  Each property
    gets a `data` sub-dict (`type`, `rank`, `isEvolvable`) plus a
    `linkedData` field naming the matching `<prop>Data` storage slot.
    Each component gains a `content.data[<linked-name>] = <data-dict>`
    entry.  Updates the module-level `linked_data_name_length_max`.

    The parent-class attribute-validation portion of the Perl original
    (`definedInParent` / cross-implementation type checking) is
    deferred — none of the currently-ported modules read those fields
    on a path our test exercises.
    """
    for component in (build.get('components') or {}).values():
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

            # Skip linked-data construction for virtual properties (or
            # ones inherited from a parent — handled when we port the
            # full validation logic).
            if attributes.get('isVirtual') or prop['definedInParent']:
                continue

            linked_data_name = prop['name'] + 'Data'
            prop['linkedData'] = linked_data_name
            component_content_data[linked_data_name] = prop['data']

            if len(linked_data_name) > _Utils.linked_data_name_length_max:
                _Utils.linked_data_name_length_max = len(linked_data_name)


# `Construct_Data` runs in the `content` phase, after Gather_Classes.
register('properties', 'content', Construct_Data)
