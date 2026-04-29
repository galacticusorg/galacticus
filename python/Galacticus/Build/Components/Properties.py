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

from Galacticus.Build.Components.Utils import register, _component_properties


def Construct_Data(build):
    """Seed each property's `data` sub-dict from its top-level
    `type`/`rank`/`isEvolvable` fields.

    Mirrors the start of `Construct_Data` (Properties.pm:260-280).
    Defers the parent-class validation, `definedInParent` flag, and
    linked-data construction to a later port stage — those features
    are not exercised by the modules ported so far.
    """
    for component in (build.get('components') or {}).values():
        for prop in _component_properties(component):
            attributes = prop.get('attributes') or {}
            prop['data'] = {
                'type':        prop.get('type'),
                'rank':        prop.get('rank'),
                'isEvolvable': attributes.get('isEvolvable'),
            }
            prop.setdefault('definedInParent', False)


# `Construct_Data` runs in the `content` phase, after Gather_Classes.
register('properties', 'content', Construct_Data)
