"""Per-class evolution methods for `floatRank0` meta-properties:
rate / scale / inactive / analytic.

Andrew Benson (ported to Python 2026)
"""

import os


from Galacticus.Build.Components.Utils import register, offset_name


def _is_debugging():
    """Return True iff `GALACTICUS_FCFLAGS` carries `-DDEBUGGING`."""
    flags = os.environ.get('GALACTICUS_FCFLAGS', '')
    for tok in flags.split():
        if tok == '-DDEBUGGING':
            return True
    return False


def Build_Meta_Rate_Functions(build, class_dict):
    """Generate `<class>FloatRank0MetaPropertyRate`."""
    name      = class_dict['name']
    type_name = 'nodeComponent' + _ucfirst(name)
    active    = name in (build.get('componentClassListActive') or [])
    debugging = _is_debugging()

    function = {
        'type':        'void',
        'name':        name + 'FloatRank0MetaPropertyRate',
        'description': (
            f"Accumulate to the rate of change of the indexed rank-0 float "
            f"meta-property of the \\mono{{{name}}} component class."
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
            {
                'intrinsic':  'double precision',
                'attributes': ['intent(in   )'],
                'variables':  ['setValue'],
            },
            {
                'intrinsic':  'logical',
                'attributes': ['optional', 'intent(inout)'],
                'variables':  ['interrupt'],
            },
            {
                'intrinsic':  'procedure',
                'type':       'interruptTask',
                'attributes': ['optional', 'intent(inout)', 'pointer'],
                'variables':  ['interruptProcedure'],
            },
            {
                'intrinsic':  'integer',
                'variables':  ['offset'],
            },
        ],
    }

    if active:
        function['modules'].append('Error')
        if debugging:
            function['modules'].append('Debugging')
            function['variables'].extend([
                {'intrinsic': 'type',      'type': 'varying_string',
                 'variables': ['message']},
                {'intrinsic': 'character', 'type': 'len=13',
                 'variables': ['label']},
            ])

        offset_all      = offset_name('all',      name, 'floatRank0MetaProperties')
        offset_active   = offset_name('active',   name, 'floatRank0MetaProperties')
        offset_inactive = offset_name('inactive', name, 'floatRank0MetaProperties')

        content = (
            "!$GLC attributes unused :: interrupt, interruptProcedure, self\n"
            f"if (.not.{name}FloatRank0MetaPropertyCreator(metaPropertyID)) "
            f"call metaPropertyNoCreator('{name}',char({name}FloatRank0"
            f"MetaPropertyLabels(metaPropertyID)),'float',0)\n"
            f"if (nodeAnalytics({offset_all}(metaPropertyID))) "
            f"call Error_Report('attempt to set rate of analytically-solved "
            "meta-property'//{introspection:location})\n"
            "if (rateComputeState == propertyTypeAll          ) then\n"
            f" offset={offset_all}(metaPropertyID)\n"
            "else if (rateComputeState == propertyTypeActive  ) then\n"
            f" if (     nodeInactives({offset_all}(metaPropertyID))) return\n"
            f" offset={offset_active}(metaPropertyID)\n"
            "else if (rateComputeState == propertyTypeInactive) then\n"
            f" if (.not.nodeInactives({offset_all}(metaPropertyID))) return\n"
            f" offset={offset_inactive}(metaPropertyID)\n"
            "else if (rateComputeState == propertyTypeNumerics) then\n"
            f" offset={offset_active}(metaPropertyID)\n"
            "else\n"
            " return\n"
            "end if\n"
            "nodeRates(offset)=nodeRates(offset)+setValue\n"
        )
        if debugging:
            content += (
                "if (isDebugging()) then\n"
                " write (label,'(e12.6)') setValue\n"
                f" message=\"   rate: (\"//getCaller()//\") {name}:floatRank0"
                "MetaProperties \"//trim(adjustl(label))\n"
                " call debugLog(message)\n"
                "end if\n"
            )
        function['content'] = content

    _bind(build, type_name, function, 'floatRank0metaPropertyRate')


def Build_Meta_Scale_Functions(build, class_dict):
    """Generate `<class>FloatRank0MetaPropertyScale`."""
    name      = class_dict['name']
    type_name = 'nodeComponent' + _ucfirst(name)
    active    = name in (build.get('componentClassListActive') or [])

    function = {
        'type':        'void',
        'name':        name + 'FloatRank0MetaPropertyScale',
        'description': (
            f"Set the absolute scale of the rank-0 float indexed "
            f"meta-property of the \\mono{{{name}}} component class."
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
            {
                'intrinsic':  'double precision',
                'attributes': ['intent(in   )'],
                'variables':  ['setValue'],
            },
        ],
    }

    if active:
        offset = offset_name('all', name, 'floatRank0MetaProperties')
        function['content'] = (
            "!$GLC attributes unused :: self\n"
            f"if (.not.{name}FloatRank0MetaPropertyCreator(metaPropertyID)) "
            f"call metaPropertyNoCreator('{name}',char({name}FloatRank0"
            f"MetaPropertyLabels(metaPropertyID)),'float',0)\n"
            f"nodeScales({offset}(metaPropertyID))=setValue\n"
        )

    _bind(build, type_name, function, 'floatRank0MetaPropertyScale')


def Build_Meta_Inactive_Functions(build, class_dict):
    """Generate `<class>FloatRank0MetaPropertyJcbnZr`."""
    name      = class_dict['name']
    type_name = 'nodeComponent' + _ucfirst(name)
    active    = name in (build.get('componentClassListActive') or [])

    function = {
        'type':        'void',
        'name':        name + 'FloatRank0MetaPropertyJcbnZr',
        'description': (
            f"Indicate that the indexed rank-0 float meta-property of the "
            f"\\mono{{{name}}} component class is inactive for differential "
            "equation solving."
        ),
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
        ],
    }

    if active:
        offset = offset_name('all', name, 'floatRank0MetaProperties')
        function['content'] = (
            "!$GLC attributes unused :: self\n"
            f"nodeInactives({offset}(metaPropertyID))=.true.\n"
        )

    _bind(build, type_name, function, 'floatRank0MetaPropertyInactive')


def Build_Meta_Analytic_Functions(build, class_dict):
    """Generate `<class>FloatRank0MetaPropertyAlytc`."""
    name      = class_dict['name']
    type_name = 'nodeComponent' + _ucfirst(name)
    active    = name in (build.get('componentClassListActive') or [])

    function = {
        'type':        'void',
        'name':        name + 'FloatRank0MetaPropertyAlytc',
        'description': (
            f"Indicate that the indexed rank-0 float meta-property of the "
            f"\\mono{{{name}}} component class is to be solved analytically "
            "for differential evolution."
        ),
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
        ],
    }

    if active:
        function['modules'] = ['Error']
        offset = offset_name('all', name, 'floatRank0MetaProperties')
        function['content'] = (
            "!$GLC attributes unused :: self\n"
            f"if (nodeAnalytics({offset}(metaPropertyID))) "
            "call Error_Report('property is already marked analytically-solvable'"
            "//{introspection:location})\n"
            f"nodeAnalytics({offset}(metaPropertyID))=.true.\n"
        )

    _bind(build, type_name, function, 'floatRank0MetaPropertyAnalytic')


def _bind(build, type_name, function, method_name):
    build.setdefault('types', {}).setdefault(type_name, {}) \
                                  .setdefault('boundFunctions', []) \
                                  .append({
        'type':       'procedure',
        'descriptor': function,
        'name':       method_name,
    })


def _ucfirst(text):
    return text[:1].upper() + text[1:] if text else text


# ---------------------------------------------------------------------------
# Hook registration.  Registration order determines the order of generated
# code — do not reorder.
# ---------------------------------------------------------------------------

register('classesEvolve', 'classIteratedFunctions', Build_Meta_Rate_Functions)
register('classesEvolve', 'classIteratedFunctions', Build_Meta_Scale_Functions)
register('classesEvolve', 'classIteratedFunctions', Build_Meta_Inactive_Functions)
register('classesEvolve', 'classIteratedFunctions', Build_Meta_Analytic_Functions)
