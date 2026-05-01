# Generic-component ODE-solver helpers.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/Components/Hierarchy/ODESolver.pm.



from Galacticus.Build.Components.Utils import register


def Component_ODE_Name_From_Index(build):
    """Generate `nodeComponentNameFromIndex` — a stub bound to the
    abstract `nodeComponent` type.

    Mirrors `Component_ODE_Name_From_Index`.  The generic parent class
    has no properties of its own, so the returned name is just `'?'`.
    """
    function = {
        'type':        'type(varying_string) => name',
        'name':        'nodeComponentNameFromIndex',
        'description': (
            r"Return the name of the property of given index for a "
            r"\mono{nodeComponent} object."
        ),
        'modules':     ['ISO_Varying_String'],
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       'nodeComponent',
                'attributes': ['intent(in   )'],
                'variables':  ['self'],
            },
            {
                'intrinsic':  'integer',
                'attributes': ['intent(inout)'],
                'variables':  ['count'],
            },
            {
                'intrinsic':  'integer',
                'attributes': ['intent(in   )'],
                'variables':  ['propertyType'],
            },
        ],
        'content': (
            "!$GLC attributes unused :: self, count, propertyType\n"
            "name='?'\n"
        ),
    }
    build.setdefault('types', {}).setdefault('nodeComponent', {}) \
                                 .setdefault('boundFunctions', []) \
                                 .append({
        'type':       'procedure',
        'descriptor': function,
        'name':       'nameFromIndex',
    })


register('hierarchyODESolver', 'functions', Component_ODE_Name_From_Index)
