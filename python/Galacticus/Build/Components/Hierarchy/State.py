"""Hierarchy initialization state variable.

Andrew Benson (ported to Python 2026)
"""



from Galacticus.Build.Components.Utils import register


def Hierarchy_State(build):
    """Insert a module-scope `hierarchyInitialized` integer counter,
    initialized to 0.
    """
    build.setdefault('variables', []).append({
        'intrinsic':  'integer',
        'variables':  ['hierarchyInitialized=0'],
    })


register('hierarchyState', 'functions', Hierarchy_State)
