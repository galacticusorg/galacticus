"""Hierarchy initialisation state variable.

Andrew Benson (ported to Python 2026)

Mirrors perl/Galacticus/Build/Components/Hierarchy/State.pm.
"""



from Galacticus.Build.Components.Utils import register


def Hierarchy_State(build):
    """Insert a module-scope `hierarchyInitialized` integer counter,
    initialised to 0.  Mirrors `Hierarchy_State`.
    """
    build.setdefault('variables', []).append({
        'intrinsic':  'integer',
        'variables':  ['hierarchyInitialized=0'],
    })


register('hierarchyState', 'functions', Hierarchy_State)
