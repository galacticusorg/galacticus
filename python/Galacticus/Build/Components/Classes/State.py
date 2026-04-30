# Per-class module-scope state variables.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/Components/Classes/State.pm.  Two
# `classIteratedFunctions` hooks: `Class_State` declares the
# class-scope variables (default component, allocator template,
# meta-property bookkeeping arrays); `Class_Size_Of` emits the
# byte-count accessor.

import os
import sys

sys.path.insert(0, os.path.join(os.environ['GALACTICUS_EXEC_PATH'], 'python'))

from Galacticus.Build.Components.Utils                   import register
from Galacticus.Build.Components.Classes.MetaProperties  import meta_property_types


def Class_State(build, class_dict):
    """Declare the module-scope variables that hold runtime state for
    a component class.

    Mirrors `Class_State`.  Always emits the `default<Class>Component`
    holder.  For active classes, also emits an allocator template plus
    six pairs of meta-property bookkeeping arrays (labels, names,
    creator flags, count) and an extra evolvable pair for `floatRank0`.
    """
    name      = class_dict['name']
    cap       = _ucfirst(name)
    type_name = 'nodeComponent' + cap
    active    = name in (build.get('componentClassListActive') or [])

    variables = build.setdefault('variables', [])

    # Default component holder — always declared, even for inactive
    # classes so callers can reference `default<Class>Component`.
    variables.append({
        'intrinsic':  'class',
        'type':       type_name,
        'attributes': ['public', 'allocatable', 'target'],
        'variables':  [f'default{cap}Component'],
    })

    if not active:
        return

    variables.append({
        'intrinsic': 'type',
        'type':      type_name,
        'variables': [f'{name}Class'],
    })

    for mpt in meta_property_types:
        cap_label = _ucfirst(mpt['label'])
        rank      = mpt['rank']
        prefix    = f"{name}{cap_label}Rank{rank}MetaProperty"
        variables.extend([
            {
                'intrinsic':  'type',
                'type':       'varying_string',
                'attributes': ['allocatable', 'dimension(:)'],
                'variables':  [f'{prefix}Labels'],
            },
            {
                'intrinsic':  'type',
                'type':       'varying_string',
                'attributes': ['allocatable', 'dimension(:)'],
                'variables':  [f'{prefix}Names'],
            },
            {
                'intrinsic':  'logical',
                'attributes': ['allocatable', 'dimension(:)'],
                'variables':  [f'{prefix}Creator'],
            },
            {
                'intrinsic':  'integer',
                'variables':  [f'{prefix}Count         =0'],
            },
        ])
        if mpt['label'] == 'float' and rank == 0:
            variables.extend([
                {
                    'intrinsic':  'logical',
                    'attributes': ['allocatable', 'dimension(:)'],
                    'variables':  [f'{prefix}Evolvable'],
                },
                {
                    'intrinsic':  'integer',
                    'variables':  [f'{prefix}EvolvableCount=0'],
                },
            ])


def Class_Size_Of(build, class_dict):
    """Generate `nodeComponent<Class>SizeOf` returning the in-memory
    size of the component instance.  Mirrors `Class_Size_Of`.
    """
    cap       = _ucfirst(class_dict['name'])
    type_name = 'nodeComponent' + cap
    function  = {
        'type':        'integer(c_size_t)',
        'name':        type_name + 'SizeOf',
        'description': f"Return the size in bytes of a {type_name} component.",
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       type_name,
                'attributes': ['intent(in   )'],
                'variables':  ['self'],
            },
        ],
        'content': f"{type_name}SizeOf=sizeof(self)\n",
    }
    build.setdefault('types', {}).setdefault(type_name, {}) \
                                  .setdefault('boundFunctions', []) \
                                  .append({
        'type':       'procedure',
        'descriptor': function,
        'name':       'sizeOf',
    })


def _ucfirst(text):
    return text[:1].upper() + text[1:] if text else text


register('classState', 'classIteratedFunctions', Class_State)
register('classState', 'classIteratedFunctions', Class_Size_Of)
