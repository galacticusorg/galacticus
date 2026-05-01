# `treeNode` HDF5-output methods.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/Components/TreeNodes/Output.pm.  Four
# hooks on the `functions` phase: `outputCount`, `outputNames`,
# `output`, and `postOutput`.  Each emits boilerplate for the
# node-level output properties (currently just `subsamplingWeight`)
# and then delegates to every active component class member.

import os
import sys


from Galacticus.Build.Components.Utils import register


# ---------------------------------------------------------------------------
# Node-level output properties.  Mirrors the small `@nodePropertiesOutputList`
# array at Output.pm:31-42.
# ---------------------------------------------------------------------------

_NODE_PROPERTIES_OUTPUT_LIST = [
    {
        'type':             'double',
        'name':             'subsamplingWeight',
        'description':      "Weight of node in the subsample.",
        'unitsInSI':        '1.0d0',
        'unitsDescription': '',
        'unitsQuantity':    '',
        'isComoving':       '.false.',
    },
]


def Tree_Node_Output_Count(build):
    """Generate `treeNodeOutputCount`.  Mirrors `Tree_Node_Output_Count`."""
    function = {
        'type':        'void',
        'name':        'treeNodeOutputCount',
        'description': (
            r"Increment the count of properties to output for a "
            r"\mono{treeNode}."
        ),
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       'treeNode',
                'attributes': ['intent(inout)'],
                'variables':  ['self'],
            },
            {
                'intrinsic':  'integer',
                'attributes': ['intent(inout)'],
                'variables':  ['integerPropertyCount', 'doublePropertyCount'],
            },
            {
                'intrinsic':  'double precision',
                'attributes': ['intent(in   )'],
                'variables':  ['time'],
            },
            {
                'intrinsic':  'integer',
                'variables':  ['i'],
            },
        ],
    }

    content = ''
    for prop in _NODE_PROPERTIES_OUTPUT_LIST:
        if prop['type'] == 'integer':
            content += "integerPropertyCount=integerPropertyCount+1\n"
        elif prop['type'] == 'double':
            content += "doublePropertyCount=doublePropertyCount+1\n"
    for class_dict in _active_classes(build):
        cap = _ucfirst(class_dict['name'])
        content += (
            f"if (allocated(self%component{cap})) then\n"
            f"  do i=1,size(self%component{cap})\n"
            f"    call self%component{cap}(i)%outputCount("
            f"integerPropertyCount,doublePropertyCount,time,instance=i)\n"
            f"  end do\n"
            f"end if\n"
        )
    function['content'] = content
    _bind(build, 'outputCount', function)


def Tree_Node_Output_Names(build):
    """Generate `treeNodeOutputNames`.  Mirrors `Tree_Node_Output_Names`."""
    function = {
        'type':        'void',
        'name':        'treeNodeOutputNames',
        'description': (
            r"Establish the names of properties to output for a "
            r"\mono{treeNode}."
        ),
        'modules':     [
            'Merger_Tree_Outputter_Buffer_Types',
            'Units_MetaData',
        ],
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       'treeNode',
                'attributes': ['intent(inout)'],
                'variables':  ['self'],
            },
            {
                'intrinsic':  'integer',
                'attributes': ['intent(inout)'],
                'variables':  ['integerProperty'],
            },
            {
                'intrinsic':  'type',
                'type':       'outputPropertyInteger',
                'attributes': ['intent(inout)', 'dimension(:)'],
                'variables':  ['integerProperties'],
            },
            {
                'intrinsic':  'integer',
                'attributes': ['intent(inout)'],
                'variables':  ['doubleProperty'],
            },
            {
                'intrinsic':  'type',
                'type':       'outputPropertyDouble',
                'attributes': ['intent(inout)', 'dimension(:)'],
                'variables':  ['doubleProperties'],
            },
            {
                'intrinsic':  'double precision',
                'attributes': ['intent(in   )'],
                'variables':  ['time'],
            },
            {
                'intrinsic':  'integer',
                'variables':  ['i'],
            },
        ],
    }

    content = ''
    for prop in _NODE_PROPERTIES_OUTPUT_LIST:
        cap_name          = _ucfirst(prop['name'])
        units_in_si       = prop['unitsInSI']
        units_description = prop.get('unitsDescription') or ''
        units_quantity    = prop.get('unitsQuantity')    or ''
        is_comoving       = prop.get('isComoving')       or '.false.'
        if prop['type'] == 'integer':
            counter, target = 'integerProperty', 'integerProperties'
        elif prop['type'] == 'double':
            counter, target = 'doubleProperty',  'doubleProperties'
        else:
            continue
        content += (
            f"{counter}={counter}+1\n"
            f"{target}({counter})%name     =\"node{cap_name}\"\n"
            f"{target}({counter})%comment  =\"{prop['description']}\"\n"
            f"{target}({counter})%units    ="
            f"unitType(unitsInSI={units_in_si},"
            f"description='{units_description}',"
            f"quantity='{units_quantity}',"
            f"isComoving={is_comoving})\n"
        )
    for class_dict in _active_classes(build):
        cap = _ucfirst(class_dict['name'])
        content += (
            f"if (allocated(self%component{cap})) then\n"
            f"  do i=1,size(self%component{cap})\n"
            f"    call self%component{cap}(i)%outputNames("
            f"integerProperty,integerProperties,"
            f"doubleProperty,doubleProperties,time,instance=i)\n"
            f"  end do\n"
            f"end if\n"
        )
    function['content'] = content
    _bind(build, 'outputNames', function)


def Tree_Node_Output(build):
    """Generate `treeNodeOutput`.  Mirrors `Tree_Node_Output`."""
    function = {
        'type':        'void',
        'name':        'treeNodeOutput',
        'description': (
            r"Populate output buffers with properties to output for a "
            r"\mono{treeNode}."
        ),
        'modules':     [
            'Multi_Counters',
            'Merger_Tree_Outputter_Buffer_Types',
        ],
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       'treeNode',
                'attributes': ['intent(inout)'],
                'variables':  ['self'],
            },
            {
                'intrinsic':  'integer',
                'attributes': ['intent(inout)'],
                'variables':  ['integerProperty', 'integerBufferCount'],
            },
            {
                'intrinsic':  'type',
                'type':       'outputPropertyInteger',
                'attributes': ['intent(inout)', 'dimension(:)'],
                'variables':  ['integerProperties'],
            },
            {
                'intrinsic':  'integer',
                'attributes': ['intent(inout)'],
                'variables':  ['doubleProperty', 'doubleBufferCount'],
            },
            {
                'intrinsic':  'type',
                'type':       'outputPropertyDouble',
                'attributes': ['intent(inout)', 'dimension(:)'],
                'variables':  ['doubleProperties'],
            },
            {
                'intrinsic':  'double precision',
                'attributes': ['intent(in   )'],
                'variables':  ['time'],
            },
            {
                'intrinsic':  'type',
                'type':       'multiCounter',
                'attributes': ['intent(in   )'],
                'variables':  ['outputInstance'],
            },
            {
                'intrinsic':  'integer',
                'variables':  ['i'],
            },
        ],
    }

    content = ''
    for prop in _NODE_PROPERTIES_OUTPUT_LIST:
        if prop['type'] == 'integer':
            content += (
                "integerProperty=integerProperty+1\n"
                f"integerProperties(integerProperty)%scalar(integerBufferCount)"
                f"=self%{prop['name']}()\n"
            )
        elif prop['type'] == 'double':
            content += (
                "doubleProperty=doubleProperty+1\n"
                f"doubleProperties(doubleProperty)%scalar(doubleBufferCount)"
                f"=self%{prop['name']}()\n"
            )
    for class_dict in _active_classes(build):
        cap = _ucfirst(class_dict['name'])
        content += (
            f"if (allocated(self%component{cap})) then\n"
            f"  do i=1,size(self%component{cap})\n"
            f"    call self%component{cap}(i)%output("
            f"integerProperty,integerBufferCount,integerProperties,"
            f"doubleProperty,doubleBufferCount,doubleProperties,"
            f"time,outputInstance,instance=i)\n"
            f"  end do\n"
            f"end if\n"
        )
    function['content'] = content
    _bind(build, 'output', function)


def Tree_Node_Post_Output(build):
    """Generate `treeNodePostOutput`.  Mirrors `Tree_Node_Post_Output`."""
    function = {
        'type':        'void',
        'name':        'treeNodePostOutput',
        'description': (
            r"Perform post-output processing of a \mono{treeNode}."
        ),
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       'treeNode',
                'attributes': ['intent(inout)'],
                'variables':  ['self'],
            },
            {
                'intrinsic':  'double precision',
                'attributes': ['intent(in   )'],
                'variables':  ['time'],
            },
            {
                'intrinsic':  'integer',
                'variables':  ['i'],
            },
        ],
    }

    content = ''
    for class_dict in _active_classes(build):
        cap = _ucfirst(class_dict['name'])
        content += (
            f"if (allocated(self%component{cap})) then\n"
            f"  do i=1,size(self%component{cap})\n"
            f"    call self%component{cap}(i)%postOutput(time)\n"
            f"  end do\n"
            f"end if\n"
        )
    function['content'] = content
    _bind(build, 'postOutput', function)


# ---------------------------------------------------------------------------
# Helpers (shared with sister modules)
# ---------------------------------------------------------------------------

def _active_classes(build):
    active = set(build.get('componentClassListActive') or [])
    for class_dict in (build.get('componentClasses') or {}).values():
        if class_dict['name'] in active:
            yield class_dict


def _bind(build, method_name, function):
    build.setdefault('types', {}).setdefault('treeNode', {}) \
                                 .setdefault('boundFunctions', []) \
                                 .append({
        'type':       'procedure',
        'descriptor': function,
        'name':       method_name,
    })


def _ucfirst(text):
    return text[:1].upper() + text[1:] if text else text


# ---------------------------------------------------------------------------
# Hook registration.  Order matches Perl Output.pm:21-26.
# ---------------------------------------------------------------------------

register('treeNodeOutput', 'functions', Tree_Node_Output_Count)
register('treeNodeOutput', 'functions', Tree_Node_Output_Names)
register('treeNodeOutput', 'functions', Tree_Node_Post_Output)
register('treeNodeOutput', 'functions', Tree_Node_Output)
