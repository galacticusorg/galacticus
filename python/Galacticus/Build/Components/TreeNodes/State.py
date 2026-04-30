# Generate `treeNodeSizeOf` — a method that reports the in-memory size of
# a `treeNode` instance.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/Components/TreeNodes/State.pm.

import os
import sys

sys.path.insert(0, os.path.join(os.environ['GALACTICUS_EXEC_PATH'], 'python'))

from Galacticus.Build.Components.Utils import register


def Tree_Node_Size_Of(build):
    """Generate the `sizeOf` type-bound method for `treeNode`.

    Mirrors `Tree_Node_Size_Of`.  The function sums `sizeof(self)`, the
    size of every component-class array currently allocated, and the
    size of every event object hanging off `self%event`.
    """
    function = {
        'type':        'integer(c_size_t)',
        'name':        'treeNodeSizeOf',
        'description': "Compute the size (in bytes) of the tree node.",
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       'treeNode',
                'attributes': ['intent(in   )'],
                'variables':  ['self'],
            },
            {
                'intrinsic':  'integer',
                'variables':  ['i'],
            },
            {
                'intrinsic':  'class',
                'type':       'nodeEvent',
                'attributes': ['pointer'],
                'variables':  ['event'],
            },
        ],
    }

    content = "treeNodeSizeOf=sizeof(self)\n"

    active_classes = set(build.get('componentClassListActive') or [])
    for class_dict in (build.get('componentClasses') or {}).values():
        if class_dict['name'] not in active_classes:
            continue
        cap = _ucfirst(class_dict['name'])
        content += (
            f"if (allocated(self%component{cap})) then\n"
            f"  do i=1,size(self%component{cap})\n"
            f"    treeNodeSizeOf=treeNodeSizeOf+self%component{cap}(i)%sizeOf()\n"
            f"  end do\n"
            f"end if\n"
        )

    content += (
        "event => self%event\n"
        "do while (associated(event))\n"
        "  treeNodeSizeOf=treeNodeSizeOf+sizeof(event)+event%nonStaticSizeOf()\n"
        "  event => event%next\n"
        "end do\n"
    )

    function['content'] = content

    build.setdefault('types', {}).setdefault('treeNode', {}) \
                                 .setdefault('boundFunctions', []) \
                                 .append({
        'type':       'procedure',
        'descriptor': function,
        'name':       'sizeOf',
    })


register('treeNodeState', 'functions', Tree_Node_Size_Of)


def _ucfirst(text):
    return text[:1].upper() + text[1:] if text else text
