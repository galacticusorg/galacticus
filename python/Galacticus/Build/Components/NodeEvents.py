# Build the `nodeEvent` class hierarchy and supporting (de)serialization.
# Andrew Benson (ported to Python 2026)
#
# Mirrors perl/Galacticus/Build/Components/NodeEvents.pm.  Seven hooks:
#
#   types      → Build_Node_Event_Class
#   interfaces → Node_Event_Task_Interface, Node_Event_Merge_Time_Set_Interface
#   functions  → Node_Event_Non_Static_Size_Of, Node_Event_Serialize_Raw,
#                Node_Event_Deserialize_Raw, Node_Event_Deserialize_Raw_Polymorphic

import os
import sys


from Galacticus.Build.Components.Utils import register


# Static table of node-event classes.  Mirrors the `@{$build->{nodeEventClasses}}`
# array assembled at NodeEvents.pm:42-166.  Index in this list is the
# `classCount` integer the polymorphic builder reads from disk.
_NODE_EVENT_CLASSES = [
    {
        'name':        'nodeEvent',
        'description': "Base class for events attached to nodes.",
        'data':        [
            {
                'intrinsic':  'integer',
                'type':       'kind=kind_int8',
                'attributes': ['public'],
                'variables':  ['ID'],
            },
            {
                'intrinsic':  'type',
                'type':       'treeNode',
                'attributes': ['pointer', 'public'],
                'variables':  ['node => null()'],
            },
            {
                'intrinsic':  'double precision',
                'attributes': ['public'],
                'variables':  ['time'],
            },
            {
                'intrinsic':  'class',
                'type':       'nodeEvent',
                'attributes': ['public', 'pointer'],
                'variables':  ['next => null()'],
            },
            {
                'intrinsic':  'procedure',
                'type':       'nodeEventTask',
                'attributes': ['public', 'pointer'],
                'variables':  ['task'],
            },
        ],
    },
    {
        'name':        'nodeEventBranchJump',
        'description': "Class for branch jump events attached to nodes.",
        'extends':     'nodeEvent',
        'data':        [],
    },
    {
        'name':        'nodeEventSubhaloPromotion',
        'description': "Class for subhalo promotion events attached to nodes.",
        'extends':     'nodeEvent',
        'data':        [],
    },
    {
        'name':        'nodeEventBranchJumpInterTree',
        'description': "Class for inter-tree branch jump events attached to nodes.",
        'extends':     'nodeEvent',
        'data':        [
            {
                'intrinsic':  'integer',
                'type':       'kind=c_size_t',
                'attributes': ['public'],
                'variables':  ['splitForestUniqueID'],
            },
            {
                'intrinsic':  'integer',
                'type':       'kind=kind_int8',
                'attributes': ['public'],
                'variables':  ['pairedNodeID'],
            },
            {
                'intrinsic':  'logical',
                'attributes': ['public'],
                'variables':  ['isPrimary', 'hasSecondary'],
            },
            {
                'intrinsic':  'class',
                'type':       '*',
                'attributes': ['pointer', 'public'],
                'variables':  ['creator'],
            },
            {
                'intrinsic':  'procedure',
                'type':       'nodeEventInterTreeMergeTimeSet',
                'attributes': ['pointer', 'nopass', 'public'],
                'variables':  ['mergeTimeSet'],
            },
        ],
    },
    {
        'name':        'nodeEventSubhaloPromotionInterTree',
        'description': "Class for inter-tree subhalo promotion events attached to nodes.",
        'extends':     'nodeEvent',
        'data':        [
            {
                'intrinsic':  'integer',
                'type':       'kind=c_size_t',
                'attributes': ['public'],
                'variables':  ['splitForestUniqueID'],
            },
            {
                'intrinsic':  'integer',
                'type':       'kind=kind_int8',
                'attributes': ['public'],
                'variables':  ['pairedNodeID'],
            },
            {
                'intrinsic':  'logical',
                'attributes': ['public'],
                'variables':  ['isPrimary'],
            },
            {
                'intrinsic':  'class',
                'type':       '*',
                'attributes': ['pointer', 'public'],
                'variables':  ['creator'],
            },
            {
                'intrinsic':  'procedure',
                'type':       'nodeEventInterTreeMergeTimeSet',
                'attributes': ['pointer', 'nopass', 'public'],
                'variables':  ['mergeTimeSet'],
            },
        ],
    },
]


def Build_Node_Event_Class(build):
    """Define one Fortran type per entry in `_NODE_EVENT_CLASSES`.
    Mirrors `Build_Node_Event_Class`.

    Stores the class table on `build['nodeEventClasses']` so the
    sister hooks can iterate it; sister hooks read this dynamic table
    rather than the static one above (they were written that way in
    the Perl original).
    """
    build['nodeEventClasses'] = _NODE_EVENT_CLASSES
    for entry in _NODE_EVENT_CLASSES:
        # Mirror Perl literal `isPublic => "true"` (string, not bool).
        type_def = {
            'name':        entry['name'],
            'comment':     entry['description'],
            'isPublic':    "true",
            'dataContent': entry['data'],
        }
        if 'extends' in entry:
            type_def['extends'] = entry['extends']
        build.setdefault('types', {})[entry['name']] = type_def


def Node_Event_Task_Interface(build):
    """Mirrors `Node_Event_Task_Interface`."""
    build.setdefault('interfaces', {})['nodeEventTask'] = {
        'name':      'nodeEventTask',
        'comment':   "Interface for node event tasks.",
        'intrinsic': 'logical',
        'data':      [
            {
                'intrinsic':  'class',
                'type':       'nodeEvent',
                'attributes': ['intent(in   )'],
                'variables':  ['thisEvent'],
            },
            {
                'intrinsic':  'type',
                'type':       'treeNode',
                'attributes': ['pointer', 'intent(inout)'],
                'variables':  ['thisNode'],
            },
            {
                'intrinsic':  'type',
                'type':       'enumerationDeadlockStatusType',
                'attributes': ['intent(inout)'],
                'variables':  ['deadlockStatus'],
            },
        ],
    }


def Node_Event_Merge_Time_Set_Interface(build):
    """Mirrors `Node_Event_Merge_Time_Set_Interface`."""
    build.setdefault('interfaces', {})['nodeEventInterTreeMergeTimeSet'] = {
        'name':      'nodeEventInterTreeMergeTimeSet',
        'comment':   "Interface for node event inter tree merge time set functions.",
        'intrinsic': 'void',
        'data':      [
            {
                'intrinsic':  'class',
                'type':       '*',
                'attributes': ['intent(inout)'],
                'variables':  ['self'],
            },
            {
                'intrinsic':  'type',
                'type':       'treeNode',
                'attributes': ['intent(inout)', 'target'],
                'variables':  ['nodeSatellite', 'nodeHost'],
            },
        ],
    }


def Node_Event_Deserialize_Raw(build):
    """Generate per-class `<class>DeserializeRaw` methods.
    Mirrors `Node_Event_Deserialize_Raw`.
    """
    for entry in build.get('nodeEventClasses') or []:
        function = {
            'type':        'void',
            'name':        entry['name'] + 'DeserializeRaw',
            'description': (
                f"Deserialize a \\mono{{{entry['name']}}} object from raw file."
            ),
            'content':     '',
            'modules':     ['ISO_C_Binding'],
            'variables':   [
                {
                    'intrinsic':  'class',
                    'type':       entry['name'],
                    'variables':  ['self'],
                    'attributes': ['intent(inout)'],
                },
                {
                    'intrinsic':  'integer',
                    'variables':  ['fileUnit'],
                    'attributes': ['intent(in   )'],
                },
            ],
        }
        if entry['name'] == 'nodeEvent':
            function['variables'].extend([
                {'intrinsic': 'integer', 'variables': ['pointerAssociated']},
                {'intrinsic': 'type',    'type': 'c_funptr',
                 'variables': ['functionLocation']},
                {'intrinsic': 'integer', 'variables': ['functionLocation_'],
                 'attributes': ['allocatable', 'dimension(:)']},
                {'intrinsic': 'integer', 'variables': ['functionLocationSize']},
            ])

        content = ''
        if 'extends' in entry:
            content += (
                "! Read the parent class.\n"
                f"call self%{entry['extends']}%deserializeRaw(fileUnit)\n"
            )
        for data in entry['data']:
            attrs = data.get('attributes') or []
            if 'pointer' in attrs:
                continue
            vars_list = ','.join(f"self%{v}" for v in data.get('variables') or [])
            content += f"read (fileUnit) {vars_list}\n"

        if entry['name'] == 'nodeEvent':
            content += (
                "read (fileUnit) pointerAssociated,functionLocationSize\n"
                "if (pointerAssociated == 1) then\n"
                "   allocate(functionLocation_(functionLocationSize))\n"
                "   read (fileUnit) functionLocation_\n"
                "   functionLocation=transfer(functionLocation_,functionLocation)\n"
                "   call c_f_ProcPointer(functionLocation,self%task)\n"
                "else\n"
                "   self%task => null()\n"
                "end if\n"
            )

        function['content'] = content
        _bind(build, entry['name'], function, 'deserializeRaw')


def Node_Event_Serialize_Raw(build):
    """Generate per-class `<class>SerializeRaw` methods.
    Mirrors `Node_Event_Serialize_Raw`.
    """
    for class_count, entry in enumerate(build.get('nodeEventClasses') or []):
        function = {
            'type':        'void',
            'name':        entry['name'] + 'SerializeRaw',
            'description': (
                f"Serialize a \\mono{{{entry['name']}}} object to raw file."
            ),
            'modules':     ['ISO_C_Binding'],
            'variables':   [
                {
                    'intrinsic':  'class',
                    'type':       entry['name'],
                    'variables':  ['self'],
                    'attributes': ['intent(in   )'],
                },
                {
                    'intrinsic':  'integer',
                    'variables':  ['fileUnit'],
                    'attributes': ['intent(in   )'],
                },
                {
                    'intrinsic':  'logical',
                    'variables':  ['includeType'],
                    'attributes': ['intent(in   ), optional'],
                },
            ],
        }
        if entry['name'] == 'nodeEvent':
            function['variables'].extend([
                {'intrinsic': 'type', 'type': 'c_funptr',
                 'variables': ['functionLocation']},
                {'intrinsic': 'integer', 'variables': ['functionLocation_'],
                 'attributes': ['allocatable', 'dimension(:)']},
            ])

        content = (
            "! Write an integer indicating the type of this event if requested.\n"
            f"if (.not.present(includeType).or.includeType) write (fileUnit) {class_count}\n"
        )
        if 'extends' in entry:
            content += (
                "! Serialize the parent class.\n"
                f"call self%{entry['extends']}%serializeRaw(fileUnit,.false.)\n"
            )
        for data in entry['data']:
            attrs = data.get('attributes') or []
            if 'pointer' in attrs:
                continue
            vars_list = ','.join(f"self%{v}" for v in data.get('variables') or [])
            content += f"write (fileUnit) {vars_list}\n"

        if entry['name'] == 'nodeEvent':
            content += (
                "if (associated(self%task)) then\n"
                "   functionLocation =c_FunLoc(self%task)\n"
                "   functionLocation_=transfer(functionLocation,functionLocation_)\n"
                "   write (fileUnit) 1,size(functionLocation_)\n"
                "   write (fileUnit) functionLocation_\n"
                "else\n"
                "   write (fileUnit) 0,0\n"
                "end if\n"
            )

        function['content'] = content
        _bind(build, entry['name'], function, 'serializeRaw')


def Node_Event_Deserialize_Raw_Polymorphic(build):
    """Generate the top-level `nodeEventBuildFromRaw` polymorphic
    deserializer.  Mirrors `Node_Event_Deserialize_Raw_Polymorphic`.
    """
    function = {
        'type':        'class(nodeEvent), pointer => event',
        'name':        'nodeEventBuildFromRaw',
        'description': (
            r"Build a \mono{nodeEvent} class object from a raw dump file."
        ),
        'modules':     ['Error'],
        'variables':   [
            {
                'intrinsic':  'integer',
                'variables':  ['fileUnit'],
                'attributes': ['intent(in   )'],
            },
            {
                'intrinsic':  'integer',
                'variables':  ['classType'],
            },
        ],
    }
    content = (
        "read (fileUnit) classType\n"
        "select case (classType)\n"
    )
    for class_count, entry in enumerate(build.get('nodeEventClasses') or []):
        content += (
            f"case ({class_count})\n"
            f"  allocate({entry['name']} :: event)\n"
        )
    content += (
        "case default\n"
        "   call Error_Report('unknown class type'//{introspection:location})\n"
        "end select\n"
        "call event%deserializeRaw(fileUnit)\n"
    )
    function['content'] = content
    build.setdefault('functions', []).append(function)


def Node_Event_Non_Static_Size_Of(build):
    """Generate `nodeEventSizeOf`.  Mirrors `Node_Event_Non_Static_Size_Of`.

    The Perl original interpolates `$code::class->{'name'}` into the
    description, but `$code::class` is never set in this hook's scope —
    it's a sister-hook iterator variable, so Perl evaluates the deref
    as undef and emits `\\mono{}`.  We fix the bug by hard-coding
    `nodeEvent` (the actual `class(...)` of `self` in the generated
    function), which is what the description was clearly meant to say.
    """
    function = {
        'type':        'integer(c_size_t)',
        'name':        'nodeEventSizeOf',
        'description': (
            r"Compute the size of the non-static parts of a \mono{nodeEvent} object."
        ),
        'variables':   [
            {
                'intrinsic':  'class',
                'type':       'nodeEvent',
                'variables':  ['self'],
                'attributes': ['intent(in   )'],
            },
        ],
        'content': "!$GLC attributes unused :: self\nnodeEventSizeOf=0_c_size_t\n",
    }
    _bind(build, 'nodeEvent', function, 'nonStaticSizeOf')


def _bind(build, type_name, function, method_name):
    build.setdefault('types', {}).setdefault(type_name, {}) \
                                  .setdefault('boundFunctions', []) \
                                  .append({
        'type':       'procedure',
        'descriptor': function,
        'name':       method_name,
    })


# ---------------------------------------------------------------------------
# Hook registration.  Order matches Perl NodeEvents.pm:19-34.
# ---------------------------------------------------------------------------

register('nodeEvents', 'types',      Build_Node_Event_Class)
register('nodeEvents', 'interfaces', Node_Event_Task_Interface)
register('nodeEvents', 'interfaces', Node_Event_Merge_Time_Set_Interface)
register('nodeEvents', 'functions',  Node_Event_Non_Static_Size_Of)
register('nodeEvents', 'functions',  Node_Event_Serialize_Raw)
register('nodeEvents', 'functions',  Node_Event_Deserialize_Raw)
register('nodeEvents', 'functions',  Node_Event_Deserialize_Raw_Polymorphic)
