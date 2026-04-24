#!/usr/bin/env python3
# Unit tests for Python utility modules ported from Perl.
# Tests the following modules:
#   - Sort.Topo
#   - List.ExtraUtils
#   - Fortran.Utils.extract_variables
#   - Galacticus.Build.SourceTree.Parse.Declarations.parse_declaration

import sys
import os
import re
import tempfile

sys.path.insert(0, os.path.join(os.environ.get('GALACTICUS_EXEC_PATH', '.'), 'python'))

from Sort.Topo import sort as topo_sort
from List.ExtraUtils import as_array, smart_push, hash_list, sorted_keys
from Fortran.Utils import extract_variables, UNIT_OPENERS
from Galacticus.Build.SourceTree import parse_file, serialize, walk_tree
from Galacticus.Build.SourceTree.Parse.Declarations import (
    parse_declaration, build_declarations, add_declarations, add_attributes,
    get_declaration, declaration_exists,
)
from Galacticus.Build.SourceTree.Parse.ModuleUses   import parse_module_uses, update_uses, add_uses
from Galacticus.Build.SourceTree.Parse.Visibilities import parse_visibilities, update_visibilities
from Galacticus.Build.SourceTree.Parse.OpenMP       import parse_openmp, update as openmp_update, copyin as openmp_copyin
from Galacticus.Build.SourceTree.Parse.Directives   import parse_directives, post_process_directives

# Import Process submodules for side-effect registration; the orchestrator is
# imported via `from ... import ...` after the submodules so they all appear
# in PROCESS_HOOKS when process_tree() is invoked.
import Galacticus.Build.SourceTree.Process.DeepCopyReset      # noqa: F401
import Galacticus.Build.SourceTree.Process.DeepCopyFinalize   # noqa: F401
import Galacticus.Build.SourceTree.Process.OptionalArgument   # noqa: F401
import Galacticus.Build.SourceTree.Process.HDF5FCInterop      # noqa: F401
import Galacticus.Build.SourceTree.Process.Dependencies       # noqa: F401
import Galacticus.Build.SourceTree.Process.Generics                  # noqa: F401
import Galacticus.Build.SourceTree.Process.NonProcessed              # noqa: F401
import Galacticus.Build.SourceTree.Process.Allocate                  # noqa: F401
import Galacticus.Build.SourceTree.Process.ForEach                   # noqa: F401
import Galacticus.Build.SourceTree.Process.SourceIntrospection       # noqa: F401
import Galacticus.Build.SourceTree.Process.ParameterMigration        # noqa: F401
import Galacticus.Build.SourceTree.Process.DebugMPI                  # noqa: F401
import Galacticus.Build.SourceTree.Process.DebugHDF5                 # noqa: F401
import Galacticus.Build.SourceTree.Process.ProfileOpenMP             # noqa: F401
import Galacticus.Build.SourceTree.Process.Constants                 # noqa: F401
import Galacticus.Build.SourceTree.Process.ConditionalCall           # noqa: F401
import Galacticus.Build.SourceTree.Process.InputParametersValidate   # noqa: F401
import Galacticus.Build.SourceTree.Process.StateStore                # noqa: F401
import Galacticus.Build.SourceTree.Process.MetaPropertyDatabase      # noqa: F401
import Galacticus.Build.SourceTree.Process.AddMetaProperty           # noqa: F401
import Galacticus.Build.SourceTree.Process.InputParameter            # noqa: F401
import Galacticus.Build.SourceTree.Process.Constructors              # noqa: F401
import Galacticus.Build.SourceTree.Process.FunctionsGlobal           # noqa: F401
import Galacticus.Build.SourceTree.Process.Enumeration               # noqa: F401
import Galacticus.Build.SourceTree.Process.DeepCopyActions           # noqa: F401
import Galacticus.Build.SourceTree.Process.StateStorable             # noqa: F401
import Galacticus.Build.SourceTree.Process.EventHooksStatic          # noqa: F401
import Galacticus.Build.SourceTree.Process.EventHooks                # noqa: F401
import Galacticus.Build.SourceTree.Process.ThreadSafeIO              # noqa: F401
import Galacticus.Build.SourceTree.Process.SourceDigest              # noqa: F401
import Galacticus.Build.SourceTree.Process.ObjectBuilder             # noqa: F401
import Galacticus.Build.SourceTree.Process.ClassDocumentation        # noqa: F401
from Galacticus.Build.SourceTree.Process.FunctionClass.Utils import (
    latex_breakable, trimlc, striplc, lctrim, strip_variable_name,
    declaration_rank, class_dependencies,
)
from Galacticus.Build.SourceTree.Process.FunctionClass.Descriptor import (
    potential_descriptor_parameters,
)
from Galacticus.Build.SourceTree.Process.FunctionClass.LinkedList import (
    linked_list_register_variable, linked_list_module,
    deep_copy_linked_list, state_store_linked_list,
    allowed_parameters_linked_list, auto_descriptor_linked_list,
    assigner_linked_list,
)
from Galacticus.Build.SourceTree.Process.FunctionClass.DeepCopy import (
    deep_copy_copied_self_block, generate_assignment_allocatable_code,
    deep_copy_declarations,
)
from Galacticus.Build.SourceTree.Process.FunctionClass.StateStore import (
    state_store_explicit_function, generate_allocatable_state_store_code,
    state_store_variables,
)
from Galacticus.Build.SourceTree.Process.FunctionClass import (
    process_function_class as _fc_process_function_class,
    _is_function_class_pointer, _init_code_content,
    _build_base_method_stubs, _build_object_type_method,
    _build_assignment_method, _build_deep_copy_methods,
    _build_state_store_methods, _build_allowed_parameters_method,
    _build_descriptor_methods, _descriptor_discover_class,
    _load_and_sort_classes,
    _generate_type_definition, _generate_constructor,
    _generate_method_functions, _generate_documentation,
    _generate_class_submodules,
    _format_variable_definitions, _source_digest_binding,
    _STATE_STORABLES_HOLDER     as _fc_state_storables_holder,
    _DEEP_COPY_ACTIONS_HOLDER   as _fc_deep_copy_actions_holder,
    _DIRECTIVE_LOCATIONS_HOLDER as _fc_directive_locations_holder,
)
from Galacticus.Build.SourceTree.Process import (
    PROCESS_HOOKS, PROCESS_DEPENDENCIES, POSTPROCESS_HOOKS,
    register_process, process_tree,
)
from Galacticus.Build.Directives import extract_directives, extract_directive
from Galacticus.Build.SourceTree import (
    insert_pre_contains, insert_post_contains, set_visibility,
)

# ============================================================================
# Test framework
# ============================================================================

PASS_COUNT = 0
FAIL_COUNT = 0


def assert_equal(actual, expected, msg):
    global PASS_COUNT, FAIL_COUNT
    if actual == expected:
        PASS_COUNT += 1
        print(f"✓ {msg}")
    else:
        FAIL_COUNT += 1
        print(f"✗ {msg}")
        print(f"  Expected: {expected}")
        print(f"  Got:      {actual}")


def assert_raises(func, exc_type, msg):
    global PASS_COUNT, FAIL_COUNT
    try:
        func()
        FAIL_COUNT += 1
        print(f"✗ {msg} (no exception raised)")
    except exc_type:
        PASS_COUNT += 1
        print(f"✓ {msg}")
    except Exception as e:
        FAIL_COUNT += 1
        print(f"✗ {msg} (wrong exception: {type(e).__name__})")


# ============================================================================
# Sort.Topo tests
# ============================================================================

def test_topo_sort():
    print("\n=== Testing Sort.Topo.sort() ===")

    # Simple DAG
    objects = ['a', 'b', 'c']
    deps = {'b': ['a'], 'c': ['b']}
    result = topo_sort(objects, deps)
    assert_equal(result, ['a', 'b', 'c'], "Simple linear dependency")

    # Independent nodes
    objects = ['x', 'y', 'z']
    deps = {}
    result = topo_sort(objects, deps)
    assert_equal(set(result), {'x', 'y', 'z'}, "Independent nodes")

    # Single node
    objects = ['only']
    deps = {}
    result = topo_sort(objects, deps)
    assert_equal(result, ['only'], "Single node")

    # Multiple dependencies
    objects = ['d', 'b', 'a', 'c']
    deps = {'d': ['b', 'c'], 'b': ['a']}
    result = topo_sort(objects, deps)
    # Check that a comes before b, and b,c come before d
    a_idx, b_idx, c_idx, d_idx = (result.index(x) for x in ['a', 'b', 'c', 'd'])
    assert_equal(a_idx < b_idx < d_idx and c_idx < d_idx, True,
                 "Multiple dependencies respected")

    # Circular dependency
    objects = ['a', 'b']
    deps = {'a': ['b'], 'b': ['a']}
    assert_raises(lambda: topo_sort(objects, deps), RuntimeError,
                  "Circular dependency raises RuntimeError")


# ============================================================================
# List.ExtraUtils tests
# ============================================================================

def test_list_extra_utils():
    print("\n=== Testing List.ExtraUtils ===")

    # as_array: None
    assert_equal(as_array(None), [], "as_array(None)")

    # as_array: scalar
    assert_equal(as_array('hello'), ['hello'], "as_array(scalar)")

    # as_array: list
    assert_equal(as_array(['a', 'b']), ['a', 'b'], "as_array(list)")

    # smart_push
    arr = []
    smart_push(arr, None)
    assert_equal(arr, [], "smart_push with None")

    arr = []
    smart_push(arr, 'x')
    assert_equal(arr, ['x'], "smart_push with scalar")

    arr = []
    smart_push(arr, ['a', 'b'])
    assert_equal(arr, ['a', 'b'], "smart_push with list")

    # hash_list
    d = {'b': {'val': 2}, 'a': {'val': 1}}
    result = hash_list(d)
    assert_equal(len(result), 2, "hash_list length")
    assert_equal([x['val'] for x in result], [1, 2], "hash_list values sorted")

    # hash_list with key_as
    d = {'x': {'id': 10}, 'y': {'id': 20}}
    result = hash_list(d, key_as='name')
    assert_equal(len(result), 2, "hash_list with key_as length")
    assert_equal(result[0]['name'], 'x', "hash_list stamped first key")
    assert_equal(result[1]['name'], 'y', "hash_list stamped second key")

    # sorted_keys
    d = {'z': 1, 'a': 2, 'm': 3}
    assert_equal(sorted_keys(d), ['a', 'm', 'z'], "sorted_keys")

    # sorted_keys with None
    assert_equal(sorted_keys(None), [], "sorted_keys(None)")

    # sorted_keys with empty dict
    assert_equal(sorted_keys({}), [], "sorted_keys({})")


# ============================================================================
# Fortran.Utils.extract_variables tests
# ============================================================================

def test_extract_variables():
    print("\n=== Testing Fortran.Utils.extract_variables() ===")

    # Simple variables
    result = extract_variables('x,y,z')
    assert_equal(result, ['x', 'y', 'z'], "Simple comma-separated list")

    # With dimensions
    result = extract_variables('x(3), y(10,20), z')
    assert_equal(result, ['x', 'y', 'z'], "Variables with dimensions (lowercased, dims removed)")

    # With qualifiers preserved
    result = extract_variables('arr(3), mat(10,20)', keep_qualifiers=True)
    assert_equal(len(result), 2, "Count with qualifiers")
    assert_equal(result[0], 'arr(3)', "First variable with dimension")
    assert_equal(result[1], 'mat(10,20)', "Second variable with dimension")

    # Upper/lower case
    result = extract_variables('MyVar,ANOTHER', lower_case=True)
    assert_equal(result, ['myvar', 'another'], "Lowercase conversion")

    result = extract_variables('MyVar,ANOTHER', lower_case=False)
    assert_equal(result, ['MyVar', 'ANOTHER'], "Original case preservation")

    # With initialization (should be stripped by default)
    result = extract_variables('x=1.0, y=2.0', keep_qualifiers=False)
    assert_equal(result, ['x', 'y'], "Initialization stripped")

    # None input
    result = extract_variables(None)
    assert_equal(result, [], "None input returns empty list")

    # Empty string
    result = extract_variables('')
    assert_equal(result, [], "Empty string returns empty list")


# ============================================================================
# Declarations.parse_declaration tests
# ============================================================================

def test_parse_declaration():
    print("\n=== Testing Declarations.parse_declaration() ===")

    # Integer declaration
    decl = parse_declaration('integer :: i')
    assert_equal(decl is not None, True, "Integer declaration recognized")
    if decl:
        assert_equal(decl['intrinsic'], 'integer', "Integer intrinsic type")
        assert_equal(decl['variables'], ['i'], "Integer variable name")

    # Real declaration
    decl = parse_declaration('real :: x, y')
    assert_equal(decl is not None, True, "Real declaration recognized")
    if decl:
        assert_equal(decl['intrinsic'], 'real', "Real intrinsic type")
        assert_equal(set(decl['variables']), {'x', 'y'}, "Real variables")

    # Double precision
    decl = parse_declaration('double precision :: pi')
    assert_equal(decl is not None, True, "Double precision recognized")
    if decl:
        assert_equal(decl['intrinsic'], 'double precision', "Double precision intrinsic")

    # Logical
    decl = parse_declaration('logical :: flag')
    assert_equal(decl is not None, True, "Logical declaration recognized")
    if decl:
        assert_equal(decl['intrinsic'], 'logical', "Logical intrinsic type")

    # With intent attribute
    decl = parse_declaration('real, intent(in) :: param')
    assert_equal(decl is not None, True, "Declaration with intent recognized")
    if decl:
        assert_equal('intent(in)' in decl['attributes'], True, "Intent attribute captured")

    # With dimension
    decl = parse_declaration('integer, dimension(:) :: arr')
    assert_equal(decl is not None, True, "Declaration with dimension recognized")
    if decl:
        assert_equal(decl['intrinsic'], 'integer', "Integer with dimension")

    # Type/class declaration
    decl = parse_declaration('type(myType) :: obj')
    assert_equal(decl is not None, True, "Type declaration recognized")
    if decl:
        assert_equal(decl['intrinsic'], 'type', "Type intrinsic")

    # Non-declaration line should return None
    decl = parse_declaration('x = y + z')
    assert_equal(decl, None, "Non-declaration returns None")

    # Comment line
    decl = parse_declaration('! this is a comment')
    assert_equal(decl, None, "Comment line returns None")

# ============================================================================
# Fortran.Utils.UNIT_OPENERS
# ============================================================================

def test_unit_openers():
    print("\n=== Testing Fortran.Utils.UNIT_OPENERS ===")

    # Interface declarations
    m = re.match(UNIT_OPENERS['interface']['regex'],'interface myInterface')
    assert_equal(m is not None, True, "Interface matched")
    if m:
        groups = m.groups()
        assert_equal(groups[UNIT_OPENERS['interface']['unit_name']], 'myInterface', "Interface name correctly extracted")

    m = re.match(UNIT_OPENERS['interface']['regex'],'abstract interface myInterface_34')
    assert_equal(m is not None, True, "Abstract interface matched")
    if m:
        groups = m.groups()
        assert_equal(groups[UNIT_OPENERS['interface']['unit_name']], 'myInterface_34', "Abstract interface name correctly extracted")

    # Function declarations
    m = re.match(UNIT_OPENERS['function']['regex'],'double precision function myFunc(a,b,c) result(d)')
    assert_equal(m is not None, True, "Function opener matched")
    if m:
        groups = m.groups()
        assert_equal(groups[UNIT_OPENERS['function']['unit_name']], 'myFunc', "Function name correctly extracted")

    # Derived type declarations
    m = re.match(UNIT_OPENERS['type']['regex'],'type :: myType')
    assert_equal(m is not None, True, "Simple derived type declaration matched")
    if m:
        groups = m.groups()
        assert_equal(groups[UNIT_OPENERS['type']['unit_name']], 'myType', "Type name correctly extracted")

    m = re.match(UNIT_OPENERS['type']['regex'],'type, abstract :: myType1')
    assert_equal(m is not None, True, "Abstract derived type declaration matched")
    if m:
        groups = m.groups()
        assert_equal(groups[UNIT_OPENERS['type']['unit_name']], 'myType1', "Type name correctly extracted")

    m = re.match(UNIT_OPENERS['type']['regex'],'type, extends(parentType) :: myType_2')
    assert_equal(m is not None, True, "Extended derived type declaration matched")
    if m:
        groups = m.groups()
        assert_equal(groups[UNIT_OPENERS['type']['unit_name']], 'myType_2', "Type name correctly extracted")


# ============================================================================
# Helper: parse inline Fortran text via a tempfile
# ============================================================================

def _parse_text(text):
    """Write text to a temp .F90 file, parse it, and return the root node."""
    fh = tempfile.NamedTemporaryFile(
        mode='w', suffix='.F90', delete=False, encoding='utf-8')
    try:
        fh.write(text)
        fh.close()
        return parse_file(fh.name)
    finally:
        os.unlink(fh.name)


def _find_node(root, node_type):
    """Return the first descendant of root with node['type'] == node_type."""
    for n in walk_tree(root):
        if n.get('type') == node_type:
            return n
    return None


def _find_nodes(root, node_type):
    return [n for n in walk_tree(root) if n.get('type') == node_type]


# ============================================================================
# Parse.ModuleUses tests
# ============================================================================

def test_parse_module_uses():
    print("\n=== Testing Parse.ModuleUses.parse_module_uses() ===")

    text = (
        "module m\n"
        "use foo\n"
        "use, intrinsic :: iso_c_binding\n"
        "!$ use omp_lib\n"
        "use bar, only : baz, qux\n"
        "x = 1\n"
        "end module m\n"
    )
    root = _parse_text(text)
    uses = _find_node(root, 'moduleUse')
    assert_equal(uses is not None, True, "moduleUse node created")
    if uses:
        order = uses.get('moduleOrder', [])
        assert_equal(order[0], 'iso_c_binding', "intrinsic module placed first")
        assert_equal(set(order), {'foo', 'iso_c_binding', 'omp_lib', 'bar'}, "all modules captured")
        mu = uses.get('moduleUse', {})
        assert_equal(mu['iso_c_binding']['intrinsic'], True, "intrinsic flag set")
        assert_equal(mu['omp_lib']['openMP'], True, "openMP flag set for !$ use")
        assert_equal(set(mu['bar'].get('only', {}).keys()), {'baz', 'qux'}, "only symbols captured")

    # Round-trip: serializing the original text should reproduce it.
    assert_equal(serialize(root), text, "moduleUses round-trip preserves source")


def test_parse_module_uses_preprocessor():
    print("\n=== Testing Parse.ModuleUses preprocessor conditions ===")

    text = (
        "module m\n"
        "#ifdef GUARDED\n"
        "use conditional_mod\n"
        "#endif\n"
        "x = 1\n"
        "end module m\n"
    )
    root = _parse_text(text)
    uses = _find_node(root, 'moduleUse')
    assert_equal(uses is not None, True, "moduleUse node created under preprocessor")
    if uses:
        conds = uses['moduleUse']['conditional_mod'].get('conditions')
        assert_equal(isinstance(conds, list) and len(conds) == 1, True,
                     "conditions list captured")
        if conds:
            assert_equal(conds[0]['name'], 'GUARDED', "condition name captured")
            assert_equal(conds[0]['invert'], False, "ifdef not inverted")


def test_update_uses():
    print("\n=== Testing Parse.ModuleUses.update_uses() ===")

    uses_node = {
        'type':        'moduleUse',
        'moduleUse':   {
            'iso_c_binding': {'openMP': False, 'intrinsic': True,  'all': True},
            'foo':           {'openMP': False, 'intrinsic': False, 'only': {'bar': True, 'baz': True}},
        },
        'moduleOrder': ['iso_c_binding', 'foo'],
        'firstChild':  {'type': 'code', 'content': '  use foo\n',
                        'parent': None, 'sibling': None, 'firstChild': None,
                        'source': 'x', 'line': 0},
        'source': 'x', 'line': 0, 'parent': None, 'sibling': None,
    }
    uses_node['firstChild']['parent'] = uses_node

    update_uses(uses_node)
    out = uses_node['firstChild']['content']
    assert_equal('use, intrinsic :: iso_c_binding' in out, True,
                 "intrinsic use line regenerated")
    assert_equal('use            :: foo' in out, True,
                 "non-intrinsic use padded to column")
    assert_equal(', only :' in out, True, "only clause emitted")
    assert_equal('bar' in out and 'baz' in out, True, "only symbols present")


def test_add_uses():
    print("\n=== Testing Parse.ModuleUses.add_uses() ===")

    root = _parse_text(
        "module m\n"
        "use foo\n"
        "x = 1\n"
        "end module m\n"
    )
    module_node = _find_node(root, 'module')
    new = {
        'moduleUse':   {'added_mod': {'openMP': False, 'intrinsic': False, 'all': True}},
        'moduleOrder': ['added_mod'],
    }
    add_uses(module_node, new)
    uses = _find_node(module_node, 'moduleUse')
    assert_equal('added_mod' in uses['moduleOrder'], True,
                 "new module appended to moduleOrder")
    assert_equal('added_mod' in uses['moduleUse'], True,
                 "new module in moduleUse dict")
    assert_equal('use :: added_mod' in serialize(root), True,
                 "new use line appears in serialized output")


# ============================================================================
# Parse.Declarations pass + helpers tests
# ============================================================================

def test_parse_declarations_pass():
    print("\n=== Testing declarations pass (via parse_file) ===")

    text = (
        "subroutine s\n"
        "implicit none\n"
        "integer :: i\n"
        "real, intent(in) :: x, y\n"
        "i = 1\n"
        "end subroutine s\n"
    )
    root = _parse_text(text)
    decl = _find_node(root, 'declaration')
    assert_equal(decl is not None, True, "declaration node created")
    if decl:
        assert_equal(decl.get('implicitNone'), True, "implicitNone captured")
        decls = decl.get('declarations', [])
        intrinsics = [d['intrinsic'] for d in decls]
        assert_equal('integer' in intrinsics and 'real' in intrinsics, True,
                     "both integer and real declarations parsed")
    assert_equal(serialize(root), text, "declaration round-trip preserves source")


def test_build_declarations():
    print("\n=== Testing Parse.Declarations.build_declarations() ===")

    decl_node = {
        'type':         'declaration',
        'implicitNone': True,
        'declarations': [
            {'intrinsic': 'integer', 'type': None, 'openMP': False,
             'attributes': ['intent(in)'], 'variables': ['i'], 'variableNames': ['i']},
            {'intrinsic': 'real', 'type': 'kind=8', 'openMP': True,
             'attributes': [], 'variables': ['x'], 'variableNames': ['x'],
             'threadprivate': True},
            {'intrinsic': 'integer', 'type': None, 'openMP': False,
             'attributes': [], 'variables': ['g'], 'variableNames': ['g'],
             'preprocessor': 'FOO'},
        ],
        'firstChild':  None,
        'source': 'x', 'line': 0, 'parent': None, 'sibling': None,
    }
    build_declarations(decl_node)
    content = decl_node['firstChild']['content']
    assert_equal(content.startswith('implicit none\n'), True, "implicit none emitted first")
    assert_equal('integer, intent(in) :: i' in content, True, "attributes rendered")
    assert_equal('!$ real(kind=8) :: x' in content, True, "openMP + type parenthesized")
    assert_equal('!$omp threadprivate(x)' in content, True, "threadprivate footer emitted")
    assert_equal('#ifdef FOO' in content and '#endif' in content, True,
                 "preprocessor guards emitted")


def test_add_declarations():
    print("\n=== Testing Parse.Declarations.add_declarations() ===")

    root = _parse_text(
        "subroutine s\n"
        "use foo\n"
        "integer :: existing\n"
        "existing = 0\n"
        "end subroutine s\n"
    )
    sub = _find_node(root, 'subroutine')
    new_decl = {
        'intrinsic': 'logical', 'type': None, 'openMP': False,
        'attributes': [], 'variables': ['flag'], 'variableNames': ['flag'],
    }
    add_declarations(sub, [new_decl])
    assert_equal(declaration_exists(sub, 'flag'), True, "new declaration visible via declaration_exists")
    out = serialize(root)
    assert_equal('logical :: flag' in out, True, "new logical declaration appears in output")


def test_add_attributes():
    print("\n=== Testing Parse.Declarations.add_attributes() ===")

    root = _parse_text(
        "subroutine s\n"
        "real :: a, b, c\n"
        "a = 1.0\n"
        "end subroutine s\n"
    )
    sub = _find_node(root, 'subroutine')
    add_attributes(sub, 'b', ['allocatable'])
    out = serialize(root)
    # 'b' should now have allocatable; a and c should remain in a non-allocatable declaration.
    assert_equal('real, allocatable :: b' in out, True,
                 "target variable 'b' split out with new attribute")
    assert_equal('real :: a, c' in out, True,
                 "other variables remain in a separate declaration")


def test_get_declaration_and_declaration_exists():
    print("\n=== Testing get_declaration() / declaration_exists() ===")

    root = _parse_text(
        "subroutine s\n"
        "integer :: i, j=5\n"
        "real :: x\n"
        "i = 0\n"
        "end subroutine s\n"
    )
    sub = _find_node(root, 'subroutine')
    assert_equal(declaration_exists(sub, 'i'), True,  "declaration_exists: i present")
    assert_equal(declaration_exists(sub, 'I'), True,  "declaration_exists is case-insensitive")
    assert_equal(declaration_exists(sub, 'k'), False, "declaration_exists: k absent")

    d = get_declaration(sub, 'x')
    assert_equal(d['intrinsic'], 'real',  "get_declaration returns correct intrinsic")
    assert_equal(d['variables'], ['x'],   "get_declaration reduces variables to target")

    # Initializer stripping: j has `=5` — still findable.
    d = get_declaration(sub, 'j')
    assert_equal(d['intrinsic'], 'integer', "get_declaration strips initializer")

    assert_raises(lambda: get_declaration(sub, 'missing'), RuntimeError,
                  "get_declaration raises on missing variable")


# ============================================================================
# Parse.Visibilities tests
# ============================================================================

def test_parse_visibilities():
    print("\n=== Testing Parse.Visibilities.parse_visibilities() ===")

    text = (
        "module m\n"
        "public :: alpha, beta\n"
        "private :: gamma\n"
        "x = 1\n"
        "end module m\n"
    )
    root = _parse_text(text)
    vis = _find_node(root, 'visibility')
    assert_equal(vis is not None, True, "visibility node created")
    if vis:
        v = vis['visibility']
        assert_equal(set(v.get('public',  {}).keys()), {'alpha', 'beta'}, "public symbols captured")
        assert_equal(set(v.get('private', {}).keys()), {'gamma'},         "private symbols captured")
    assert_equal(serialize(root), text, "visibility round-trip preserves source")


def test_update_visibilities():
    print("\n=== Testing Parse.Visibilities.update_visibilities() ===")

    vis_node = {
        'type':       'visibility',
        'visibility': {'public': {'beta': True, 'alpha': True}, 'private': {'z': True}},
        'firstChild': None,
        'source': 'x', 'line': 0, 'parent': None, 'sibling': None,
    }
    update_visibilities(vis_node)
    content = vis_node['firstChild']['content']
    assert_equal('! Galacticus::Build::SourceTree::Parse::Visibilities(): updated' in content,
                 True, "comment header emitted")
    assert_equal('public :: alpha, beta' in content,  True, "public list sorted")
    assert_equal('private :: z' in content,           True, "private list emitted")


# ============================================================================
# Parse.OpenMP tests
# ============================================================================

def test_parse_openmp():
    print("\n=== Testing Parse.OpenMP.parse_openmp() ===")

    text = (
        "subroutine s\n"
        "!$omp parallel do schedule(static) private(i)\n"
        "do i=1,10\n"
        "end do\n"
        "!$omp end parallel do\n"
        "end subroutine s\n"
    )
    root = _parse_text(text)
    omp_nodes = _find_nodes(root, 'openMP')
    assert_equal(len(omp_nodes), 4,
                 "four openMP nodes (parallel + do + end do + end parallel)")
    names = [n['name'] for n in omp_nodes]
    assert_equal('parallel' in names and 'do' in names, True,
                 "parallel and composite 'do' both present")

    # schedule should attach to the composite 'do' node, private to parallel.
    for n in omp_nodes:
        opt_names = [o['name'] for o in n.get('options', [])]
        if n['name'] == 'do' and not n.get('isCloser'):
            assert_equal('schedule' in opt_names, True,
                         "schedule attached to composite 'do'")
        if n['name'] == 'parallel' and not n.get('isCloser'):
            assert_equal('private' in opt_names, True,
                         "private attached to 'parallel'")


def test_openmp_update():
    print("\n=== Testing Parse.OpenMP.update() ===")

    node = {
        'type':       'openMP',
        'name':       'parallel',
        'isCloser':   False,
        'options':    [{'name': 'private', 'content': '(i,j)'}],
        'firstChild': None,
        'source': 'x', 'line': 0, 'parent': None, 'sibling': None,
    }
    openmp_update(node)
    assert_equal(node['firstChild']['content'], "!$omp parallel private(i,j)\n",
                 "update() reconstructs !$omp line")

    node['isCloser'] = True
    node['options'] = []
    openmp_update(node)
    assert_equal(node['firstChild']['content'], "!$omp end parallel\n",
                 "update() handles isCloser with no options")


def test_openmp_copyin():
    print("\n=== Testing Parse.OpenMP.copyin() ===")

    node = {
        'type':       'openMP',
        'name':       'parallel',
        'isCloser':   False,
        'options':    [],
        'firstChild': None,
        'source': 'x', 'line': 0, 'parent': None, 'sibling': None,
    }
    openmp_copyin(node, ['a', 'b'])
    openmp_copyin(node, ['b', 'c'])  # 'b' should not duplicate

    copyin_opts = [o for o in node['options'] if o['name'] == 'copyin']
    assert_equal(len(copyin_opts), 1, "single copyin option after repeated calls")
    assert_equal(copyin_opts[0]['content'], '(a,b,c)',
                 "copyin content deduplicated and comma-joined")
    assert_equal(node['firstChild']['content'], "!$omp parallel copyin(a,b,c)\n",
                 "code content regenerated by copyin()")


# ============================================================================
# Parse.Directives tests
# ============================================================================

def test_parse_directives():
    print("\n=== Testing Parse.Directives.parse_directives() ===")

    text = (
        "subroutine s\n"
        "  !![\n"
        "  <foo bar=\"baz\"/>\n"
        "  !!]\n"
        "x = 1\n"
        "end subroutine s\n"
    )
    root = _parse_text(text)
    foo = _find_node(root, 'foo')
    assert_equal(foo is not None, True, "directive node with directive name as type")
    if foo:
        assert_equal(foo['directive'].get('bar'), 'baz',
                     "directive attribute parsed into dict")
    assert_equal(serialize(root), text, "directive round-trip preserves raw text")


def test_post_process_directives():
    print("\n=== Testing Parse.Directives.post_process_directives() ===")

    text = (
        "subroutine s\n"
        "  !![\n"
        "  <foo bar=\"baz\"/>\n"
        "  !!]\n"
        "end subroutine s\n"
    )
    root = _parse_text(text)
    foo = _find_node(root, 'foo')
    assert_equal(foo is not None, True, "directive node located")

    # Unprocessed → should raise.
    assert_raises(lambda: post_process_directives(root), RuntimeError,
                  "post_process_directives raises on unprocessed directive")

    # Mark processed → should pass silently.
    foo['directive']['processed'] = True
    raised = False
    try:
        post_process_directives(root)
    except Exception:
        raised = True
    assert_equal(raised, False,
                 "post_process_directives passes when all directives are processed")


# ============================================================================
# Process.* tests
# ============================================================================

def _make_directive_node(directive_type, directive, parent):
    """Build a minimal directive node and link it as `parent`'s only child."""
    node = {
        'type':       directive_type,
        'directive':  dict(directive),
        'parent':     parent,
        'firstChild': None,
        'sibling':    None,
        'source':     'test',
        'line':       1,
    }
    node['firstChild'] = {
        'type':       'code',
        'content':    '',
        'parent':     node,
        'firstChild': None,
        'sibling':    None,
        'source':     'test',
        'line':       1,
    }
    # Link into parent's child chain as the first child.
    existing = parent.get('firstChild')
    parent['firstChild'] = node
    node['sibling'] = existing
    return node


def test_process_tree_orchestrator():
    print("\n=== Testing Process.process_tree() orchestrator ===")

    # Snapshot current registrations and install two throw-away hooks with
    # a dependency edge, then verify execution order respects it.
    saved_hooks = dict(PROCESS_HOOKS)
    saved_deps  = dict(PROCESS_DEPENDENCIES)
    saved_post  = dict(POSTPROCESS_HOOKS)
    PROCESS_HOOKS.clear()
    PROCESS_DEPENDENCIES.clear()
    POSTPROCESS_HOOKS.clear()
    try:
        call_order = []
        def hook_a(tree, options): call_order.append('a')
        def hook_b(tree, options): call_order.append('b')
        # Perl convention: `before=['b']` means b runs AFTER a.
        register_process('a', hook_a, before=['b'])
        register_process('b', hook_b)
        dummy_tree = {'type': 'file', 'firstChild': None, 'parent': None, 'sibling': None}
        process_tree(dummy_tree)
        assert_equal(call_order, ['a', 'b'],
                     "process_tree respects `before` dependency")
    finally:
        PROCESS_HOOKS.clear();        PROCESS_HOOKS.update(saved_hooks)
        PROCESS_DEPENDENCIES.clear(); PROCESS_DEPENDENCIES.update(saved_deps)
        POSTPROCESS_HOOKS.clear();    POSTPROCESS_HOOKS.update(saved_post)


def test_process_deep_copy_reset():
    print("\n=== Testing Process.DeepCopyReset ===")

    root = {'type': 'file', 'firstChild': None, 'parent': None, 'sibling': None,
            'source': 'test', 'line': 0}
    sub = {'type': 'subroutine', 'name': 's', 'parent': root,
           'firstChild': None, 'sibling': None, 'source': 'test', 'line': 0}
    root['firstChild'] = sub

    directive = _make_directive_node(
        'deepCopyReset', {'variables': 'objA objB', 'processed': False}, sub)

    from Galacticus.Build.SourceTree.Process.DeepCopyReset import process_deep_copy_reset
    process_deep_copy_reset(root, {})

    assert_equal(directive['directive']['processed'], True,
                 "deepCopyReset directive marked processed")
    emitted = directive.get('sibling')
    assert_equal(emitted is not None and emitted.get('type') == 'code', True,
                 "code node inserted after deepCopyReset directive")
    assert_equal(emitted['content'],
                 "call objA%deepCopyReset()\ncall objB%deepCopyReset()\n",
                 "emitted code calls deepCopyReset on each named object")


def test_process_deep_copy_finalize():
    print("\n=== Testing Process.DeepCopyFinalize ===")

    root = {'type': 'file', 'firstChild': None, 'parent': None, 'sibling': None,
            'source': 'test', 'line': 0}
    sub = {'type': 'subroutine', 'name': 's', 'parent': root,
           'firstChild': None, 'sibling': None, 'source': 'test', 'line': 0}
    root['firstChild'] = sub

    directive = _make_directive_node(
        'deepCopyFinalize', {'variables': 'only_one', 'processed': False}, sub)

    from Galacticus.Build.SourceTree.Process.DeepCopyFinalize import process_deep_copy_finalize
    process_deep_copy_finalize(root, {})

    assert_equal(directive['directive']['processed'], True,
                 "deepCopyFinalize directive marked processed")
    emitted = directive.get('sibling')
    assert_equal(emitted['content'], "call only_one%deepCopyFinalize()\n",
                 "emitted code calls deepCopyFinalize on the named object")


def test_process_optional_argument():
    print("\n=== Testing Process.OptionalArgument ===")

    # Parse a subroutine with a real, optional, intent(in) declaration,
    # then hand-inject an optionalArgument directive as the first child.
    root = _parse_text(
        "subroutine s(foo)\n"
        "real, optional, intent(in) :: foo\n"
        "foo = 0.0\n"
        "end subroutine s\n"
    )
    sub = _find_node(root, 'subroutine')
    directive = _make_directive_node(
        'optionalArgument',
        {'name': 'foo', 'defaultsTo': '1.5', 'processed': False},
        sub,
    )

    from Galacticus.Build.SourceTree.Process.OptionalArgument import process_optional_arguments
    process_optional_arguments(root, {})

    assert_equal(directive['directive']['processed'], True,
                 "optionalArgument directive marked processed")
    assert_equal(declaration_exists(sub, 'foo_'), True,
                 "new `_` declaration added for the optional argument")
    # The new declaration should have neither `optional` nor any `intent(...)` attribute.
    new_decl = get_declaration(sub, 'foo_')
    assert_equal(any(a == 'optional' or a.startswith('intent(')
                     for a in new_decl.get('attributes') or []), False,
                 "new declaration strips `optional` and `intent(...)` attributes")

    out = serialize(root)
    assert_equal("foo_=1.5" in out, True,
                 "setter assigns defaultsTo value to the underscore variable")
    assert_equal("if (present(foo)) foo_=foo" in out, True,
                 "setter copies caller-supplied value when present")


def test_process_hdf5_fc_interop():
    print("\n=== Testing Process.HDF5FCInterop ===")

    tmpdir = tempfile.mkdtemp()
    try:
        # Minimal hdf5FCInterop.dat mapping one HDF5 kind to a C kind.
        with open(os.path.join(tmpdir, 'hdf5FCInterop.dat'), 'w') as fh:
            fh.write("h5s_t = c_size_t\n")
        saved_build = os.environ.get('BUILDPATH')
        os.environ['BUILDPATH'] = tmpdir
        try:
            root = _parse_text(
                "module m\n"
                "use foo\n"
                "integer(kind=h5s_t) :: n\n"
                "end module m\n"
            )
            from Galacticus.Build.SourceTree.Process.HDF5FCInterop import process_hdf5_fc_interop
            process_hdf5_fc_interop(root, {})

            out = serialize(root)
            assert_equal('c_size_t' in out,
                         True, "declaration kind rewritten to C-interop kind")
            assert_equal('h5s_t' not in out,
                         True, "original HDF5 kind replaced")
            assert_equal('ISO_C_Binding' in out,
                         True, "ISO_C_Binding use statement injected")
            assert_equal('c_size_t' in out.split('only :')[-1] if 'only :' in out else False,
                         True, "imported symbol listed in `only :` clause")
        finally:
            if saved_build is None:
                del os.environ['BUILDPATH']
            else:
                os.environ['BUILDPATH'] = saved_build
    finally:
        import shutil
        shutil.rmtree(tmpdir)


def test_process_dependencies():
    print("\n=== Testing Process.Dependencies ===")

    tmpdir = tempfile.mkdtemp()
    try:
        aux = os.path.join(tmpdir, 'aux')
        os.makedirs(aux)
        with open(os.path.join(aux, 'dependencies.yml'), 'w') as fh:
            fh.write("libfoo: 1.2.3\n")
            fh.write("libbar: 4.5\n")

        saved_exec = os.environ.get('GALACTICUS_EXEC_PATH')
        os.environ['GALACTICUS_EXEC_PATH'] = tmpdir
        try:
            root = {'type': 'file', 'firstChild': None, 'parent': None,
                    'sibling': None, 'source': 'test', 'line': 0}
            directive = _make_directive_node(
                'dependenciesInitialize', {'processed': False}, root)

            from Galacticus.Build.SourceTree.Process.Dependencies import process_dependencies
            process_dependencies(root, {})

            assert_equal(directive['directive']['processed'], True,
                         "dependenciesInitialize directive marked processed")
            emitted = directive.get('sibling')
            assert_equal(emitted is not None and emitted.get('type') == 'code', True,
                         "code node inserted after directive")
            lines = emitted['content'].strip().split('\n')
            assert_equal(len(lines), 2, "one line per dependency")
            # Dependencies are emitted in sorted order.
            assert_equal("call dependencies_%set(var_str('libbar'),var_str('4.5'))"     in lines[0], True,
                         "first line initializes alphabetically-first dependency")
            assert_equal("call dependencies_%set(var_str('libfoo'),var_str('1.2.3'))"   in lines[1], True,
                         "second line initializes alphabetically-second dependency")
        finally:
            if saved_exec is None:
                del os.environ['GALACTICUS_EXEC_PATH']
            else:
                os.environ['GALACTICUS_EXEC_PATH'] = saved_exec
    finally:
        import shutil
        shutil.rmtree(tmpdir)


# ============================================================================
# Process.Generics tests (and its Generics-dependent siblings)
# ============================================================================

def test_process_generics_subtree_expansion():
    print("\n=== Testing Process.Generics subtree expansion ===")

    src = (
        "module m\n"
        "  !![\n"
        "  <generic identifier=\"cType\">\n"
        "   <instance label=\"Double\" intrinsic=\"real(kind=c_double)\" />\n"
        "   <instance label=\"Long\"   intrinsic=\"integer(kind=c_long)\" />\n"
        "  </generic>\n"
        "  !!]\n"
        "  subroutine worker_{cType¦label}(x)\n"
        "    {cType¦intrinsic}, intent(in) :: x\n"
        "    print *, x, \" as {cType¦label}\"\n"
        "  end subroutine worker_{cType¦label}\n"
        "end module m\n"
    )
    root = _parse_text(src)

    from Galacticus.Build.SourceTree.Process.Generics import process_generics
    process_generics(root, {})

    out = serialize(root)
    assert_equal('subroutine worker_Double' in out, True,
                 "Double-labeled subroutine generated")
    assert_equal('subroutine worker_Long'   in out, True,
                 "Long-labeled subroutine generated")
    assert_equal('{cType¦label}' not in out, True,
                 "{cType¦label} placeholders fully substituted")
    assert_equal('real(kind=c_double), intent(in) :: x'    in out, True,
                 "Double instance's `intrinsic` inlined into the type")
    assert_equal('integer(kind=c_long), intent(in) :: x'   in out, True,
                 "Long instance's `intrinsic` inlined into the type")
    assert_equal('\" as Double\"' in out and '\" as Long\"' in out, True,
                 "string-literal placeholders substituted per instance")


def test_process_generics_regex_form():
    print("\n=== Testing Process.Generics regEx¦FROM¦TO¦ form ===")

    src = (
        "module m\n"
        "  !![\n"
        "  <generic identifier=\"Type\">\n"
        "   <instance label=\"Logical\" outputConverter=\"regEx¦(.*)¦char($1)¦\"/>\n"
        "   <instance label=\"Integer\" outputConverter=\"regEx¦(.*)¦$1¦\"/>\n"
        "  </generic>\n"
        "  !!]\n"
        "  subroutine s\n"
        "    print *, {Type¦outputConverter¦myFlag}\n"
        "  end subroutine s\n"
        "end module m\n"
    )
    root = _parse_text(src)
    from Galacticus.Build.SourceTree.Process.Generics import process_generics
    process_generics(root, {})
    out = serialize(root)
    assert_equal('print *, char(myFlag)' in out, True,
                 "Logical instance wraps the tag body in char(...)")
    assert_equal('print *, myFlag' in out,       True,
                 "Integer instance emits the tag body verbatim")


def test_process_generics_conditional():
    print("\n=== Testing Process.Generics match-conditional form ===")

    src = (
        "module m\n"
        "  !![\n"
        "  <generic identifier=\"Type\">\n"
        "   <instance label=\"Array\"  />\n"
        "   <instance label=\"Scalar\" />\n"
        "  </generic>\n"
        "  !!]\n"
        "  subroutine s_{Type¦label}\n"
        "    print *, {Type¦match¦Array¦is array¦is scalar}\n"
        "  end subroutine s_{Type¦label}\n"
        "end module m\n"
    )
    root = _parse_text(src)
    from Galacticus.Build.SourceTree.Process.Generics import process_generics
    process_generics(root, {})
    out = serialize(root)
    assert_equal('is array'  in out, True, "Array instance picked MATCH branch")
    assert_equal('is scalar' in out, True, "Scalar instance picked NO-MATCH branch")


def test_process_non_processed():
    print("\n=== Testing Process.NonProcessed ===")

    root = {'type': 'file', 'firstChild': None, 'parent': None, 'sibling': None,
            'source': 'test', 'line': 0}
    # Build a few directive-like nodes directly — NonProcessed walks the raw tree.
    d1 = _make_directive_node('methods',        {'processed': False}, root)
    d2 = _make_directive_node('someCustomTask', {'processed': False}, root)
    d3 = _make_directive_node('notListed',      {'processed': False}, root)

    from Galacticus.Build.SourceTree.Process.NonProcessed import process_non_processed
    process_non_processed(root, {})

    assert_equal(d1['directive']['processed'], True,
                 "listed directive ('methods') marked processed")
    assert_equal(d2['directive']['processed'], True,
                 "*Task-suffixed directive marked processed")
    assert_equal(d3['directive']['processed'], False,
                 "unrelated directive left alone")


def test_process_allocate():
    print("\n=== Testing Process.Allocate ===")

    # Rank-2 array: declaration provides `dimension(:,:)`, directive uses `shape`.
    root = _parse_text(
        "subroutine s\n"
        "real, allocatable, dimension(:,:) :: arr\n"
        "end subroutine s\n"
    )
    sub = _find_node(root, 'subroutine')
    directive = _make_directive_node(
        'allocate', {'variable': 'arr', 'shape': 'myShape', 'processed': False}, sub)

    from Galacticus.Build.SourceTree.Process.Allocate import process_allocate
    process_allocate(root, {})

    out = serialize(root)
    assert_equal('allocate(arr(myShape(1),myShape(2)))' in out, True,
                 "rank-2 allocation with shape() references emitted")
    assert_equal(directive['directive']['processed'], True,
                 "allocate directive marked processed")

    # Rank-1 via `size` directive attribute.
    root = _parse_text(
        "subroutine s\n"
        "real, allocatable, dimension(:) :: v\n"
        "end subroutine s\n"
    )
    sub = _find_node(root, 'subroutine')
    _make_directive_node(
        'allocate', {'variable': 'v', 'size': 'src', 'processed': False}, sub)
    process_allocate(root, {})
    assert_equal('allocate(v(size(src,dim=1)))' in serialize(root), True,
                 "rank-1 allocation via size(src,dim=i) emitted")

    # Scalar (rank=0): no indices emitted.
    root = _parse_text(
        "subroutine s\n"
        "integer :: scalar_var\n"
        "end subroutine s\n"
    )
    sub = _find_node(root, 'subroutine')
    _make_directive_node(
        'allocate', {'variable': 'scalar_var', 'shape': 'anything', 'processed': False}, sub)
    process_allocate(root, {})
    assert_equal('allocate(scalar_var)' in serialize(root), True,
                 "rank-0 allocation has no index list")


def test_process_for_each_rank1():
    print("\n=== Testing Process.ForEach (rank-1) ===")

    root = _parse_text(
        "subroutine s\n"
        "real, dimension(:) :: v\n"
        "end subroutine s\n"
    )
    sub = _find_node(root, 'subroutine')
    directive = _make_directive_node(
        'forEach',
        {'variable': 'v', 'content': '  v{index}=0.0\n', 'processed': False},
        sub,
    )

    from Galacticus.Build.SourceTree.Process.ForEach import process_for_each
    process_for_each(root, {})

    assert_equal(directive['directive']['processed'], True,
                 "forEach directive marked processed")
    out = serialize(root)
    assert_equal('do foreach__1=1,size(v,dim=1)' in out, True,
                 "do loop generated with size(v,dim=1)")
    assert_equal('v(foreach__1)=0.0' in out, True,
                 "{index} expanded to parenthesized index list")
    assert_equal('end do' in out, True, "loop closer emitted")
    assert_equal(declaration_exists(sub, 'foreach__1'), True,
                 "foreach__1 index variable declared")


def test_process_for_each_rank0():
    print("\n=== Testing Process.ForEach (scalar / rank-0) ===")

    root = _parse_text(
        "subroutine s\n"
        "real :: scalar_var\n"
        "end subroutine s\n"
    )
    sub = _find_node(root, 'subroutine')
    _make_directive_node(
        'forEach',
        {'variable': 'scalar_var',
         'content': "  write(*,'(a1)') '.'\n",
         'processed': False},
        sub,
    )
    from Galacticus.Build.SourceTree.Process.ForEach import process_for_each
    process_for_each(root, {})
    out = serialize(root)
    # Rank-0 should emit no `do`/`end do` statements.
    assert_equal('do foreach__' in out, False, "no do loop for rank-0 variable")
    assert_equal(declaration_exists(sub, 'foreach__1'), False,
                 "no index variable declared for rank-0")


# ============================================================================
# Process.SourceIntrospection tests
# ============================================================================

def test_source_introspection_instrument():
    print("\n=== Testing Process.SourceIntrospection.instrument() ===")

    from Galacticus.Build.SourceTree.Process.SourceIntrospection import instrument

    src = (
        "line 1 without placeholder\n"
        "line 2 has {introspection:location} here\n"
        "line 3 with compact {introspection:location:compact}\n"
        "line 4 without\n"
    )
    out = instrument(src)
    assert_equal('{introspection:location:2}' in out, True,
                 "plain placeholder tagged with line 2")
    assert_equal('{introspection:location:compact:3}' in out, True,
                 "compact placeholder tagged with line 3")
    # Lines that had no placeholder must be preserved byte-for-byte.
    assert_equal(out.splitlines(keepends=True)[0], "line 1 without placeholder\n",
                 "lines without placeholder unchanged")


def test_source_introspection_location():
    print("\n=== Testing Process.SourceIntrospection.location() ===")

    from Galacticus.Build.SourceTree.Process.SourceIntrospection import location

    # Build a nested file > module > subroutine chain.
    root = {'type': 'file', 'name': 'thing.F90', 'firstChild': None,
            'parent': None, 'sibling': None}
    mod  = {'type': 'module', 'name': 'mThing', 'firstChild': None,
            'parent': root, 'sibling': None}
    sub  = {'type': 'subroutine', 'name': 'doIt', 'firstChild': None,
            'parent': mod, 'sibling': None}
    root['firstChild'] = mod
    mod['firstChild']  = sub

    expr = location(sub, 42)
    # Non-compact form names every ancestor, ends with line info.
    assert_equal("subroutine:doIt"   in expr, True, "subroutine name in location")
    assert_equal("module:mThing"     in expr, True, "module name in location")
    assert_equal("file:thing.F90"    in expr, True, "file name in location")
    assert_equal("[line 42]"         in expr, True, "line number in location")

    expr = location(sub, 7, compact=True)
    assert_equal("subroutine(doIt)"  in expr, True, "compact form uses parenthesized name")
    assert_equal("':7'"              in expr, True, "compact form appends :line")


def test_source_introspection_process():
    print("\n=== Testing Process.SourceIntrospection.process_source_introspection() ===")

    # Synthetic tree with a code node containing a pre-tagged placeholder.
    root = {'type': 'file', 'name': 'x.F90', 'firstChild': None,
            'parent': None, 'sibling': None, 'source': 'x.F90', 'line': 0}
    sub  = {'type': 'subroutine', 'name': 's', 'parent': root,
            'firstChild': None, 'sibling': None, 'source': 'x.F90', 'line': 0}
    root['firstChild'] = sub
    code = {'type': 'code', 'content': 'err("x"//{introspection:location:5})\n',
            'parent': sub, 'firstChild': None, 'sibling': None,
            'source': 'x.F90', 'line': 0}
    sub['firstChild'] = code

    from Galacticus.Build.SourceTree.Process.SourceIntrospection import \
        process_source_introspection
    process_source_introspection(root, {})

    assert_equal('{introspection:location' not in code['content'], True,
                 "placeholder fully expanded")
    assert_equal("file:x.F90"       in code['content'], True, "expansion names file")
    assert_equal("subroutine:s"     in code['content'], True, "expansion names subroutine")
    assert_equal("[line 5]"         in code['content'], True, "expansion includes line number")


# ============================================================================
# Process.ParameterMigration tests
# ============================================================================

def test_process_parameter_migration():
    print("\n=== Testing Process.ParameterMigration ===")

    tmpdir = tempfile.mkdtemp()
    try:
        aux = os.path.join(tmpdir, 'scripts', 'aux')
        os.makedirs(aux)
        with open(os.path.join(aux, 'migrations.xml'), 'w') as fh:
            fh.write('<migrations>\n'
                     '  <migration commit="aaaa"/>\n'
                     '  <migration commit="bbbb"/>\n'
                     '</migrations>\n')
        saved = os.environ.get('GALACTICUS_EXEC_PATH')
        os.environ['GALACTICUS_EXEC_PATH'] = tmpdir
        try:
            root = _parse_text(
                "subroutine s\n"
                "integer :: placeholder\n"
                "end subroutine s\n"
            )
            sub = _find_node(root, 'subroutine')
            directive = _make_directive_node(
                'parameterMigration', {'processed': False}, sub)

            from Galacticus.Build.SourceTree.Process.ParameterMigration import \
                process_parameter_migration
            process_parameter_migration(root, {})

            assert_equal(directive['directive']['processed'], True,
                         "parameterMigration directive marked processed")
            out = serialize(root)
            assert_equal('commitHash(1)="aaaa"//c_null_char' in out, True,
                         "first commit emitted")
            assert_equal('commitHash(2)="bbbb"//c_null_char' in out, True,
                         "second commit emitted")
            assert_equal(declaration_exists(sub, 'commitHash'), True,
                         "commitHash array declared")
        finally:
            if saved is None:
                del os.environ['GALACTICUS_EXEC_PATH']
            else:
                os.environ['GALACTICUS_EXEC_PATH'] = saved
    finally:
        import shutil
        shutil.rmtree(tmpdir)


# ============================================================================
# Process.DebugMPI tests
# ============================================================================

def test_process_debug_mpi():
    print("\n=== Testing Process.DebugMPI ===")

    saved = os.environ.get('GALACTICUS_FCFLAGS')
    os.environ['GALACTICUS_FCFLAGS'] = '-DDEBUGMPI -O2'
    try:
        root = _parse_text(
            "subroutine s\n"
            "call mpiSelf%broadcast(x)\n"
            "call mpiBarrier()\n"
            "end subroutine s\n"
        )
        from Galacticus.Build.SourceTree.Process.DebugMPI import process_debug_mpi
        process_debug_mpi(root, {})

        out = serialize(root)
        assert_equal("mpiSelf call to method \"broadcast\"" in out, True,
                     "debug write emitted for mpiSelf call")
        assert_equal("']: mpiBarrier'" in out, True,
                     "debug write emitted for mpiBarrier")
        # Location expression fragments should appear too.
        assert_equal("subroutine:s" in out, True,
                     "location expression embedded in debug writes")
    finally:
        if saved is None:
            del os.environ['GALACTICUS_FCFLAGS']
        else:
            os.environ['GALACTICUS_FCFLAGS'] = saved


def test_process_debug_mpi_disabled():
    print("\n=== Testing Process.DebugMPI (disabled, no flag) ===")

    saved = os.environ.pop('GALACTICUS_FCFLAGS', None)
    try:
        root = _parse_text(
            "subroutine s\n"
            "call mpiSelf%broadcast(x)\n"
            "end subroutine s\n"
        )
        original_out = serialize(root)
        from Galacticus.Build.SourceTree.Process.DebugMPI import process_debug_mpi
        process_debug_mpi(root, {})
        assert_equal(serialize(root), original_out,
                     "DebugMPI is a no-op when -DDEBUGMPI is not set")
    finally:
        if saved is not None:
            os.environ['GALACTICUS_FCFLAGS'] = saved


# ============================================================================
# Process.DebugHDF5 tests
# ============================================================================

def test_process_debug_hdf5():
    print("\n=== Testing Process.DebugHDF5 ===")

    saved = os.environ.get('GALACTICUS_FCFLAGS')
    os.environ['GALACTICUS_FCFLAGS'] = '-DDEBUGHDF5'
    try:
        root = _parse_text(
            "module someModule\n"
            "subroutine s\n"
            "use foo\n"
            "call hdf5Access%set()\n"
            "call hdf5Access%unset()\n"
            "end subroutine s\n"
            "end module someModule\n"
        )
        from Galacticus.Build.SourceTree.Process.DebugHDF5 import process_debug_hdf5
        process_debug_hdf5(root, {})

        out = serialize(root)
        assert_equal('call IO_HDF5_Start_Locked()' in out, True,
                     "hdf5Access%set() replaced with IO_HDF5_Start_Locked()")
        assert_equal('call IO_HDF5_End_Locked()' in out, True,
                     "hdf5Access%unset() replaced with IO_HDF5_End_Locked()")
        assert_equal('IO_HDF5_Start_Locked' in out and 'IO_HDF5_End_Locked' in out, True,
                     "use IO_HDF5 with only clause injected")

        # Inside the IO_HDF5 module itself, the substitution must be suppressed.
        root = _parse_text(
            "module IO_HDF5\n"
            "subroutine s\n"
            "call hdf5Access%set()\n"
            "end subroutine s\n"
            "end module IO_HDF5\n"
        )
        process_debug_hdf5(root, {})
        out = serialize(root)
        assert_equal('call hdf5Access%set()' in out, True,
                     "inside IO_HDF5 module the original call is preserved")
    finally:
        if saved is None:
            del os.environ['GALACTICUS_FCFLAGS']
        else:
            os.environ['GALACTICUS_FCFLAGS'] = saved


# ============================================================================
# Process.ProfileOpenMP tests
# ============================================================================

def test_process_profile_openmp():
    print("\n=== Testing Process.ProfileOpenMP ===")

    tmpdir = tempfile.mkdtemp()
    try:
        with open(os.path.join(tmpdir, 'openMPCriticalSections.xml'), 'w') as fh:
            fh.write('<sections>\n'
                     '  <critical name="my_section" id="7"/>\n'
                     '</sections>\n')
        saved_flags = os.environ.get('GALACTICUS_FCFLAGS')
        saved_build = os.environ.get('BUILDPATH')
        os.environ['GALACTICUS_FCFLAGS'] = '-DOMPPROFILE'
        os.environ['BUILDPATH']           = tmpdir
        try:
            root = _parse_text(
                "subroutine s\n"
                "  !$omp critical(my_section)\n"
                "  x = 1\n"
                "  !$omp end critical\n"
                "end subroutine s\n"
            )
            from Galacticus.Build.SourceTree.Process.ProfileOpenMP import process_profile_openmp
            process_profile_openmp(root, {})

            out = serialize(root)
            assert_equal('ompProfileTimeWaitStart=OMP_Get_WTime()' in out, True,
                         "wait-start timing emitted")
            assert_equal('ompProfileTimeWaitEnd=OMP_Get_WTime()'   in out, True,
                         "wait-end timing emitted")
            assert_equal('criticalSectionWaitTime(7)'              in out, True,
                         "accumulator uses id from XML")
        finally:
            if saved_flags is None:
                del os.environ['GALACTICUS_FCFLAGS']
            else:
                os.environ['GALACTICUS_FCFLAGS'] = saved_flags
            if saved_build is None:
                del os.environ['BUILDPATH']
            else:
                os.environ['BUILDPATH'] = saved_build
    finally:
        import shutil
        shutil.rmtree(tmpdir)


# ============================================================================
# Process.Constants tests (inline-value path only — GSL/kernel paths need a C compiler)
# ============================================================================

def test_process_constant_inline():
    print("\n=== Testing Process.Constants (inline value) ===")

    root = _parse_text(
        "module m\n"
        "  integer :: placeholder\n"
        "end module m\n"
    )
    mod = _find_node(root, 'module')
    directive = _make_directive_node(
        'constant',
        {
            'type':     'double precision',
            'variable': 'my_pi',
            'value':    '3.1415926535897932d0',
            'processed': False,
        },
        mod,
    )

    from Galacticus.Build.SourceTree.Process.Constants import process_constant
    process_constant(root, {})

    out = serialize(root)
    assert_equal(directive['directive']['processed'], True,
                 "constant directive marked processed")
    assert_equal('double precision, parameter, public :: my_pi=3.1415926535897932d0' in out,
                 True, "inline value rendered as a `parameter` declaration")


# ============================================================================
# Process.ConditionalCall tests
# ============================================================================

def test_process_conditional_call():
    print("\n=== Testing Process.ConditionalCall ===")

    root = _parse_text(
        "subroutine s\n"
        "integer :: x\n"
        "end subroutine s\n"
    )
    sub = _find_node(root, 'subroutine')
    directive = _make_directive_node(
        'conditionalCall',
        {
            'argument': [
                {'name': 'a', 'value': '1.0', 'condition': 'haveA'},
                {'name': 'b', 'value': '2.0', 'condition': 'haveB'},
            ],
            'call': 'call target({conditions})\n',
            'processed': False,
        },
        sub,
    )

    from Galacticus.Build.SourceTree.Process.ConditionalCall import process_conditional_call
    process_conditional_call(root, {})

    assert_equal(directive['directive']['processed'], True,
                 "conditionalCall directive marked processed")
    out = serialize(root)
    # 2^2 = 4 if-blocks, one per boolean combination.
    assert_equal(out.count('end if'), 4, "four if-then blocks emitted (2^N, N=2)")
    # Condition variables should be declared.
    assert_equal(declaration_exists(sub, 'condition1__'), True,
                 "condition1__ declared")
    assert_equal(declaration_exists(sub, 'condition2__'), True,
                 "condition2__ declared")
    # Both-present state: call target(a=1.0,b=2.0)
    assert_equal('call target(a=1.0,b=2.0)' in out, True,
                 "both-present state emits both arguments")
    # Both-absent state: call target()
    assert_equal('call target()' in out, True,
                 "both-absent state emits the empty-argument call")


def test_process_conditional_call_parameter_present():
    print("\n=== Testing Process.ConditionalCall (parameterPresent form) ===")

    root = _parse_text(
        "subroutine s\n"
        "integer :: x\n"
        "end subroutine s\n"
    )
    sub = _find_node(root, 'subroutine')
    _make_directive_node(
        'conditionalCall',
        {
            'argument': {'name': 'foo', 'value': '3.14',
                         'parameterPresent': 'parameters',
                         'parameterName': 'fooParam'},
            'call': 'call doStuff(base,{conditions})\n',
            'processed': False,
        },
        sub,
    )
    from Galacticus.Build.SourceTree.Process.ConditionalCall import process_conditional_call
    process_conditional_call(root, {})
    out = serialize(root)
    assert_equal("condition1__=parameters%isPresent('fooParam')" in out, True,
                 "parameterPresent form emits isPresent() call")


# ============================================================================
# Process.InputParametersValidate tests
# ============================================================================

def test_process_input_parameters_validate():
    print("\n=== Testing Process.InputParametersValidate ===")

    tmpdir = tempfile.mkdtemp()
    try:
        # Provide a stateStorables.xml with one functionClass whose name is
        # "myThingClass" — the walk will then tag nodes of type "myThing" as
        # the enclosing function class.
        with open(os.path.join(tmpdir, 'stateStorables.xml'), 'w') as fh:
            fh.write('<stateStorables>\n'
                     '  <functionClasses>\n'
                     '    <functionClass name="myThingClass"/>\n'
                     '  </functionClasses>\n'
                     '</stateStorables>\n')
        saved = os.environ.get('BUILDPATH')
        os.environ['BUILDPATH'] = tmpdir
        try:
            # Parse a function inside a module, and hand-wrap the function in
            # a `myThing` type node (as the Process pipeline would after
            # FunctionClass-expansion) so the walk finds the class marker
            # before the validate directive.
            root = _parse_text(
                "function f(params)\n"
                "use dummy\n"
                "end function f\n"
            )
            func = _find_node(root, 'function')
            # Artificially wrap the function in a "myThing" type node.
            mything = {'type': 'myThing', 'name': 'instance',
                       'parent': root, 'firstChild': func,
                       'sibling': None, 'source': 'x', 'line': 0}
            func['parent'] = mything
            root['firstChild'] = mything

            directive = _make_directive_node(
                'inputParametersValidate', {'source': 'params', 'processed': False},
                func,
            )

            from Galacticus.Build.SourceTree.Process.InputParametersValidate import \
                process_input_parameters_validate
            process_input_parameters_validate(root, {})

            assert_equal(directive['directive']['processed'], True,
                         "inputParametersValidate directive marked processed")
            out = serialize(root)
            assert_equal('allocatable' in out and 'allowedParameterNames_' in out, True,
                         "allowedParameterNames_ allocatable declaration added")
            assert_equal("call f%allowedParameters(allowedParameterNames_,'params',.false.)" in out, True,
                         "allowedParameters call emitted with result name")
            assert_equal("myThingDsblVldtn == 0" in out, True,
                         "function-class name prefix used in the disable check")
            assert_equal("params%checkParameters(allowedParameterNames=allowedParameterNames_)" in out,
                         True, "checkParameters invoked on the named source")
        finally:
            if saved is None:
                del os.environ['BUILDPATH']
            else:
                os.environ['BUILDPATH'] = saved
    finally:
        import shutil
        shutil.rmtree(tmpdir)


# ============================================================================
# Galacticus.Build.Directives.extract_directives tests
# ============================================================================

def _write_temp_f90(text):
    fh = tempfile.NamedTemporaryFile(
        mode='w', suffix='.F90', delete=False, encoding='utf-8')
    fh.write(text)
    fh.close()
    return fh.name


def test_extract_directives_basic():
    print("\n=== Testing Galacticus.Build.Directives.extract_directives (basic) ===")

    path = _write_temp_f90(
        "module m\n"
        "  !![\n"
        "  <addMetaProperty name=\"foo\" type=\"float\" isCreator=\"yes\"/>\n"
        "  !!]\n"
        "  integer :: x\n"
        "  !![\n"
        "  <addMetaProperty name=\"bar\" type=\"int\" isCreator=\"no\"/>\n"
        "  !!]\n"
        "end module m\n"
    )
    try:
        found = extract_directives(path, 'addMetaProperty')
        assert_equal(len(found), 2, "both addMetaProperty directives returned")
        names = sorted(d.get('name') for d in found)
        assert_equal(names, ['bar', 'foo'], "directive names captured")
        first = extract_directive(path, 'addMetaProperty')
        assert_equal(first.get('name'), 'foo',
                     "extract_directive returns only the first match")
    finally:
        os.unlink(path)


def test_extract_directives_multiline_and_wildcard():
    print("\n=== Testing extract_directives (multi-line + wildcard + rootElementType) ===")

    path = _write_temp_f90(
        "  !![\n"
        "  <functionClass>\n"
        "    <name>fooImpl</name>\n"
        "  </functionClass>\n"
        "  !!]\n"
        "  !![\n"
        "  <addMetaProperty name=\"x\" isCreator=\"yes\"/>\n"
        "  !!]\n"
    )
    try:
        all_ = extract_directives(path, '*', set_root_element_type=True)
        assert_equal(len(all_), 2, "two directives across wildcards found")
        roots = sorted(d['rootElementType'] for d in all_)
        assert_equal(roots, ['addMetaProperty', 'functionClass'],
                     "rootElementType tagged on every returned dict")
        # Multi-line <functionClass>...</functionClass> inner <name> captured.
        fc = next(d for d in all_ if d['rootElementType'] == 'functionClass')
        assert_equal(fc.get('name'), 'fooImpl',
                     "multi-line directive body parsed into fields")
    finally:
        os.unlink(path)


def test_extract_directives_conditions():
    print("\n=== Testing extract_directives (conditions filter) ===")

    path = _write_temp_f90(
        "  !![\n"
        "  <addMetaProperty name=\"a\" isCreator=\"yes\"/>\n"
        "  !!]\n"
        "  !![\n"
        "  <addMetaProperty name=\"b\" isCreator=\"no\"/>\n"
        "  !!]\n"
        "  !![\n"
        "  <addMetaProperty name=\"c\" isCreator=\"yes\"/>\n"
        "  !!]\n"
    )
    try:
        creators = extract_directives(
            path, 'addMetaProperty', conditions={'isCreator': 'yes'})
        names = sorted(d.get('name') for d in creators)
        assert_equal(names, ['a', 'c'],
                     "conditions={isCreator:yes} filters non-creators out")
    finally:
        os.unlink(path)


def test_extract_directives_missing_file():
    print("\n=== Testing extract_directives (missing file) ===")

    result = extract_directives('/nonexistent/path/does.not.exist.F90', 'foo')
    assert_equal(result, [], "missing file yields an empty list, not an error")


# ============================================================================
# Process.StateStore tests (stateStore + stateRestore)
# ============================================================================

def test_process_state_store():
    print("\n=== Testing Process.StateStore (stateStore directive) ===")

    root = {'type': 'file', 'name': 'x.F90', 'firstChild': None,
            'parent': None, 'sibling': None, 'source': 'x', 'line': 0}
    sub  = {'type': 'subroutine', 'name': 's', 'parent': root,
            'firstChild': None, 'sibling': None, 'source': 'x', 'line': 0}
    root['firstChild'] = sub
    directive = _make_directive_node(
        'stateStore', {'variables': 'alpha beta', 'processed': False}, sub)

    from Galacticus.Build.SourceTree.Process.StateStore import process_state_store
    process_state_store(root, {})

    assert_equal(directive['directive']['processed'], True,
                 "stateStore directive marked processed")
    emitted = directive.get('sibling')
    assert_equal(emitted is not None and emitted.get('type') == 'code', True,
                 "code node inserted after stateStore directive")
    content = emitted['content']
    assert_equal('write (stateFile) associated(alpha)' in content, True,
                 "associated() write emitted for alpha")
    assert_equal('call alpha%stateStore(stateFile,gslStateFile,stateOperationID)'
                 in content, True,
                 "conditional stateStore call emitted for alpha")
    assert_equal('write (stateFile) associated(beta)' in content, True,
                 "associated() write emitted for beta")


def test_process_state_restore():
    print("\n=== Testing Process.StateStore (stateRestore directive) ===")

    root = {'type': 'file', 'name': 'x.F90', 'firstChild': None,
            'parent': None, 'sibling': None, 'source': 'x', 'line': 0}
    sub  = {'type': 'subroutine', 'name': 's', 'parent': root,
            'firstChild': None, 'sibling': None, 'source': 'x', 'line': 0}
    root['firstChild'] = sub
    directive = _make_directive_node(
        'stateRestore', {'variables': 'alpha', 'processed': False}, sub)

    from Galacticus.Build.SourceTree.Process.StateStore import process_state_store
    process_state_store(root, {})

    out = serialize(root)
    assert_equal(directive['directive']['processed'], True,
                 "stateRestore directive marked processed")
    assert_equal('read (stateFile) wasAllocated_' in out, True,
                 "read wasAllocated_ emitted")
    assert_equal("'alpha' was stored, but is now not allocated" in out, True,
                 "error message references the variable name")
    assert_equal('call alpha%stateRestore(stateFile,gslStateFile,stateOperationID)'
                 in out, True, "stateRestore call emitted")
    assert_equal(declaration_exists(sub, 'wasAllocated_'), True,
                 "wasAllocated_ added as a logical declaration")
    assert_equal('use :: Error' in out or 'use Error' in out or 'Error' in out, True,
                 "Error module imported for Error_Report")


# ============================================================================
# Process.MetaPropertyDatabase tests
# ============================================================================

def test_process_add_meta_property():
    print("\n=== Testing Process.AddMetaProperty ===")
    from Galacticus.Build.SourceTree.Process.AddMetaProperty import (
        process_add_meta_property,
    )

    # Happy path: float rank-0, isCreator=yes, isEvolvable=yes.
    root = {'type': 'file', 'name': 'x.F90', 'firstChild': None,
            'parent': None, 'sibling': None, 'source': 'x', 'line': 0}
    sub  = {'type': 'subroutine', 'name': 's', 'parent': root,
            'firstChild': None, 'sibling': None, 'source': 'x', 'line': 0}
    root['firstChild'] = sub
    directive = _make_directive_node(
        'addMetaProperty',
        {
            'id':        'mpID',
            'component': 'darkMatter',
            'name':      'spin',
            'type':      'float',
            'rank':      0,
            'isCreator': 'yes',
            'isEvolvable': 'yes',
        },
        sub,
    )
    process_add_meta_property(root, {})
    assert_equal(directive['directive']['processed'], True,
                 "addMetaProperty directive marked processed")
    emitted = directive.get('sibling')
    assert_equal(emitted is not None and emitted.get('type') == 'code', True,
                 "code node inserted after directive")
    content = emitted['content']
    assert_equal(
        "mpID=defaultDarkMatterComponent%addFloatRank0MetaProperty"
        "(var_str('spin'),'darkMatter:spin',isCreator=.true.,"
        "isEvolvable=.true.)" in content, True,
        "float rank-0 emission includes isEvolvable and prefixed name")

    # Happy path: integer rank-1 (omits isEvolvable; only float rank-0 carries it).
    root2 = {'type': 'file', 'name': 'y.F90', 'firstChild': None,
             'parent': None, 'sibling': None, 'source': 'y', 'line': 0}
    sub2  = {'type': 'subroutine', 'name': 's', 'parent': root2,
             'firstChild': None, 'sibling': None, 'source': 'y', 'line': 0}
    root2['firstChild'] = sub2
    d2 = _make_directive_node(
        'addMetaProperty',
        {
            'id':        'mpID2',
            'component': 'basic',
            'name':      'indices',
            'type':      'integer',
            'rank':      1,
            'isCreator': 'no',
        },
        sub2,
    )
    process_add_meta_property(root2, {})
    c2 = d2['sibling']['content']
    assert_equal(
        "mpID2=defaultBasicComponent%addIntegerRank1MetaProperty"
        "(var_str('indices'),'basic:indices',isCreator=.false.)" in c2, True,
        "integer rank-1 emission omits isEvolvable, uses .false. isCreator")

    # Raw-string name starting with a quote is passed through untouched.
    root3 = {'type': 'file', 'name': 'z.F90', 'firstChild': None,
             'parent': None, 'sibling': None, 'source': 'z', 'line': 0}
    sub3  = {'type': 'subroutine', 'name': 's', 'parent': root3,
             'firstChild': None, 'sibling': None, 'source': 'z', 'line': 0}
    root3['firstChild'] = sub3
    d3 = _make_directive_node(
        'addMetaProperty',
        {
            'id':        'mpID3',
            'component': 'basic',
            'name':      "'dynamicName'",
            'type':      'float',
            'rank':      0,
            'isCreator': 'no',
            'isEvolvable': 'no',
        },
        sub3,
    )
    process_add_meta_property(root3, {})
    c3 = d3['sibling']['content']
    assert_equal(
        "var_str('dynamicName')" in c3 and
        "'basic:'//'dynamicName'" in c3, True,
        "quoted name is passed through var_str(...) with // concatenation")

    # Validation: rank > 1 is rejected.
    root4 = {'type': 'file', 'name': 'w.F90', 'firstChild': None,
             'parent': None, 'sibling': None, 'source': 'w', 'line': 0}
    sub4  = {'type': 'subroutine', 'name': 's', 'parent': root4,
             'firstChild': None, 'sibling': None, 'source': 'w', 'line': 0}
    root4['firstChild'] = sub4
    _make_directive_node(
        'addMetaProperty',
        {
            'id':        'mpID4',
            'component': 'basic',
            'name':      'tensor',
            'type':      'float',
            'rank':      2,
        },
        sub4,
    )
    try:
        process_add_meta_property(root4, {})
        raised = False
    except RuntimeError as e:
        raised = 'rank > 1' in str(e)
    assert_equal(raised, True,
                 "rank > 1 meta-property raises a RuntimeError")


def test_process_meta_property_database():
    print("\n=== Testing Process.MetaPropertyDatabase ===")

    tmpdir = tempfile.mkdtemp()
    try:
        # One source file containing a functionClass definition plus a creator
        # addMetaProperty directive.  The creator's name is the implementation
        # name with the functionClass prefix; the stripped, lcfirst'd remainder
        # is what the generated code records as implementationName.
        src_file = os.path.join(tmpdir, 'dummy.F90')
        with open(src_file, 'w') as fh:
            fh.write(
                "  !![\n"
                "  <myFunc>\n"
                "    <name>myFuncImpl1</name>\n"
                "  </myFunc>\n"
                "  !!]\n"
                "  !![\n"
                '  <addMetaProperty component="node" name="mass" type="float" rank="0" isCreator="yes"/>\n'
                "  !!]\n"
            )
        with open(os.path.join(tmpdir, 'directiveLocations.xml'), 'w') as fh:
            fh.write(
                "<directiveLocations>\n"
                "  <addMetaProperty>\n"
                f"    <file>{src_file}</file>\n"
                "  </addMetaProperty>\n"
                "</directiveLocations>\n"
            )
        with open(os.path.join(tmpdir, 'stateStorables.xml'), 'w') as fh:
            fh.write(
                "<stateStorables>\n"
                "  <functionClasses>\n"
                "    <functionClass name=\"myFuncClass\"/>\n"
                "  </functionClasses>\n"
                "</stateStorables>\n"
            )
        saved = os.environ.get('BUILDPATH')
        os.environ['BUILDPATH'] = tmpdir
        try:
            root = {'type': 'file', 'name': 'driver.F90', 'firstChild': None,
                    'parent': None, 'sibling': None, 'source': 'driver.F90', 'line': 0}
            mod = {'type': 'module', 'name': 'driver', 'parent': root,
                   'firstChild': None, 'sibling': None,
                   'source': 'driver.F90', 'line': 0}
            root['firstChild'] = mod
            directive = _make_directive_node(
                'metaPropertyDatabase', {'processed': False}, mod)

            from Galacticus.Build.SourceTree.Process.MetaPropertyDatabase \
                import process_meta_property_database
            process_meta_property_database(root, {})

            assert_equal(directive['directive']['processed'], True,
                         "metaPropertyDatabase directive marked processed")
            emitted = directive.get('sibling')
            assert_equal(emitted is not None and emitted.get('type') == 'code',
                         True, "code node inserted after directive")
            content = emitted['content']
            assert_equal(
                'subroutine metaPropertyNoCreator' in content,
                True, "metaPropertyNoCreator subroutine synthesized",
            )
            assert_equal(
                "if (component_ == 'node'"
                " .and. name_ == 'mass'"
                " .and. type_ == 'float'"
                " .and. rank_ == 0) then"
                in content,
                True, "if-branch emitted for the registered creator",
            )
            assert_equal("className='myFunc'" in content, True,
                         "className populated with functionClass name")
            assert_equal("implementationName='impl1'" in content, True,
                         "implementationName is the lcfirst-stripped suffix")
        finally:
            if saved is None:
                del os.environ['BUILDPATH']
            else:
                os.environ['BUILDPATH'] = saved
    finally:
        import shutil
        shutil.rmtree(tmpdir)


# ============================================================================
# Process.InputParameter tests
# ============================================================================

def _write_input_parameter_build_dir():
    """Create a minimal BUILDPATH with stub Makefile_All_Execs and
    directiveLocations.xml for InputParameter tests.  Returns the path.
    """
    tmpdir = tempfile.mkdtemp()
    with open(os.path.join(tmpdir, 'Makefile_All_Execs'), 'w') as fh:
        fh.write("all_exes = Galacticus.exe tests.foo.exe\n")
    with open(os.path.join(tmpdir, 'directiveLocations.xml'), 'w') as fh:
        fh.write("<directiveLocations></directiveLocations>\n")
    return tmpdir


def test_process_input_parameter_simple():
    print("\n=== Testing Process.InputParameter (simple form) ===")

    tmpdir = _write_input_parameter_build_dir()
    saved = os.environ.get('BUILDPATH')
    os.environ['BUILDPATH'] = tmpdir
    try:
        root = _parse_text(
            "subroutine s\n"
            "integer :: counter\n"
            "end subroutine s\n"
        )
        sub = _find_node(root, 'subroutine')
        _make_directive_node(
            'inputParameter',
            {'source': 'parameters', 'name': 'counter', 'processed': False},
            sub,
        )

        from Galacticus.Build.SourceTree.Process.InputParameter \
            import process_input_parameters
        process_input_parameters(root, {})

        out = serialize(root)
        assert_equal("call parameters%value('counter',counter)" in out, True,
                     "value() call emitted with quoted parameter name")
        assert_equal('use :: Input_Parameters' in out or 'use Input_Parameters' in out,
                     True, "Input_Parameters module imported")
    finally:
        if saved is None:
            del os.environ['BUILDPATH']
        else:
            os.environ['BUILDPATH'] = saved
        import shutil
        shutil.rmtree(tmpdir)


def test_process_input_parameter_with_default_and_no_output():
    print("\n=== Testing Process.InputParameter (defaultValue + writeOutput=no) ===")

    tmpdir = _write_input_parameter_build_dir()
    saved = os.environ.get('BUILDPATH')
    os.environ['BUILDPATH'] = tmpdir
    try:
        root = _parse_text(
            "subroutine s\n"
            "real :: tolerance\n"
            "end subroutine s\n"
        )
        sub = _find_node(root, 'subroutine')
        _make_directive_node(
            'inputParameter',
            {
                'source':        'parameters',
                'name':          'tolerance',
                'variable':      'tolerance',
                'defaultValue':  '1.0d-6',
                'writeOutput':   'no',
                'processed':     False,
            },
            sub,
        )
        from Galacticus.Build.SourceTree.Process.InputParameter \
            import process_input_parameters
        process_input_parameters(root, {})
        out = serialize(root)
        assert_equal('defaultValue=1.0d-6' in out, True, "defaultValue passed through")
        assert_equal('writeOutput=.false.' in out, True,
                     "writeOutput='no' renders as .false.")
    finally:
        if saved is None:
            del os.environ['BUILDPATH']
        else:
            os.environ['BUILDPATH'] = saved
        import shutil
        shutil.rmtree(tmpdir)


# ============================================================================
# SourceTree insert_pre_contains / insert_post_contains tests
# ============================================================================

def test_insert_pre_post_contains():
    print("\n=== Testing insert_pre_contains / insert_post_contains ===")

    # Module that already has a `contains` marker.
    root = _parse_text(
        "module m\n"
        "integer :: x\n"
        "contains\n"
        "subroutine s\n"
        "end subroutine s\n"
        "end module m\n"
    )
    mod = _find_node(root, 'module')
    insert_pre_contains(mod, [{
        'type': 'code', 'content': 'real :: injected\n',
        'parent': None, 'firstChild': None, 'sibling': None,
        'source': 'test', 'line': 0,
    }])
    insert_post_contains(mod, [{
        'type': 'code', 'content': '! injected post\n',
        'parent': None, 'firstChild': None, 'sibling': None,
        'source': 'test', 'line': 0,
    }])
    out = serialize(root)
    assert_equal('real :: injected' in out.split('contains')[0], True,
                 "pre-contains insertion lands before the contains marker")
    assert_equal('! injected post' in out.split('contains')[1], True,
                 "post-contains insertion lands after the contains marker")

    # Module without a contains marker: insert_post_contains should create one.
    root = _parse_text(
        "module m\n"
        "integer :: x\n"
        "end module m\n"
    )
    mod = _find_node(root, 'module')
    insert_post_contains(mod, [{
        'type': 'code', 'content': '! added sub\n',
        'parent': None, 'firstChild': None, 'sibling': None,
        'source': 'test', 'line': 0,
    }])
    out = serialize(root)
    assert_equal('contains' in out, True,
                 "contains marker auto-created when absent")
    assert_equal('! added sub' in out.split('contains')[-1], True,
                 "new code appears after the auto-created contains")


# ============================================================================
# Process.Constructors tests
# ============================================================================

def _write_state_storables(dir_path,
                           function_classes=(),
                           function_class_instances=()):
    """Write a minimal stateStorables.xml exposing the given class / instance names."""
    lines = ["<stateStorables>"]
    if function_classes:
        lines.append("  <functionClasses>")
        for name in function_classes:
            lines.append(f"    <functionClass name=\"{name}\"/>")
        lines.append("  </functionClasses>")
    for inst in function_class_instances:
        lines.append(f"  <functionClassInstances>{inst}</functionClassInstances>")
    lines.append("</stateStorables>\n")
    with open(os.path.join(dir_path, 'stateStorables.xml'), 'w') as fh:
        fh.write("\n".join(lines))


def test_process_constructors_basic_assignment():
    print("\n=== Testing Process.Constructors (scalar + optional) ===")

    tmpdir = tempfile.mkdtemp()
    try:
        _write_state_storables(tmpdir)
        # Invalidate module cache so the test uses our fresh stateStorables.xml.
        import Galacticus.Build.SourceTree.Process.Constructors as C
        C._STATE_STORABLES = None
        saved = os.environ.get('BUILDPATH')
        os.environ['BUILDPATH'] = tmpdir
        try:
            root = _parse_text(
                "function s(a, b) result(self)\n"
                "integer, intent(in) :: a\n"
                "real, intent(in), optional :: b\n"
                "end function s\n"
            )
            func = _find_node(root, 'function')
            _make_directive_node(
                'constructorAssign',
                {'variables': 'a, b', 'processed': False},
                func,
            )
            from Galacticus.Build.SourceTree.Process.Constructors import \
                process_constructors
            process_constructors(root, {})

            out = serialize(root)
            assert_equal('self%a=a' in out, True,
                         "mandatory scalar argument assigned to self")
            assert_equal('if (present(b)) self%b=b' in out, True,
                         "optional argument wrapped in present() guard")
        finally:
            if saved is None:
                del os.environ['BUILDPATH']
            else:
                os.environ['BUILDPATH'] = saved
            C._STATE_STORABLES = None
    finally:
        import shutil; shutil.rmtree(tmpdir)


def test_process_constructors_pointer_refcount():
    print("\n=== Testing Process.Constructors (pointer + refcount) ===")

    tmpdir = tempfile.mkdtemp()
    try:
        _write_state_storables(tmpdir, function_classes=['myfuncclass'])
        import Galacticus.Build.SourceTree.Process.Constructors as C
        C._STATE_STORABLES = None
        saved = os.environ.get('BUILDPATH')
        os.environ['BUILDPATH'] = tmpdir
        try:
            root = _parse_text(
                "function s(obj) result(self)\n"
                "type(myFuncClass), intent(in), target :: obj\n"
                "end function s\n"
            )
            func = _find_node(root, 'function')
            _make_directive_node(
                'constructorAssign',
                {'variables': '*obj', 'processed': False},
                func,
            )
            from Galacticus.Build.SourceTree.Process.Constructors import \
                process_constructors
            process_constructors(root, {})

            out = serialize(root)
            assert_equal('self%obj => obj' in out, True,
                         "pointer assignment uses => operator")
            assert_equal('call self%obj%referenceCountIncrement()' in out, True,
                         "functionClass pointer triggers referenceCountIncrement call")
        finally:
            if saved is None:
                del os.environ['BUILDPATH']
            else:
                os.environ['BUILDPATH'] = saved
            C._STATE_STORABLES = None
    finally:
        import shutil; shutil.rmtree(tmpdir)


def test_process_constructors_wildcard_type():
    print("\n=== Testing Process.Constructors (class(*) polymorphic pointer) ===")

    tmpdir = tempfile.mkdtemp()
    try:
        _write_state_storables(tmpdir)
        import Galacticus.Build.SourceTree.Process.Constructors as C
        C._STATE_STORABLES = None
        saved = os.environ.get('BUILDPATH')
        os.environ['BUILDPATH'] = tmpdir
        try:
            root = _parse_text(
                "function s(obj) result(self)\n"
                "class(*), intent(in), target :: obj\n"
                "end function s\n"
            )
            func = _find_node(root, 'function')
            _make_directive_node(
                'constructorAssign',
                {'variables': '*obj', 'processed': False},
                func,
            )
            from Galacticus.Build.SourceTree.Process.Constructors import \
                process_constructors
            process_constructors(root, {})
            out = serialize(root)
            assert_equal('select type(s__ => self%obj)' in out, True,
                         "wildcard class pointer emits select-type block")
            assert_equal('class is (functionClass)' in out, True,
                         "select-type narrows to functionClass branch")
            assert_equal('call s__%referenceCountIncrement()' in out, True,
                         "referenceCountIncrement called via the type-guarded alias")
        finally:
            if saved is None:
                del os.environ['BUILDPATH']
            else:
                os.environ['BUILDPATH'] = saved
            C._STATE_STORABLES = None
    finally:
        import shutil; shutil.rmtree(tmpdir)


def test_process_constructors_slash_suppresses_count():
    print("\n=== Testing Process.Constructors (*/ suppresses referenceCountIncrement) ===")

    tmpdir = tempfile.mkdtemp()
    try:
        _write_state_storables(tmpdir, function_classes=['myfuncclass'])
        import Galacticus.Build.SourceTree.Process.Constructors as C
        C._STATE_STORABLES = None
        saved = os.environ.get('BUILDPATH')
        os.environ['BUILDPATH'] = tmpdir
        try:
            root = _parse_text(
                "function s(obj) result(self)\n"
                "type(myFuncClass), intent(in), target :: obj\n"
                "end function s\n"
            )
            func = _find_node(root, 'function')
            _make_directive_node(
                'constructorAssign',
                {'variables': '*/obj', 'processed': False},
                func,
            )
            from Galacticus.Build.SourceTree.Process.Constructors import \
                process_constructors
            process_constructors(root, {})
            out = serialize(root)
            assert_equal('self%obj => obj' in out, True,
                         "pointer assignment still emitted for `*/obj`")
            assert_equal('referenceCountIncrement' not in out, True,
                         "*/prefix suppresses referenceCountIncrement call")
        finally:
            if saved is None:
                del os.environ['BUILDPATH']
            else:
                os.environ['BUILDPATH'] = saved
            C._STATE_STORABLES = None
    finally:
        import shutil; shutil.rmtree(tmpdir)


# ============================================================================
# Process.FunctionsGlobal tests
# ============================================================================

def _write_functions_global_build_dir(source_files):
    """Set up a BUILDPATH tempdir whose directiveLocations.xml points at the
    given source files (each already containing <functionGlobal> directives).
    Returns the tempdir path.
    """
    tmpdir = tempfile.mkdtemp()
    file_lines = "\n".join(f"    <file>{f}</file>" for f in source_files)
    with open(os.path.join(tmpdir, 'directiveLocations.xml'), 'w') as fh:
        fh.write(
            "<directiveLocations>\n"
            "  <functionGlobal>\n"
            f"{file_lines}\n"
            "  </functionGlobal>\n"
            "</directiveLocations>\n"
        )
    with open(os.path.join(tmpdir, 'stateStorables.xml'), 'w') as fh:
        fh.write("<stateStorables></stateStorables>\n")
    # FunctionsGlobal's pointers path reparses the generated code via
    # process_tree(); that invocation runs every registered Process hook, so
    # provide empty versions of the artifacts the other hooks read eagerly
    # (HDF5FCInterop, InputParameter, …) to keep the test hermetic.
    open(os.path.join(tmpdir, 'hdf5FCInterop.dat'), 'w').close()
    with open(os.path.join(tmpdir, 'Makefile_All_Execs'), 'w') as fh:
        fh.write("all_exes =\n")
    return tmpdir


def test_process_functions_global_pointers():
    print("\n=== Testing Process.FunctionsGlobal (type='pointers') ===")

    # Write a source file with two functionGlobal directives.
    src_file = _write_temp_f90(
        "  !![\n"
        "  <functionGlobal>\n"
        "   <unitName>alpha</unitName>\n"
        "   <type>double precision</type>\n"
        "   <arguments>integer, intent(in) :: i</arguments>\n"
        "  </functionGlobal>\n"
        "  !!]\n"
        "  !![\n"
        "  <functionGlobal>\n"
        "   <unitName>beta</unitName>\n"
        "   <type>void</type>\n"
        "   <arguments>integer, intent(in) :: j</arguments>\n"
        "  </functionGlobal>\n"
        "  !!]\n"
    )
    tmpdir = _write_functions_global_build_dir([src_file])

    # Invalidate module caches.
    import Galacticus.Build.SourceTree.Process.FunctionsGlobal as FG
    FG._DIRECTIVE_LOCATIONS = None
    FG._STATE_STORABLES     = None

    saved = os.environ.get('BUILDPATH')
    os.environ['BUILDPATH'] = tmpdir
    try:
        root = _parse_text(
            "module m\n"
            "public\n"
            "end module m\n"
        )
        mod = _find_node(root, 'module')
        _make_directive_node(
            'functionsGlobal', {'type': 'pointers', 'processed': False}, mod,
        )
        from Galacticus.Build.SourceTree.Process.FunctionsGlobal import \
            process_functions_global
        process_functions_global(root, {})

        out = serialize(root)
        assert_equal('procedure(alpha_Null), pointer :: alpha_ => alpha_Null' in out,
                     True, "pointer declaration emitted for alpha")
        assert_equal('procedure(beta_Null), pointer :: beta_ => beta_Null' in out,
                     True, "pointer declaration emitted for beta")
        assert_equal('function alpha_Null' in out, True,
                     "alpha_Null function definition emitted")
        assert_equal('subroutine beta_Null' in out, True,
                     "beta_Null subroutine emitted for void return type")
        assert_equal('alpha_Null=0.0d0' in out, True,
                     "double-precision default return initialized to 0.0d0")
        assert_equal('contains' in out, True,
                     "contains marker present (auto-inserted when missing)")
    finally:
        os.unlink(src_file)
        if saved is None:
            del os.environ['BUILDPATH']
        else:
            os.environ['BUILDPATH'] = saved
        FG._DIRECTIVE_LOCATIONS = None
        FG._STATE_STORABLES     = None
        import shutil; shutil.rmtree(tmpdir)


def test_process_functions_global_establish():
    print("\n=== Testing Process.FunctionsGlobal (type='establish') ===")

    src_file = _write_temp_f90(
        "module provider\n"
        "contains\n"
        "  subroutine alpha\n"
        "  end subroutine alpha\n"
        "end module provider\n"
        "  !![\n"
        "  <functionGlobal>\n"
        "   <unitName>alpha</unitName>\n"
        "   <type>void</type>\n"
        "  </functionGlobal>\n"
        "  !!]\n"
    )
    tmpdir = _write_functions_global_build_dir([src_file])

    import Galacticus.Build.SourceTree.Process.FunctionsGlobal as FG
    FG._DIRECTIVE_LOCATIONS = None
    FG._STATE_STORABLES     = None

    saved = os.environ.get('BUILDPATH')
    os.environ['BUILDPATH'] = tmpdir
    try:
        root = _parse_text(
            "subroutine wireup\n"
            "end subroutine wireup\n"
        )
        sub = _find_node(root, 'subroutine')
        _make_directive_node(
            'functionsGlobal', {'type': 'establish', 'processed': False}, sub,
        )
        from Galacticus.Build.SourceTree.Process.FunctionsGlobal import \
            process_functions_global
        process_functions_global(root, {})

        out = serialize(root)
        assert_equal('alpha_ => alpha' in out, True,
                     "establish path emits the pointer assignment")
        assert_equal('use' in out and 'provider' in out and 'alpha' in out, True,
                     "use statement for the providing module injected")
    finally:
        os.unlink(src_file)
        if saved is None:
            del os.environ['BUILDPATH']
        else:
            os.environ['BUILDPATH'] = saved
        FG._DIRECTIVE_LOCATIONS = None
        FG._STATE_STORABLES     = None
        import shutil; shutil.rmtree(tmpdir)


# ============================================================================
# set_visibility tests
# ============================================================================

def test_set_visibility():
    print("\n=== Testing set_visibility ===")

    root = _parse_text(
        "module m\n"
        "use foo\n"
        "integer :: x\n"
        "end module m\n"
    )
    mod = _find_node(root, 'module')
    set_visibility(mod, 'doStuff', 'public')
    set_visibility(mod, 'secretVar', 'private')
    set_visibility(mod, 'alsoDo',  'public')
    out = serialize(root)
    assert_equal('public :: alsoDo, doStuff' in out, True,
                 "public list sorted and joined")
    assert_equal('private :: secretVar' in out, True,
                 "private list emitted")
    # Visibility block must come AFTER the `use foo` (matches Perl ordering).
    use_pos = out.index('use foo')
    vis_pos = out.index('public')
    assert_equal(vis_pos > use_pos, True,
                 "visibility block placed after moduleUse sibling")


# ============================================================================
# Process.Enumeration tests
# ============================================================================

def _build_enumeration_tree(directive):
    """Build a minimal module→enumeration-directive tree and return (root, mod, dir)."""
    root = _parse_text(
        "module m\n"
        "end module m\n"
    )
    mod = _find_node(root, 'module')
    dir_node = _make_directive_node('enumeration', directive, mod)
    return root, mod, dir_node


def test_process_enumeration_basic():
    print("\n=== Testing Process.Enumeration (basic type + equality) ===")

    root, mod, _ = _build_enumeration_tree({
        'name':  'color',
        'entry': [
            {'label': 'red',   'description': 'the red member'},
            {'label': 'green', 'description': 'the green member'},
        ],
        'processed': False,
    })
    from Galacticus.Build.SourceTree.Process.Enumeration import \
        process_enumerations
    process_enumerations(root, {})

    out = serialize(root)
    assert_equal('type, extends(enumerationType) :: enumerationcolorType' in out,
                 True, "enumeration type declared")
    assert_equal('colorRed=enumerationcolorType(0)' in out, True,
                 "first member uses indexing=0 by default")
    assert_equal('colorGreen=enumerationcolorType(1)' in out, True,
                 "second member increments")
    assert_equal('enumerationColorIsEqual' in out, True,
                 "equality function emitted (name ucfirst'd per Perl)")
    assert_equal('enumerationColorDescribe' in out, True,
                 "describe function always emitted")


def test_process_enumeration_validator():
    print("\n=== Testing Process.Enumeration (validator=yes + Min/Max/Count) ===")

    root, mod, _ = _build_enumeration_tree({
        'name':      'mode',
        'validator': 'yes',
        'indexing':  '1',
        'entry': [
            {'label': 'slow'},
            {'label': 'fast'},
        ],
        'processed': False,
    })
    from Galacticus.Build.SourceTree.Process.Enumeration import \
        process_enumerations
    process_enumerations(root, {})
    out = serialize(root)
    assert_equal('modeMin  =1' in out, True, "Min param emitted")
    assert_equal('modeMax  =2' in out, True, "Max param emitted")
    assert_equal('modeCount=2' in out, True, "Count param emitted")
    assert_equal('enumerationModeIsValid' in out, True,
                 "validator function emitted")


def test_process_enumeration_encode_decode():
    print("\n=== Testing Process.Enumeration (encode/decode function families) ===")

    root, mod, _ = _build_enumeration_tree({
        'name':           'kind',
        'encodeFunction': 'yes',
        'decodeFunction': 'yes',
        'entry': [
            {'label': 'flat', 'description': 'flat geometry'},
            {'label': 'open', 'description': 'open geometry'},
        ],
        'processed': False,
    })
    from Galacticus.Build.SourceTree.Process.Enumeration import \
        process_enumerations
    process_enumerations(root, {})
    out = serialize(root)
    assert_equal('interface enumerationKindEncode' in out, True,
                 "encode interface emitted")
    assert_equal('enumerationKindEncodeCharchar' not in out, True,
                 "function names not accidentally concatenated")
    assert_equal('interface enumerationKindDecode' in out, True,
                 "decode interface emitted")
    assert_equal('interface enumerationKindDescription' in out, True,
                 "description interface emitted alongside decode")
    assert_equal("case ('kindFlat')" in out, True,
                 "encode case emitted with prefix")
    assert_equal("case ('flat')" in out, True,
                 "encode case emitted without prefix")


# ============================================================================
# Process.DeepCopyActions tests
# ============================================================================

def test_process_deep_copy_actions_setto():
    print("\n=== Testing Process.DeepCopyActions (setTo action) ===")

    # A base class + one derived non-abstract class, with a directive
    # requesting a setTo on the derived class's `counter` variable.
    root = _parse_text(
        "module m\n"
        "type, abstract :: baseT\n"
        "end type baseT\n"
        "type, extends(baseT) :: subT\n"
        "integer :: counter\n"
        "end type subT\n"
        "end module m\n"
    )
    mod = _find_node(root, 'module')
    _make_directive_node(
        'deepCopyActions',
        {
            'class':    'baseT',
            'subT':     {'setTo': {'variables': 'counter', 'state': '0'}},
            'processed': False,
        },
        mod,
    )

    # deepCopyActions reads $BUILDPATH/deepCopyActions.xml eagerly in Perl;
    # the Python port delays the read until it's needed, but process_tree()
    # fires every registered hook, so provide empty stub XML for the others.
    tmpdir = tempfile.mkdtemp()
    open(os.path.join(tmpdir, 'deepCopyActions.xml'), 'w').close()
    with open(os.path.join(tmpdir, 'stateStorables.xml'), 'w') as fh:
        fh.write("<stateStorables></stateStorables>\n")
    with open(os.path.join(tmpdir, 'directiveLocations.xml'), 'w') as fh:
        fh.write("<directiveLocations></directiveLocations>\n")
    open(os.path.join(tmpdir, 'hdf5FCInterop.dat'), 'w').close()
    with open(os.path.join(tmpdir, 'Makefile_All_Execs'), 'w') as fh:
        fh.write("all_exes =\n")

    saved = os.environ.get('BUILDPATH')
    os.environ['BUILDPATH'] = tmpdir
    try:
        from Galacticus.Build.SourceTree.Process.DeepCopyActions import \
            process_deep_copy_actions
        process_deep_copy_actions(root, {})

        out = serialize(root)
        assert_equal('subroutine baseTDeepCopyActions(self)' in out, True,
                     "generated subroutine named after target class")
        assert_equal('type is (subT)' in out, True,
                     "select-type branch emitted for non-abstract descendant")
        assert_equal('self%counter=0' in out, True,
                     "setTo action emits self%var=state")
        assert_equal('procedure :: deepCopyActions => baseTDeepCopyActions' in out,
                     True, "type-binding inserted into base class")
    finally:
        if saved is None:
            del os.environ['BUILDPATH']
        else:
            os.environ['BUILDPATH'] = saved
        import shutil
        shutil.rmtree(tmpdir)


def test_process_deep_copy_actions_method_call():
    print("\n=== Testing Process.DeepCopyActions (methodCall action) ===")

    root = _parse_text(
        "module m\n"
        "type, abstract :: baseT\n"
        "end type baseT\n"
        "type, extends(baseT) :: subT\n"
        "integer :: counter\n"
        "end type subT\n"
        "end module m\n"
    )
    mod = _find_node(root, 'module')
    _make_directive_node(
        'deepCopyActions',
        {
            'class':    'baseT',
            'subT':     {'methodCall': {'method': 'reset', 'arguments': '42'}},
            'processed': False,
        },
        mod,
    )

    tmpdir = tempfile.mkdtemp()
    for f in ('deepCopyActions.xml', 'hdf5FCInterop.dat'):
        open(os.path.join(tmpdir, f), 'w').close()
    with open(os.path.join(tmpdir, 'stateStorables.xml'), 'w') as fh:
        fh.write("<stateStorables></stateStorables>\n")
    with open(os.path.join(tmpdir, 'directiveLocations.xml'), 'w') as fh:
        fh.write("<directiveLocations></directiveLocations>\n")
    with open(os.path.join(tmpdir, 'Makefile_All_Execs'), 'w') as fh:
        fh.write("all_exes =\n")

    saved = os.environ.get('BUILDPATH')
    os.environ['BUILDPATH'] = tmpdir
    try:
        from Galacticus.Build.SourceTree.Process.DeepCopyActions import \
            process_deep_copy_actions
        process_deep_copy_actions(root, {})
        out = serialize(root)
        assert_equal('call self%reset(42)' in out, True,
                     "methodCall action emits `call self%method(args)`")
    finally:
        if saved is None:
            del os.environ['BUILDPATH']
        else:
            os.environ['BUILDPATH'] = saved
        import shutil
        shutil.rmtree(tmpdir)


# ============================================================================
# Process.StateStorable tests
# ============================================================================

def _write_state_storable_build_dir(storable_type_names):
    """Create a BUILDPATH tempdir with a stateStorables.xml listing the given
    type names and the other stub files every hook eagerly opens.  Returns
    the tempdir.
    """
    tmpdir = tempfile.mkdtemp()
    lines = ['<stateStorables>']
    for name in storable_type_names:
        lines.append(f'  <stateStorables type="{name}"/>')
    lines.append('</stateStorables>\n')
    with open(os.path.join(tmpdir, 'stateStorables.xml'), 'w') as fh:
        fh.write("\n".join(lines))
    # Stubs used by OTHER hooks' eager loads if anything ever triggers
    # process_tree() alongside this test.
    with open(os.path.join(tmpdir, 'directiveLocations.xml'), 'w') as fh:
        fh.write("<directiveLocations></directiveLocations>\n")
    open(os.path.join(tmpdir, 'hdf5FCInterop.dat'), 'w').close()
    open(os.path.join(tmpdir, 'deepCopyActions.xml'), 'w').close()
    with open(os.path.join(tmpdir, 'Makefile_All_Execs'), 'w') as fh:
        fh.write("all_exes =\n")
    return tmpdir


def test_process_state_storable_basic():
    print("\n=== Testing Process.StateStorable (scalar stateStorable class) ===")

    tmpdir = _write_state_storable_build_dir(['subT'])
    saved = os.environ.get('BUILDPATH')
    os.environ['BUILDPATH'] = tmpdir
    try:
        root = _parse_text(
            "module m\n"
            "type, abstract :: baseT\n"
            "end type baseT\n"
            "type, extends(baseT) :: subT\n"
            "integer :: counter\n"
            "real    :: weight\n"
            "end type subT\n"
            "end module m\n"
        )
        mod = _find_node(root, 'module')
        _make_directive_node(
            'stateStorable',
            {'class': 'baseT', 'processed': False},
            mod,
        )
        from Galacticus.Build.SourceTree.Process.StateStorable import \
            process_state_storable
        process_state_storable(root, {})

        out = serialize(root)
        assert_equal('subroutine baseTStateStore' in out, True,
                     "StateStore subroutine synthesized")
        assert_equal('subroutine baseTStateRestore' in out, True,
                     "StateRestore subroutine synthesized")
        assert_equal('type is (subT)' in out, True,
                     "non-abstract descendant gets a select-type branch")
        assert_equal('procedure :: stateStore   => baseTStateStore' in out, True,
                     "stateStore type-binding emitted on base class")
        assert_equal('procedure :: stateRestore => baseTStateRestore' in out, True,
                     "stateRestore type-binding emitted on base class")
        assert_equal('write (stateFile) self%counter, &' in out
                     or 'write (stateFile) self%counter,' in out, True,
                     "static intrinsic members batched into a single write")
    finally:
        if saved is None:
            del os.environ['BUILDPATH']
        else:
            os.environ['BUILDPATH'] = saved
        import shutil
        shutil.rmtree(tmpdir)


def test_process_state_storable_excludes_and_restore_to():
    print("\n=== Testing Process.StateStorable (exclude + restoreTo) ===")

    tmpdir = _write_state_storable_build_dir([])
    saved = os.environ.get('BUILDPATH')
    os.environ['BUILDPATH'] = tmpdir
    try:
        root = _parse_text(
            "module m\n"
            "type :: baseT\n"
            "integer :: keep\n"
            "integer :: skip_me\n"
            "integer :: reset_me\n"
            "end type baseT\n"
            "end module m\n"
        )
        mod = _find_node(root, 'module')
        _make_directive_node(
            'stateStorable',
            {
                'class': 'baseT',
                'baseT': {
                    'exclude':   {'variables': 'skip_me'},
                    'restoreTo': {'variables': 'reset_me', 'state': '42'},
                },
                'processed': False,
            },
            mod,
        )
        from Galacticus.Build.SourceTree.Process.StateStorable import \
            process_state_storable
        process_state_storable(root, {})
        out = serialize(root)
        assert_equal('self%keep' in out, True,
                     "non-excluded variable present in generated code")
        assert_equal('self%skip_me' not in out, True,
                     "excluded variable absent from generated code")
        assert_equal('self%reset_me=42' in out, True,
                     "restoreTo directive emits a state assignment rather than a read")
    finally:
        if saved is None:
            del os.environ['BUILDPATH']
        else:
            os.environ['BUILDPATH'] = saved
        import shutil
        shutil.rmtree(tmpdir)


def test_process_state_storable_class_restore():
    print("\n=== Testing Process.StateStorable (class-restore dispatchers) ===")

    tmpdir = _write_state_storable_build_dir([])
    saved = os.environ.get('BUILDPATH')
    os.environ['BUILDPATH'] = tmpdir
    try:
        root = _parse_text(
            "module m\n"
            "type, abstract :: baseT\n"
            "end type baseT\n"
            "type, extends(baseT) :: subT\n"
            "integer :: n\n"
            "end type subT\n"
            "end module m\n"
        )
        mod = _find_node(root, 'module')
        _make_directive_node(
            'stateStorable',
            {'class': 'baseT', 'processed': False},
            mod,
        )
        from Galacticus.Build.SourceTree.Process.StateStorable import \
            process_state_storable
        process_state_storable(root, {})

        out = serialize(root)
        assert_equal('subroutine baseTClassRestore(self,stateFile)' in out, True,
                     "rank-0 class-restore dispatcher emitted")
        assert_equal('subroutine baseTClassRestore1D(self,stateFile,storedShape)'
                     in out, True,
                     "rank-1 class-restore dispatcher emitted")
        assert_equal('allocate(subT :: self' in out, True,
                     "dispatcher allocates the concrete subclass")
        assert_equal('baseTClassRestore' in out and 'public' in out, True,
                     "ClassRestore symbol made public via set_visibility")
    finally:
        if saved is None:
            del os.environ['BUILDPATH']
        else:
            os.environ['BUILDPATH'] = saved
        import shutil
        shutil.rmtree(tmpdir)


# ============================================================================
# Galacticus.Build.Dependencies.dependency_sort tests
# ============================================================================

def test_dependency_sort_after_and_before():
    print("\n=== Testing Galacticus.Build.Dependencies.dependency_sort ===")

    from Galacticus.Build.Dependencies import dependency_sort

    # `X.after = Y` means Perl puts X before Y (inverted English reading —
    # this is the semantics the callers rely on).
    order = dependency_sort({'X': {'after': 'Y'}, 'Y': {}})
    assert_equal(order, ['X', 'Y'],
                 "`after` edge places the tagged name before its reference")

    # `X.before = Y` means Perl puts Y before X.
    order = dependency_sort({'X': {'before': 'Y'}, 'Y': {}})
    assert_equal(order, ['Y', 'X'],
                 "`before` edge places the tagged name after its reference")

    # No constraints → alphabetical.
    order = dependency_sort({'c': {}, 'a': {}, 'b': {}})
    assert_equal(order, ['a', 'b', 'c'],
                 "constraint-free tasks fall back to alphabetical order")


# ============================================================================
# Process.EventHooksStatic tests
# ============================================================================

def test_process_event_hooks_static_single_hook():
    print("\n=== Testing Process.EventHooksStatic (single hook) ===")

    tmpdir = tempfile.mkdtemp()
    # File containing the hooked function (sub `foo_hook` inside module `mFoo`).
    hooked_file = os.path.join(tmpdir, 'hooked.F90')
    with open(hooked_file, 'w') as fh:
        fh.write(
            "module mFoo\n"
            "contains\n"
            "  subroutine foo_hook\n"
            "  end subroutine foo_hook\n"
            "  !![\n"
            "  <myHook function=\"foo_hook\"/>\n"
            "  !!]\n"
            "end module mFoo\n"
        )
    # directiveLocations.xml announces both `myHook` locations and
    # eventHookStatic locations.
    with open(os.path.join(tmpdir, 'directiveLocations.xml'), 'w') as fh:
        fh.write(
            "<directiveLocations>\n"
            f"  <myHook><file>{hooked_file}</file></myHook>\n"
            f"  <eventHookStatic><file>{hooked_file}</file></eventHookStatic>\n"
            "</directiveLocations>\n"
        )
    with open(os.path.join(tmpdir, 'stateStorables.xml'), 'w') as fh:
        fh.write("<stateStorables></stateStorables>\n")

    # Invalidate caches so the test sees our temp XMLs.
    import Galacticus.Build.SourceTree.Process.EventHooksStatic as EHS
    EHS._DIRECTIVE_LOCATIONS = None
    EHS._STATE_STORABLES     = None
    EHS._EVENT_HOOK_NAMES    = None

    saved = os.environ.get('BUILDPATH')
    os.environ['BUILDPATH'] = tmpdir
    try:
        root = _parse_text(
            "subroutine driver\n"
            "end subroutine driver\n"
        )
        sub = _find_node(root, 'subroutine')
        _make_directive_node(
            'eventHookStatic',
            {'name': 'myHook', 'processed': False},
            sub,
        )
        from Galacticus.Build.SourceTree.Process.EventHooksStatic import \
            process_event_hooks_static
        process_event_hooks_static(root, {})

        out = serialize(root)
        assert_equal('call foo_hook()' in out, True,
                     "call to hooked function emitted at directive site")
        assert_equal('mFoo' in out and 'foo_hook' in out, True,
                     "providing module imported alongside the hooked function")
    finally:
        if saved is None:
            del os.environ['BUILDPATH']
        else:
            os.environ['BUILDPATH'] = saved
        EHS._DIRECTIVE_LOCATIONS = None
        EHS._STATE_STORABLES     = None
        EHS._EVENT_HOOK_NAMES    = None
        import shutil
        shutil.rmtree(tmpdir)


def test_process_event_hooks_static_ordering():
    print("\n=== Testing Process.EventHooksStatic (after/before ordering) ===")

    tmpdir = tempfile.mkdtemp()
    first_file  = os.path.join(tmpdir, 'first.F90')
    second_file = os.path.join(tmpdir, 'second.F90')
    with open(first_file, 'w') as fh:
        fh.write(
            "module mFirst\n"
            "contains\n"
            "  subroutine first_hook\n"
            "  end subroutine first_hook\n"
            "  !![\n"
            "  <myHook function=\"first_hook\" after=\"second_hook\"/>\n"
            "  !!]\n"
            "end module mFirst\n"
        )
    with open(second_file, 'w') as fh:
        fh.write(
            "module mSecond\n"
            "contains\n"
            "  subroutine second_hook\n"
            "  end subroutine second_hook\n"
            "  !![\n"
            "  <myHook function=\"second_hook\"/>\n"
            "  !!]\n"
            "end module mSecond\n"
        )
    with open(os.path.join(tmpdir, 'directiveLocations.xml'), 'w') as fh:
        fh.write(
            "<directiveLocations>\n"
            f"  <myHook><file>{first_file}</file><file>{second_file}</file></myHook>\n"
            f"  <eventHookStatic><file>{first_file}</file><file>{second_file}</file></eventHookStatic>\n"
            "</directiveLocations>\n"
        )
    with open(os.path.join(tmpdir, 'stateStorables.xml'), 'w') as fh:
        fh.write("<stateStorables></stateStorables>\n")

    import Galacticus.Build.SourceTree.Process.EventHooksStatic as EHS
    EHS._DIRECTIVE_LOCATIONS = None
    EHS._STATE_STORABLES     = None
    EHS._EVENT_HOOK_NAMES    = None

    saved = os.environ.get('BUILDPATH')
    os.environ['BUILDPATH'] = tmpdir
    try:
        root = _parse_text(
            "subroutine driver\n"
            "end subroutine driver\n"
        )
        sub = _find_node(root, 'subroutine')
        _make_directive_node(
            'eventHookStatic',
            {'name': 'myHook', 'processed': False},
            sub,
        )
        from Galacticus.Build.SourceTree.Process.EventHooksStatic import \
            process_event_hooks_static
        process_event_hooks_static(root, {})
        out = serialize(root)

        # Perl's `after` inversion: first_hook (which has after=second_hook)
        # is emitted BEFORE second_hook.
        first_pos  = out.index('call first_hook')
        second_pos = out.index('call second_hook')
        assert_equal(first_pos < second_pos, True,
                     "after=other places the tagged call before `other`")
    finally:
        if saved is None:
            del os.environ['BUILDPATH']
        else:
            os.environ['BUILDPATH'] = saved
        EHS._DIRECTIVE_LOCATIONS = None
        EHS._STATE_STORABLES     = None
        EHS._EVENT_HOOK_NAMES    = None
        import shutil
        shutil.rmtree(tmpdir)


# ============================================================================
# Process.EventHooks tests
# ============================================================================

def _write_event_hooks_build_dir(source_files=()):
    """Minimal BUILDPATH for EventHooks: directiveLocations.xml pointing at
    the given source files (each containing `<eventHook>` directives), plus
    all the other stubs eager loaders read."""
    tmpdir = tempfile.mkdtemp()
    files_xml = "\n".join(f"    <file>{f}</file>" for f in source_files)
    with open(os.path.join(tmpdir, 'directiveLocations.xml'), 'w') as fh:
        fh.write(
            "<directiveLocations>\n"
            "  <eventHook>\n"
            f"{files_xml}\n"
            "  </eventHook>\n"
            "</directiveLocations>\n"
        )
    with open(os.path.join(tmpdir, 'stateStorables.xml'), 'w') as fh:
        fh.write("<stateStorables></stateStorables>\n")
    open(os.path.join(tmpdir, 'hdf5FCInterop.dat'), 'w').close()
    open(os.path.join(tmpdir, 'deepCopyActions.xml'), 'w').close()
    with open(os.path.join(tmpdir, 'Makefile_All_Execs'), 'w') as fh:
        fh.write("all_exes =\n")
    return tmpdir


def test_event_hooks_interface_type_get():
    print("\n=== Testing EventHooks interface_type_get ===")
    from Galacticus.Build.SourceTree.Process.EventHooks import _interface_type_get
    assert_equal(_interface_type_get({'name': 'foo'}), 'Unspecified',
                 "hook without <interface> → 'Unspecified'")
    # md5("foo") in hex.
    import hashlib
    assert_equal(
        _interface_type_get({'name': 'foo', 'interface': 'x'}),
        hashlib.md5(b'foo').hexdigest(),
        "hook with <interface> → md5_hex of the hook name")


def test_process_event_hook_call_site():
    print("\n=== Testing Process.EventHooks (call-site dispatch) ===")

    tmpdir = _write_event_hooks_build_dir([])
    import Galacticus.Build.SourceTree.Process.EventHooks as EH
    EH._DIRECTIVE_LOCATIONS = None
    saved = os.environ.get('BUILDPATH')
    os.environ['BUILDPATH'] = tmpdir
    try:
        root = _parse_text(
            "subroutine driver\n"
            "end subroutine driver\n"
        )
        sub = _find_node(root, 'subroutine')
        directive = _make_directive_node(
            'eventHook',
            {'name': 'tick', 'callWith': 'a,b', 'processed': False},
            sub,
        )
        from Galacticus.Build.SourceTree.Process.EventHooks import \
            process_event_hooks
        process_event_hooks(root, {})
        out = serialize(root)
        assert_equal('tickEventGlobal%count()' in out, True,
                     "global dispatch emitted for hook 'tick'")
        assert_equal('call hook_%function_(hook_%object_,a,b)' in out, True,
                     "hook callsite uses the callWith argument list")
        assert_equal(declaration_exists(sub, 'tickIterator'), True,
                     "iterator declaration added to enclosing scope")
    finally:
        if saved is None:
            del os.environ['BUILDPATH']
        else:
            os.environ['BUILDPATH'] = saved
        EH._DIRECTIVE_LOCATIONS = None
        import shutil
        shutil.rmtree(tmpdir)


def test_process_event_hook_openmp_wrappers():
    print("\n=== Testing Process.EventHooks (OpenMP parallel wrappers) ===")

    tmpdir = _write_event_hooks_build_dir([])
    import Galacticus.Build.SourceTree.Process.EventHooks as EH
    EH._DIRECTIVE_LOCATIONS = None
    saved = os.environ.get('BUILDPATH')
    os.environ['BUILDPATH'] = tmpdir
    try:
        root = _parse_text(
            "subroutine s\n"
            "!$omp parallel\n"
            "x = 1\n"
            "!$omp end parallel\n"
            "end subroutine s\n"
        )
        from Galacticus.Build.SourceTree.Process.EventHooks import \
            process_event_hooks
        process_event_hooks(root, {})
        out = serialize(root)
        assert_equal('call eventsHooksFilterCopyOut()' in out, True,
                     "copy-out injected before !$omp parallel")
        assert_equal('call eventsHooksFilterCopyIn()' in out, True,
                     "copy-in injected after !$omp parallel")
        assert_equal('call eventsHooksFilterRestore()' in out, True,
                     "restore injected before !$omp end parallel")
    finally:
        if saved is None:
            del os.environ['BUILDPATH']
        else:
            os.environ['BUILDPATH'] = saved
        EH._DIRECTIVE_LOCATIONS = None
        import shutil
        shutil.rmtree(tmpdir)


def test_process_event_hook_manager():
    print("\n=== Testing Process.EventHooks (manager directive) ===")

    # Build the BUILDPATH first, then drop the hook file in it and list it.
    tmpdir = _write_event_hooks_build_dir([])
    hook_file = os.path.join(tmpdir, 'hook.F90')
    with open(hook_file, 'w') as fh:
        fh.write(
            "  !![\n"
            "  <eventHook name=\"tick\">\n"
            "   <interface>\n"
            "    integer, intent(in) :: i\n"
            "   </interface>\n"
            "  </eventHook>\n"
            "  !!]\n"
        )
    with open(os.path.join(tmpdir, 'directiveLocations.xml'), 'w') as fh:
        fh.write(
            "<directiveLocations>\n"
            f"  <eventHook><file>{hook_file}</file></eventHook>\n"
            "</directiveLocations>\n"
        )

    import shutil
    import Galacticus.Build.SourceTree.Process.EventHooks as EH
    EH._DIRECTIVE_LOCATIONS = None

    saved = os.environ.get('BUILDPATH')
    os.environ['BUILDPATH'] = tmpdir
    try:
        root = _parse_text(
            "module m\n"
            "end module m\n"
        )
        mod = _find_node(root, 'module')
        _make_directive_node(
            'eventHookManager', {'processed': False}, mod,
        )
        from Galacticus.Build.SourceTree.Process.EventHooks import \
            process_event_hooks
        process_event_hooks(root, {})
        out = serialize(root)
        import hashlib
        iface = hashlib.md5(b'tick').hexdigest()
        assert_equal(f'type, extends(hook) :: hook{iface}' in out, True,
                     "per-interface hook type emitted")
        assert_equal(f'type, extends(eventHook) :: eventHook{iface}' in out, True,
                     "per-interface eventHook type emitted")
        assert_equal(f'subroutine eventHook{iface}Attach' in out, True,
                     "attach subroutine emitted for the interface")
        assert_equal(f'subroutine eventHook{iface}Detach' in out, True,
                     "detach subroutine emitted for the interface")
        assert_equal(f'function eventHook{iface}IsAttached' in out, True,
                     "isAttached function emitted for the interface")
        assert_equal('subroutine eventsHooksInitialize' in out, True,
                     "initializer subroutine emitted")
        assert_equal('subroutine eventsHooksFilterCopyOut_' in out, True,
                     "copy-out subroutine emitted")
        assert_equal('subroutine eventsHooksFilterRestore_' in out, True,
                     "restore subroutine emitted")
        assert_equal('subroutine eventsHooksWaitTimes' in out, True,
                     "wait-times subroutine emitted (OMPPROFILE-gated)")
        assert_equal(f'type(eventHook{iface}), public  :: tickEvent' in out, True,
                     "module-level event variable declared")
    finally:
        if saved is None:
            del os.environ['BUILDPATH']
        else:
            os.environ['BUILDPATH'] = saved
        EH._DIRECTIVE_LOCATIONS = None
        shutil.rmtree(tmpdir)


# ============================================================================
# Process.ThreadSafeIO tests
# ============================================================================

def test_process_thread_safe_io_enabled():
    print("\n=== Testing Process.ThreadSafeIO (enabled, wraps I/O) ===")

    saved = os.environ.get('GALACTICUS_FCFLAGS')
    os.environ['GALACTICUS_FCFLAGS'] = '-DTHREADSAFEIO'
    try:
        root = _parse_text(
            "subroutine s\n"
            "write (unit,*) 'hello'\n"
            "x = 1\n"
            "read (unit,*) y\n"
            "end subroutine s\n"
        )
        from Galacticus.Build.SourceTree.Process.ThreadSafeIO import \
            process_thread_safe_io
        process_thread_safe_io(root, {})
        out = serialize(root)
        assert_equal('!$omp critical(gfortranInternalIO_)' in out, True,
                     "I/O block opened with gfortranInternalIO_ critical")
        assert_equal('!$omp end critical(gfortranInternalIO_)' in out, True,
                     "I/O block closed with matching end critical")
        # The write/read lines must still appear verbatim.
        assert_equal("write (unit,*) 'hello'" in out, True,
                     "original write statement preserved")
    finally:
        if saved is None:
            del os.environ['GALACTICUS_FCFLAGS']
        else:
            os.environ['GALACTICUS_FCFLAGS'] = saved


def test_process_thread_safe_io_disabled():
    print("\n=== Testing Process.ThreadSafeIO (disabled) ===")
    saved = os.environ.pop('GALACTICUS_FCFLAGS', None)
    try:
        root = _parse_text(
            "subroutine s\n"
            "write (unit,*) 'hi'\n"
            "end subroutine s\n"
        )
        before = serialize(root)
        from Galacticus.Build.SourceTree.Process.ThreadSafeIO import \
            process_thread_safe_io
        process_thread_safe_io(root, {})
        assert_equal(serialize(root), before,
                     "no changes when -DTHREADSAFEIO is absent")
    finally:
        if saved is not None:
            os.environ['GALACTICUS_FCFLAGS'] = saved


def test_process_thread_safe_io_skips_error_module():
    print("\n=== Testing Process.ThreadSafeIO (skips Error module) ===")
    saved = os.environ.get('GALACTICUS_FCFLAGS')
    os.environ['GALACTICUS_FCFLAGS'] = '-DTHREADSAFEIO'
    try:
        root = _parse_text(
            "module Error\n"
            "subroutine r\n"
            "write (unit,*) 'err'\n"
            "end subroutine r\n"
            "end module Error\n"
        )
        before = serialize(root)
        from Galacticus.Build.SourceTree.Process.ThreadSafeIO import \
            process_thread_safe_io
        process_thread_safe_io(root, {})
        assert_equal(serialize(root), before,
                     "Error module is left untouched (no deadlock risk)")
    finally:
        if saved is None:
            del os.environ['GALACTICUS_FCFLAGS']
        else:
            os.environ['GALACTICUS_FCFLAGS'] = saved


# ============================================================================
# Process.SourceDigest tests
# ============================================================================

def test_process_source_digest():
    print("\n=== Testing Process.SourceDigest ===")

    root = _parse_text(
        "module m\n"
        "end module m\n"
    )
    mod = _find_node(root, 'module')
    _make_directive_node(
        'sourceDigest',
        {'name': 'someFile', 'processed': False},
        mod,
    )
    from Galacticus.Build.SourceTree.Process.SourceDigest import \
        process_source_digests
    process_source_digests(root, {})
    out = serialize(root)
    assert_equal('character(C_Char)' in out and 'dimension(23)' in out, True,
                 "MD5 character array declaration emitted")
    assert_equal('bind(C, name="someFileMD5")' in out, True,
                 "C binding tag matches Perl-generated symbol name")
    assert_equal('use' in out and 'ISO_C_Binding' in out and 'C_Char' in out,
                 True, "ISO_C_Binding, only : C_Char imported")


# ============================================================================
# Process.ObjectBuilder tests (lighter directives only; objectBuilder proper
# is exercised by real Galacticus sources and requires stateStorables.xml)
# ============================================================================

def test_process_object_destructor():
    print("\n=== Testing Process.ObjectBuilder (objectDestructor) ===")

    root = _parse_text(
        "subroutine s\n"
        "end subroutine s\n"
    )
    sub = _find_node(root, 'subroutine')
    _make_directive_node(
        'objectDestructor',
        {'name': 'obj_', 'processed': False},
        sub,
    )
    from Galacticus.Build.SourceTree.Process.ObjectBuilder import \
        process_object_builder

    # Provide a minimal BUILDPATH with an empty stateStorables.xml so the
    # lazy load inside process_object_builder doesn't raise.
    tmpdir = tempfile.mkdtemp()
    with open(os.path.join(tmpdir, 'stateStorables.xml'), 'w') as fh:
        fh.write("<stateStorables></stateStorables>\n")
    import Galacticus.Build.SourceTree.Process.ObjectBuilder as OB
    OB._STATE_STORABLES = None
    saved = os.environ.get('BUILDPATH')
    os.environ['BUILDPATH'] = tmpdir
    try:
        process_object_builder(root, {})
        out = serialize(root)
        assert_equal('if (associated(obj_)) then' in out, True,
                     "associated() guard emitted")
        assert_equal('referenceCount_=obj_%referenceCountDecrement()' in out, True,
                     "reference count decrement call emitted")
        assert_equal('deallocate(obj_)' in out, True,
                     "deallocate emitted on zero-count branch")
        assert_equal('nullify(obj_)' in out, True,
                     "nullify emitted on nonzero-count branch")
        assert_equal(declaration_exists(sub, 'referenceCount_'), True,
                     "referenceCount_ local declared")
    finally:
        if saved is None:
            del os.environ['BUILDPATH']
        else:
            os.environ['BUILDPATH'] = saved
        OB._STATE_STORABLES = None
        import shutil; shutil.rmtree(tmpdir)


def test_process_reference_increment_acquire_construct():
    print("\n=== Testing Process.ObjectBuilder (reference* and deepCopy) ===")

    root = _parse_text("subroutine s\nend subroutine s\n")
    sub = _find_node(root, 'subroutine')
    _make_directive_node(
        'referenceCountIncrement',
        {'object': 'ptr_', 'owner': 'self', 'processed': False},
        sub,
    )
    _make_directive_node(
        'referenceAcquire',
        {'owner': 'self', 'target': 'dst_', 'source': 'src_', 'processed': False},
        sub,
    )
    _make_directive_node(
        'referenceConstruct',
        {'object': 'thing_', 'constructor': 'thingType(args)', 'processed': False},
        sub,
    )
    _make_directive_node(
        'deepCopy',
        {'source': 'a_', 'destination': 'b_', 'processed': False},
        sub,
    )

    tmpdir = tempfile.mkdtemp()
    with open(os.path.join(tmpdir, 'stateStorables.xml'), 'w') as fh:
        fh.write("<stateStorables></stateStorables>\n")
    import Galacticus.Build.SourceTree.Process.ObjectBuilder as OB
    OB._STATE_STORABLES = None
    saved = os.environ.get('BUILDPATH')
    os.environ['BUILDPATH'] = tmpdir
    try:
        from Galacticus.Build.SourceTree.Process.ObjectBuilder import \
            process_object_builder
        process_object_builder(root, {})
        out = serialize(root)
        assert_equal('call self%ptr_%referenceCountIncrement()' in out, True,
                     "referenceCountIncrement emits owner-prefixed call")
        assert_equal('self%dst_ => src_' in out
                     and 'call self%dst_%referenceCountIncrement()' in out,
                     True,
                     "referenceAcquire emits pointer assignment + increment")
        assert_equal('thing_=thingType(args)' in out
                     and 'call thing_%autoHook()' in out,
                     True,
                     "referenceConstruct emits assignment + autoHook")
        assert_equal('call a_%deepCopy(b_)' in out
                     and 'a_%copiedSelf => b_' in out
                     and 'call b_%autoHook()' in out,
                     True,
                     "deepCopy emits deepCopy + copiedSelf + autoHook lines")
    finally:
        if saved is None:
            del os.environ['BUILDPATH']
        else:
            os.environ['BUILDPATH'] = saved
        OB._STATE_STORABLES = None
        import shutil; shutil.rmtree(tmpdir)


# ============================================================================
# Process.ClassDocumentation tests
# ============================================================================

def test_process_class_documentation_disabled():
    print("\n=== Testing Process.ClassDocumentation (disabled by default) ===")

    saved = os.environ.pop('GALACTICUS_BUILD_DOCS', None)
    try:
        root = _parse_text(
            "module m\n"
            "type :: aT\n"
            "end type aT\n"
            "end module m\n"
        )
        before = serialize(root)
        from Galacticus.Build.SourceTree.Process.ClassDocumentation import \
            process_class_documentation
        process_class_documentation(root, {})
        assert_equal(serialize(root), before,
                     "no changes when GALACTICUS_BUILD_DOCS is not 'yes'")
    finally:
        if saved is not None:
            os.environ['GALACTICUS_BUILD_DOCS'] = saved


def test_process_class_documentation_enabled():
    print("\n=== Testing Process.ClassDocumentation (enabled, writes XML) ===")

    tmpdir = tempfile.mkdtemp()
    saved_build = os.environ.get('BUILDPATH')
    saved_docs  = os.environ.get('GALACTICUS_BUILD_DOCS')
    os.environ['BUILDPATH']             = tmpdir
    os.environ['GALACTICUS_BUILD_DOCS'] = 'yes'

    # Write the source to a real file so `tree.get('name')` matches the
    # `<base>.F90` shape that triggers the XML dump.
    src_path = os.path.join(tmpdir, 'driver.F90')
    with open(src_path, 'w') as fh:
        fh.write(
            "module m\n"
            "type :: aT\n"
            "contains\n"
            "  procedure :: doIt => aT_doIt\n"
            "end type aT\n"
            "contains\n"
            "subroutine aT_doIt(self, x)\n"
            "class(aT), intent(inout) :: self\n"
            "real,       intent(in   ) :: x\n"
            "end subroutine aT_doIt\n"
            "end module m\n"
        )
    try:
        # Reset module-level output-previous tracker so we get a fresh write.
        import Galacticus.Build.SourceTree.Process.ClassDocumentation as CD
        CD._OUTPUT_PREVIOUS.clear()

        tree = parse_file(src_path)
        from Galacticus.Build.SourceTree.Process.ClassDocumentation import \
            process_class_documentation
        process_class_documentation(tree, {})

        out_xml = os.path.join(tmpdir, 'driver.classes.xml')
        assert_equal(os.path.exists(out_xml), True,
                     "<basename>.classes.xml written to BUILDPATH")
        with open(out_xml, 'r') as fh:
            payload = fh.read()
        assert_equal('<classes>' in payload, True, "root element is <classes>")
        assert_equal('<aT>' in payload or '<aT/>' in payload, True,
                     "class 'aT' recorded in the XML")
    finally:
        if saved_build is None:
            del os.environ['BUILDPATH']
        else:
            os.environ['BUILDPATH'] = saved_build
        if saved_docs is None:
            del os.environ['GALACTICUS_BUILD_DOCS']
        else:
            os.environ['GALACTICUS_BUILD_DOCS'] = saved_docs
        import shutil
        shutil.rmtree(tmpdir)


# ============================================================================
# Process.FunctionClass.Utils tests
# ============================================================================

def test_functionclass_utils_small_helpers():
    print("\n=== Testing FunctionClass.Utils (small helpers) ===")
    assert_equal(latex_breakable('fooBarBaz'), r'foo\-Bar\-Baz',
                 "latex_breakable inserts `\\-` on lower→upper transitions")
    assert_equal(trimlc('  Hello  '), 'hello',
                 "trimlc lowercases + trims leading/trailing whitespace")
    assert_equal(striplc('Hello World'), 'helloworld',
                 "striplc lowercases + removes ALL whitespace")
    assert_equal(lctrim('Hello   '), 'hello',
                 "lctrim strips trailing whitespace + lowercases")
    assert_equal(strip_variable_name('foo(1,2)=5'), 'foo',
                 "strip_variable_name returns leading identifier only")
    assert_equal(declaration_rank({'attributes': ['intent(in)']}), 0,
                 "declaration_rank returns 0 for scalars")
    assert_equal(declaration_rank({'attributes': ['dimension(:,:,:)']}), 3,
                 "declaration_rank returns dim-count for arrays")


def test_functionclass_utils_class_dependencies():
    print("\n=== Testing FunctionClass.Utils.class_dependencies ===")
    root = _parse_text(
        "module m\n"
        "  type, extends(myFoo) :: myFooBar\n"
        "    class(myFooBaz), pointer :: related\n"
        "    class(myFooClass), pointer :: abstract_class_ref\n"
        "  end type myFooBar\n"
        "end module m\n"
    )
    mod = _find_node(root, 'module')
    # Inject a directive node of type `myFoo` on the module's child chain
    # so class_dependencies has a directive to capture attributes from.
    # The walk starts at the directive itself (or before the type node) and
    # finds the first matching `type myFoo<suffix>` node.
    type_node = _find_node(root, 'type')
    directive = _make_directive_node(
        'myFoo',
        {'name': 'myFoo', 'processed': False,
         'node': 'foo', 'foo_attr': 'bar'},
        mod,
    )
    # Run class_dependencies starting at the directive.
    class_record, deps = class_dependencies(directive, 'myFoo')
    assert_equal(class_record.get('extends'), 'myFoo',
                 "extends parent name captured")
    assert_equal(class_record.get('type'), 'myFooBar',
                 "class type = directive + suffix captured")
    assert_equal('myFoo' in deps, True,
                 "parent type recorded as dependency")
    assert_equal('myFooBaz' in deps, True,
                 "cross-referenced non-Class member type recorded")
    # ...but the `myFooClass` abstract class reference is excluded because
    # its name ends with 'Class'.
    assert_equal('myFooClass' not in deps, True,
                 "abstract `…Class` references excluded from dependencies")


# ============================================================================
# Process.FunctionClass.Descriptor tests
# ============================================================================

def test_functionclass_descriptor_classification():
    print("\n=== Testing FunctionClass.Descriptor.potential_descriptor_parameters ===")
    state_storables = {
        'functionClasses': {'functionClass': [
            {'name': 'myFooClass'}, {'name': 'myBarClass'},
        ]},
        'functionClassInstances': ['myFoo', 'myBar'],
    }
    declarations = [
        # pointer to functionClass → object
        {'intrinsic': 'class', 'type': 'myFooClass',
         'attributes': ['pointer'], 'variables': ['obj']},
        # stateful type
        {'intrinsic': 'type', 'type': 'statefulDouble',
         'attributes': [], 'variables': ['sDbl']},
        # enumeration type
        {'intrinsic': 'type', 'type': 'enumerationKindType',
         'attributes': [], 'variables': ['ek']},
        # intrinsic integer → parameter
        {'intrinsic': 'integer', 'type': None,
         'attributes': [], 'variables': ['counter']},
        # varying_string type → parameter
        {'intrinsic': 'type', 'type': 'varying_string',
         'attributes': [], 'variables': ['name']},
        # custom descriptor procedure
        {'intrinsic': 'procedure', 'type': 'customD',
         'attributes': [], 'variables': ['descriptor=>customD']},
    ]
    non_abstract_class = {}
    names = {}
    potential_descriptor_parameters(
        declarations, non_abstract_class,
        {'linkedList': {'object': 'alpha beta'}},
        state_storables, names)
    assert_equal(names.get('objects'), ['obj'],
                 "pointer-to-functionClass classified as object")
    assert_equal(len(names.get('statefulTypes') or []), 1,
                 "stateful type captured")
    assert_equal(len(names.get('enumerations') or []), 1,
                 "enumeration type captured")
    assert_equal(len(names.get('parameters') or []), 2,
                 "intrinsic + varying_string captured as parameters")
    assert_equal(names.get('linkedListObjects'), ['alpha', 'beta'],
                 "linkedList objects extracted from class metadata")
    assert_equal(non_abstract_class.get('hasCustomDescriptor'), True,
                 "procedure :: descriptor=>... flips hasCustomDescriptor")


# ============================================================================
# Process.FunctionClass.LinkedList tests
# ============================================================================

def _linked_list_spec():
    return {
        'type':       'myItem',
        'variable':   'first',
        'next':       'next',
        'object':     'foo bar',
        'objectType': 'fooType barType',
        'module':     'myMod',
    }


def test_functionclass_linkedlist_register_and_module():
    print("\n=== Testing FunctionClass.LinkedList (register + module helpers) ===")
    ll = _linked_list_spec()
    assert_equal(linked_list_module(ll), 'myMod',
                 "linked_list_module returns the `module` attr")
    vars_ = []
    linked_list_register_variable(ll, vars_, 'x')
    linked_list_register_variable(ll, vars_, 'x')
    assert_equal(len(vars_), 1,
                 "linked_list_register_variable is idempotent on identical type")
    assert_equal(vars_[0]['attributes'], ['pointer'],
                 "registered variable carries `pointer` attribute")


def test_functionclass_linkedlist_deep_copy():
    print("\n=== Testing FunctionClass.LinkedList.deep_copy_linked_list ===")
    ll = _linked_list_spec()
    dc, dcr, dcf, mod = deep_copy_linked_list(
        {'linkedList': ll, 'node': {'line': 1}}, {}, [], [], [])
    assert_equal('do while (associated(myItemitem))' in dc, True,
                 "deepCopy body iterates over myItemitem")
    assert_equal('deepCopy(myItemitemNew%foo)' in dc, True,
                 "deepCopy body invokes %deepCopy() on each object")
    assert_equal('deepCopyReset' in dcr, True,
                 "reset code calls deepCopyReset")
    assert_equal('deepCopyFinalize' in dcf, True,
                 "finalize code calls deepCopyFinalize")
    assert_equal(mod, 'myMod',
                 "module attribute returned as the 4th tuple element")


def test_functionclass_linkedlist_state_store():
    print("\n=== Testing FunctionClass.LinkedList.state_store_linked_list ===")
    ll = _linked_list_spec()
    inp, out, mod = state_store_linked_list(
        {'linkedList': ll}, {}, [])
    assert_equal('stateRestore(stateFile,gslStateFile,stateOperationID)' in inp,
                 True, "input (restore) code calls %stateRestore(…)")
    assert_equal('stateStore(stateFile,gslStateFile,stateOperationID)' in out,
                 True, "output (store) code calls %stateStore(…)")
    assert_equal(mod, 'myMod', "module passed through")


def test_functionclass_linkedlist_allowed_and_descriptor_and_assigner():
    print("\n=== Testing FunctionClass.LinkedList (allowed / descriptor / assigner) ===")
    ll = _linked_list_spec()
    iter_allowed, _ = allowed_parameters_linked_list(
        {'linkedList': ll}, [], 'mySrc')
    assert_equal("'mySrc'" in iter_allowed, True,
                 "allowed_parameters passes source string literal")
    iter_desc, _ = auto_descriptor_linked_list(ll, [])
    assert_equal('descriptor(parameters)' in iter_desc, True,
                 "auto_descriptor calls %descriptor(parameters)")
    iter_asn, _ = assigner_linked_list(ll, [])
    assert_equal('referenceCountIncrement()' in iter_asn, True,
                 "assigner bumps reference count per member")
    assert_equal('nullify(self%first)' in iter_asn, True,
                 "assigner nullifies the destination list head first")


# ============================================================================
# Process.FunctionClass.DeepCopy tests
# ============================================================================

def test_functionclass_deepcopy_copied_self_block():
    print("\n=== Testing FunctionClass.DeepCopy.deep_copy_copied_self_block ===")
    dc = {}
    deep_copy_copied_self_block(
        dc, 'obj',
        {'intrinsic': 'class', 'type': 'myType'},
        {'node': {'line': 5}},
        indent='  ',
    )
    code = dc['assignments']
    assert_equal('  if (associated(self%obj%copiedSelf)) then' in code, True,
                 "indent applied to all lines")
    assert_equal('class is (myType)' in code, True,
                 "intrinsic + type spliced into class-is clause")
    assert_equal('call self%obj%copiedSelf%referenceCountIncrement()' in code,
                 True, "reference-count increment call emitted")
    assert_equal('allocate(destination%obj,mold=self%obj)' in code, True,
                 "mold-based allocate on else branch emitted")


def test_functionclass_deepcopy_generate_assignment_allocatable():
    print("\n=== Testing FunctionClass.DeepCopy.generate_assignment_allocatable_code ===")
    asn = {}
    rank_max = [0]
    generate_assignment_allocatable_code(
        asn,
        {'intrinsic': 'integer', 'attributes': []},
        'v', 'allocated', rank_max)
    assert_equal('if (allocated(self%v)) deallocate(self%v)' in asn['code'], True,
                 "allocated self guarded deallocate")
    assert_equal('self%v=from%v' in asn['code'], True,
                 "simple scalar assignment for non-type")

    asn = {}
    rank_max = [0]
    generate_assignment_allocatable_code(
        asn,
        {'intrinsic': 'type', 'type': 'thing',
         'attributes': ['dimension(:,:)']},
        'v', 'allocated', rank_max)
    assert_equal('do i1__=lbound(from%v,dim=1)' in asn['code'], True,
                 "rank-2 type assignment uses lbound/ubound loop")
    assert_equal('self%v(i1__,i2__)=from%v(i1__,i2__)' in asn['code'], True,
                 "element-wise assignment body emitted")
    assert_equal(rank_max[0], 2,
                 "rank-maximum-ref updated when rank exceeds previous max")


def test_functionclass_deepcopy_declarations_functionclass_pointer():
    print("\n=== Testing FunctionClass.DeepCopy.deep_copy_declarations (functionClass pointer) ===")
    state_storables = {
        'functionClasses': {'functionClass': [{'name': 'myType'}]},
        'functionClassInstances': [],
    }
    declarations = [
        {'intrinsic': 'class', 'type': 'myType',
         'attributes': ['pointer'],
         'variables': ['obj'], 'variableNames': ['obj']},
    ]
    deep_copy = {}
    found = []
    deep_copy_declarations(
        {},                           # class_record
        {'node': {'line': 1}},       # non_abstract_class
        None,                        # node (unused by this path)
        declarations, [], 1,
        deep_copy, found,
        state_storables, {},
    )
    assert_equal(deep_copy.get('needReferenceCount'), 1,
                 "needReferenceCount flag set for functionClass pointer")
    assert_equal('call self%obj%deepCopy(destination%obj)' in deep_copy['assignments'],
                 True, "deepCopy call emitted for the pointer member")
    assert_equal('call self%obj%deepCopyReset   ()' in deep_copy['resetCode'],
                 True, "reset code emitted")
    assert_equal('call self%obj%deepCopyFinalize()' in deep_copy['finalizeCode'],
                 True, "finalize code emitted")


def test_functionclass_deepcopy_declarations_explicit_actions():
    print("\n=== Testing FunctionClass.DeepCopy.deep_copy_declarations (explicit deepCopyActions) ===")
    actions = {'deepCopyActions': [{'type': 'myThing'}]}
    declarations = [
        {'intrinsic': 'type', 'type': 'myThing',
         'attributes': ['dimension(:)', 'allocatable'],
         'variables': ['items'], 'variableNames': ['items']},
    ]
    deep_copy = {}
    deep_copy_declarations(
        {}, {'node': {'line': 1}}, None,
        declarations, [], 1, deep_copy, [], {}, actions,
    )
    code = deep_copy['assignments']
    assert_equal('if (allocated(self%items)) then' in code, True,
                 "allocatable guard emitted for type with deepCopyActions")
    assert_equal('do i1=lbound(self%items,dim=1)' in code, True,
                 "rank-1 do-loop emitted")
    assert_equal('call destination%items(i1)%deepCopyActions()' in code, True,
                 "deepCopyActions call emitted per element")


def test_functionclass_deepcopy_declarations_hdf5_and_setto():
    print("\n=== Testing FunctionClass.DeepCopy.deep_copy_declarations (HDF5 + setTo) ===")
    state_storables = {'functionClasses': {}, 'functionClassInstances': []}
    declarations = [
        {'intrinsic': 'type', 'type': 'hdf5Object',
         'attributes': [],
         'variables': ['obj'], 'variableNames': ['obj']},
        {'intrinsic': 'integer', 'type': None, 'attributes': [],
         'variables': ['host'], 'variableNames': ['host']},
    ]
    class_record = {
        'deepCopy': {
            'setTo': {'variables': 'host%counter', 'value': '0'},
        },
    }
    deep_copy = {}
    deep_copy_declarations(
        class_record, {'node': {'line': 1}}, None,
        declarations, [], 1, deep_copy, [], state_storables, {},
    )
    code = deep_copy['assignments']
    assert_equal('HDF5_Access' in deep_copy['modules'], True,
                 "HDF5_Access module recorded for hdf5object member")
    assert_equal('call self%obj%deepCopy(destination%obj)' in code, True,
                 "hdf5 member deep-copy call emitted")
    assert_equal('destination%host%counter=0' in code, True,
                 "setTo directive emits destination assignment")


# ============================================================================
# Process.FunctionClass.StateStore tests
# ============================================================================

def test_functionclass_statestore_explicit_function():
    print("\n=== Testing FunctionClass.StateStore.state_store_explicit_function ===")
    inp, out, mods = state_store_explicit_function({
        'stateStore': {'stateStore': {
            'variables': 'p q',
            'store':     'myStore',
            'restore':   'myRestore',
            'module':    'myMod',
        }},
    })
    assert_equal('call myStore(self%p,stateFile,gslStateFile,stateOperationID)'
                 in out, True, "output code uses explicit store function")
    assert_equal('call myRestore(self%p,stateFile,gslStateFile,stateOperationID)'
                 in inp, True, "input code uses explicit restore function")
    assert_equal('nullify(self%q)' in inp, True,
                 "input code nullifies on not-associated branch")
    assert_equal('myMod' in mods, True, "module recorded in modules dict")


def test_functionclass_statestore_generate_allocatable():
    print("\n=== Testing FunctionClass.StateStore.generate_allocatable_state_store_code ===")
    ss = {}
    sss = {}
    generate_allocatable_state_store_code(
        ss, sss, 'arr', 3, 'self%arr', 'self%arr')
    assert_equal(sss['allocatablesFound'], True, "allocatablesFound flipped")
    assert_equal(sss['dimensionalsFound'], True, "dimensionalsFound flipped")
    assert_equal(sss['stateFileUsed'],     True, "stateFileUsed flipped")
    assert_equal(sss['labelUsed'],         True, "labelUsed flipped")
    assert_equal('shape(self%arr,kind=c_size_t)' in ss['outputCode'], True,
                 "shape() line emitted")
    assert_equal('allocate(storedShape(3))' in ss['inputCode'], True,
                 "rank-3 storedShape allocation emitted")
    assert_equal('allocate(self%arr(storedShape(1),storedShape(2),storedShape(3)))'
                 in ss['inputCode'], True,
                 "rank-3 shaped allocate emitted on restore")


def test_functionclass_statestore_variables_static_and_restore():
    print("\n=== Testing FunctionClass.StateStore.state_store_variables (static + restoreTo) ===")
    declarations = [
        {'intrinsic': 'integer', 'type': None, 'attributes': [],
         'variables': ['a'], 'variableNames': ['a']},
        {'intrinsic': 'real', 'type': None, 'attributes': [],
         'variables': ['b'], 'variableNames': ['b']},
    ]
    class_record = {
        'stateStorable': {
            'restoreTo': {'variables': 'b', 'state': '0.0'},
        },
    }
    state_store = {}
    state_stores = {}
    state_store_variables(state_stores, state_store, class_record,
                          declarations, [], {})
    assert_equal(state_store['staticVariables'], ['a'],
                 "non-restoreTo static variables collected in order")
    assert_equal('self%b=0.0' in state_store.get('inputCode', ''), True,
                 "restoreTo directive writes state into inputCode")


def test_functionclass_statestore_custom_hooks_flag():
    print("\n=== Testing FunctionClass.StateStore.state_store_variables (custom hooks) ===")
    declarations = [
        {'intrinsic': 'procedure',
         'variables': ['stateStore=>myStore']},
        {'intrinsic': 'procedure',
         'variables': ['stateRestore=>myRestore']},
    ]
    state_store = {}
    state_store_variables({}, state_store, {}, declarations, [], {})
    assert_equal(state_store.get('hasCustomStateStore'),   True,
                 "hasCustomStateStore flag set on stateStore=>… binding")
    assert_equal(state_store.get('hasCustomStateRestore'), True,
                 "hasCustomStateRestore flag set on stateRestore=>… binding")


def test_functionclass_statestore_variables_allocatable_intrinsic():
    print("\n=== Testing FunctionClass.StateStore.state_store_variables (allocatable intrinsic) ===")
    declarations = [
        {'intrinsic': 'integer', 'type': None,
         'attributes': ['allocatable', 'dimension(:)'],
         'variables': ['items'], 'variableNames': ['items']},
    ]
    ss = {}
    sss = {}
    state_store_variables(sss, ss, None, declarations, [], {})
    assert_equal('allocate(storedShape(1))' in ss['inputCode'], True,
                 "allocatable intrinsic triggers storedShape(1) allocation")
    assert_equal(sss.get('allocatablesFound'), True,
                 "allocatablesFound flag set")


# ============================================================================
# Process.FunctionClass scaffolding tests (PR D.7.3a)
# ============================================================================

def _reset_fc_caches():
    _fc_state_storables_holder['value']     = None
    _fc_deep_copy_actions_holder['value']   = None
    _fc_directive_locations_holder['value'] = None


def test_functionclass_is_function_class_pointer():
    print("\n=== Testing FunctionClass._is_function_class_pointer ===")
    state_storables = {
        'functionClasses': {'functionClass': [{'name': 'myFooClass'}]},
        'functionClassInstances': ['myFooCore'],
    }
    # `class(myFooClass), pointer` → True
    assert_equal(_is_function_class_pointer(
        {'intrinsic': 'class', 'type': 'myFooClass', 'attributes': ['pointer']},
        'myFooClass', state_storables,
    ), True, "class + pointer + registered → True")
    # instance match
    assert_equal(_is_function_class_pointer(
        {'intrinsic': 'type', 'type': 'myFooCore', 'attributes': ['pointer']},
        'myFooCore', state_storables,
    ), True, "type + pointer + functionClassInstance → True")
    # missing pointer attribute
    assert_equal(_is_function_class_pointer(
        {'intrinsic': 'class', 'type': 'myFooClass', 'attributes': []},
        'myFooClass', state_storables,
    ), False, "no pointer attribute → False")
    # wrong intrinsic
    assert_equal(_is_function_class_pointer(
        {'intrinsic': 'integer', 'type': 'myFooClass',
         'attributes': ['pointer']},
        'myFooClass', state_storables,
    ), False, "non-class/type intrinsic → False")
    # type unknown
    assert_equal(_is_function_class_pointer(
        {'intrinsic': 'class', 'type': 'elsewhere',
         'attributes': ['pointer']},
        'elsewhere', state_storables,
    ), False, "type not in state storables → False")


def test_functionclass_init_code_content():
    print("\n=== Testing FunctionClass._init_code_content ===")
    code, pre, post = _init_code_content()
    assert_equal(isinstance(code['module']['preContains'],  list),  True,
                 "preContains is a list")
    assert_equal(isinstance(code['module']['postContains'], list),  True,
                 "postContains is a list")
    assert_equal(code['module']['preContains'][0]  is pre,           True,
                 "first preContains entry is the returned pre-node")
    assert_equal(code['module']['postContains'][0] is post,          True,
                 "first postContains entry is the returned post-node")
    assert_equal(code['submodule'], {},
                 "submodule container starts empty")


def test_functionclass_build_base_method_stubs():
    print("\n=== Testing FunctionClass._build_base_method_stubs ===")
    methods = {}
    # Directive with autoHook override + destructor.
    directive = {
        'autoHook': {
            'modules': [{'name': 'MyMod', 'only': 'foo, bar'}],
            'code':    '! custom autoHook body\n',
        },
        'destructor': {
            'modules': [{'name': 'DMod', 'only': 'baz'}],
            'code':    '! custom destructor body\n',
        },
    }
    _build_base_method_stubs(directive, methods)
    assert_equal('stateStore' in methods and 'stateRestore' in methods, True,
                 "stateStore + stateRestore stubs populated")
    assert_equal(methods['stateStore']['modules'][0]['name'], 'ISO_C_Binding',
                 "stateStore imports ISO_C_Binding, only : c_ptr")
    assert_equal(methods['autoHook']['code'], '! custom autoHook body\n',
                 "directive autoHook code overrides default stub")
    mod = next(m for m in methods['autoHook']['modules']
               if m['name'] == 'MyMod')
    assert_equal(mod['only'], ['foo', 'bar'],
                 "autoHook `only` list comma-split correctly")
    assert_equal('destructor' in methods, True,
                 "destructor added when directive requests it")
    assert_equal(methods['destructor']['code'], '! custom destructor body\n',
                 "destructor carries the directive-supplied code")


def test_functionclass_build_object_type_method():
    print("\n=== Testing FunctionClass._build_object_type_method ===")
    methods = {}
    _build_object_type_method(
        {'name': 'myFoo'},
        [
            {'name': 'myFooCore',    'shortName': 'core'},
            {'name': 'myFooAdvanced','shortName': 'advanced'},
            {'name': 'myFooXYZ',     'shortName': 'XYZ'},
        ],
        methods,
    )
    code = methods['objectType']['code']
    assert_equal('type is (myFooCore)' in code, True,
                 "select-type branch emitted per non-abstract class")
    assert_equal("myFooObjectType='core'" in code, True,
                 "short-name branch uses lcfirst for typical names")
    assert_equal("myFooObjectType='XYZ'" in code, True,
                 "short-name keeps acronyms with 2+ caps intact")
    assert_equal("myFooObjectType='myFooAdvanced'" in code, True,
                 "long-name branch returns full class name")
    assert_equal(methods['objectType']['modules'], 'ISO_Varying_String',
                 "objectType method imports ISO_Varying_String")


def test_functionclass_load_and_sort_classes(tmp_source_write=None):
    print("\n=== Testing FunctionClass._load_and_sort_classes ===")
    tmpdir = tempfile.mkdtemp()
    try:
        # Two classes: base myFooCore extends myFoo; derived myFooAdv extends
        # myFooCore.  Alphabetical sort would be [Adv, Core] but the topo
        # sort must place Core first.  We write two source files.
        core_path = os.path.join(tmpdir, 'core.F90')
        adv_path  = os.path.join(tmpdir, 'adv.F90')
        with open(core_path, 'w') as fh:
            fh.write(
                "module core_mod\n"
                "  !![\n"
                "  <myFoo>\n"
                "   <name>myFooCore</name>\n"
                "  </myFoo>\n"
                "  !!]\n"
                "  type, extends(myFooClass) :: myFooCore\n"
                "  end type myFooCore\n"
                "end module core_mod\n"
            )
        with open(adv_path, 'w') as fh:
            fh.write(
                "module adv_mod\n"
                "  !![\n"
                "  <myFoo>\n"
                "   <name>myFooAdv</name>\n"
                "  </myFoo>\n"
                "  !!]\n"
                "  type, extends(myFooCore) :: myFooAdv\n"
                "  end type myFooAdv\n"
                "end module adv_mod\n"
            )
        # stateStorables + deepCopyActions stubs so process_tree doesn't fail.
        # The stateStorables entry for `myFooClass` causes the functionClass
        # hook's early pass to mark each `<myFoo>` directive as processed
        # (mirroring what happens in a real Galacticus build where
        # stateStorables.xml is populated before FunctionClass runs).
        with open(os.path.join(tmpdir, 'stateStorables.xml'), 'w') as fh:
            fh.write(
                "<stateStorables>\n"
                "  <functionClasses>\n"
                "    <functionClass name=\"myFooClass\"/>\n"
                "  </functionClasses>\n"
                "</stateStorables>\n"
            )
        with open(os.path.join(tmpdir, 'directiveLocations.xml'), 'w') as fh:
            fh.write("<directiveLocations></directiveLocations>\n")
        open(os.path.join(tmpdir, 'deepCopyActions.xml'), 'w').close()
        open(os.path.join(tmpdir, 'hdf5FCInterop.dat'), 'w').close()
        with open(os.path.join(tmpdir, 'Makefile_All_Execs'), 'w') as fh:
            fh.write("all_exes =\n")

        saved = os.environ.get('BUILDPATH')
        os.environ['BUILDPATH'] = tmpdir
        try:
            _reset_fc_caches()
            directive = {'name': 'myFoo'}
            directive_locations = {
                'myFoo': {'file': [adv_path, core_path]},
            }
            classes, ordered, non_abstract = _load_and_sort_classes(
                directive, directive_locations)
            assert_equal(set(classes.keys()), {'myFooCore', 'myFooAdv'},
                         "all classes discovered")
            # Perl's loadAndSortClasses builds `dep[parent] += [child]` and
            # feeds that to Sort::Topo, whose docstring says `dep[X]` is
            # "things X depends on".  The emitted order therefore puts the
            # extending (child) class BEFORE its parent — a deliberate
            # (and Perl-matching) consequence of that edge direction.
            assert_equal([c['type'] for c in ordered],
                         ['myFooAdv', 'myFooCore'],
                         "topo-sort matches Perl: Adv before Core")
            assert_equal(
                [c['shortName'] for c in non_abstract],
                ['adv', 'core'],
                "short names derived with lcfirst, same order as classes",
            )
        finally:
            if saved is None:
                del os.environ['BUILDPATH']
            else:
                os.environ['BUILDPATH'] = saved
            _reset_fc_caches()
    finally:
        import shutil
        shutil.rmtree(tmpdir)


def test_functionclass_load_and_sort_classes_default_validation():
    print("\n=== Testing FunctionClass._load_and_sort_classes (invalid default) ===")
    tmpdir = tempfile.mkdtemp()
    core_path = os.path.join(tmpdir, 'core.F90')
    with open(core_path, 'w') as fh:
        fh.write(
            "module core_mod\n"
            "  !![\n"
            "  <myFoo><name>myFooCore</name></myFoo>\n"
            "  !!]\n"
            "  type, extends(myFooClass) :: myFooCore\n"
            "  end type myFooCore\n"
            "end module core_mod\n"
        )
    for f in ('stateStorables.xml', 'directiveLocations.xml',
              'deepCopyActions.xml', 'hdf5FCInterop.dat'):
        open(os.path.join(tmpdir, f), 'w').close()
    with open(os.path.join(tmpdir, 'stateStorables.xml'), 'w') as fh:
        fh.write("<stateStorables></stateStorables>\n")
    with open(os.path.join(tmpdir, 'directiveLocations.xml'), 'w') as fh:
        fh.write("<directiveLocations></directiveLocations>\n")
    with open(os.path.join(tmpdir, 'Makefile_All_Execs'), 'w') as fh:
        fh.write("all_exes =\n")

    saved = os.environ.get('BUILDPATH')
    os.environ['BUILDPATH'] = tmpdir
    try:
        _reset_fc_caches()
        directive_locations = {'myFoo': {'file': core_path}}
        assert_raises(
            lambda: _load_and_sort_classes(
                {'name': 'myFoo', 'default': 'nonexistent'},
                directive_locations,
            ),
            RuntimeError,
            "invalid default value raises RuntimeError")
    finally:
        if saved is None:
            del os.environ['BUILDPATH']
        else:
            os.environ['BUILDPATH'] = saved
        _reset_fc_caches()
        import shutil
        shutil.rmtree(tmpdir)


def test_functionclass_process_hook_is_no_op_when_no_directive():
    print("\n=== Testing process_function_class is a no-op on trees without functionClass ===")
    root = _parse_text(
        "subroutine s\n"
        "integer :: counter\n"
        "end subroutine s\n"
    )
    before = serialize(root)
    tmpdir = tempfile.mkdtemp()
    with open(os.path.join(tmpdir, 'stateStorables.xml'), 'w') as fh:
        fh.write("<stateStorables></stateStorables>\n")
    saved = os.environ.get('BUILDPATH')
    os.environ['BUILDPATH'] = tmpdir
    try:
        _reset_fc_caches()
        _fc_process_function_class(root, {})
        assert_equal(serialize(root), before,
                     "tree without <functionClass> unchanged by the hook")
    finally:
        if saved is None:
            del os.environ['BUILDPATH']
        else:
            os.environ['BUILDPATH'] = saved
        _reset_fc_caches()
        import shutil
        shutil.rmtree(tmpdir)


def test_functionclass_process_hook_marks_class_nodes():
    print("\n=== Testing process_function_class early pass marks <class> directive nodes ===")
    # Synthesise a tree with a node of type 'myFoo' whose parent is a
    # module, AND make sure `myFooClass` is registered as a functionClass
    # in stateStorables.xml.  The hook's early pass should flip the
    # directive's `processed` flag without requiring a <functionClass>
    # directive elsewhere in the tree.
    root = _parse_text("module m\nend module m\n")
    mod = _find_node(root, 'module')
    directive_node = _make_directive_node(
        'myFoo', {'processed': False}, mod,
    )
    tmpdir = tempfile.mkdtemp()
    with open(os.path.join(tmpdir, 'stateStorables.xml'), 'w') as fh:
        fh.write(
            "<stateStorables>\n"
            "  <functionClasses>\n"
            "    <functionClass name=\"myFooClass\"/>\n"
            "  </functionClasses>\n"
            "</stateStorables>\n"
        )
    saved = os.environ.get('BUILDPATH')
    os.environ['BUILDPATH'] = tmpdir
    try:
        _reset_fc_caches()
        _fc_process_function_class(root, {})
        assert_equal(directive_node['directive']['processed'], True,
                     "early pass flips processed=True on class-instance nodes")
    finally:
        if saved is None:
            del os.environ['BUILDPATH']
        else:
            os.environ['BUILDPATH'] = saved
        _reset_fc_caches()
        import shutil
        shutil.rmtree(tmpdir)


# ============================================================================
# Process.FunctionClass body-emitter tests (PR D.7.3b)
# ============================================================================

def test_functionclass_format_variable_definitions():
    print("\n=== Testing FunctionClass._format_variable_definitions ===")
    text = _format_variable_definitions([
        {'intrinsic': 'integer', 'type': None,
         'attributes': ['intent(in)'],
         'variables': ['n']},
        {'intrinsic': 'real', 'type': 'kind=8',
         'attributes': ['intent(inout)', 'dimension(:)'],
         'variables': ['x', 'y']},
    ])
    assert_equal('integer, intent(in) :: n' in text, True,
                 "simple integer declaration rendered")
    assert_equal(
        'real(kind=8), intent(inout), dimension(:) :: x, y' in text,
        True, "real with kind + multi-attr + multi-var rendered")


def test_functionclass_source_digest_binding():
    print("\n=== Testing FunctionClass._source_digest_binding ===")
    out = _source_digest_binding('myClass')
    assert_equal(
        'character(C_Char), dimension(23), bind(C, name="myClassMD5") '
        ':: myClass5' in out, True,
        "C-binding declaration matches Perl SourceDigest::Binding format")


def test_functionclass_generate_type_definition():
    print("\n=== Testing FunctionClass._generate_type_definition ===")
    root = _parse_text("module m\nend module m\n")
    mod = _find_node(root, 'module')
    # Create a fake node+pre; we don't go through process_function_class()
    # here, we just invoke the generator directly.
    pre = {'type': 'code', 'content': '', 'parent': None,
           'firstChild': None, 'sibling': None, 'source': 'x', 'line': 1}
    fc_directive = {
        'name':        'myFoo',
        'description': 'Just a test class.',
        'data': [
            {'scope': 'self',   'content': '    integer :: counter=0'},
            {'scope': 'module', 'content': '    integer :: globalFlag=0',
             'threadprivate': 'yes'},
        ],
        'generic': [
            {'name': 'evaluate', 'method': ['evaluateRank0']},
        ],
    }
    methods = {
        'evaluateRank0': {
            'type':        'double precision',
            'description': 'Evaluate at a single point.',
            'argument':    ['double precision, intent(in   ) :: x'],
        },
        'destructor': {'type': 'void'},
    }
    directive_node = _make_directive_node('myFoo', fc_directive, mod)
    _generate_type_definition(fc_directive, methods, pre, directive_node)
    text = pre['content']
    assert_equal('type, extends(functionClass) :: myFooClass' in text, True,
                 "type declaration emitted with default functionClass parent")
    assert_equal("stateOperationID=0" in text, True,
                 "stateOperationID declaration emitted")
    assert_equal('integer :: counter=0' in text, True,
                 "self-scope <data> entry inlined inside type")
    assert_equal('integer :: globalFlag=0' in text, True,
                 "module-scope <data> entry emitted after end type")
    assert_equal('!$omp threadprivate(globalFlag)' in text, True,
                 "threadprivate entry emitted for module-scope <data>")
    assert_equal('procedure :: evaluateRank0 => myFooEvaluateRank0__' in text,
                 True, "method binding line emitted with __ (no code given)")
    assert_equal('generic :: evaluate => evaluateRank0' in text, True,
                 "generic interface binding emitted")
    assert_equal('final :: myFooDestructor' in text, True,
                 "final binding emitted when destructor is present")
    assert_equal('integer :: myFooDsblVldtn=0' in text, True,
                 "disable-validation state variable emitted")


def test_functionclass_generate_constructor_with_default_and_recursion():
    print("\n=== Testing FunctionClass._generate_constructor (default + recursion) ===")
    root = _parse_text("module m\nend module m\n")
    mod = _find_node(root, 'module')
    pre = {'type': 'code', 'content': '', 'parent': None,
           'firstChild': None, 'sibling': None, 'source': 'x', 'line': 1}
    post = {'type': 'code', 'content': '', 'parent': None,
            'firstChild': None, 'sibling': None, 'source': 'x', 'line': 1}
    fc_directive = {'name': 'myFoo', 'default': 'core'}
    classes_ordered = [
        {'name': 'myFooCore', 'shortName': 'core', 'recursive': 'yes'},
        {'name': 'myFooAdv',  'shortName': 'adv'},
    ]
    non_abstract = classes_ordered
    directive_node = _make_directive_node('myFoo', fc_directive, mod)
    # Top-level tree is the fake root; pretend it's a file-typed node.
    fake_tree = {'type': 'file', 'name': 'test.F90'}
    tmpdir = tempfile.mkdtemp()
    saved = os.environ.get('BUILDPATH')
    os.environ['BUILDPATH'] = tmpdir
    try:
        _generate_constructor(fc_directive, classes_ordered, non_abstract,
                              pre, post, directive_node, fake_tree)
    finally:
        if saved is None:
            del os.environ['BUILDPATH']
        else:
            os.environ['BUILDPATH'] = saved
        import shutil
        shutil.rmtree(tmpdir)

    assert_equal('interface myFoo' in pre['content'], True,
                 "constructor interface declared")
    assert_equal('RecursiveBuildNode' in pre['content'], True,
                 "recursive-build state variable declared when any class is recursive")
    body = post['content']
    assert_equal('function myFooCnstrctrPrmtrs(parameters,copyInstance,'
                 'parameterName) result(self)' in body, True,
                 "constructor function declared with all expected parameters")
    assert_equal("parameterName_='myFoo'" in body, True,
                 "default parameterName falls back to directive name")
    assert_equal("case ('core')" in body, True,
                 "select case branch emitted per non-abstract class")
    assert_equal('self=myFooCore(subParameters)' in body, True,
                 "per-class self-construction emitted")
    assert_equal('message=message//char(10)' in body, True,
                 "case-default error message lists alternatives")


def test_functionclass_generate_method_functions():
    print("\n=== Testing FunctionClass._generate_method_functions ===")
    post = {'type': 'code', 'content': '', 'parent': None,
            'firstChild': None, 'sibling': None, 'source': 'x', 'line': 1}
    node = {'type': 'functionClass', 'line': 0, 'source': 't.F90',
            'parent': None}
    fc_directive = {'name': 'myFoo'}
    methods = {
        # void method with explicit code
        'tick': {
            'type': 'void',
            'pass': 'yes',
            'code': 'call advance(self)',
            'argument': ['integer, intent(in) :: steps'],
        },
        # function stub — no code provided; method returns integer
        'count': {
            'type': 'integer',
            'pass': 'yes',
        },
        # function stub returning class pointer
        'child': {
            'type':        'class(myFooClass)',
            'pass':        'yes',
            'description': 'Return a child pointer.',
        },
        # skip — has explicit function override
        'overridden': {'type': 'void', 'function': 'myImpl'},
    }
    _generate_method_functions(fc_directive, methods, post, node)
    body = post['content']
    assert_equal('subroutine myFooTick(self,steps)' in body, True,
                 "subroutine stub emitted for void method with code")
    assert_equal('call advance(self)' in body, True,
                 "explicit <code> block spliced into body")
    assert_equal('integer function myFooCount__' in body, True,
                 "function stub uses __ suffix when no code is given")
    assert_equal('myFooCount__=0' in body, True,
                 "default initial return value emitted for integer")
    assert_equal('function myFooChild__' in body, True,
                 "function stub for class-pointer return type emitted")
    assert_equal('myFooChild__=> null()' in body, True,
                 "class-pointer default is => null()")
    assert_equal('overridden' not in body, True,
                 "methods with explicit function= are skipped")


def test_functionclass_generate_documentation_basic():
    print("\n=== Testing FunctionClass._generate_documentation (basic write) ===")
    # Set up a minimal tree for one class that has an interface block naming
    # a module-procedure constructor — the walker needs both to emit the
    # class's `\\subsection{}` header without parameter details.
    class_src = (
        "module m\n"
        "  interface myFooCore\n"
        "    module procedure myFooCoreConstructor\n"
        "  end interface myFooCore\n"
        "end module m\n"
    )
    tree = _parse_text(class_src)
    classes = {
        'myFooCore': {
            'name':        'myFooCore',
            'description': 'The core implementation.',
            'tree':        tree,
            'extends':     'myFoo',
        },
    }
    non_abstract_classes = [classes['myFooCore']]
    directive = {
        'name':            'myFoo',
        'descriptiveName': 'My Foo Example',
        'description':     'A toy functionClass for unit tests.',
    }

    cwd = os.getcwd()
    tmpdir = tempfile.mkdtemp()
    try:
        os.chdir(tmpdir)
        _generate_documentation(directive, classes, non_abstract_classes)
        out_path = os.path.join(
            tmpdir, 'doc', 'physics', 'my_foo_example.tex')
        assert_equal(os.path.exists(out_path), True,
                     "documentation file written at expected path")
        with open(out_path, 'r') as fh:
            payload = fh.read()
        assert_equal('\\section{My Foo Example}' in payload, True,
                     "\\section header with descriptive name")
        assert_equal('\\refClass{myFooCore}' in payload, True,
                     "\\refClass reference to the implementation emitted")
        assert_equal('The core implementation.' in payload, True,
                     "class description spliced into subsection body")
    finally:
        os.chdir(cwd)
        import shutil
        shutil.rmtree(tmpdir)


# ============================================================================
# Main
# ============================================================================

# ---------------------------------------------------------------------------
# D.7.3c build_* method tests
# ---------------------------------------------------------------------------

def _make_simple_class(name, extends, member_declarations=None,
                       constructor_children=None, opener=None):
    """Helper: build a minimal class_record tree with a type block and
    an optional interface+constructor function.

    `member_declarations` populates the type body; `constructor_children`
    are siblings attached after a type(inputParameters) :: parameters
    declaration inside a function named <name>Constructor.
    """
    type_body_sibling = None
    if member_declarations:
        type_body_sibling = {
            'type': 'declaration',
            'declarations': member_declarations,
            'sibling':    None,
            'firstChild': None,
        }
    type_node = {
        'type':       'type',
        'name':       name,
        'opener':     opener or f'type, extends({extends}) :: {name}',
        'firstChild': type_body_sibling,
        'sibling':    None,
    }
    if constructor_children is not None:
        first = {
            'type':         'declaration',
            'declarations': [{
                'intrinsic':  'type',
                'type':       'inputParameters',
                'variables':  ['parameters'],
                'attributes': [],
            }],
            'sibling':   None,
            'firstChild': None,
        }
        prev = first
        for child in constructor_children:
            prev['sibling'] = child
            prev = child
            prev['sibling'] = None
        func = {
            'type':       'function',
            'name':       name + 'Constructor',
            'opener':     f'  function {name}Constructor(parameters)',
            'firstChild': first,
            'sibling':    None,
        }
        iface_body = {
            'type':       'moduleProcedure',
            'names':      [name + 'Constructor'],
            'firstChild': None,
            'sibling':    None,
        }
        iface = {
            'type':       'interface',
            'name':       name,
            'firstChild': iface_body,
            'sibling':    func,
        }
        type_node['sibling'] = iface
    return {
        'name':    name,
        'extends': extends,
        'node':    {'line': 1, 'source': f'{name}.F90'},
        'tree':    {'firstChild': type_node},
    }


def test_functionclass_build_assignment_method():
    print("\n=== Testing FunctionClass._build_assignment_method ===")
    directive = {'name': 'testFoo', 'data': []}
    cls = _make_simple_class(
        'testFooSimple', 'testFoo',
        member_declarations=[
            {'intrinsic':  'integer', 'type': None, 'attributes': [],
             'variables': ['x']},
            {'intrinsic':  'class', 'type': 'testFoo',
             'attributes': ['pointer'], 'variables': ['obj']},
        ])
    methods = {}
    state_storables = {
        'functionClasses':        {'functionClass': [{'name': 'testFoo'}]},
        'functionClassInstances': [],
    }
    _build_assignment_method(
        directive, [cls], {'testFooSimple': cls}, methods, state_storables)
    code = methods['assignment(=)']['code']
    assert_equal('select type (self)' in code, True,
                 "outer select-type frame emitted")
    assert_equal('type is (testFooSimple)' in code, True,
                 "per-class type guard emitted")
    assert_equal('self%x=from%x' in code, True,
                 "scalar member direct-assigned with =")
    assert_equal('self%obj=>from%obj' in code, True,
                 "pointer member associated with =>")
    assert_equal(
        'if (associated(self%obj)) call self%obj%referenceCountIncrement()'
        in code, True,
        "reference-count increment guarded by associated()")
    assert_equal("class default" in code and
                 "call Error_Report('self and from types do not match'" in code,
                 True, "class-default branch emits Error_Report")
    assert_equal('self%isDefaultOfClass=from%isDefaultOfClass' in code, True,
                 "functionClass base-class fields assigned")
    assert_equal(methods['assignment(=)']['selfIntent'], 'out',
                 "assignment operator uses selfIntent=out")


def test_functionclass_build_deep_copy_methods():
    print("\n=== Testing FunctionClass._build_deep_copy_methods ===")
    directive = {'name': 'testFoo', 'data': []}
    cls = _make_simple_class(
        'testFooSimple', 'testFoo',
        member_declarations=[
            {'intrinsic':  'integer', 'type': None, 'attributes': [],
             'variables': ['x']},
            {'intrinsic':  'class', 'type': 'testFoo',
             'attributes': ['pointer'], 'variables': ['obj']},
        ])
    methods = {}
    state_storables = {
        'functionClasses':        {'functionClass': [{'name': 'testFoo'}]},
        'functionClassInstances': [],
    }
    _build_deep_copy_methods(
        directive, [cls], {'testFooSimple': cls}, 42, methods,
        state_storables, {})
    assert_equal(methods['deepCopy']['code'],
                 'call self%deepCopy_(destination)',
                 "deepCopy wrapper delegates to deepCopy_")
    assert_equal('integer :: referenceCount__' in methods['deepCopy_']['code'],
                 True, "referenceCount__ declared when pointer member present")
    assert_equal('nullify(destination%obj)'
                 in methods['deepCopy_']['code'], True,
                 "deep-copy nullifies destination pointer member")
    assert_equal('call self%obj%deepCopy(destination%obj)'
                 in methods['deepCopy_']['code'], True,
                 "deep-copy invokes inner deepCopy on pointer member")
    assert_equal('self%copiedSelf => null()'
                 in methods['deepCopyReset']['code'], True,
                 "deepCopyReset starts with copiedSelf nullification")
    assert_equal(
        'call self%obj%deepCopyReset   ()' in methods['deepCopyReset']['code'],
        True, "deepCopyReset delegates to pointer member")
    assert_equal(
        'call self%obj%deepCopyFinalize()'
        in methods['deepCopyFinalize']['code'], True,
        "deepCopyFinalize delegates to pointer member")
    assert_equal('Error' in methods['deepCopy_']['modules'].split(), True,
                 "deepCopy_ modules include Error")


def test_functionclass_build_state_store_methods():
    print("\n=== Testing FunctionClass._build_state_store_methods ===")
    directive = {'name': 'testFoo', 'data': []}
    cls = _make_simple_class(
        'testFooSimple', 'testFoo',
        member_declarations=[
            {'intrinsic':  'integer', 'type': None, 'attributes': [],
             'variables': ['x'], 'variableNames': ['x']},
            {'intrinsic':  'double precision', 'type': None,
             'attributes': [], 'variables': ['y'], 'variableNames': ['y']},
        ])
    methods = {}
    state_storables = {
        'functionClasses':        {'functionClass': [{'name': 'testFoo'}]},
        'functionClassInstances': [],
        'stateStorables':         [],
    }
    _build_state_store_methods(
        directive, [cls], {'testFooSimple': cls}, methods, state_storables)
    assert_equal(methods['stateStore']['code'],
                 'call self%stateStore_(stateFile,gslStateFile,'
                 'stateOperationID)',
                 "stateStore wrapper delegates to stateStore_")
    ss_code = methods['stateStore_']['code']
    assert_equal('position=FTell(stateFile)' in ss_code, True,
                 "stateStore_ tracks stream position via FTell")
    assert_equal("call displayIndent(var_str('storing state for " in ss_code,
                 True, "stateStore_ emits displayIndent log line")
    assert_equal('if (self%stateOperationID == stateOperationID)' in ss_code,
                 True, "stateStore_ guards on stateOperationID dedup")
    assert_equal('write (stateFile) self%x, &\n  & self%y' in ss_code, True,
                 "stateStore_ batches static variables in one write")
    assert_equal('!$GLC attributes unused :: gslStateFile' in ss_code, True,
                 "stateStore_ emits unused-attributes for gslStateFile")
    sr_code = methods['stateRestore_']['code']
    assert_equal('read (stateFile) self%x, &\n  & self%y' in sr_code, True,
                 "stateRestore_ batches reads symmetrically")
    assert_equal(
        "call displayMessage('restoring \"x\"',"
        "verbosity=verbosityLevelWorking)" in sr_code, True,
        "stateRestore_ emits a displayMessage per static variable")


def test_functionclass_build_allowed_parameters_method():
    print("\n=== Testing FunctionClass._build_allowed_parameters_method ===")
    directive = {'name': 'testFoo'}

    # Positive path: constructor declares `type(inputParameters) ::
    # parameters` and has two inputParameter directives.
    ip_a = {'type': 'inputParameter',
            'directive': {'source': 'parameters', 'name': 'paramA'},
            'firstChild': None, 'sibling': None}
    ip_b = {'type': 'inputParameter',
            'directive': {'source': 'parameters', 'name': 'paramB'},
            'firstChild': None, 'sibling': None}
    cls = _make_simple_class(
        'testFooSimple', 'testFoo',
        constructor_children=[ip_a, ip_b])
    methods = {}
    _build_allowed_parameters_method(directive, [cls], methods)
    code = methods['allowedParameters']['code']
    assert_equal('logical                                            '
                 ':: isNew' in code, True,
                 "positive path prepends logical isNew decl")
    assert_equal("if (sourceName == 'parameters')" in code, True,
                 "source-guard conditional emitted")
    assert_equal("allowedParameters(1)='paramA'" in code, True,
                 "paramA listed in the else-branch allocate")
    assert_equal("allowedParameters(2)='paramB'" in code, True,
                 "paramB listed second in the else-branch allocate")
    assert_equal('call move_alloc(allowedParameters,allowedParametersTmp)'
                 in code, True,
                 "reallocation uses move_alloc (gfortran PR 37336 workaround)")
    assert_equal('recursive' in methods['allowedParameters'], True,
                 "method flagged recursive")

    # Unused-attributes fallback: no constructor found -> no params found.
    empty_cls = _make_simple_class('testFooEmpty', 'testFoo')
    methods2 = {}
    _build_allowed_parameters_method(directive, [empty_cls], methods2)
    assert_equal(methods2['allowedParameters']['code'].strip(),
                 '!$GLC attributes unused :: self, allowedParameters, '
                 'sourceName',
                 "empty path collapses to GLC unused-attributes stub")


def test_functionclass_build_descriptor_methods():
    print("\n=== Testing FunctionClass._build_descriptor_methods ===")
    directive = {'name': 'testFoo', 'data': []}

    # Single class with one scalar double-precision parameter reached via
    # the parameters constructor -> inputParameter directive.
    ip = {'type': 'inputParameter',
          'directive': {'source': 'parameters', 'name': 'paramA'},
          'firstChild': None, 'sibling': None}
    cls = _make_simple_class(
        'testFooSimple', 'testFoo',
        member_declarations=[{
            'intrinsic':     'double precision', 'type': None,
            'attributes':    [],
            'variables':     ['parama'],
            'variableNames': ['paramA'],
        }],
        constructor_children=[ip])
    methods = {}
    state_storables = {
        'functionClasses':        {'functionClass': [{'name': 'testFoo'}]},
        'functionClassInstances': [],
    }
    _build_descriptor_methods(
        directive, [cls], {'testFooSimple': cls}, methods,
        {'source': 'mod.F90'}, state_storables)
    d_code = methods['descriptor']['code']
    assert_equal('logical :: includeFileModificationTimes_' in d_code, True,
                 "descriptor emits includeFileModificationTimes_ decl")
    assert_equal("if (includeClass_) call descriptor%addParameter"
                 "('testFoo','simple')" in d_code, True,
                 "descriptor class header emitted under includeClass_")
    assert_equal("write (parameterLabel,'(e17.10)') self%paramA" in d_code,
                 True, "double-precision param written with e17.10")
    assert_equal("call parameters%addParameter('paramA',"
                 "trim(adjustl(parameterLabel)))" in d_code, True,
                 "addParameter call follows the write for scalar numeric")

    # Discovery-only probe: parameter classified into 'parameters' bucket.
    pn, dp, sp, dm, pcu, eo, fm, su = _descriptor_discover_class(
        cls, directive, {'testFooSimple': cls}, state_storables)
    assert_equal(dm, True, "discovery sets declaration_matches on inputParameters")
    assert_equal(eo, 'testFoo', "discovery extracts extension_of from type opener")
    assert_equal(su, 1, "discovery marks supported=1 on clean path")
    assert_equal(len(dp.get('parameters') or []), 1,
                 "one descriptor_parameters entry captured")

    # hashedDescriptor companion.
    h_code = methods['hashedDescriptor']['code']
    assert_equal('descriptor=inputParameters()' in h_code, True,
                 "hashedDescriptor builds a fresh inputParameters tree")
    assert_equal('call self%descriptor(descriptor,includeClass=.true.,'
                 'includeFileModificationTimes=includeFileModificationTimes)'
                 in h_code, True,
                 "hashedDescriptor delegates to descriptor method")
    assert_equal('type is (testFooSimple)' in h_code, True,
                 "hashedDescriptor select-type arm per class")
    assert_equal('String_C_To_Fortran(testFooSimple5)' in h_code, True,
                 "hashedDescriptor references per-class MD5 symbol <name>5")
    assert_equal(
        'testFooHashedDescriptor=Hash_MD5(descriptorString)' in h_code, True,
        "hashedDescriptor finalises via Hash_MD5 assignment")
    assert_equal(methods['hashedDescriptor']['type'], 'type(varying_string)',
                 "hashedDescriptor returns a varying_string")


def main():
    print("=" * 70)
    print("Python Utils Unit Tests")
    print("=" * 70)

    test_topo_sort()
    test_list_extra_utils()
    test_extract_variables()
    test_parse_declaration()
    test_unit_openers()

    # Parser module tests
    test_parse_module_uses()
    test_parse_module_uses_preprocessor()
    test_update_uses()
    test_add_uses()
    test_parse_declarations_pass()
    test_build_declarations()
    test_add_declarations()
    test_add_attributes()
    test_get_declaration_and_declaration_exists()
    test_parse_visibilities()
    test_update_visibilities()
    test_parse_openmp()
    test_openmp_update()
    test_openmp_copyin()
    test_parse_directives()
    test_post_process_directives()

    # Process module tests
    test_process_tree_orchestrator()
    test_process_deep_copy_reset()
    test_process_deep_copy_finalize()
    test_process_optional_argument()
    test_process_hdf5_fc_interop()
    test_process_dependencies()
    test_process_generics_subtree_expansion()
    test_process_generics_regex_form()
    test_process_generics_conditional()
    test_process_non_processed()
    test_process_allocate()
    test_process_for_each_rank1()
    test_process_for_each_rank0()
    test_source_introspection_instrument()
    test_source_introspection_location()
    test_source_introspection_process()
    test_process_parameter_migration()
    test_process_debug_mpi()
    test_process_debug_mpi_disabled()
    test_process_debug_hdf5()
    test_process_profile_openmp()
    test_process_constant_inline()
    test_process_conditional_call()
    test_process_conditional_call_parameter_present()
    test_process_input_parameters_validate()
    test_extract_directives_basic()
    test_extract_directives_multiline_and_wildcard()
    test_extract_directives_conditions()
    test_extract_directives_missing_file()
    test_process_state_store()
    test_process_state_restore()
    test_process_add_meta_property()
    test_process_meta_property_database()
    test_process_input_parameter_simple()
    test_process_input_parameter_with_default_and_no_output()
    test_insert_pre_post_contains()
    test_process_constructors_basic_assignment()
    test_process_constructors_pointer_refcount()
    test_process_constructors_wildcard_type()
    test_process_constructors_slash_suppresses_count()
    test_process_functions_global_pointers()
    test_process_functions_global_establish()
    test_set_visibility()
    test_process_enumeration_basic()
    test_process_enumeration_validator()
    test_process_enumeration_encode_decode()
    test_process_deep_copy_actions_setto()
    test_process_deep_copy_actions_method_call()
    test_process_state_storable_basic()
    test_process_state_storable_excludes_and_restore_to()
    test_process_state_storable_class_restore()
    test_dependency_sort_after_and_before()
    test_process_event_hooks_static_single_hook()
    test_process_event_hooks_static_ordering()
    test_event_hooks_interface_type_get()
    test_process_event_hook_call_site()
    test_process_event_hook_openmp_wrappers()
    test_process_event_hook_manager()
    test_process_thread_safe_io_enabled()
    test_process_thread_safe_io_disabled()
    test_process_thread_safe_io_skips_error_module()
    test_process_source_digest()
    test_process_object_destructor()
    test_process_reference_increment_acquire_construct()
    test_process_class_documentation_disabled()
    test_process_class_documentation_enabled()
    test_functionclass_utils_small_helpers()
    test_functionclass_utils_class_dependencies()
    test_functionclass_descriptor_classification()
    test_functionclass_linkedlist_register_and_module()
    test_functionclass_linkedlist_deep_copy()
    test_functionclass_linkedlist_state_store()
    test_functionclass_linkedlist_allowed_and_descriptor_and_assigner()
    test_functionclass_deepcopy_copied_self_block()
    test_functionclass_deepcopy_generate_assignment_allocatable()
    test_functionclass_deepcopy_declarations_functionclass_pointer()
    test_functionclass_deepcopy_declarations_explicit_actions()
    test_functionclass_deepcopy_declarations_hdf5_and_setto()
    test_functionclass_statestore_explicit_function()
    test_functionclass_statestore_generate_allocatable()
    test_functionclass_statestore_variables_static_and_restore()
    test_functionclass_statestore_custom_hooks_flag()
    test_functionclass_statestore_variables_allocatable_intrinsic()
    test_functionclass_is_function_class_pointer()
    test_functionclass_init_code_content()
    test_functionclass_build_base_method_stubs()
    test_functionclass_build_object_type_method()
    test_functionclass_load_and_sort_classes()
    test_functionclass_load_and_sort_classes_default_validation()
    test_functionclass_process_hook_is_no_op_when_no_directive()
    test_functionclass_process_hook_marks_class_nodes()
    test_functionclass_format_variable_definitions()
    test_functionclass_source_digest_binding()
    test_functionclass_generate_type_definition()
    test_functionclass_generate_constructor_with_default_and_recursion()
    test_functionclass_generate_method_functions()
    test_functionclass_generate_documentation_basic()
    test_functionclass_build_assignment_method()
    test_functionclass_build_deep_copy_methods()
    test_functionclass_build_state_store_methods()
    test_functionclass_build_allowed_parameters_method()
    test_functionclass_build_descriptor_methods()

    print("\n" + "=" * 70)
    print(f"Results: {PASS_COUNT} passed, {FAIL_COUNT} failed")
    print("=" * 70)

    if FAIL_COUNT > 0:
        print("FAIL: Some tests failed")
        sys.exit(1)
    else:
        print("PASS: All tests passed")
        sys.exit(0)


if __name__ == '__main__':
    main()
