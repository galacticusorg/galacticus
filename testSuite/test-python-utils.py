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
import Galacticus.Build.SourceTree.Process.InputParameter            # noqa: F401
from Galacticus.Build.SourceTree.Process import (
    PROCESS_HOOKS, PROCESS_DEPENDENCIES, POSTPROCESS_HOOKS,
    register_process, process_tree,
)
from Galacticus.Build.Directives import extract_directives, extract_directive

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
# Main
# ============================================================================

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
    test_process_meta_property_database()
    test_process_input_parameter_simple()
    test_process_input_parameter_with_default_and_no_output()

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
