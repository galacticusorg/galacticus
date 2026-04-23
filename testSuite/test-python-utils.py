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
from Galacticus.Build.SourceTree.Process import (
    PROCESS_HOOKS, PROCESS_DEPENDENCIES, POSTPROCESS_HOOKS,
    register_process, process_tree,
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
