#!/usr/bin/env python3
# Unit tests for Python utility modules ported from Perl.
# Tests the following modules:
#   - Sort.Topo
#   - List.ExtraUtils
#   - Fortran.Utils.extract_variables
#   - Galacticus.Build.SourceTree.Parse.Declarations.parse_declaration

import sys
import os

sys.path.insert(0, os.path.join(os.environ.get('GALACTICUS_EXEC_PATH', '.'), 'python'))

from Sort.Topo import sort as topo_sort
from List.ExtraUtils import as_array, smart_push, hash_list, sorted_keys
from Fortran.Utils import extract_variables
from Galacticus.Build.SourceTree.Parse.Declarations import parse_declaration

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
