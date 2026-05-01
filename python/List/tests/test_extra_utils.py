# Tests for `List.ExtraUtils` (port of Perl List::ExtraUtils).
#
# These functions are small, pure, and used widely across the codebase but
# previously had zero coverage.  The intent of these tests is to lock down
# the documented behaviour against the Perl reference (None handling,
# scalar-vs-list dispatch, dict-keyed-vs-list-of-dicts coercion) so any
# accidental change is caught.

from List.ExtraUtils import smart_push, as_array, hash_list, sorted_keys


# ---------------------------------------------------------------------------
# smart_push
# ---------------------------------------------------------------------------

def test_smart_push_none_is_noop():
    """A None item must not modify the array (mirrors Perl `defined($item)`)."""
    arr = [1, 2]
    smart_push(arr, None)
    assert arr == [1, 2]


def test_smart_push_scalar_appends_one_element():
    arr = []
    smart_push(arr, 'x')
    smart_push(arr, 7)
    assert arr == ['x', 7]


def test_smart_push_list_extends_in_place():
    """A list item is unpacked: each element is appended individually."""
    arr = ['a']
    smart_push(arr, ['b', 'c', 'd'])
    assert arr == ['a', 'b', 'c', 'd']


def test_smart_push_empty_list_is_noop():
    arr = [1]
    smart_push(arr, [])
    assert arr == [1]


# ---------------------------------------------------------------------------
# as_array
# ---------------------------------------------------------------------------

def test_as_array_none_returns_empty_list():
    assert as_array(None) == []


def test_as_array_scalar_returns_singleton():
    assert as_array('x') == ['x']
    assert as_array(0)   == [0]


def test_as_array_list_returns_independent_copy():
    """Mutating the returned list must not affect the input — important
    because callers commonly accumulate further into the result."""
    src = [1, 2, 3]
    out = as_array(src)
    out.append(4)
    assert src == [1, 2, 3]
    assert out == [1, 2, 3, 4]


# ---------------------------------------------------------------------------
# hash_list
# ---------------------------------------------------------------------------

def test_hash_list_none_returns_empty():
    assert hash_list(None) == []


def test_hash_list_empty_dict_returns_empty():
    assert hash_list({}) == []


def test_hash_list_returns_values_sorted_by_key():
    """Returned values are emitted in `sorted(keys())` order — the Perl
    reference always sorted lexicographically, and downstream code relies on
    deterministic ordering for code generation."""
    d = {'banana': 2, 'apple': 1, 'cherry': 3}
    assert hash_list(d) == [1, 2, 3]


def test_hash_list_with_key_as_stamps_key_into_each_value_dict():
    """When `key_as` is given, the dict key is stamped into each value
    dict under that field name (mutating in place, like Perl)."""
    d = {'A': {'value': 1}, 'B': {'value': 2}}
    out = hash_list(d, key_as='name')
    assert out == [{'value': 1, 'name': 'A'}, {'value': 2, 'name': 'B'}]
    # In-place mutation is part of the contract:
    assert d['A'] == {'value': 1, 'name': 'A'}


def test_hash_list_with_key_as_promotes_non_dict_value():
    """Non-dict values get replaced by a fresh `{key_as: key}` dict so
    callers can rely on dict-shaped results regardless of input."""
    d = {'A': 'scalar', 'B': 42}
    out = hash_list(d, key_as='name')
    assert out == [{'name': 'A'}, {'name': 'B'}]


# ---------------------------------------------------------------------------
# sorted_keys
# ---------------------------------------------------------------------------

def test_sorted_keys_none_returns_empty():
    assert sorted_keys(None) == []


def test_sorted_keys_empty_returns_empty():
    assert sorted_keys({}) == []


def test_sorted_keys_returns_lexicographic_order():
    assert sorted_keys({'b': 1, 'a': 1, 'c': 1}) == ['a', 'b', 'c']
