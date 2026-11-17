"""Regression test for enumeration members in `stateStorable` types.

Bug (#1786): `StateStorable` passed every `type(...)` member to
`_process_derived_declaration`, which emits nothing for a type that is not
itself registered as state-storable. Enumeration members
(`type(enumeration...Type)`) were therefore silently left at their default
values on restore. For `rangeLattice` this left `scheme` undefined, so a
restored lattice reported itself as undefined. Enumeration members are now
stored by their integer `ID`, as the functionClass state store already does.
"""

import pytest

from Galacticus.Build.SourceTree import parse_code, walk_tree
from Galacticus.Build.SourceTree.Process.StateStorable import (
    _emit_store_restore, _parse_type_opener,
)


def _emit(source, target):
    tree    = parse_code(source, name='<test>', instrument=False)
    classes = {}
    for node in walk_tree(tree):
        if node.get('type') == 'type':
            name, extends, abstract = _parse_type_opener(node['opener'])
            classes[name] = {'node': node, 'extends': extends, 'abstract': abstract}
    code, _, _ = _emit_store_restore({'class': target}, classes, set())
    store, restore = code.split('subroutine ' + target + 'StateRestore', 1)
    return store, restore


_SOURCE = """module m
  type, public :: lattice
     type   (enumerationSchemeType) :: scheme=schemeUndefined
     integer                        :: count =0
  end type lattice
end module m
"""


def test_enumeration_member_is_stored_and_restored_by_id():
    store, restore = _emit(_SOURCE, 'lattice')
    assert 'write (stateFile) self%scheme%ID' in store
    assert 'read  (stateFile) self%scheme%ID' in restore
    # The intrinsic member is still stored alongside it.
    assert 'self%count' in store
    assert 'self%count' in restore


def test_allocatable_enumeration_member_is_rejected():
    source = _SOURCE.replace(
        'type   (enumerationSchemeType) :: scheme=schemeUndefined',
        'type   (enumerationSchemeType), allocatable, dimension(:) :: scheme')
    with pytest.raises(RuntimeError, match='enumeration'):
        _emit(source, 'lattice')


def test_inline_array_members_are_stored_whole():
    # A member declared with an inline array specification must be stored as
    # a whole array: keeping the specification in the name (as `self%s(2)`)
    # would store only its last element.
    source = """module m
  type, public :: holder
     type   (enumerationSchemeType) :: schemes(2)
     double precision               :: values (3)=0.0d0
  end type holder
end module m
"""
    store, restore = _emit(source, 'holder')
    for code in (store, restore):
        assert 'self%schemes%ID' in code
        assert 'self%values' in code
        assert 'self%schemes(' not in code
        assert 'self%values(' not in code
