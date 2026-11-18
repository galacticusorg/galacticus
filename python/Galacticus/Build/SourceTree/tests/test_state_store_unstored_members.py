"""Tests for how the state store generators treat members that are not of a
registered state-storable type (issue #1786).

Such members used to be skipped silently, so after a state restore they kept
whatever value the freshly constructed object had - which broke
reproducibility whenever a member that *was* restored (an "initialized" flag,
or a memoization key) referred to one that was not. Now:

* `tabulationInverse` members are stored through that type's own
  `stateStore(stateFile)`/`stateRestore(stateFile)` procedures;
* type names are matched case-insensitively, as Fortran does (a member
  declared `type(table1dGeneric)` was skipped because its type is registered
  as `table1DGeneric`);
* any other owned member that is skipped is reported by a build-time warning,
  unless it is explicitly excluded, or its type is listed in `UNSTORED_TYPES`
  (which includes `varying_string`: strings are configuration, some of which,
  such as output file names, may legitimately differ in the restoring run).
"""

import contextlib
import io

from Galacticus.Build.StateStorables import warn_unstored_member
from Galacticus.Build.SourceTree import parse_code, walk_tree
from Galacticus.Build.SourceTree.Parse.Declarations import parse_declaration
from Galacticus.Build.SourceTree.Process.FunctionClass.StateStore import (
    state_store_variables,
)
from Galacticus.Build.SourceTree.Process.FunctionClass.Utils import (
    variable_rank,
)
from Galacticus.Build.SourceTree.Process.StateStorable import (
    _emit_store_restore, _parse_type_opener,
)


# ---------------------------------------------------------------------------
# functionClass state store
# ---------------------------------------------------------------------------

def _function_class_walk(lines, excludes=(), state_storables=None):
    state_store  = {'excludes': list(excludes), 'owner': 'widgetFancy',
                    'ownerSource': 'widget.F90'}
    state_stores = {}
    declarations = [parse_declaration(line) for line in lines]
    state_store_variables(state_stores, state_store, None, declarations, [],
                          state_storables or {})
    return state_store, state_stores


def test_tabulation_inverse_member_is_stored_by_its_own_procedures():
    state_store, state_stores = _function_class_walk(
        ['type(tabulationInverse) :: massScaleFree_'])
    assert ' call self%massscalefree_%stateStore  (stateFile)' in state_store['outputCode']
    assert ' call self%massscalefree_%stateRestore(stateFile)' in state_store['inputCode']
    assert not state_store.get('unstored')
    # It needs no GSL state file.
    assert not state_stores.get('gslStateFileUsed')


def test_self_storing_array_with_inline_dimension_is_stored_element_wise():
    state_store, _ = _function_class_walk(
        ['type(tabulationInverse) :: tables(3)'])
    assert 'do i1=lbound(self%tables,dim=1),ubound(self%tables,dim=1)' in state_store['outputCode']
    assert 'call self%tables(i1)%stateStore  (stateFile)' in state_store['outputCode']
    assert 'call self%tables(i1)%stateRestore(stateFile)' in state_store['inputCode']


def test_varying_string_member_is_neither_stored_nor_reported():
    # Output file names, for example, may differ in the restoring run, so must
    # not be overwritten by the values from the storing run.
    state_store, _ = _function_class_walk(['type(varying_string) :: fileName'])
    assert 'filename' not in state_store['outputCode']
    assert 'filename' not in state_store['inputCode']
    assert not state_store.get('unstored')


def test_variable_rank_honors_inline_dimensions():
    declaration = parse_declaration('type(tabulationInverse) :: a(3), b(0:2,4), c')
    variables   = declaration['variables']
    assert [variable_rank(declaration, v) for v in variables] == [1, 2, 0]


def test_unregistered_member_is_reported_with_its_owner_and_original_case():
    state_store, _ = _function_class_walk(['type(fooType) :: widgetState'])
    assert state_store['unstored'] == [
        ('widgetFancy', 'widgetState', 'fooType', 'widget.F90')]
    assert 'widgetstate' not in state_store['outputCode']


def test_pointer_excluded_and_unstored_type_members_are_not_reported():
    state_store, _ = _function_class_walk(
        ['type(fooType), pointer :: borrowed',
         'type(fooType) :: excludedMember',
         'type(ompLock) :: lock',
         'type(rootFinder) :: finder'],
        excludes=['excludedMember'])
    assert not state_store.get('unstored')


def test_registered_type_is_matched_case_insensitively():
    state_storables = {'stateStorables': [
        {'type': 'table1DGeneric', 'module': 'Tables'}]}
    state_store, _ = _function_class_walk(
        ['type(table1dGeneric) :: perturbations'],
        state_storables=state_storables)
    assert 'call self%perturbations%stateStore' in state_store['outputCode']
    assert not state_store.get('unstored')


def test_warning_names_member_type_and_remedy():
    stderr = io.StringIO()
    with contextlib.redirect_stderr(stderr):
        warn_unstored_member('widgetFancy', 'widgetState', 'fooType', 'widget.F90')
    message = stderr.getvalue()
    # Starts with "Warning:", which the CI build checks detect.
    assert message.startswith("Warning: state store for 'widgetFancy' (widget.F90) skips member 'widgetState' of type 'fooType'")
    assert '<exclude variables="widgetState"/>' in message


# ---------------------------------------------------------------------------
# stateStorable types
# ---------------------------------------------------------------------------

def _state_storable_emit(source, target, storable_types=()):
    tree    = parse_code(source, name='<test>', instrument=False)
    classes = {}
    for node in walk_tree(tree):
        if node.get('type') == 'type':
            name, extends, abstract = _parse_type_opener(node['opener'])
            classes[name] = {'node': node, 'extends': extends, 'abstract': abstract}
    stderr = io.StringIO()
    with contextlib.redirect_stderr(stderr):
        code, _, _ = _emit_store_restore(
            {'class': target}, classes, {t.lower() for t in storable_types})
    store, restore = code.split('subroutine ' + target + 'StateRestore', 1)
    return store, restore, stderr.getvalue()


_SOURCE = """module m
  type, public :: holder
     type   (tabulationInverse)                         :: table
     type   (fooType       )                            :: gadget
     type   (fooType       ), pointer                   :: borrowed => null()
     type   (barType       ), allocatable, dimension(:) :: items
  end type holder
end module m
"""


def test_state_storable_tabulation_inverse_member_is_stored():
    store, restore, _ = _state_storable_emit(_SOURCE, 'holder', ['barType'])
    assert 'call self%table%stateStore  (stateFile)' in store
    assert 'call self%table%stateRestore(stateFile)' in restore


def test_state_storable_unregistered_member_is_reported():
    _, _, warnings = _state_storable_emit(_SOURCE, 'holder', ['barType'])
    assert "skips member 'gadget' of type 'fooType'" in warnings
    # Pointers are not reported, and neither are registered types.
    assert 'borrowed' not in warnings
    assert 'items' not in warnings


def test_state_storable_array_elements_are_looped_over_their_bounds():
    # Arrays may have non-unit lower bounds, which are now restored, so the
    # element loops must run over the actual bounds, not 1:size().
    store, restore, _ = _state_storable_emit(_SOURCE, 'holder', ['barType'])
    for code in (store, restore):
        assert 'do i1=lbound(self%items,dim=1),ubound(self%items,dim=1)' in code
        assert 'size(self%items' not in code
