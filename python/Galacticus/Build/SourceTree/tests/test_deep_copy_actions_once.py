"""Regression test for repeated actions in generated `deepCopyActions` procedures.

Bug: `DeepCopyActions` emitted a level's `setTo` and `methodCall` actions once for every declaration node in the type body (each
run of declarations broken by a comment or directive), rather than once for the level. So `interpolator`, whose body has three
such nodes, called `GSLReallocate` three times on every deep copy (and `interpolator2D` twice). The calls are idempotent, but each
does a redundant allocation and initialization of the GSL interpolator. Each action must be emitted exactly once for each level
of the type hierarchy - including a level whose type body has no declarations.

Andrew Benson, Claude (19-November-2026).
"""

from Galacticus.Build.SourceTree import parse_code, walk_tree
from Galacticus.Build.SourceTree.Process.DeepCopyActions import (
    _emit_deep_copy_action, _parse_type_opener,
)


def _emit(source, target, directive):
    tree    = parse_code(source, name='<test>', instrument=False)
    classes = {}
    for node in walk_tree(tree):
        if node.get('type') == 'type':
            name, extends, abstract = _parse_type_opener(node['opener'])
            classes[name] = {'node': node, 'extends': extends, 'abstract': abstract}
    return _emit_deep_copy_action(target, classes, directive)


# As in `interpolator`, a comment separates the declarations, so the body holds more than one declaration node.
_SOURCE = """module m
  type, public :: holder
     !!{RST
     A type holding some values.
     !!}
     private
     double precision                            :: x
     integer                                     :: n
     ! A comment between declarations.
     logical                                     :: initialized
     double precision, allocatable, dimension(:) :: values
  end type holder
  type, public, extends(holder) :: holderChild
     integer :: m
  end type holderChild
  type, public, extends(holder) :: holderEmpty
  end type holderEmpty
end module m
"""


def _branches(code):
    # Split the generated code into its `type is` branches, keyed by type name.
    branches = {}
    for branch in code.split('type is (')[1:]:
        branches[branch.split(')')[0]] = branch
    # One branch for each concrete type: holder, holderChild, holderEmpty.
    assert sorted(branches) == ['holder', 'holderChild', 'holderEmpty'], code
    return branches


def test_method_call_is_emitted_once_per_level():
    directive = {'class': 'holder', 'holder': {'methodCall': {'method': 'reallocate'}}}
    code = _emit(_SOURCE, 'holder', directive)
    for branch in _branches(code).values():
        # Each concrete type inherits the one action declared for the `holder` level, once.
        assert branch.count('call self%reallocate()') == 1, code


def test_set_to_is_emitted_once_per_level():
    directive = {'class': 'holder', 'holder': {'setTo': {'variables': 'initialized', 'state': '.false.'}}}
    code = _emit(_SOURCE, 'holder', directive)
    for branch in _branches(code).values():
        assert branch.count('self%initialized=.false.') == 1, code


def test_child_level_actions_are_emitted_once_and_only_for_the_child():
    directive = {
        'class'      : 'holder',
        'holder'     : {'methodCall': {'method': 'reallocate'}},
        'holderChild': {'methodCall': {'method': 'resetChild'}, 'setTo': {'variables': 'm', 'state': '0'}},
        'holderEmpty': {'methodCall': {'method': 'resetEmpty'}},
    }
    branches = _emit_branches(directive)
    assert branches['holderChild'].count('call self%resetChild()') == 1
    assert branches['holderChild'].count('self%m=0'              ) == 1
    assert branches['holderChild'].count('call self%reallocate()') == 1
    assert 'resetChild' not in branches['holder'     ]
    assert 'resetChild' not in branches['holderEmpty']
    # A level whose type body has no declarations still has its actions emitted.
    assert branches['holderEmpty'].count('call self%resetEmpty()') == 1


def _emit_branches(directive):
    return _branches(_emit(_SOURCE, 'holder', directive))
