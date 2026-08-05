"""Tests for `Galacticus.Build.SourceTree.Parse.Signatures`.

The module renders the notation in which method arguments and return types are
documented, and recovers that notation from a parsed procedure.  Both the code
generators (which embed the notation in the `<methods>` metadata they emit) and
the documentation extractor (which renders it onto the class pages) depend on
it, so the notation is pinned here.
"""

from Galacticus.Build.SourceTree import parse_code, walk_tree
from Galacticus.Build.SourceTree.Parse.Signatures import (
    arguments_to_rst, is_declaration, normalize_argument, procedure_signature,
    return_type_to_rst, split_arguments,
)


def test_derived_type_argument_keeps_its_parentheses():
    """A `type(…)` spec must not collapse to the meaningless `typetreeNode`."""
    assert arguments_to_rst(
        ['type(treeNode), intent(inout), target :: node']
    ) == ['``type(treeNode) node`` [inout]']


def test_argument_name_keeps_its_original_case():
    assert arguments_to_rst(
        ['double precision, intent(in) :: massSplitting']
    ) == ['``double precision massSplitting`` [in]']


def test_every_variable_of_a_declaration_is_rendered():
    assert arguments_to_rst(
        ['double precision, intent(in), optional :: massSplitting, velocityKick']
    ) == ['``double precision massSplitting`` (optional) [in]',
          '``double precision velocityKick`` (optional) [in]']


def test_optional_is_spelled_out_not_bracketed():
    """`[name]` would read as a third direction marker alongside
    `[in]`/`[out]`/`[inout]`."""
    rendered, = arguments_to_rst(
        ['integer, intent(in), optional :: weightIndex'])
    assert rendered == '``integer weightIndex`` (optional) [in]'
    assert '[weightIndex]' not in rendered


def test_rank_belongs_to_the_type_not_the_name():
    """Declared either way — on the variable or as an attribute."""
    assert arguments_to_rst(
        ['double precision, intent(in) :: position(3)']
    ) == ['``double precision(3) position`` [in]']
    assert arguments_to_rst(
        ['double precision, intent(out), allocatable, dimension(:) :: array']
    ) == ['``double precision(:) array`` [out]']


def test_trailing_underscore_names_need_no_escaping():
    """A trailing underscore — which Galacticus uses for object variables —
    would read as a reference in text, but the name is inside the literal."""
    assert arguments_to_rst(
        ['class(cosmologyFunctionsClass), intent(in) :: cosmologyFunctions_']
    ) == ['``class(cosmologyFunctionsClass) cosmologyFunctions_`` [in]']


def test_unparseable_declaration_is_skipped():
    assert arguments_to_rst(['this is not a declaration']) == []


def test_void_has_no_return_type():
    """`void` is how a directive spells "subroutine": there is nothing to
    render."""
    assert return_type_to_rst('void')       is None
    assert return_type_to_rst('``void``')   is None
    assert return_type_to_rst('')           is None
    assert return_type_to_rst(None)         is None


def test_return_type_accepts_fortran_or_notation():
    assert return_type_to_rst('double precision')     == '``double precision``'
    assert return_type_to_rst('``double precision``') == '``double precision``'


def test_return_type_rank_and_pointer_markers():
    assert return_type_to_rst('double precision, dimension(3)') \
        == '``double precision(3)``'
    # A leading `*` is how the component metadata marks a pointer result.
    assert return_type_to_rst('``*type(treeNode)``') \
        == '``type(treeNode)`` (pointer)'


def test_arguments_split_at_argument_boundaries_not_every_comma():
    assert split_arguments(
        '``double precision(:,:)`` matrix [out], ``integer`` count [in]'
    ) == ['``double precision(:,:) matrix`` [out]', '``integer count`` [in]']
    assert split_arguments('') == []
    assert split_arguments(None) == []


def test_declarations_are_told_apart_from_the_notation():
    """Both reach the extractor: hand-written `<argument>` elements carry
    declarations, the generators emit the notation."""
    assert is_declaration('type(treeNode), intent(in) :: node')
    assert not is_declaration('``type(treeNode) node`` [in]')


def test_a_hand_written_argument_gets_its_name_moved_into_the_literal():
    """Metadata is authored with the two apart; rendering wants them together,
    and an argument already in that form is left alone."""
    assert normalize_argument('``double precision`` timeStep [in]') \
        == '``double precision timeStep`` [in]'
    assert normalize_argument('``double precision timeStep`` [in]') \
        == '``double precision timeStep`` [in]'
    assert normalize_argument('``void``') == '``void``'


_FIXTURE = """\
module Test_Fixture
  type :: testType
   contains
     procedure :: radiusStripped => testRadiusStripped
     procedure :: reset          => testReset
  end type testType
contains

  double precision function testRadiusStripped(self,node,factor)
    implicit none
    class           (testType), intent(inout)           :: self
    type            (treeNode), intent(inout), target   :: node
    double precision          , intent(in   ), optional :: factor

    testRadiusStripped=0.0d0
    return
  end function testRadiusStripped

  subroutine testReset(self)
    implicit none
    class(testType), intent(inout) :: self

    return
  end subroutine testReset

  function testCount(self) result(countTotal)
    implicit none
    integer        , intent(  out) :: countTotal
    class(testType), intent(inout) :: self

    countTotal=0
    return
  end function testCount
end module Test_Fixture
"""


def _procedure(name):
    tree = parse_code(_FIXTURE, name='testFixture.F90', instrument=False)
    for node in walk_tree(tree):
        if node['type'] in ('function', 'subroutine') and node.get('name') == name:
            return node
    raise AssertionError(f'{name} not found in the fixture')


def test_procedure_signature_reads_type_and_arguments():
    return_type, arguments = procedure_signature(_procedure('testRadiusStripped'))
    assert return_type == '``double precision``'
    assert arguments   == ['``type(treeNode) node`` [inout]',
                           '``double precision factor`` (optional) [in]']


def test_procedure_signature_drops_the_passed_object():
    """A method's own object is not one of its arguments — unless the binding
    is `nopass`, in which case the first dummy is a real argument."""
    _, arguments = procedure_signature(_procedure('testRadiusStripped'))
    assert not any('self' in argument for argument in arguments)
    _, arguments = procedure_signature(_procedure('testRadiusStripped'),
                                       skip_passed_object=False)
    assert arguments[0] == '``class(testType) self`` [inout]'


def test_procedure_signature_of_a_subroutine_has_no_return_type():
    assert procedure_signature(_procedure('testReset')) == (None, [])


def test_procedure_signature_reads_a_result_variable_declaration():
    """A function with no type prefix declares its result in the body."""
    return_type, _ = procedure_signature(_procedure('testCount'))
    assert return_type == '``integer``'
