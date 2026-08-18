"""Tests for the variable-declaration formatter.

The emitter in `Parse.Declarations` reproduces the column layout described in
the developer guide: the intrinsic type, its parenthetical, the attributes, the
`::` and the variables each occupy a column whose width is shared across the
whole alignment group.

Covered here:

  1. The columns themselves, including the parenthetical field, the positional
     attribute columns, and `=` / `=>` alignment down a variable column.
  2. Canonical attribute ordering and the padded `intent` field, which is what
     makes "the same attribute in the same column" achievable at all.
  3. Alignment groups that span comments and `#ifdef` — `_pass_declarations`
     ends a declaration node at any line it does not recognize, so without
     rejoining, alignment would silently reset at every comment.
  4. Statements the formatter must *not* touch: derived-type definitions that
     slipped past the unit parser, and declarations whose hand-wrapped
     initializer no layout can fit.
  5. Idempotence, which two separate bugs broke — the OpenMP `!$` lines were
     skipped when measuring the block indent (so blocks marched three columns
     right per pass), and `final` bindings were dropped entirely.
"""

import os
import re

import pytest

from Galacticus.Build.SourceTree import (
    parse_code_for_format,
    serialize_for_format,
    walk_tree,
)
from Galacticus.Build.SourceTree.Parse.Declarations import (
    ATTRIBUTE_ORDER,
    LINE_LENGTH_MAX,
    format_declarations_tree,
)


def _format(source):
    tree = parse_code_for_format(source, name='<test>')
    format_declarations_tree(tree)
    return serialize_for_format(tree['firstChild'])


def _column_of(text, needle):
    """Column at which `needle` appears on each line that contains it."""
    return {line.index(needle) for line in text.splitlines() if needle in line}


# ---------------------------------------------------------------------------
# 1. Columns
# ---------------------------------------------------------------------------

def test_double_colon_aligns_across_group():
    source = (
        "  subroutine s\n"
        "    implicit none\n"
        "    integer :: i\n"
        "    double precision :: x\n"
        "    type(varying_string) :: name\n"
        "  end subroutine s\n"
    )
    out = _format(source)
    assert len(_column_of(out, '::')) == 1, out


def test_parenthetical_field_pads_declarations_without_one():
    """`double precision` has no type-spec, but must still leave room for the
    parenthetical so that the `::` of the group lines up."""
    source = (
        "  subroutine s\n"
        "    implicit none\n"
        "    type(someVeryLongTypeName) :: a\n"
        "    double precision :: b\n"
        "  end subroutine s\n"
    )
    out = _format(source)
    assert len(_column_of(out, '::')) == 1, out
    assert '(someVeryLongTypeName)' in out, out


def test_intrinsic_is_padded_before_the_parenthetical():
    """The house style pads the intrinsic and then opens the bracket —
    `class(x)` and `type (x)` — rather than padding inside the brackets."""
    source = (
        "  subroutine s\n"
        "    implicit none\n"
        "    type(a) :: p\n"
        "    class(b) :: q\n"
        "  end subroutine s\n"
    )
    out = _format(source)
    assert len(_column_of(out, '(')) == 1, out


def test_initializers_align_down_a_column():
    source = (
        "  subroutine s\n"
        "    implicit none\n"
        "    logical :: shortName=.false., otherFlag=.true.\n"
        "    logical :: aMuchLongerName=.false., b=.true.\n"
        "  end subroutine s\n"
    )
    out = _format(source)
    equals = [line.index('=') for line in out.splitlines() if '=' in line]
    assert len(set(equals)) == 1, out


def test_pointer_arrow_keeps_its_spaces_on_the_longest_name():
    """` => ` is spaced on both sides.  Because the padding of the name column
    is what separates shorter names from the arrow, the *widest* name in a
    column would otherwise come out as `widestName=> target`."""
    source = (
        "  type t\n"
        "     class(fooClass), pointer :: shortOne => null()\n"
        "     class(fooClass), pointer :: aDecidedlyLongerName => null()\n"
        "  end type t\n"
    )
    out = _format(source)
    assert 'aDecidedlyLongerName => null()' in out, out
    assert len(_column_of(out, '=>')) == 1, out


# ---------------------------------------------------------------------------
# 2. Attributes
# ---------------------------------------------------------------------------

def test_attributes_are_canonically_ordered():
    source = (
        "  subroutine s(x)\n"
        "    implicit none\n"
        "    double precision, dimension(:), allocatable, intent(inout) :: x\n"
        "  end subroutine s\n"
    )
    out = _format(source)
    line = next(l for l in out.splitlines() if '::' in l and 'x' in l)
    positions = [line.index(name) for name in
                 ('intent', 'allocatable', 'dimension') if name in line]
    assert positions == sorted(positions), line
    # And that is the order the canonical list prescribes.
    assert (ATTRIBUTE_ORDER.index('intent')
            < ATTRIBUTE_ORDER.index('allocatable')
            < ATTRIBUTE_ORDER.index('dimension'))


def test_intent_field_is_padded_to_a_single_column():
    """`in`, `out` and `inout` share one column: `in` left-aligned within it,
    `out` right-aligned, so all three `)` line up."""
    source = (
        "  subroutine s(a, b, c)\n"
        "    implicit none\n"
        "    integer, intent(in) :: a\n"
        "    integer, intent(out) :: b\n"
        "    integer, intent(inout) :: c\n"
        "  end subroutine s\n"
    )
    out = _format(source)
    assert 'intent(in   )' in out, out
    assert 'intent(  out)' in out, out
    assert 'intent(inout)' in out, out


def test_differing_attributes_share_a_positional_column():
    """A column holds whichever attribute comes first on each row — `pointer`
    on one, `intent(inout)` on another — sized to the wider of the two."""
    source = (
        "  subroutine s(p, q)\n"
        "    implicit none\n"
        "    class(fooClass), pointer :: p\n"
        "    type(barType), intent(inout) :: q\n"
        "  end subroutine s\n"
    )
    out = _format(source)
    assert len(_column_of(out, '::')) == 1, out


def test_attribute_case_is_preserved():
    """`dimension(nBins)` must not come back as `dimension(nbins)`."""
    source = (
        "  subroutine s(x)\n"
        "    implicit none\n"
        "    double precision, dimension(nBinsTotal) :: x\n"
        "  end subroutine s\n"
    )
    out = _format(source)
    assert 'dimension(nBinsTotal)' in out, out


def test_variable_case_is_preserved():
    source = (
        "  subroutine s\n"
        "    implicit none\n"
        "    double precision :: massStellarDisk\n"
        "  end subroutine s\n"
    )
    assert 'massStellarDisk' in _format(source)


# ---------------------------------------------------------------------------
# 3. Alignment groups
# ---------------------------------------------------------------------------

def test_alignment_spans_a_comment():
    source = (
        "  subroutine s\n"
        "    implicit none\n"
        "    integer :: i\n"
        "    ! A dividing comment.\n"
        "    double precision :: someLongerVariable\n"
        "  end subroutine s\n"
    )
    out = _format(source)
    assert '! A dividing comment.' in out, out
    assert len(_column_of(out, '::')) == 1, out


def test_alignment_spans_a_preprocessor_guard():
    source = (
        "  subroutine s\n"
        "    implicit none\n"
        "    integer :: i\n"
        "#ifdef USEMPI\n"
        "    type(mpiCommunicator) :: communicator\n"
        "#endif\n"
        "    double precision :: x\n"
        "  end subroutine s\n"
    )
    out = _format(source)
    assert out.count('#ifdef USEMPI') == 1, out
    assert out.count('#endif') == 1, out
    assert len(_column_of(out, '::')) == 1, out


def test_real_code_ends_an_alignment_group():
    """A statement between two runs of declarations is a genuine boundary, so
    the runs are sized independently."""
    source = (
        "  subroutine s\n"
        "    implicit none\n"
        "    integer :: i\n"
        "    call somethingEntirely()\n"
        "    double precision :: aVeryMuchLongerVariableName\n"
        "  end subroutine s\n"
    )
    out = _format(source)
    assert len(_column_of(out, '::')) == 2, out


# ---------------------------------------------------------------------------
# 4. What must not be touched
# ---------------------------------------------------------------------------

def test_derived_type_definition_is_left_alone():
    """`type, extends(...) :: name` opens a derived type.  A few in the tree
    carry enough padding to slip past the unit parser and reach the declaration
    parser instead; tidying that padding away would change the tree shape on
    the next parse, so such groups are skipped."""
    source = (
        "  type, extends(hdf5Group             ) :: hdf5File\n"
        "     integer :: identifier\n"
        "  end type hdf5File\n"
    )
    out = _format(source)
    assert 'extends(hdf5Group             )' in out, out


def test_unfittable_initializer_keeps_its_hand_wrapping():
    """A single variable with a large initializer cannot be broken by a
    formatter that only breaks between variables, and its continuation lines
    have already been joined away by the parser.  Re-emitting it from the
    parsed fields would collapse it onto one enormous line, so the original
    text is written back instead."""
    values = ', '.join(f'{n}.0d0' for n in range(60))
    source = (
        "  subroutine s\n"
        "    implicit none\n"
        f"    double precision, dimension(60), parameter :: table=[ &\n"
        f"         & {values} &\n"
        "         & ]\n"
        "  end subroutine s\n"
    )
    out = _format(source)
    # The statement comes back exactly as written — continuation lines, spacing
    # and all — rather than rebuilt from the parsed fields.
    for line in source.splitlines():
        if 'table=' in line or line.strip().startswith('&'):
            assert line in out.splitlines(), (line, out)


def test_generic_placeholder_statement_is_never_wrapped():
    """`Process.Generics` expands `{Type¦label}` placeholders one *physical
    line* at a time, emitting a line carrying one once per instance. A
    statement split across a continuation would expand into N copies of its
    first line followed by N copies of its second, and the trailing `&` of each
    copy of the first line would glue it to the *next* copy — fusing N
    statements into one. Such statements must come back exactly as written."""
    source = (
        "  type t\n"
        "   contains\n"
        "     generic :: readAttribute => IO_HDF5_Read_Attribute_{Type\u00a6label}_Scalar, "
        "IO_HDF5_Read_Attribute_{Type\u00a6label}_1D_Array_Allocatable\n"
        "     procedure :: shortOne => fooShortOne\n"
        "  end type t\n"
    )
    out = _format(source)
    generic = [l for l in out.splitlines() if 'generic' in l]
    assert len(generic) == 1, out
    assert not generic[0].rstrip().endswith('&'), generic[0]
    # Both specific procedures remain on that one line.
    assert 'IO_HDF5_Read_Attribute_{Type\u00a6label}_Scalar' in generic[0], generic[0]
    assert 'IO_HDF5_Read_Attribute_{Type\u00a6label}_1D_Array_Allocatable' in generic[0], generic[0]


def test_oversized_declaration_does_not_pad_its_neighbours():
    """A declaration left verbatim takes no part in sizing the columns — one
    huge data table used to pad every declaration beside it out to its width."""
    values = ', '.join(f'{n}.0d0' for n in range(60))
    source = (
        "  subroutine s\n"
        "    implicit none\n"
        f"    double precision, dimension(60), parameter :: table=[ &\n"
        f"         & {values} &\n"
        "         & ]\n"
        "    integer :: i\n"
        "  end subroutine s\n"
    )
    out = _format(source)
    short = next(l for l in out.splitlines() if l.strip().startswith('integer'))
    assert len(short) <= LINE_LENGTH_MAX, repr(short)


# ---------------------------------------------------------------------------
# 5. Idempotence
# ---------------------------------------------------------------------------

def test_openmp_declaration_indent_is_stable():
    """Unconditional declarations in a group that also has `!$` ones are
    indented three further columns.  The `!$` lines start with `!`, and
    skipping them as comments when measuring the block indent left only the
    padded lines to measure — so the block moved right on every pass."""
    source = (
        "  subroutine s\n"
        "    implicit none\n"
        "    integer :: i\n"
        "    !$ integer :: threadCount\n"
        "  end subroutine s\n"
    )
    once  = _format(source)
    twice = _format(once)
    assert twice == once, (once, twice)
    assert '\n    !$ integer' in once, once


def test_final_binding_survives():
    """A `final` binding names no target and is placed in the *target* column.
    Emitting the value only when the column carried an operator dropped these
    bindings from the output entirely."""
    source = (
        "  type t\n"
        "   contains\n"
        "     final :: fooDestructor\n"
        "     procedure :: value => fooValue\n"
        "  end type t\n"
    )
    out = _format(source)
    assert 'fooDestructor' in out, out
    assert _format(out) == out


def test_formatting_is_idempotent_over_a_mixed_block():
    source = (
        "  subroutine s(a, b)\n"
        "    implicit none\n"
        "    double precision, intent(in), dimension(:), allocatable :: a\n"
        "    class(fooClass), pointer :: b => null()\n"
        "    ! A comment inside the block.\n"
        "#ifdef USEMPI\n"
        "    integer :: rank\n"
        "#endif\n"
        "    logical :: flagOne=.false., flagTwo=.true.\n"
        "  end subroutine s\n"
    )
    once  = _format(source)
    twice = _format(once)
    assert twice == once, (once, twice)


# ---------------------------------------------------------------------------
# 6. The whole source tree
# ---------------------------------------------------------------------------

def _source_files():
    root = os.path.join(os.environ.get('GALACTICUS_EXEC_PATH', '.'), 'source')
    if not os.path.isdir(root):
        return []
    return sorted(
        os.path.join(directory, name)
        for directory, _, names in os.walk(root)
        for name in names
        if name.endswith(('.F90', '.Inc'))
    )


def _declaration_summary(tree):
    summary = []
    for node in walk_tree(tree):
        if node.get('type') != 'declaration':
            continue
        for declaration in node.get('declarations', []):
            type_spec = declaration.get('type')
            summary.append((
                declaration['intrinsic'],
                re.sub(r'\s+', '', type_spec).lower() if type_spec else '',
                frozenset(re.sub(r'\s+', '', a).lower()
                          for a in (declaration.get('attributes') or [])),
                tuple(re.sub(r'\s+', '', v).lower()
                      for v in declaration.get('variables', [])),
                bool(declaration.get('openMP')),
            ))
    return sorted(summary)


@pytest.mark.slow
def test_declaration_formatting_preserves_semantics_tree_wide():
    """Reformatting must never change what is declared.

    Attribute *order* is deliberately excluded from the comparison — the
    formatter canonicalizes it — but the intrinsic, type-spec, attribute set,
    variables and initializers must all survive untouched.
    """
    files = _source_files()
    assert files, 'no source files found — is GALACTICUS_EXEC_PATH set?'

    failures = []
    for path in files:
        with open(path, errors='replace') as fh:
            original = fh.read()
        tree   = parse_code_for_format(original, source=path)
        before = _declaration_summary(tree)
        format_declarations_tree(tree)
        formatted = serialize_for_format(tree['firstChild'])
        after = _declaration_summary(
            parse_code_for_format(formatted, source=path))
        if after != before:
            failures.append(path)
    assert not failures, '\n'.join(failures[:20])


@pytest.mark.slow
def test_declaration_formatting_is_idempotent_tree_wide():
    files = _source_files()
    assert files, 'no source files found — is GALACTICUS_EXEC_PATH set?'

    def once(text, path):
        tree = parse_code_for_format(text, source=path)
        format_declarations_tree(tree)
        return serialize_for_format(tree['firstChild'])

    failures = []
    for path in files:
        with open(path, errors='replace') as fh:
            original = fh.read()
        first = once(original, path)
        if once(first, path) != first:
            failures.append(path)
    assert not failures, '\n'.join(failures[:20])
