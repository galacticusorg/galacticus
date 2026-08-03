"""Tests that the build parser loudly rejects parameterized derived types (PDTs).

Galacticus does not support PDTs: the parser and the generators that consume
its output (state storage, deep copy, digests, module dependencies) have no
representation for `kind`/`len` type parameters, and would silently drop a PDT
from generated code.  `check_no_parameterized_derived_type` turns that silent
misparse into a clear error, and the central `parse_code` loop calls it on
every logical line.  See issue #114.
"""

import pytest

from Fortran.Utils import check_no_parameterized_derived_type
from Galacticus.Build.SourceTree import parse_code


# --- Lines that ARE parameterized derived types (must raise). ---------------

PDT_LINES = [
    # Definition openers, with and without `::` and attributes.
    "type :: foo(k, n)",
    "type foo(k, n)",
    "  type, public :: foo(k)",
    "type, extends(base) :: foo(kind, len)",
    "type, abstract, extends(base) :: foo(k)",
    "TYPE :: FOO(K)",                     # case-insensitive
    # Parameterized declarations.
    "type(foo(real64, 3)) :: x",
    "  class(foo(real64)) :: y",
    "type(matrix(real64, n, m)), pointer :: p",
]


@pytest.mark.parametrize("line", PDT_LINES)
def test_pdt_lines_are_rejected(line):
    with pytest.raises(ValueError, match="parameterized derived types"):
        check_no_parameterized_derived_type(line)


# --- Lines that are NOT PDTs (must pass through unchanged). ------------------

NON_PDT_LINES = [
    # Ordinary derived-type definitions.
    "type :: foo",
    "type foo",
    "  type, public :: foo",
    "type, extends(base) :: foo",
    "type, abstract :: foo",
    # Ordinary declarations (single bare type name inside the parens).
    "type(foo) :: x",
    "type(foo), pointer :: p",
    "class(foo) :: y",
    "class(*) :: z",
    "type(varying_string) :: s",
    # Select-type / associate constructs that superficially resemble a PDT.
    "type is (foo)",
    "  type is (integer)",
    "class is (foo)",
    "class default",
    # Unrelated statements.
    "integer :: i",
    "double precision :: c(6)",
    "call foo(a, b)",
]


@pytest.mark.parametrize("line", NON_PDT_LINES)
def test_non_pdt_lines_pass_through(line):
    # Returns the line unchanged and does not raise.
    assert check_no_parameterized_derived_type(line) == line


# --- End-to-end through the parser. -----------------------------------------

def test_parse_code_rejects_pdt_definition():
    source = (
        "module foo\n"
        "  implicit none\n"
        "  type :: matrix(k, n)\n"
        "     integer, kind :: k\n"
        "     integer, len  :: n\n"
        "     real(k) :: values(n, n)\n"
        "  end type matrix\n"
        "end module foo\n"
    )
    with pytest.raises(ValueError, match="parameterized derived types"):
        parse_code(source, name="<test>", instrument=False)


def test_parse_code_rejects_pdt_declaration():
    source = (
        "subroutine s()\n"
        "  type(matrix(real64, 3)) :: m\n"
        "end subroutine s\n"
    )
    with pytest.raises(ValueError, match="parameterized derived types"):
        parse_code(source, name="<test>", instrument=False)


def test_parse_code_accepts_ordinary_derived_type():
    source = (
        "module foo\n"
        "  implicit none\n"
        "  type, public :: point\n"
        "     double precision :: position(3)\n"
        "  end type point\n"
        "  type(point) :: origin\n"
        "end module foo\n"
    )
    # Must not raise.
    tree = parse_code(source, name="<test>", instrument=False)
    assert tree is not None
