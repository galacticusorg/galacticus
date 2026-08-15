"""Tests for reformatting Fortran source in place.

Reformatting differs from the build's code generation in one demanding
respect: whatever the formatter does not deliberately change, it must leave
byte-identical.  `parse_code()` does not offer that — it is free to rewrite
anything the compiler will not notice, and does — so formatters go through
`parse_code_for_format()` / `serialize_for_format()` instead.

Covered here:

  1. `uncomment_embedded()` inverts `_comment_embedded()`, which prefixes the
     body of every `!!{ … !!}` and `!![ … !!]` block with `"!< "`.  Without an
     inverse, serializing a parsed tree back to source rewrote every directive
     and documentation block in the file.

  2. `parse_code_for_format()` suppresses the other two rewrites: the
     `{introspection:location}` line-number instrumentation, and the
     directives pass, which re-emits a `!![`/`!!]` marker pair around each
     individual directive and so split one multi-directive block into several.

  3. The `use`-statement line-length ruler: `update_uses()` emits
     `SYMBOLS_PER_ROW_MAX` symbols per row, but drops to fewer when that would
     push a row past `LINE_LENGTH_MAX`.

  4. Idempotence.  Two bugs broke it, and both are easy to reintroduce: the
     indent was read off the first line of the block, picking up the OpenMP
     padding the previous pass had added; and intrinsic modules were hoisted
     to index zero, which reversed them relative to each other on every pass.
"""

import os
import re

import pytest

from Galacticus.Build.SourceTree import (
    _comment_embedded,
    uncomment_embedded,
    parse_code_for_format,
    serialize_for_format,
    walk_tree,
)
from Galacticus.Build.SourceTree.Parse.ModuleUses import (
    LINE_LENGTH_MAX,
    SYMBOLS_PER_ROW_MAX,
    update_uses,
)


def _format(source):
    """Apply the `use`-block formatter to `source` and return the result."""
    tree = parse_code_for_format(source, name='<test>')
    for node in walk_tree(tree):
        if node.get('type') == 'moduleUse':
            update_uses(node)
    return serialize_for_format(tree['firstChild'])


def _use_lines(text):
    return [line for line in text.splitlines()
            if re.match(r'^\s*(!\$\s*)?(&|use(?![a-zA-Z0-9_]))', line)]


def _symbols_in_row(line):
    """Count the `only` symbols on one physical line of a `use` statement.

    Everything up to and including `, only : ` is the statement's own prefix,
    whose comma must not be counted; continuation rows carry only symbols.  A
    trailing `&` marks a row that continues, so the last symbol on it is
    followed by a comma like the others.
    """
    body = line.split(", only : ", 1)[-1].lstrip("& ").rstrip()
    return len([symbol for symbol in body.rstrip("&").split(",") if symbol.strip()])


# ---------------------------------------------------------------------------
# 1. uncomment_embedded
# ---------------------------------------------------------------------------

def test_uncomment_embedded_inverts_comment_embedded():
    source = (
        "module foo\n"
        "  !!{RST\n"
        "  Documentation body.\n"
        "  !!}\n"
        "  !![\n"
        "  <directive name=\"a\"/>\n"
        "  <directive name=\"b\"/>\n"
        "  !!]\n"
        "  implicit none\n"
        "end module foo\n"
    )
    commented = _comment_embedded(source)
    # Sanity: the transform really did something to the block bodies.
    assert "!<   Documentation body.\n" in commented, commented
    assert uncomment_embedded(commented) == source


def test_uncomment_embedded_leaves_plain_source_alone():
    """Lines outside a `!!{`/`!![` block never carry the marker, so a file
    with no embedded blocks must pass through untouched."""
    source = (
        "subroutine bar\n"
        "  ! An ordinary comment.\n"
        "  integer :: i\n"
        "end subroutine bar\n"
    )
    assert uncomment_embedded(source) == source


def test_uncomment_embedded_ignores_nested_markers():
    """`_comment_embedded()` prefixes a `!!{` that appears *inside* a block
    rather than treating it as opening a second one.  The inverse has to make
    the same call, or its notion of where the block ends drifts out of step
    and it strips markers off lines that never had them."""
    source = (
        "  !![\n"
        "  <description>\n"
        "  !!{\n"
        "  </description>\n"
        "  !!]\n"
        "  integer :: i\n"
    )
    assert uncomment_embedded(_comment_embedded(source)) == source


# ---------------------------------------------------------------------------
# 2. Round-trip fidelity of the format-oriented parse
# ---------------------------------------------------------------------------

def test_directive_block_survives_round_trip():
    """The directives pass emits one `!![ … !!]` pair per directive, so a
    block holding three of them came back as three blocks.  The
    format-oriented parse must not run that pass."""
    source = (
        "  function f()\n"
        "    !![\n"
        "    <objectDestructor name=\"a_\"/>\n"
        "    <objectDestructor name=\"b_\"/>\n"
        "    <objectDestructor name=\"c_\"/>\n"
        "    !!]\n"
        "    return\n"
        "  end function f\n"
    )
    out = serialize_for_format(parse_code_for_format(source)['firstChild'])
    assert out == source
    assert out.count("!![") == 1, out


def test_introspection_location_not_instrumented():
    """`parse_code(instrument=True)` tags each placeholder with its line
    number, turning `…:location}` into `…:location:3}`."""
    source = (
        "  subroutine s\n"
        "    implicit none\n"
        "    call Error_Report('bad'//{introspection:location})\n"
        "  end subroutine s\n"
    )
    out = serialize_for_format(parse_code_for_format(source)['firstChild'])
    assert out == source
    assert "{introspection:location}" in out, out


def test_round_trip_preserves_unformatted_use_block():
    """Nothing may change until a formatter asks for it — including a `use`
    block whose layout does not match the house style."""
    source = (
        "module foo\n"
        "  use::Error,only:Error_Report\n"
        "  use :: Display, only : displayMessage\n"
        "  implicit none\n"
        "end module foo\n"
    )
    out = serialize_for_format(parse_code_for_format(source)['firstChild'])
    assert out == source


# ---------------------------------------------------------------------------
# 3. The line-length ruler
# ---------------------------------------------------------------------------

def test_short_symbols_use_full_row():
    """A block whose rows fit comfortably keeps the full complement of
    symbols per row."""
    source = (
        "module foo\n"
        "  use :: M, only : a, b, c, d, e\n"
        "  implicit none\n"
        "end module foo\n"
    )
    out = _format(source)
    first = _use_lines(out)[0]
    assert _symbols_in_row(first) == SYMBOLS_PER_ROW_MAX, out
    assert first.rstrip().endswith("&"), out


def test_long_symbols_reduce_symbols_per_row():
    """Four symbols per row would run past the limit here, so the formatter
    must fall back to fewer — and every emitted row must then fit."""
    symbols = [f"aVeryLongSymbolNameIndeed{i:02d}" * 2 for i in range(8)]
    source = (
        "module foo\n"
        f"  use :: M, only : {', '.join(symbols)}\n"
        "  implicit none\n"
        "end module foo\n"
    )
    out = _format(source)
    lines = _use_lines(out)
    assert all(len(line) <= LINE_LENGTH_MAX for line in lines), lines
    # Fewer than the maximum symbols per row, but still more than one.
    assert 1 < max(_symbols_in_row(line) for line in lines) < SYMBOLS_PER_ROW_MAX, lines
    # No symbol was lost.
    for symbol in symbols:
        assert symbol in out, symbol


def test_single_oversized_symbol_overflows_rather_than_breaking():
    """A symbol that cannot fit on a line of its own is emitted anyway —
    there is nowhere to break it, and the build passes
    `-ffree-line-length-none`, so the limit is a style rule."""
    symbol = "x" * (LINE_LENGTH_MAX + 20)
    source = (
        "module foo\n"
        f"  use :: M, only : {symbol}, y\n"
        "  implicit none\n"
        "end module foo\n"
    )
    out = _format(source)
    assert symbol in out, out
    # One symbol per row: the long one overflows, the short one does not.
    assert max(len(line) for line in _use_lines(out)) > LINE_LENGTH_MAX


def test_ruler_keeps_columns_aligned_across_statements():
    """Symbol columns are shared across the whole block, so the count per row
    is chosen once for the block rather than per statement."""
    source = (
        "module foo\n"
        "  use :: Alpha, only : one, two, three, four, five\n"
        "  use :: Beta , only : six, seven\n"
        "  implicit none\n"
        "end module foo\n"
    )
    out = _format(source)
    lines = [line for line in out.splitlines() if ", only : " in line]
    columns = {line.index(", only : ") for line in lines}
    assert len(columns) == 1, lines


# ---------------------------------------------------------------------------
# 4. Idempotence
# ---------------------------------------------------------------------------

def test_openmp_block_indent_is_stable():
    """Unconditional uses in a block that also has `!$ use` lines are indented
    three columns further, so that `use` lines up across the sentinel.  Reading
    the block indent off its first line picked that padding up and added three
    more columns on every pass."""
    source = (
        "  subroutine s\n"
        "    use :: Error, only : Error_Report\n"
        "    !$ use :: OMP_Lib, only : OMP_Get_Max_Threads\n"
        "    implicit none\n"
        "  end subroutine s\n"
    )
    once  = _format(source)
    twice = _format(once)
    assert twice == once, (once, twice)
    # And the sentinel line still sits at the original indent.
    assert "\n    !$ use :: OMP_Lib" in once, once


def test_multiple_intrinsic_modules_keep_relative_order():
    """Intrinsic modules are hoisted ahead of the rest but must keep their own
    order.  Inserting each at index zero reversed them, so reformatting
    oscillated between two layouts and never reached a fixed point."""
    source = (
        "module foo\n"
        "  use, intrinsic :: ISO_C_Binding   , only : c_size_t\n"
        "  use, intrinsic :: ISO_Fortran_Env , only : output_unit\n"
        "  use            :: Display         , only : displayMessage\n"
        "  implicit none\n"
        "end module foo\n"
    )
    once  = _format(source)
    twice = _format(once)
    assert twice == once, (once, twice)
    assert once.index("ISO_C_Binding") < once.index("ISO_Fortran_Env"), once
    # Intrinsics still come first.
    assert once.index("ISO_Fortran_Env") < once.index("Display"), once


def test_formatting_is_idempotent_over_mixed_block():
    source = (
        "  subroutine s\n"
        "    use, intrinsic :: ISO_C_Binding, only : c_size_t, c_int\n"
        "    use :: Error, only : Error_Report\n"
        "#ifdef USEMPI\n"
        "    use :: MPI_Utilities, only : mpiSelf, mpiBarrier\n"
        "#endif\n"
        "    !$ use :: OMP_Lib, only : OMP_Get_Max_Threads\n"
        "    implicit none\n"
        "  end subroutine s\n"
    )
    once  = _format(source)
    twice = _format(once)
    assert twice == once, (once, twice)
    # The preprocessor guard survives, exactly once.
    assert once.count("#ifdef USEMPI") == 1, once
    assert once.count("#endif") == 1, once


# ---------------------------------------------------------------------------
# 5. The whole source tree
# ---------------------------------------------------------------------------

def _source_files():
    """Every Fortran source file in the repository's `source/` tree."""
    root = os.path.join(os.environ.get('GALACTICUS_EXEC_PATH', '.'), 'source')
    if not os.path.isdir(root):
        return []
    return sorted(
        os.path.join(directory, name)
        for directory, _, names in os.walk(root)
        for name in names
        if name.endswith(('.F90', '.Inc'))
    )


@pytest.mark.slow
def test_every_source_file_round_trips_unchanged():
    """Parsing and reserializing any file in the tree must return it byte for
    byte.  Fixtures cannot stand in for this: the failures worth catching are
    the constructs nobody thought to write a fixture for.

    A file that fails here is not merely formatted oddly — it means the parser
    cannot represent something the file contains, so reformatting it would
    silently drop or mangle that construct.
    """
    files = _source_files()
    assert files, 'no source files found — is GALACTICUS_EXEC_PATH set?'

    failures = []
    for path in files:
        with open(path, errors='replace') as fh:
            original = fh.read()
        try:
            out = serialize_for_format(
                parse_code_for_format(original, source=path)['firstChild'])
        except Exception as error:                      # noqa: BLE001
            failures.append(f'{path}: raised {error!r}')
            continue
        if out != original:
            failures.append(f'{path}: round-trip differs')
    assert not failures, '\n'.join(failures[:20])


@pytest.mark.slow
def test_formatting_every_source_file_is_idempotent():
    """Reformatting a file that is already formatted must be a no-op.

    Without this, two separate bugs went unnoticed — the OpenMP indent grew by
    three columns per pass, and pairs of intrinsic modules swapped places on
    every pass — because both need a *second* run to become visible.
    """
    files = _source_files()
    assert files, 'no source files found — is GALACTICUS_EXEC_PATH set?'

    failures = []
    for path in files:
        with open(path, errors='replace') as fh:
            original = fh.read()
        once = _format(original)
        if _format(once) != once:
            failures.append(path)
    assert not failures, '\n'.join(failures[:20])
