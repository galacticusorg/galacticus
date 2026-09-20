"""Tests for the run-time enforcement of the bounds declared by `inputParameter`
directives.

For each `<minimum>` or `<maximum>` element, `InputParameter.py` must emit a
call to `inputParameterBoundCheck` immediately after reading the parameter,
passing the bound as a double precision literal together with its declared
text, and must import `Input_Parameters_Bounds` into the enclosing
subprogram. A parameter which declares no bound must be left unchanged, and a
bound which is not a number must stop the build.
"""

import pytest

from Galacticus.Build.SourceTree import parse_code, serialize


def _process(source, tmp_path, monkeypatch):
    """Run the inputParameter process on `source`, returning the serialized result."""
    from Galacticus.Build.SourceTree.Process.InputParameter import process_input_parameters
    (tmp_path / 'Makefile_All_Execs').write_text("all_exes = Galacticus.exe\n")
    monkeypatch.setenv('BUILDPATH', str(tmp_path))
    tree = parse_code(source, name='test.F90')
    process_input_parameters(tree, {})
    return serialize(tree)


def _subroutine(elements):
    return (
        "subroutine foo(parameters)\n"
        "  implicit none\n"
        "  type(inputParameters), intent(inout) :: parameters\n"
        "  double precision :: fraction\n"
        "  !![\n"
        "  <inputParameter>\n"
        "    <name>fraction</name>\n"
        "    <source>parameters</source>\n"
        "    <defaultValue>0.5d0</defaultValue>\n"
        f"{elements}"
        "    <description>A fraction.</description>\n"
        "  </inputParameter>\n"
        "  !!]\n"
        "end subroutine foo\n"
    )


def test_bounds_emit_checks(tmp_path, monkeypatch):
    output = _process(
        _subroutine(
            "    <minimum inclusive=\"false\">0.0</minimum>\n"
            "    <maximum>1</maximum>\n"
        ),
        tmp_path, monkeypatch,
    )
    lines = output.splitlines()
    read  = next(i for i, line in enumerate(lines) if "call parameters%value('fraction'" in line)
    assert lines[read+1].strip() == (
        "call inputParameterBoundCheck(parameters,'fraction',fraction,0.0d0,'0.0',isMinimum=.true.,isInclusive=.false.)"
    )
    assert lines[read+2].strip() == (
        "call inputParameterBoundCheck(parameters,'fraction',fraction,1d0,'1',isMinimum=.false.,isInclusive=.true.)"
    )
    assert any(
        line.strip().startswith('use') and 'Input_Parameters_Bounds' in line and 'inputParameterBoundCheck' in line
        for line in lines
    )


def test_exponent_bound_becomes_double_literal(tmp_path, monkeypatch):
    output = _process(_subroutine("    <maximum>1.5e3</maximum>\n"), tmp_path, monkeypatch)
    assert ",1.5d3,'1.5e3',isMinimum=.false.," in output


def test_no_bounds_no_check(tmp_path, monkeypatch):
    output = _process(_subroutine(""), tmp_path, monkeypatch)
    assert 'inputParameterBoundCheck' not in output
    assert 'Input_Parameters_Bounds' not in output


def test_non_numeric_bound_is_rejected(tmp_path, monkeypatch):
    with pytest.raises(RuntimeError, match="is not a number"):
        _process(_subroutine("    <minimum>fractionMinimum</minimum>\n"), tmp_path, monkeypatch)
