"""Tests for scripts/build/consolidateUses.py, which moves procedure-level USE statements of a preprocessed module or
submodule to its specification part where that cannot change the meaning of any name.

The structural tests check what is (and is not) moved. Where gfortran is available, the compile tests check that the
consolidated source compiles, and that a program using it behaves exactly as the original does.
"""

import os
import shutil
import subprocess
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir))
import consolidateUses  # noqa: E402
from consolidateUses import Unsupported, consolidate  # noqa: E402

SCRIPT = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, 'consolidateUses.py')


def run(source, parent=None):
    lines = source.strip('\n').split('\n')
    output, removed, moved = consolidate(lines, parent)
    assert len(output) == len(lines)
    return '\n'.join(output), removed, moved


def procedureUses(text):
    """USE statements left after the CONTAINS of the unit."""
    body = text.lower().split('\ncontains', 1)[1]
    return [line.strip() for line in body.split('\n') if line.strip().startswith('use')]


MODULE = """
module user
  implicit none
  private
  public :: a, b
contains
  integer function a()
    use :: lib, only : f, t
    type(t) :: x
    a=f(x)
  end function a
  integer function b()
    use :: lib, only : f
    b=f()
  end function b
end module user
"""


def test_moves_shared_import():
    text, removed, moved = run(MODULE)
    assert text.split('\n')[0] == 'module user ; use :: lib, only : f, t'
    assert procedureUses(text) == []
    assert (removed, moved) == (2, 2)


def test_line_count_and_numbers_preserved():
    original = MODULE.strip('\n').split('\n')
    text, _, _ = run(MODULE)
    lines = text.split('\n')
    for i, line in enumerate(original):
        if 'use ::' not in line and i > 0:
            assert lines[i] == line


def test_blocked_when_another_procedure_mentions_name_without_import():
    # In b, "size" is the intrinsic: moving lib's "size" to the module would change its meaning there.
    source = MODULE.replace("use :: lib, only : f, t", "use :: lib, only : f, t, size").replace(
        "    use :: lib, only : f\n    b=f()\n", "    b=size([1,2])\n")
    text, _, _ = run(source)
    assert procedureUses(text) == ['use :: lib, only : size']
    assert text.split('\n')[0] == 'module user ; use :: lib, only : f, t'


def test_blocked_when_another_procedure_declares_name():
    # b declares its own "f". (Moving would in fact be safe, as the local declaration hides the moved name, but a
    # declaration cannot always be told from a reference, so this is treated in the same way.)
    source = MODULE.replace("    use :: lib, only : f\n    b=f()\n", "    integer :: f\n    f=1\n    b=f\n")
    text, _, _ = run(source)
    assert procedureUses(text) == ['use :: lib, only : f']


def test_blocked_by_different_import_of_same_name():
    source = MODULE.replace("use :: lib, only : f\n", "use :: other, only : f\n")
    text, _, _ = run(source)
    assert 'use :: lib, only : f' in procedureUses(text)
    assert 'use :: other, only : f' in procedureUses(text)


def test_blocked_by_specification_part():
    source = MODULE.replace("  public :: a, b\n", "  public :: a, b, f\n")
    text, _, _ = run(source)
    assert 'use :: lib, only : f' in ' '.join(procedureUses(text))


def test_import_of_a_declared_name_is_left_for_the_compiler_to_report():
    # Importing a name and declaring it in the same procedure is an error: moving the import would hide it.
    source = MODULE.replace("    type(t) :: x\n", "    type(t) :: x\n    integer :: y, f(3)\n")
    text, _, _ = run(source)
    assert any('only : f' in use for use in procedureUses(text))


def test_internal_procedure_sees_host_import():
    source = """
module user
  implicit none
  private
contains
  subroutine a()
    use :: lib, only : f
    call g()
  contains
    subroutine g()
      print *, f()
    end subroutine g
  end subroutine a
  subroutine b()
    use :: lib, only : f
    print *, f()
  end subroutine b
end module user
"""
    text, removed, moved = run(source)
    assert procedureUses(text) == []
    assert moved == 1


def test_assignment_and_intrinsic_operators_not_moved():
    source = """
module user
  implicit none
  private
contains
  subroutine a()
    use :: iso_varying_string, only : varying_string, assignment(=), operator(//), operator(.cat.)
    type(varying_string) :: s
    s='a'//'b'
    s=s.cat.s
  end subroutine a
end module user
"""
    text, _, _ = run(source)
    assert procedureUses(text) == ['use :: iso_varying_string, only : assignment(=), operator(//)']
    assert 'only : operator(.cat.), varying_string' in text.split('\n')[0]


def test_conditional_sentinel_interface_and_unrestricted_uses_not_moved():
    source = """
module user
  implicit none
  private
contains
  subroutine a()
#ifdef USEMPI
    use :: mpi, only : mpi_comm_world
#endif
    !$ use :: omp_lib, only : omp_get_thread_num
    use :: error
    interface
      subroutine c(x)
        use :: lib, only : t
        type(t) :: x
      end subroutine c
    end interface
  end subroutine a
end module user
"""
    text, removed, moved = run(source)
    assert (removed, moved) == (0, 0)
    assert text == source.strip('\n')


def test_public_module_gets_private_statement():
    source = MODULE.replace("  private\n", "")
    text, _, _ = run(source)
    assert 'private :: f, t ; contains' in text.split('\n')


def test_names_provided_by_unit_are_dropped_not_moved():
    source = MODULE.replace("  implicit none\n", "  use :: lib, only : f\n  implicit none\n")
    text, removed, moved = run(source)
    assert text.split('\n')[0] == 'module user ; use :: lib, only : t'
    assert procedureUses(text) == []


SUBMODULE = """
submodule (parent) child
  implicit none
  type(t), pointer :: p
contains
  module procedure a
    use :: lib, only : g, t
    a=g(p)
  end procedure a
  module procedure b
    use :: lib, only : g
    b=g(p)
  end procedure b
end submodule child
"""


def test_submodule_drops_names_the_parent_provides():
    text, removed, moved = run(SUBMODULE, parent={('lib', 'g', 'g'), ('lib', 't', 't')})
    assert procedureUses(text) == []
    assert moved == 0 and removed == 2
    assert text.split('\n')[0] == 'submodule (parent) child'


def test_submodule_drops_names_of_parents_unrestricted_use():
    text, _, moved = run(SUBMODULE, parent={('lib', '*', '*')})
    assert procedureUses(text) == [] and moved == 0


def test_submodule_type_name_in_specification_blocks_move():
    # The specification part refers to "t", which (not being provided by the parent as lib's "t") could be a different
    # entity: moving lib's "t" there could change it.
    text, _, _ = run(SUBMODULE, parent={('lib', 'g', 'g')})
    assert 'use :: lib, only : t' in procedureUses(text)


def test_submodule_without_parent_is_unsupported():
    with pytest.raises(Unsupported):
        consolidate(SUBMODULE.strip('\n').split('\n'), None)


def test_type_private_statement_is_not_default_private():
    source = MODULE.replace("  private\n  public :: a, b\n", "  type :: u\n    private\n    integer :: i\n  end type u\n")
    text, _, _ = run(source)
    assert any(line.strip().startswith('private :: f, t') for line in text.split('\n'))


def test_programs_are_unsupported():
    with pytest.raises(Unsupported):
        consolidate(['program p', 'contains', 'end program p'])


def test_mismatched_procedure_end_is_unsupported():
    with pytest.raises(Unsupported):
        consolidate(MODULE.replace('end function b', 'end function c').strip('\n').split('\n'))


def test_script_leaves_unsupported_files_unchanged(tmp_path):
    source = tmp_path / 'p.p.F90'
    source.write_text('program p\n  use :: lib\nend program p\n')
    output = tmp_path / 'p.pc.F90'
    result = subprocess.run([sys.executable, SCRIPT, str(source), str(output), '--report'], capture_output=True,
                            text=True, check=True)
    assert 'unchanged' in result.stdout
    assert output.read_text() == source.read_text()


def test_find_parent_source(tmp_path):
    parentSource = tmp_path / 'build' / 'dir' / 'parent.p.F90'
    parentSource.parent.mkdir(parents=True)
    parentSource.write_text('module parent\nend module parent\n')
    dependencies = tmp_path / 'Makefile_Module_Dependencies'
    dependencies.write_text(f'{tmp_path}/build/parent.mod: {tmp_path}/build/dir/parent.o\n')
    assert consolidateUses.findParentSource('Parent', str(dependencies)) == str(parentSource)
    assert consolidateUses.findParentSource('other', str(dependencies)) is None
    # The parent's consolidated source, if any, is what was compiled, and its specification part includes the names
    # moved there.
    consolidatedSource = tmp_path / 'build' / 'dir' / 'parent.pc.F90'
    consolidatedSource.write_text('module parent ; use :: lib, only : g\nend module parent\n')
    assert consolidateUses.findParentSource('Parent', str(dependencies)) == str(consolidatedSource)
    assert ('lib', 'g', 'g') in consolidateUses.specificationImports(consolidatedSource.read_text().split('\n'))


# Compile tests: a library module, and a user module whose procedures each USE it, exercising names that are moved, a
# name blocked by a local variable of the same name, an internal procedure, a defined operator, and an intrinsic module.
LIBRARY = """
module lib
  implicit none
  type :: t
     integer :: i=3
  end type t
  interface operator(.twice.)
     module procedure twice
  end interface operator(.twice.)
contains
  integer function f(x)
    type(t), intent(in) :: x
    f=x%i
  end function f
  integer function twice(i)
    integer, intent(in) :: i
    twice=2*i
  end function twice
end module lib
"""

USER = """
module user
  implicit none
  public
contains
  integer function a()
    use :: lib, only : f, t
    type(t) :: x
    a=f(x)
  end function a
  integer function b()
    use :: lib, only : t, operator(.twice.)
    type(t) :: x
    b=.twice.x%i
  end function b
  integer function c()
    integer :: f
    f=7
    c=f
  end function c
  integer function d()
    use :: lib, only : f, t
    d=g()
  contains
    integer function g()
      type(t) :: y
      g=f(y)+1
    end function g
  end function d
  integer function e()
    use, intrinsic :: iso_fortran_env, only : int64
    integer(int64) :: k
    k=5
    e=int(k)
  end function e
end module user
"""

PROGRAM = """
program main
  use :: user
  implicit none
  print '(5i4)', a(), b(), c(), d(), e()
end program main
"""


@pytest.mark.skipif(shutil.which('gfortran') is None, reason='gfortran is not available')
def test_consolidated_source_compiles_and_behaves_identically(tmp_path):
    results = {}
    for variant in ('original', 'consolidated'):
        directory = tmp_path / variant
        directory.mkdir()
        (directory / 'lib.F90').write_text(LIBRARY)
        (directory / 'main.F90').write_text(PROGRAM)
        user = directory / 'user.F90'
        if variant == 'original':
            user.write_text(USER)
        else:
            (directory / 'user.p.F90').write_text(USER)
            subprocess.run([sys.executable, SCRIPT, str(directory / 'user.p.F90'), str(user)], check=True)
            text = user.read_text()
            assert 'module user ; use' in text
            assert 'private ::' in text
            assert 'use :: lib, only : f' in text.split('contains', 1)[1].lower() or 'integer :: f' in text
        for name in ('lib', 'user', 'main'):
            subprocess.run(['gfortran', '-c', '-Wall', '-Werror', f'{name}.F90'], cwd=directory, check=True)
        subprocess.run(['gfortran', 'lib.o', 'user.o', 'main.o', '-o', 'main'], cwd=directory, check=True)
        results[variant] = subprocess.run(['./main'], cwd=directory, capture_output=True, text=True,
                                          check=True).stdout
    assert results['original'] == results['consolidated'] == '   3   6   7   4   5\n'
