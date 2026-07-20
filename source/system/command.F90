!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
!!    Andrew Benson <abenson@carnegiescience.edu>
!!
!! This file is part of Galacticus.
!!
!!    Galacticus is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License as published by
!!    the Free Software Foundation, either version 3 of the License, or
!!    (at your option) any later version.
!!
!!    Galacticus is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with Galacticus.  If not, see <http://www.gnu.org/licenses/>.

!!{RST
Contains a module which executes system commands.
!!}

module System_Command
  !!{RST
  Executes system commands.
  !!}
  implicit none
  private
  public :: System_Command_Do, shellEscape

  interface System_Command_Do
     module procedure System_Command_Char
     module procedure System_Command_VarStr
  end interface System_Command_Do

  interface shellEscape
     module procedure shellEscapeChar
     module procedure shellEscapeVarStr
  end interface shellEscape

contains

  function shellEscapeVarStr(token) result(escaped)
    !!{RST
    Return the string ``token`` wrapped in single quotes, with any embedded single quotes escaped, such that it can be used
    safely as a single token (e.g. a file system path) in a shell command even if it contains spaces or other characters
    that are special to the shell. This should be applied to any path substituted into a command passed to
    ``System_Command_Do``.

    Note that the result of ``shellEscape`` must first be assigned to a local variable, which is then used in building the
    command string. The result of ``shellEscape`` must {\normalfont \bfseries not} be used directly as an actual argument (for
    example inline within a ``System_Command_Do`` argument, or as an operand of a concatenation that is passed as an actual
    argument), as gfortran leaks the function result in that case.

    For the same reason, the {\normalfont \ttfamily token} argument passed to ``shellEscape`` must not contain a call to a
    function that returns a {\normalfont \ttfamily varying\_string} (such as {\normalfont \ttfamily inputPath}), whether as the
    bare argument ({\normalfont \ttfamily shellEscape(inputPath(...))}) or nested within a concatenation
    ({\normalfont \ttfamily shellEscape(inputPath(...)//"...")}), since gfortran leaks that nested function result. Assign the
    argument to a local variable first, then escape the variable, e.g. {\normalfont \ttfamily localVar=inputPath(...)//"..."}
    followed by {\normalfont \ttfamily localVar=shellEscape(localVar)}. Building the value in an assignment is safe; only a
    {\normalfont \ttfamily varying\_string}-returning function call that appears inside another call's argument list leaks.
    !!}
    use :: ISO_Varying_String, only : varying_string, replace, operator(//)
    implicit none
    type(varying_string), intent(in   ) :: token
    type(varying_string)                :: escaped

    ! Wrap in single quotes, replacing each embedded single quote with the sequence '\'' (close quote, escaped quote, reopen
    ! quote), which is the standard POSIX-shell-safe encoding.
    escaped="'"//replace(token,"'","'\''",every=.true.)//"'"
    return
  end function shellEscapeVarStr

  function shellEscapeChar(token) result(escaped)
    !!{RST
    Return the string ``token`` wrapped in single quotes, with any embedded single quotes escaped. See {\normalfont \ttfamily
    shellEscapeVarStr}.
    !!}
    use :: ISO_Varying_String, only : varying_string, var_str
    implicit none
    character(len=*), intent(in   ) :: token
    type(varying_string)            :: escaped

    escaped=shellEscapeVarStr(var_str(token))
    return
  end function shellEscapeChar

  subroutine System_Command_VarStr(command,iStatus)
    !!{RST
    Executes the system command ``command``, optionally returning the resulting status in ``iStatus``.
    !!}
    use :: ISO_Varying_String, only : char, varying_string
    implicit none
    type   (varying_string), intent(in   )           :: command
    integer                , intent(  out), optional :: iStatus

    call System_Command_Char(char(command),iStatus)
    return
  end subroutine System_Command_VarStr

  subroutine System_Command_Char(command,iStatus)
    !!{RST
    Executes the system command ``command``, optionally returning the resulting status in ``iStatus``.
    !!}
    use :: Error, only : Error_Report
    implicit none
    character(len=*), intent(in   )           :: command
    integer         , intent(  out), optional :: iStatus
    integer                                   :: iStatusActual

    call System(command,iStatusActual)
    if (present(iStatus)) then
       iStatus=iStatusActual
    else
       if (iStatusActual /= 0) call Error_Report('failed to execute system command:'//char(10)//' --> '//command//{introspection:location})
    end if
    return
  end subroutine System_Command_Char

end module System_Command
