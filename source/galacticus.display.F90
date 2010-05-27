!! Copyright 2009, Andrew Benson <abenson@caltech.edu>
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


!% Contains a module which implements outputting of formatted, indented messages at various vebosity levels from \glc.

module Galacticus_Display
  !% Implements outputting of formatted, indented messages at various vebosity levels from \glc.
  !$ use OMP_Lib
  use ISO_Varying_String
  private
  public :: Galacticus_Display_Message,Galacticus_Display_Indent,Galacticus_Display_Unindent,Galacticus_Verbosity_Level

  integer                                      :: maxThreads
  integer,           allocatable, dimension(:) :: indentationLevel
  character(len=10), allocatable, dimension(:) :: indentationFormat
  character(len=10), allocatable, dimension(:) :: indentationFormatNoNewLine

  logical                                      :: displayInitialized=.false.
  integer                                      :: verbosityLevel
  integer,           parameter, public         :: verbosityWarn=2, verbosityInfo=3, verbosityDebug=4

  interface Galacticus_Display_Message
     module procedure Galacticus_Display_Message_Char
     module procedure Galacticus_Display_Message_VarStr
  end interface

  interface Galacticus_Display_Indent
     module procedure Galacticus_Display_Indent_Char
     module procedure Galacticus_Display_Indent_VarStr
  end interface

contains

  integer function Galacticus_Verbosity_Level()
    !% Returns the verbositly level in \glc.
    implicit none
    
    Galacticus_Verbosity_Level=verbosityLevel
    return
  end function Galacticus_Verbosity_Level

  subroutine Initialize_Display
    !% Initialize the module by determining the requested verbosity level.
    use Input_Parameters
    implicit none

    !$omp critical (Initialize_Display)
    if (.not.displayInitialized) then
       ! Get the verbosity level parameter.
       !@ <inputParameter>
       !@   <name>verbosityLevel</name>
       !@   <defaultValue>1</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The level of verbosity for \glc\ (higher values give more verbosity).
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('verbosityLevel',verbosityLevel,1)

       ! For OpenMP runs, create an array of indentation levels.
       maxThreads=1
       !$ maxThreads=omp_get_max_threads()
       ! Do not use Alloc_Array() routines here as they may trigger recursive calls to this module which can lead to unbreakable
       ! parallel locks.
       allocate(indentationLevel          (maxThreads))
       allocate(indentationFormat         (maxThreads))
       allocate(indentationFormatNoNewLine(maxThreads))
       indentationLevel=0
       indentationFormat='(a)'
       indentationFormatNoNewLine='(a,$)'

       displayInitialized=.true.
    end if
    !$omp end critical (Initialize_Display)
    return
  end subroutine Initialize_Display

  subroutine Galacticus_Display_Indent_VarStr(message,verbosity)
    !% Increase the indentation level and display a message.
    implicit none
    type(varying_string), intent(in)           :: message
    integer,              intent(in), optional :: verbosity

    if (present(verbosity)) then
       call Galacticus_Display_Indent_Char(char(message),verbosity)
    else
       call Galacticus_Display_Indent_Char(char(message))
    end if
    return
  end subroutine Galacticus_Display_Indent_VarStr

  subroutine Galacticus_Display_Indent_Char(message,verbosity)
    !% Increase the indentation level and display a message.
    implicit none
    character(len=*), intent(in)           :: message
    integer,          intent(in), optional :: verbosity
    logical                                :: showMessage
    integer                                :: threadNumber

    !$omp critical(Galacticus_Message_Lock)
    call Initialize_Display
    if (present(verbosity)) then
       showMessage=(verbosity<=verbosityLevel)
    else
       showMessage=.true.
    end if
    if (showMessage) then
       !$ if (omp_in_parallel()) then
       !$    write (0,'(i2,a2,$)') omp_get_thread_num(),": "
       !$ else
       !$    write (0,'(a2,a2,$)') "MM",": "
       !$ end if
       threadNumber=1
       !$ if (omp_in_parallel()) threadNumber=omp_get_thread_num()+1
       write (0,indentationFormatNoNewLine(threadNumber)) '-> '
       write (0,'(a)') trim(message)
       !$ if (omp_in_parallel()) then
       !$    indentationLevel(omp_get_thread_num()+1)=indentationLevel(omp_get_thread_num()+1)+1
       !$ else
             indentationLevel=indentationLevel+1
       !$ end if
       call Create_Indentation_Format
    end if
    !$omp end critical(Galacticus_Message_Lock)
    return
  end subroutine Galacticus_Display_Indent_Char

  subroutine Galacticus_Display_Unindent(message,verbosity)
    !% Decrease the indentation level and display a message.
    implicit none
    character(len=*), intent(in)           :: message
    integer,          intent(in), optional :: verbosity
    integer                                :: threadNumber
    logical                                :: showMessage

    !$omp critical(Galacticus_Message_Lock)
    call Initialize_Display
    if (present(verbosity)) then
       showMessage=(verbosity<=verbosityLevel)
    else
       showMessage=.true.
    end if
    if (showMessage) then
       !$ if (omp_in_parallel()) then
       !$    indentationLevel(omp_get_thread_num()+1)=max(indentationLevel(omp_get_thread_num()+1)-1,0)
       !$ else
             indentationLevel=max(indentationLevel-1,0)
       !$ end if
       call Create_Indentation_Format
       !$ if (omp_in_parallel()) then
       !$    write (0,'(i2,a2,$)') omp_get_thread_num(),": "
       !$ else
       !$    write (0,'(a2,a2,$)') "MM",": "
       !$ end if
       threadNumber=1
       !$ if (omp_in_parallel()) threadNumber=omp_get_thread_num()+1
       write (0,indentationFormatNoNewLine(threadNumber)) '<- '
       write (0,'(a)') trim(message)
    end if
    !$omp end critical(Galacticus_Message_Lock)
    return
  end subroutine Galacticus_Display_Unindent

  subroutine Galacticus_Display_Message_Char(message,verbosity)
    !% Display a message (input as a {\tt character} variable).
    implicit none
    character(len=*), intent(in)           :: message
    integer,          intent(in), optional :: verbosity
    integer                                :: threadNumber
    logical                                :: showMessage

    !$omp critical(Galacticus_Message_Lock)
    call Initialize_Display
    if (present(verbosity)) then
       showMessage=(verbosity<=verbosityLevel)
    else
       showMessage=.true.
    end if
    if (showMessage) then
       !$ if (omp_in_parallel()) then
       !$    write (0,'(i2,a2,$)') omp_get_thread_num(),": "
       !$ else
       !$    write (0,'(a2,a2,$)') "MM",": "
       !$ end if
       threadNumber=1
       !$ if (omp_in_parallel()) threadNumber=omp_get_thread_num()+1
       write (0,indentationFormat(threadNumber)) trim(message)
    end if
    !$omp end critical(Galacticus_Message_Lock)
    return
  end subroutine Galacticus_Display_Message_Char

  subroutine Galacticus_Display_Message_VarStr(message,verbosity)
    !% Display a message (input as a {\tt varying\_string} variable).
    implicit none
    type(varying_string), intent(in)           :: message
    integer,              intent(in), optional :: verbosity
    integer                                    :: threadNumber
    logical                                    :: showMessage

    !$omp critical(Galacticus_Message_Lock)
    call Initialize_Display
    if (present(verbosity)) then
       showMessage=(verbosity<=verbosityLevel)
    else
       showMessage=.true.
    end if
    if (showMessage) then
       !$ if (omp_in_parallel()) then
       !$    write (0,'(i2,a2,$)') omp_get_thread_num(),": "
       !$ else
       !$    write (0,'(a2,a2,$)') "MM",": "
       !$ end if
       threadNumber=1
       !$ if (omp_in_parallel()) threadNumber=omp_get_thread_num()+1
       write (0,indentationFormat(threadNumber)) char(message)
    end if
    !$omp end critical(Galacticus_Message_Lock)
    return
  end subroutine Galacticus_Display_Message_VarStr

  subroutine Create_Indentation_Format
    !% Create a format for indentation.
    implicit none
    integer,    parameter :: indentSpaces=4
    integer               :: threadNumber

    threadNumber=1
    !$ if (omp_in_parallel()) threadNumber=omp_get_thread_num()+1
    select case (indentationLevel(threadNumber)*indentSpaces)
    case (0)
       write (indentationFormat(threadNumber),'(a)') '(a)'
       write (indentationFormatNoNewLine(threadNumber),'(a)') '(a,$)'
    case (1:9)
       write (indentationFormat(threadNumber),'(a,i1,a)') '(',indentationLevel(threadNumber)*indentSpaces,'x,a)'
       write (indentationFormatNoNewLine(threadNumber),'(a,i1,a)') '(',indentationLevel(threadNumber)*indentSpaces,'x,a,$)'
    case (10:99)
       write (indentationFormat(threadNumber),'(a,i2,a)') '(',indentationLevel(threadNumber)*indentSpaces,'x,a)'
       write (indentationFormatNoNewLine(threadNumber),'(a,i2,a)') '(',indentationLevel(threadNumber)*indentSpaces,'x,a,$)'
    end select
    !$ if (.not.omp_in_parallel()) then
    !$    indentationLevel          =indentationLevel(1)
    !$    indentationFormat         =indentationFormat(1)
    !$    indentationFormatNoNewLine=indentationFormatNoNewLine(1)
    !$ end if
    return
  end subroutine Create_Indentation_Format

end module Galacticus_Display
