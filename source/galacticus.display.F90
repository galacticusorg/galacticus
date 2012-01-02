!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


!% Contains a module which implements outputting of formatted, indented messages at various vebosity levels from \glc.

module Galacticus_Display
  !% Implements outputting of formatted, indented messages at various vebosity levels from \glc.
  !$ use OMP_Lib
  use ISO_Varying_String
  implicit none
  private
  public :: Galacticus_Display_Message,Galacticus_Display_Indent,Galacticus_Display_Unindent,Galacticus_Verbosity_Level&
       &,Galacticus_Verbosity_Level_Set, Galacticus_Display_Counter, Galacticus_Display_Counter_Clear

  integer                                      :: maxThreads
  integer,           allocatable, dimension(:) :: indentationLevel
  character(len=10), allocatable, dimension(:) :: indentationFormat
  character(len=10), allocatable, dimension(:) :: indentationFormatNoNewLine

  logical                                      :: displayInitialized=.false.
  integer                                      :: verbosityLevel  =1
  integer,           parameter, public         :: verbosityWorking=2, &
       &                                          verbosityWarn   =3, &
       &                                          verbosityInfo   =4, &
       &                                          verbosityDebug  =5

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

  subroutine Galacticus_Verbosity_Level_Set(verbosityLevelNew)
    !% Set the verbosity level.
    implicit none
    integer, intent(in) :: verbosityLevelNew

    verbosityLevel=verbosityLevelNew
    return
  end subroutine Galacticus_Verbosity_Level_Set

  subroutine Initialize_Display
    !% Initialize the module by determining the requested verbosity level.
    implicit none

    !$omp critical (Initialize_Display)
    if (.not.displayInitialized) then
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

  subroutine Galacticus_Display_Counter(percentageComplete,isNew,verbosity)
    !% Displays a percentage counter and bar to show progress.
    implicit none
    integer,          intent(in)           :: percentageComplete
    logical,          intent(in)           :: isNew
    integer,          intent(in), optional :: verbosity
    character(len=50)                      :: bar
    integer                                :: percentage,minorCount,majorCount
    logical                                :: showMessage
    
    !$omp critical(Galacticus_Message_Lock)
    call Initialize_Display
    if (present(verbosity)) then
       showMessage=(verbosity<=verbosityLevel)
    else
       showMessage=.true.
    end if
    if (showMessage) then
       if (.not.isNew) call Galacticus_Display_Counter_Clear_Lockless()
       percentage=max(0,min(percentageComplete,100))
       majorCount=percentage/2
       minorCount=percentage-majorCount*2
       bar=repeat("=",majorCount)//repeat("-",minorCount)//repeat(" ",50-majorCount-minorCount)
       write (0,'(1x,i3,"% [",a50,"]",$)') percentage,bar
    end if
    !$omp end critical(Galacticus_Message_Lock)
    return
  end subroutine Galacticus_Display_Counter

  subroutine Galacticus_Display_Counter_Clear(verbosity)
    !% Clears a percentage counter.
    implicit none
    integer, intent(in), optional :: verbosity
    logical                       :: showMessage
    
    !$omp critical(Galacticus_Message_Lock)
    call Galacticus_Display_Counter_Clear_Lockless(verbosity)
    !$omp end critical(Galacticus_Message_Lock)
    return
  end subroutine Galacticus_Display_Counter_Clear

  subroutine Galacticus_Display_Counter_Clear_Lockless(verbosity)
    !% Clears a percentage counter.
    implicit none
    integer, intent(in), optional :: verbosity
    logical                       :: showMessage
    
    call Initialize_Display
    if (present(verbosity)) then
       showMessage=(verbosity<=verbosityLevel)
    else
       showMessage=.true.
    end if
    if (showMessage) then
       write (0,'(a58,$)') repeat(char(8),58)
       write (0,'(a58,$)') repeat(" "    ,58)
       write (0,'(a58,$)') repeat(char(8),58)
    end if
    return
  end subroutine Galacticus_Display_Counter_Clear_Lockless

end module Galacticus_Display
