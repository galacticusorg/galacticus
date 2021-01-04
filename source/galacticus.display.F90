!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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

!% Contains a module which implements outputting of formatted, indented messages at various vebosity levels from \glc.

! Specify an explicit dependence on the C interface files.
!: $(BUILDPATH)/isatty.o

module Galacticus_Display
  !% Implements outputting of formatted, indented messages at various vebosity levels from \glc.
  use, intrinsic :: ISO_C_Binding, only : c_int
  implicit none
  private
  public :: Galacticus_Display_Message,Galacticus_Display_Indent,Galacticus_Display_Unindent,Galacticus_Verbosity_Level&
       &,Galacticus_Verbosity_Level_Set, Galacticus_Display_Counter, Galacticus_Display_Counter_Clear

  integer                                      :: maxThreads
  integer          , allocatable, dimension(:) :: indentationLevel
  character(len=10), allocatable, dimension(:) :: indentationFormat
  character(len=10), allocatable, dimension(:) :: indentationFormatNoNewLine

  character(len=20)                            :: threadFormat                              , masterFormat

  logical                                      :: displayInitialized        =.false.
  integer          , parameter  , public       :: verbosityDebug            =5              , verbosityInfo   =4, &
       &                                          verbosityWarn             =3              , verbosityWorking=2, &
       &                                          verbosityStandard         =1              , verbositySilent =0
  integer                                      :: verbosityLevel            =verbositySilent

  ! Progress bar state.
  logical                                      :: barVisible                =.false.
  integer                                      :: barPercentage             =0

  ! Output type.
  logical                                      :: stdOutIsFile
  
  interface Galacticus_Display_Message
     module procedure Galacticus_Display_Message_Char
     module procedure Galacticus_Display_Message_VarStr
  end interface Galacticus_Display_Message

  interface Galacticus_Display_Indent
     module procedure Galacticus_Display_Indent_Char
     module procedure Galacticus_Display_Indent_VarStr
  end interface Galacticus_Display_Indent

  interface Galacticus_Display_Unindent
     module procedure Galacticus_Display_Unindent_Char
     module procedure Galacticus_Display_Unindent_VarStr
  end interface Galacticus_Display_Unindent

  interface
     function stdOutIsATTY() bind(c,name='stdOutIsATTY')
       !% Template for a C function that determines if stdout is a TTY.
       import
       integer(c_int) :: stdOutIsATTY
     end function stdOutIsATTY
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
    integer, intent(in   ) :: verbosityLevelNew

    verbosityLevel=verbosityLevelNew
    return
  end subroutine Galacticus_Verbosity_Level_Set

  subroutine Initialize_Display
    !% Initialize the module by determining the requested verbosity level.
#ifdef USEMPI
    use    :: MPI    , only : MPI_Comm_Size      , MPI_Comm_Rank, MPI_Comm_World
#endif
    !$ use :: OMP_Lib, only : OMP_Get_Max_Threads
    implicit none
    integer           :: ompDigitsMaximum
#ifdef USEMPI
    integer           :: mpiDigitsMaximum , iError           , &
         &               mpiRank          , mpiCount
    character(len=32) :: masterHyperFormat, threadHyperFormat
#endif

    if (.not.displayInitialized) then
       !$omp critical (Initialize_Display)
       if (.not.displayInitialized) then
          ! For OpenMP runs, create an array of indentation levels, and formats.
          maxThreads=1
          !$ maxThreads=omp_get_max_threads()
          ! Do not use allocateArray() routines here as they may trigger recursive calls to this module which can lead to unbreakable
          ! parallel locks.
          allocate(indentationLevel          (maxThreads))
          allocate(indentationFormat         (maxThreads))
          allocate(indentationFormatNoNewLine(maxThreads))
          indentationLevel          =0
          indentationFormat         ='(a)'
          indentationFormatNoNewLine='(a,$)'
          ompDigitsMaximum             =int(log10(float(maxThreads)))+1
#ifdef USEMPI
          call MPI_Comm_Size(MPI_Comm_World,mpiCount,iError)
          if (iError /= 0) mpiCount=1
          call MPI_Comm_Rank(MPI_Comm_World,mpiRank ,iError)
          if (iError /= 0) mpiRank =0
          mpiDigitsMaximum=int(log10(float(mpiCount  )))+1
          !# <workaround type="gfortran" PR="92836" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=92836">
          !#  <description>Internal file I/O in gfortran can be non-thread safe.</description>
          !# </workaround>
#ifdef THREADSAFEIO
          !$omp critical(gfortranInternalIO)
#endif
          write (threadHyperFormat,'(a,i1,a,i1,a)'  ) '(a,i',mpiDigitsMaximum,'.',mpiDigitsMaximum,',a,i1,a)'
          write (masterHyperFormat,'(a,i1,a,i1,a)'  ) '(a,i',mpiDigitsMaximum,'.',mpiDigitsMaximum,',a,a ,a)'
          write (threadFormat     ,threadHyperFormat) '("',mpiRank,':",i'        ,ompDigitsMaximum ,' ,a2,$)'
          write (masterFormat     ,masterHyperFormat) '("',mpiRank,':',repeat("M",ompDigitsMaximum),'",a2,$)'
#ifdef THREADSAFEIO
          !$omp end critical(gfortranInternalIO)
#endif
#else
          !# <workaround type="gfortran" PR="92836" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=92836">
          !#  <description>Internal file I/O in gfortran can be non-thread safe.</description>
          !# </workaround>
#ifdef THREADSAFEIO
          !$omp critical(gfortranInternalIO)
#endif
          write (threadFormat     ,'(a,i1,a)'       ) '(i'           ,ompDigitsMaximum ,' ,a2,$)'
          write (masterFormat     ,'(a,a ,a)'       ) '("',repeat("M",ompDigitsMaximum),'",a2,$)'
#ifdef THREADSAFEIO
          !$omp end critical(gfortranInternalIO)
#endif
#endif
          ! Determine if stdout is connected to a file or pipe.
          stdOutIsFile=stdOutIsATTY() /= 1
          ! Mark display as initialized.
          displayInitialized=.true.
       end if
       !$omp end critical (Initialize_Display)
    end if
    return
  end subroutine Initialize_Display

  subroutine Galacticus_Display_Indent_VarStr(message,verbosity)
    !% Increase the indentation level and display a message.
    use :: ISO_Varying_String, only : varying_string, char
    implicit none
    type   (varying_string), intent(in   )           :: message
    integer                , intent(in   ), optional :: verbosity

    call Galacticus_Display_Indent_Char(char(message),verbosity)
    return
  end subroutine Galacticus_Display_Indent_VarStr

  subroutine Galacticus_Display_Indent_Char(message,verbosity)
    !% Increase the indentation level and display a message.
    !$ use :: OMP_Lib, only : OMP_In_Parallel, OMP_Get_Thread_Num
    implicit none
    character(len=*), intent(in   )           :: message
    integer         , intent(in   ), optional :: verbosity
    logical                                   :: showMessage
    integer                                   :: threadNumber

    !$omp critical(Galacticus_Message_Lock)
    call Initialize_Display
    if (present(verbosity)) then
       showMessage=(verbosity       <= verbosityLevel)
    else
       showMessage=(verbositySilent <  verbosityLevel)
    end if
    if (showMessage) then
       !# <workaround type="gfortran" PR="92836" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=92836">
       !#  <description>Internal file I/O in gfortran can be non-thread safe.</description>
       !# </workaround>
#ifdef THREADSAFEIO
       !$omp critical(gfortranInternalIO)
#endif
       !$ if (omp_in_parallel()) then
       !$    write (0,threadFormat) omp_get_thread_num(),": "
       !$ else
       !$    write (0,masterFormat)                      ": "
       !$ end if
       threadNumber=1
       !$ if (omp_in_parallel()) threadNumber=omp_get_thread_num()+1
       write (0,indentationFormatNoNewLine(threadNumber)) '-> '
       write (0,'(a)') trim(message)
#ifdef THREADSAFEIO
       !$omp end critical(gfortranInternalIO)
#endif
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

  subroutine Galacticus_Display_Unindent_VarStr(message,verbosity)
    !% Decrease the indentation level and display a message.
    use :: ISO_Varying_String, only : varying_string, char
    implicit none
    type   (varying_string), intent(in   )           :: message
    integer                , intent(in   ), optional :: verbosity

    call Galacticus_Display_Unindent_Char(char(message),verbosity)
    return
  end subroutine Galacticus_Display_Unindent_VarStr

  subroutine Galacticus_Display_Unindent_Char(message,verbosity)
    !% Decrease the indentation level and display a message.
    !$ use :: OMP_Lib, only : OMP_In_Parallel, OMP_Get_Thread_Num
    implicit none
    character(len=*), intent(in   )           :: message
    integer         , intent(in   ), optional :: verbosity
    integer                                   :: threadNumber
    logical                                   :: showMessage

    !$omp critical(Galacticus_Message_Lock)
    call Initialize_Display
    if (present(verbosity)) then
       showMessage=(verbosity       <= verbosityLevel)
    else
       showMessage=(verbositySilent <  verbosityLevel)
    end if
    if (showMessage) then
       !$ if (omp_in_parallel()) then
       !$    indentationLevel(omp_get_thread_num()+1)=max(indentationLevel(omp_get_thread_num()+1)-1,0)
       !$ else
       indentationLevel=max(indentationLevel-1,0)
       !$ end if
       call Create_Indentation_Format
       !# <workaround type="gfortran" PR="92836" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=92836">
       !#  <description>Internal file I/O in gfortran can be non-thread safe.</description>
       !# </workaround>
#ifdef THREADSAFEIO
       !$omp critical(gfortranInternalIO)
#endif
       !$ if (omp_in_parallel()) then
       !$    write (0,threadFormat) omp_get_thread_num(),": "
       !$ else
       !$    write (0,masterFormat)                      ": "
       !$ end if
       threadNumber=1
       !$ if (omp_in_parallel()) threadNumber=omp_get_thread_num()+1
       write (0,indentationFormatNoNewLine(threadNumber)) '<- '
       write (0,'(a)') trim(message)
#ifdef THREADSAFEIO
       !$omp end critical(gfortranInternalIO)
#endif
    end if
    !$omp end critical(Galacticus_Message_Lock)
    return
  end subroutine Galacticus_Display_Unindent_Char

  subroutine Galacticus_Display_Message_Char(message,verbosity)
    !% Display a message (input as a {\normalfont \ttfamily character} variable).
    !$ use :: OMP_Lib, only : OMP_In_Parallel, OMP_Get_Thread_Num
    implicit none
    character(len=*), intent(in   )           :: message
    integer         , intent(in   ), optional :: verbosity
    integer                                   :: threadNumber
    logical                                   :: showMessage

    !$omp critical(Galacticus_Message_Lock)
    call Initialize_Display
    if (present(verbosity)) then
       showMessage=(verbosity       <= verbosityLevel)
    else
       showMessage=(verbositySilent <  verbosityLevel)
    end if
    if (showMessage) then
       if (barVisible) call Galacticus_Display_Counter_Clear_Lockless()
       !# <workaround type="gfortran" PR="92836" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=92836">
       !#  <description>Internal file I/O in gfortran can be non-thread safe.</description>
       !# </workaround>
#ifdef THREADSAFEIO
       !$omp critical(gfortranInternalIO)
#endif
       !$ if (omp_in_parallel()) then
       !$    write (0,threadFormat) omp_get_thread_num(),": "
       !$ else
       !$    write (0,masterFormat)                      ": "
       !$ end if
       threadNumber=1
       !$ if (omp_in_parallel()) threadNumber=omp_get_thread_num()+1
       write (0,indentationFormat(threadNumber)) trim(message)
#ifdef THREADSAFEIO
       !$omp end critical(gfortranInternalIO)
#endif
       if (barVisible) call Galacticus_Display_Counter_Lockless(barPercentage,.true.)
    end if
    !$omp end critical(Galacticus_Message_Lock)
    return
  end subroutine Galacticus_Display_Message_Char

  subroutine Galacticus_Display_Message_VarStr(message,verbosity)
    !% Display a message (input as a {\normalfont \ttfamily varying\_string} variable).
    !$ use :: OMP_Lib           , only : OMP_In_Parallel, OMP_Get_Thread_Num
    use    :: ISO_Varying_String, only : varying_string , char
    implicit none
    type   (varying_string), intent(in   )           :: message
    integer                , intent(in   ), optional :: verbosity
    integer                                          :: threadNumber
    logical                                          :: showMessage

    !$omp critical(Galacticus_Message_Lock)
    call Initialize_Display
    if (present(verbosity)) then
       showMessage=(verbosity       <= verbosityLevel)
    else
       showMessage=(verbositySilent <  verbosityLevel)
    end if
    if (showMessage) then
       if (barVisible) call Galacticus_Display_Counter_Clear_Lockless()
       !# <workaround type="gfortran" PR="92836" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=92836">
       !#  <description>Internal file I/O in gfortran can be non-thread safe.</description>
       !# </workaround>
#ifdef THREADSAFEIO
       !$omp critical(gfortranInternalIO)
#endif
       !$ if (omp_in_parallel()) then
       !$    write (0,threadFormat) omp_get_thread_num(),": "
       !$ else
       !$    write (0,masterFormat)                      ": "
       !$ end if
       threadNumber=1
       !$ if (omp_in_parallel()) threadNumber=omp_get_thread_num()+1
       write (0,indentationFormat(threadNumber)) char(message)
#ifdef THREADSAFEIO
       !$omp end critical(gfortranInternalIO)
#endif
       if (barVisible) call Galacticus_Display_Counter_Lockless(barPercentage,.true.)
    end if
    !$omp end critical(Galacticus_Message_Lock)
    return
  end subroutine Galacticus_Display_Message_VarStr

  subroutine Create_Indentation_Format
    !% Create a format for indentation.
    !$ use :: OMP_Lib, only : OMP_In_Parallel, OMP_Get_Thread_Num
    implicit none
    integer, parameter :: indentSpaces=4
    integer            :: threadNumber

    threadNumber=1
    !$ if (omp_in_parallel()) threadNumber=omp_get_thread_num()+1
    !# <workaround type="gfortran" PR="92836" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=92836">
    !#  <description>Internal file I/O in gfortran can be non-thread safe.</description>
    !# </workaround>
#ifdef THREADSAFEIO
    !$omp critical(gfortranInternalIO)
#endif
    select case (indentationLevel(threadNumber)*indentSpaces)
    case (0)
       write (indentationFormat         (threadNumber),'(a)'     ) '(a)'
       write (indentationFormatNoNewLine(threadNumber),'(a)'     ) '(a,$)'
    case (1:9)
       write (indentationFormat         (threadNumber),'(a,i1,a)') '(',indentationLevel(threadNumber)*indentSpaces,'x,a)'
       write (indentationFormatNoNewLine(threadNumber),'(a,i1,a)') '(',indentationLevel(threadNumber)*indentSpaces,'x,a,$)'
    case (10:99)
       write (indentationFormat         (threadNumber),'(a,i2,a)') '(',indentationLevel(threadNumber)*indentSpaces,'x,a)'
       write (indentationFormatNoNewLine(threadNumber),'(a,i2,a)') '(',indentationLevel(threadNumber)*indentSpaces,'x,a,$)'
    end select
#ifdef THREADSAFEIO
    !$omp end critical(gfortranInternalIO)
#endif
    !$ if (.not.omp_in_parallel()) then
    !$    indentationLevel          =indentationLevel          (1)
    !$    indentationFormat         =indentationFormat         (1)
    !$    indentationFormatNoNewLine=indentationFormatNoNewLine(1)
    !$ end if
    return
  end subroutine Create_Indentation_Format

  subroutine Galacticus_Display_Counter(percentageComplete,isNew,verbosity)
    !% Displays a percentage counter and bar to show progress.
    implicit none
    integer, intent(in   )           :: percentageComplete
    logical, intent(in   )           :: isNew
    integer, intent(in   ), optional :: verbosity

    !$omp critical(Galacticus_Message_Lock)
    call Galacticus_Display_Counter_Lockless(percentageComplete,isNew,verbosity)
    !$omp end critical(Galacticus_Message_Lock)
    return
  end subroutine Galacticus_Display_Counter

  subroutine Galacticus_Display_Counter_Lockless(percentageComplete,isNew,verbosity)
    !% Displays a percentage counter and bar to show progress.
    implicit none
    integer          , intent(in   )           :: percentageComplete
    logical          , intent(in   )           :: isNew
    integer          , intent(in   ), optional :: verbosity
    character(len=50)                          :: bar
    integer                                    :: majorCount        , minorCount, percentage
    logical                                    :: showMessage

    call Initialize_Display
    if (present(verbosity)) then
       showMessage=(verbosity       <= verbosityLevel)
    else
       showMessage=(verbositySilent <  verbosityLevel)
    end if
    if (showMessage) then
       if (percentageComplete == barPercentage .and. .not.isNew) return
       if (.not.isNew) call Galacticus_Display_Counter_Clear_Lockless()
       percentage=max(0,min(percentageComplete,100))
       majorCount=percentage/2
       minorCount=percentage-majorCount*2
       bar=repeat("=",majorCount)//repeat("-",minorCount)//repeat(" ",50-majorCount-minorCount)
       !# <workaround type="gfortran" PR="92836" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=92836">
       !#  <description>Internal file I/O in gfortran can be non-thread safe.</description>
       !# </workaround>
#ifdef THREADSAFEIO
       !$omp critical(gfortranInternalIO)
#endif
       write (0,'(1x,i3,"% [",a50,"]",$)') percentage,bar
       ! For output to a file, add a newline, since we will not be deleting the bar.
       if (stdOutIsFile) write (0,*)
#ifdef THREADSAFEIO
       !$omp end critical(gfortranInternalIO)
#endif
       barVisible   =.true.
       barPercentage=percentageComplete
    end if
    return
  end subroutine Galacticus_Display_Counter_Lockless

  subroutine Galacticus_Display_Counter_Clear(verbosity)
    !% Clears a percentage counter.
    implicit none
    integer, intent(in   ), optional :: verbosity

    !$omp critical(Galacticus_Message_Lock)
    call Galacticus_Display_Counter_Clear_Lockless(verbosity)
    barVisible   =.false.
    barPercentage=0
    !$omp end critical(Galacticus_Message_Lock)
    return
  end subroutine Galacticus_Display_Counter_Clear

  subroutine Galacticus_Display_Counter_Clear_Lockless(verbosity)
    !% Clears a percentage counter.
    implicit none
    integer, intent(in   ), optional :: verbosity
    logical                          :: showMessage

    call Initialize_Display()
    ! If output is to a file we do not attempt to clear the bar (which is useful only on a TTY).
    if (stdOutIsFile) return
    if (present(verbosity)) then
       showMessage=(verbosity       <= verbosityLevel)
    else
       showMessage=(verbositySilent <  verbosityLevel)
    end if
    if (showMessage) then
       !# <workaround type="gfortran" PR="92836" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=92836">
       !#  <description>Internal file I/O in gfortran can be non-thread safe.</description>
       !# </workaround>
#ifdef THREADSAFEIO
       !$omp critical(gfortranInternalIO)
#endif
       write (0,'(a58,$)') repeat(char(8),58)
       write (0,'(a58,$)') repeat(" "    ,58)
       write (0,'(a58,$)') repeat(char(8),58)
#ifdef THREADSAFEIO
       !$omp end critical(gfortranInternalIO)
#endif
    end if
    return
  end subroutine Galacticus_Display_Counter_Clear_Lockless

end module Galacticus_Display
