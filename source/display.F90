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

!!{
Contains a module which implements outputting of formatted, indented messages at various verbosity levels from \glc.
!!}

module Display
  !!{
  Implements outputting of formatted, indented messages at various verbosity levels from \glc.
  !!}
  use, intrinsic :: ISO_C_Binding, only : c_int
  implicit none
  private
  public :: displayMessage     , displayIndent , displayUnindent    , displayVerbosity, &
       &    displayVerbositySet, displayCounter, displayCounterClear, displayRed      , &
       &    displayMagenta     , displayGreen  , displayBlue        , displayYellow   , &
       &    displayBold        , displayReset

  !![
  <enumeration>
   <name>verbosityLevel</name>
   <description>Verbosity levels for message display.</description>
   <indexing>0</indexing>
   <visibility>public</visibility>
   <encodeFunction>yes</encodeFunction>
   <decodeFunction>yes</decodeFunction>
   <errorValue>-1</errorValue>
   <validator>yes</validator>
   <entry label="silent"  />
   <entry label="standard"/>
   <entry label="working" />
   <entry label="warn"    />
   <entry label="info"    />
   <entry label="debug"   />
  </enumeration>
  !!]

  integer                                                             :: maxThreads
  integer                                 , allocatable, dimension(:) :: indentationLevel
  character(len=10                       ), allocatable, dimension(:) :: indentationFormat
  character(len=10                       ), allocatable, dimension(:) :: indentationFormatNoNewLine

  character(len=20                       )                            :: threadFormat                                     , masterFormat

  logical                                                             :: displayInitialized        =.false.               , verbositySet=.false.
  type     (enumerationVerbosityLevelType)                            :: verbosityLevel            =verbosityLevelStandard

  ! Progress bar state.
  logical                                                             :: barVisible                =.false.
  integer                                                             :: barPercentage             =0

  ! Output type.
  logical                                                             :: stdOutIsFile
  
  ! ANSI codes.
  character(len=*                        ), parameter                 :: ESC                       =achar(27)
   
  interface displayMessage
     module procedure displayMessageChar
     module procedure displayMessageVarStr
  end interface displayMessage

  interface displayIndent
     module procedure displayIndentChar
     module procedure displayIndentVarStr
  end interface displayIndent

  interface displayUnindent
     module procedure displayUnindentChar
     module procedure displayUnindentVarStr
  end interface displayUnindent

contains

  function displayVerbosity()
    !!{
    Returns the verbosity level in \glc.
    !!}
    implicit none
    type(enumerationVerbosityLevelType) :: displayVerbosity

    displayVerbosity=verbosityLevel
    return
  end function displayVerbosity

  subroutine displayVerbositySet(verbosityLevelNew)
    !!{
    Set the verbosity level.
    !!}
    implicit none
    type(enumerationVerbosityLevelType), intent(in   ) :: verbosityLevelNew

    if (enumerationVerbosityLevelIsValid(verbosityLevelNew)) then
       ! Requested value is valid, set it.
       verbosityLevel=verbosityLevelNew
       verbositySet  =.true.
    else if (.not.verbositySet) then
       ! Verbosity has not been set, and the requested value is invalid. Use standard verbosity.
       verbosityLevel=verbosityLevelStandard
    end if
    return
  end subroutine displayVerbositySet

  subroutine initialize()
    !!{
    Initialize the module by determining the requested verbosity level.
    !!}
#ifdef USEMPI
    use    :: MPI_F08      , only : MPI_Comm_Size      , MPI_Comm_Rank, MPI_Comm_World
#endif
    !$ use :: OMP_Lib      , only : OMP_Get_Max_Threads
    use    :: System_Output, only : stdOutIsATTY
    implicit none
    integer           :: ompDigitsMaximum
#ifdef USEMPI
    integer           :: mpiDigitsMaximum , iError           , &
         &               mpiRank          , mpiCount
    character(len=32) :: masterHyperFormat, threadHyperFormat
#endif

    if (.not.displayInitialized) then
       !$omp critical (displayInitialize)
       if (.not.displayInitialized) then
          ! For OpenMP runs, create an array of indentation levels, and formats.
          maxThreads=1
          !$ maxThreads=omp_get_max_threads()
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
          write (threadHyperFormat,'(a,i1,a,i1,a)'  ) '(a,i',mpiDigitsMaximum,'.',mpiDigitsMaximum                      ,',a,i1,a1,i1,a)'
          write (masterHyperFormat,'(a,i1,a,i1,a)'  ) '(a,i',mpiDigitsMaximum,'.',mpiDigitsMaximum                      ,',a,a       ,a)'
          write (threadFormat     ,threadHyperFormat) '("',mpiRank,':",i'        ,ompDigitsMaximum ,'.',ompDigitsMaximum,' ,a2,$)'
          write (masterFormat     ,masterHyperFormat) '("',mpiRank,':',repeat("M",ompDigitsMaximum)                     ,'",a2,$)'
#else
          write (threadFormat     ,'(a,i1,a1,i1,a)'       ) '(i'           ,ompDigitsMaximum ,'.',ompDigitsMaximum,' ,a2,$)'
          write (masterFormat     ,'(a,a       ,a)'       ) '("',repeat("M",ompDigitsMaximum)                     ,'",a2,$)'
#endif
          ! Determine if stdout is connected to a file or pipe.
          stdOutIsFile=.not.stdOutIsATTY()
          ! Mark display as initialized.
          displayInitialized=.true.
       end if
       !$omp end critical (displayInitialize)
    end if
    return
  end subroutine initialize

  subroutine displayIndentVarStr(message,verbosity)
    !!{
    Increase the indentation level and display a message.
    !!}
    use :: ISO_Varying_String, only : varying_string, char
    implicit none
    type(varying_string               ), intent(in   )           :: message
    type(enumerationVerbosityLevelType), intent(in   ), optional :: verbosity

    call displayIndentChar(char(message),verbosity)
    return
  end subroutine displayIndentVarStr

  subroutine displayIndentChar(message,verbosity)
    !!{
    Increase the indentation level and display a message.
    !!}
    use   , intrinsic :: ISO_Fortran_Env, only : output_unit
    !$ use            :: OMP_Lib        , only : OMP_In_Parallel, OMP_Get_Thread_Num
    implicit none
    character(len=*                        ), intent(in   )           :: message
    type     (enumerationVerbosityLevelType), intent(in   ), optional :: verbosity
    integer                                                           :: threadNumber

    !$omp critical(Display_Lock)
    call initialize()
    if (showMessage(verbosity)) then
       !$ if (omp_in_parallel()) then
       !$    write (output_unit,threadFormat) omp_get_thread_num(),": "
       !$ else
       !$    write (output_unit,masterFormat)                      ": "
       !$ end if
       threadNumber=1
       !$ if (omp_in_parallel()) threadNumber=omp_get_thread_num()+1
       write (output_unit,indentationFormatNoNewLine(threadNumber)) '-> '
       write (output_unit,'(a)') trim(message)
       !$ if (omp_in_parallel()) then
       !$    indentationLevel(omp_get_thread_num()+1)=indentationLevel(omp_get_thread_num()+1)+1
       !$ else
             indentationLevel=indentationLevel+1
       !$ end if
       call formatIndentationCreate()
    end if
    !$omp end critical(Display_Lock)
    return
  end subroutine displayIndentChar

  subroutine displayUnindentVarStr(message,verbosity)
    !!{
    Decrease the indentation level and display a message.
    !!}
    use :: ISO_Varying_String, only : varying_string, char
    implicit none
    type(varying_string               ), intent(in   )           :: message
    type(enumerationVerbosityLevelType), intent(in   ), optional :: verbosity

    call displayUnindentChar(char(message),verbosity)
    return
  end subroutine displayUnindentVarStr

  subroutine displayUnindentChar(message,verbosity)
    !!{
    Decrease the indentation level and display a message.
    !!}
    use   , intrinsic :: ISO_Fortran_Env, only : output_unit
    !$ use            :: OMP_Lib        , only : OMP_In_Parallel, OMP_Get_Thread_Num
    implicit none
    character(len=*                        ), intent(in   )           :: message
    type     (enumerationVerbosityLevelType), intent(in   ), optional :: verbosity
    integer                                                           :: threadNumber

    !$omp critical(Display_Lock)
    call initialize()
    if (showMessage(verbosity)) then
       !$ if (omp_in_parallel()) then
       !$    indentationLevel(omp_get_thread_num()+1)=max(indentationLevel(omp_get_thread_num()+1)-1,0)
       !$ else
       indentationLevel=max(indentationLevel-1,0)
       !$ end if
       call formatIndentationCreate()
       !$ if (omp_in_parallel()) then
       !$    write (output_unit,threadFormat) omp_get_thread_num(),": "
       !$ else
       !$    write (output_unit,masterFormat)                      ": "
       !$ end if
       threadNumber=1
       !$ if (omp_in_parallel()) threadNumber=omp_get_thread_num()+1
       write (output_unit,indentationFormatNoNewLine(threadNumber)) '<- '
       write (output_unit,'(a)') trim(message)
    end if
    !$omp end critical(Display_Lock)
    return
  end subroutine displayUnindentChar

  subroutine displayMessageChar(message,verbosity)
    !!{
    Display a message (input as a {\normalfont \ttfamily character} variable).
    !!}
    use   , intrinsic :: ISO_Fortran_Env, only : output_unit
    !$ use            :: OMP_Lib        , only : OMP_In_Parallel, OMP_Get_Thread_Num
    implicit none
    character(len=*                        ), intent(in   )           :: message
    type     (enumerationVerbosityLevelType), intent(in   ), optional :: verbosity
    integer                                                           :: threadNumber

    !$omp critical(Display_Lock)
    call initialize()
    if (showMessage(verbosity)) then
       if (barVisible) call counterClearLockless()
       !$ if (omp_in_parallel()) then
       !$    write (output_unit,threadFormat) omp_get_thread_num(),": "
       !$ else
       !$    write (output_unit,masterFormat)                      ": "
       !$ end if
       threadNumber=1
       !$ if (omp_in_parallel()) threadNumber=omp_get_thread_num()+1
       write (output_unit,indentationFormat(threadNumber)) trim(message)
       if (barVisible) call displayCounterLockless(barPercentage,.true.)
    end if
    !$omp end critical(Display_Lock)
    return
  end subroutine displayMessageChar

  subroutine displayMessageVarStr(message,verbosity)
    !!{
    Display a message (input as a {\normalfont \ttfamily varying\_string} variable).
    !!}
    use   , intrinsic :: ISO_Fortran_Env   , only : output_unit
    !$ use            :: OMP_Lib           , only : OMP_In_Parallel, OMP_Get_Thread_Num
    use               :: ISO_Varying_String, only : varying_string , char
    implicit none
    type   (varying_string               ), intent(in   )           :: message
    type   (enumerationVerbosityLevelType), intent(in   ), optional :: verbosity
    integer                                                         :: threadNumber

    !$omp critical(Display_Lock)
    call initialize()
    if (showMessage(verbosity)) then
       if (barVisible) call counterClearLockless()
       !$ if (omp_in_parallel()) then
       !$    write (output_unit,threadFormat) omp_get_thread_num(),": "
       !$ else
       !$    write (output_unit,masterFormat)                      ": "
       !$ end if
       threadNumber=1
       !$ if (omp_in_parallel()) threadNumber=omp_get_thread_num()+1
       write (output_unit,indentationFormat(threadNumber)) char(message)
       if (barVisible) call displayCounterLockless(barPercentage,.true.)
    end if
    !$omp end critical(Display_Lock)
    return
  end subroutine displayMessageVarStr

  subroutine formatIndentationCreate()
    !!{
    Create a format for indentation.
    !!}
    !$ use :: OMP_Lib, only : OMP_In_Parallel, OMP_Get_Thread_Num
    implicit none
    integer, parameter :: indentSpaces=4
    integer            :: threadNumber

    threadNumber=1
    !$ if (omp_in_parallel()) threadNumber=omp_get_thread_num()+1
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
    !$ if (.not.omp_in_parallel()) then
    !$    indentationLevel          =indentationLevel          (1)
    !$    indentationFormat         =indentationFormat         (1)
    !$    indentationFormatNoNewLine=indentationFormatNoNewLine(1)
    !$ end if
    return
  end subroutine formatIndentationCreate

  subroutine displayCounter(percentageComplete,isNew,verbosity)
    !!{
    Displays a percentage counter and bar to show progress.
    !!}
    implicit none
    integer                               , intent(in   )           :: percentageComplete
    logical                               , intent(in   )           :: isNew
    type   (enumerationVerbosityLevelType), intent(in   ), optional :: verbosity

    !$omp critical(Display_Lock)
    call displayCounterLockless(percentageComplete,isNew,verbosity)
    !$omp end critical(Display_Lock)
    return
  end subroutine displayCounter

  subroutine displayCounterLockless(percentageComplete,isNew,verbosity)
    !!{
    Displays a percentage counter and bar to show progress.
    !!}
    use, intrinsic :: ISO_Fortran_Env, only : output_unit
    implicit none
    integer                                 , intent(in   )           :: percentageComplete
    logical                                 , intent(in   )           :: isNew
    type     (enumerationVerbosityLevelType), intent(in   ), optional :: verbosity
    character(len=50                       )                          :: bar
    integer                                                           :: majorCount        , minorCount, percentage

    call initialize()
    if (showMessage(verbosity)) then
       if (percentageComplete == barPercentage .and. .not.isNew) return
       if (.not.isNew) call counterClearLockless()
       percentage=max(0,min(percentageComplete,100))
       majorCount=percentage/2
       minorCount=percentage-majorCount*2
       bar=repeat("=",majorCount)//repeat("-",minorCount)//repeat(" ",50-majorCount-minorCount)
       write (output_unit,'(1x,i3,"% [",a50,"]",$)') percentage,bar
       ! For output to a file, add a newline, since we will not be deleting the bar.
       if (stdOutIsFile) write (output_unit,*)
       barVisible   =.true.
       barPercentage=percentageComplete
    end if
    return
  end subroutine displayCounterLockless

  subroutine displayCounterClear(verbosity)
    !!{
    Clears a percentage counter.
    !!}
    implicit none
    type(enumerationVerbosityLevelType), intent(in   ), optional :: verbosity

    !$omp critical(Display_Lock)
    call counterClearLockless(verbosity)
    barVisible   =.false.
    barPercentage=0
    !$omp end critical(Display_Lock)
    return
  end subroutine displayCounterClear

  subroutine counterClearLockless(verbosity)
    !!{
    Clears a percentage counter.
    !!}
    use, intrinsic :: ISO_Fortran_Env, only : output_unit
    implicit none
    type(enumerationVerbosityLevelType), intent(in   ), optional :: verbosity

    call initialize()
    ! If output is to a file we do not attempt to clear the bar (which is useful only on a TTY).
    if (stdOutIsFile) return
    if (showMessage(verbosity)) then
       write (output_unit,'(a58,$)') repeat(char(8),58)
       write (output_unit,'(a58,$)') repeat(" "    ,58)
       write (output_unit,'(a58,$)') repeat(char(8),58)
    end if
    return
  end subroutine counterClearLockless

  logical function showMessage(verbosity)
    !!{
    Return true if the message should be displayed at the current verbosity level.
    !!}
    implicit none
    type(enumerationVerbosityLevelType), intent(in   ), optional :: verbosity

    if (present(verbosity)) then
       showMessage=(verbosity            <= verbosityLevel)
    else
       showMessage=(verbosityLevelSilent <  verbosityLevel)
    end if
    return
  end function showMessage

  function displayRed()
    !!{
    Return the ANSI escape code for red text.
    !!}
#ifndef USEMPI
    use :: System_Output, only : stdOutIsATTY
#endif
    implicit none
    character(len=:), allocatable :: displayRed

#ifdef USEMPI
    displayRed=""
#else
    if (stdOutIsATTY()) then
       displayRed=ESC//"[31m"
    else
       displayRed=""
    end if
#endif
    return
  end function displayRed

  function displayBlue()
    !!{
    Return the ANSI escape code for blue text.
    !!}
#ifndef USEMPI
    use :: System_Output, only : stdOutIsATTY
#endif
    implicit none
    character(len=:), allocatable :: displayBlue

#ifdef USEMPI
    displayBlue=""
#else
    if (stdOutIsATTY()) then
       displayBlue=ESC//"[34m"
    else
       displayBlue=""
    end if
#endif
    return
  end function displayBlue
  
  function displayYellow()
    !!{
    Return the ANSI escape code for yellow text.
    !!}
#ifndef USEMPI
    use :: System_Output, only : stdOutIsATTY
#endif
    implicit none
    character(len=:), allocatable :: displayYellow

#ifdef USEMPI
    displayYellow=""
#else
    if (stdOutIsATTY()) then
       displayYellow=ESC//"[33m"
    else
       displayYellow=""
    end if
#endif
    return
  end function displayYellow
  
  function displayMagenta()
    !!{
    Return the ANSI escape code for magenta text.
    !!}
#ifndef USEMPI
    use :: System_Output, only : stdOutIsATTY
#endif
    implicit none
    character(len=:), allocatable :: displayMagenta

#ifdef USEMPI
    displayMagenta=""
#else
    if (stdOutIsATTY()) then
       displayMagenta=ESC//"[35m"
    else
       displayMagenta=""
    end if
#endif
    return
  end function displayMagenta
  
  function displayGreen()
    !!{
    Return the ANSI escape code for green text.
    !!}
#ifndef USEMPI
    use :: System_Output, only : stdOutIsATTY
#endif
    implicit none
    character(len=:), allocatable :: displayGreen

#ifdef USEMPI
    displayGreen=""
#else
    if (stdOutIsATTY()) then
       displayGreen=ESC//"[32m"
    else
       displayGreen=""
    end if
#endif
    return
  end function displayGreen
  
  function displayBold()
    !!{
    Return the ANSI escape code for bold text.
    !!}
#ifndef USEMPI
    use :: System_Output, only : stdOutIsATTY
#endif
    implicit none
    character(len=:), allocatable :: displayBold

#ifdef USEMPI
    displayBold=""
#else
    if (stdOutIsATTY()) then
       displayBold=ESC//"[1m"
    else
       displayBold=""
    end if
#endif
    return
  end function displayBold
  
  function displayReset()
    !!{
    Return the ANSI escape code to reset text.
    !!}
#ifndef USEMPI
    use :: System_Output, only : stdOutIsATTY
#endif
    implicit none
    character(len=:), allocatable :: displayReset

#ifdef USEMPI
    displayReset=""
#else
    if (stdOutIsATTY()) then
       displayReset=ESC//"[0m"
    else
       displayReset=""
    end if
#endif
    return
  end function displayReset
  
end module Display
