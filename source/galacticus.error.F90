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

!% Contains a module which implements error reporting for the {\normalfont \scshape Galacticus} package.

module Galacticus_Error
  !% Implements error reporting for the {\normalfont \scshape Galacticus} package.
  use, intrinsic :: ISO_C_Binding     , only : c_int
  use            :: ISO_Varying_String, only : varying_string
  use            :: Interface_GSL     , only : GSL_Failure   , GSL_Success , GSL_eDom, GSL_eRange, &
          &                                    GSL_eUndrFlw  , GSL_eZeroDiv
  implicit none
  private
  public :: Galacticus_Error_Report               , Galacticus_Error_Handler_Register    , &
       &    Galacticus_Component_List             , Galacticus_GSL_Error_Handler_Abort_On, &
       &    Galacticus_GSL_Error_Handler_Abort_Off, Galacticus_GSL_Error_Status          , &
       &    Galacticus_Warn                       , Galacticus_Error_Wait_Set

  interface Galacticus_Error_Report
     module procedure Galacticus_Error_Report_Char
     module procedure Galacticus_Error_Report_VarStr
  end interface Galacticus_Error_Report

  interface Galacticus_Warn
     module procedure Galacticus_Warn_Char
     module procedure Galacticus_Warn_VarStr
  end interface Galacticus_Warn

  ! Specify an explicit dependence on the hdf5_cFuncs.o object file.
  !: $(BUILDPATH)/hdf5_cFuncs.o
  interface
     subroutine H5Close_C() bind(c,name='H5Close_C')
     end subroutine H5Close_C
  end interface

  ! Public error codes. Where relevant these copy GSL error codes, otherwise values above 1024
  ! are used so as not to conflict with GSL error codes.
  integer, parameter, public :: errorStatusSuccess     =GSL_Success  ! Success.
  integer, parameter, public :: errorStatusFail        =GSL_Failure  ! Generic failure.
  integer, parameter, public :: errorStatusInputDomain =GSL_eDom     ! Input domain error.
  integer, parameter, public :: errorStatusOutOfRange  =GSL_eRange   ! Output range error.
  integer, parameter, public :: errorStatusDivideByZero=GSL_eZeroDiv ! Divide by zero.
  integer, parameter, public :: errorStatusUnderflow   =GSL_eUndrFlw ! Floating point underflow.
  integer, parameter, public :: errorStatusXCPU        =1025         ! CPU time limit exceeded.

  ! Time to wait after errors under MPI.
  integer                    :: errorWaitTime          =86400

  ! GSL error status.
  logical             :: abortOnErrorGSL=.true.
  integer(kind=c_int) :: errorStatusGSL
  !$omp threadprivate(abortOnErrorGSL,errorStatusGSL)

  ! Type used to accumulate warning messages.
  type :: warning
     type(varying_string)          :: message
     type(warning       ), pointer :: next
  end type warning

  ! Record of warnings.
  type   (warning), pointer :: warningList
  logical                   :: warningsFound=.false.

contains

  subroutine Galacticus_Error_Report_VarStr(message)
    !% Display an error message.
    use :: ISO_Varying_String, only : char
    implicit none
    type(varying_string), intent(in   ) :: message

    call Galacticus_Error_Report_Char(char(message))
    return
  end subroutine Galacticus_Error_Report_VarStr

  subroutine Galacticus_Error_Report_Char(message)
    !% Display an error message.
#ifndef UNCLEANEXIT
    use    :: HDF5   , only : H5Close_F
#endif
    !$ use :: OMP_Lib, only : OMP_Get_Thread_Num, OMP_In_Parallel
    implicit none
    character(len=*), intent(in   ) :: message
    integer                         :: error

    write (0,'(a)') 'Fatal error:'
    write (0,'(a)') trim(message)
    !$ if (omp_in_parallel()) then
    !$    write (0,*) " => Error occurred in thread ",omp_get_thread_num()
    !$ else
    !$    write (0,*) " => Error occurred in master thread"
    !$ end if
    call BackTrace             (     )
    call Galacticus_Warn_Review(     )
    call Flush                 (    0)
#ifdef UNCLEANEXIT
    call Exit(1)
#else
    call H5Close_F             (error)
    call H5Close_C             (     )
    call Abort                 (     )
#endif
    return
  end subroutine Galacticus_Error_Report_Char

  subroutine Galacticus_Warn_VarStr(message)
    !% Display a warning message
    use :: ISO_Varying_String, only : char
    implicit none
    type(varying_string), intent(in   ) :: message

    call Galacticus_Warn_Char(char(message))
    return
  end subroutine Galacticus_Warn_VarStr

  subroutine Galacticus_Warn_Char(message)
    !% Display a warning message.
    use :: Display           , only : displayMessage, displayVerbosity, verbosityLevelWarn
    use :: ISO_Varying_String, only : assignment(=)
    implicit none
    character(len=*  ), intent(in   ) :: message
    type     (warning), pointer       :: newWarning

    ! Display the message.
    call displayMessage(message,verbosity=verbosityLevelWarn)
    ! Add this warning message to the list of warnings in case we need to display them on an
    ! error condition.
    !$omp critical (Galacticus_Warn)
    if (displayVerbosity() < verbosityLevelWarn) then
       if (.not.warningsFound) then
          allocate(warningList)
          newWarning => warningList
       else
          newWarning => warningList
          do while (associated(newWarning%next))
             newWarning => newWarning%next
          end do
          allocate(newWarning%next)
          newWarning => newWarning%next
       end if
       newWarning%next    => null   ()
       newWarning%message =  message
    end if
    warningsFound=.true.
    !$omp end critical (Galacticus_Warn)
    return
  end subroutine Galacticus_Warn_Char

  subroutine Galacticus_Warn_Review()
    !% Review any warning messages emitted during the run.
    use :: ISO_Varying_String, only : char
    implicit none
    type(warning), pointer :: warning_

    !$omp critical (Galacticus_Warn)
    if (warningsFound) then
       write (0,*) " => The following warnings were issued:"
       warning_ => warningList
       do while (associated(warning_))
          write (0,*) char(warning_%message)
          warning_ => warning_%next
       end do
    end if
    !$omp end critical (Galacticus_Warn)
    return
  end subroutine Galacticus_Warn_Review

  subroutine Galacticus_Error_Handler_Register()
    !% Register signal handlers.
    use, intrinsic :: ISO_C_Binding, only : c_funptr
    use            :: Interface_GSL, only : gslSetErrorHandler
    implicit none
    type(c_funptr) :: standardGslErrorHandler

    call Signal( 2,Galacticus_Signal_Handler_SIGINT )
    call Signal( 4,Galacticus_Signal_Handler_SIGILL )
    call Signal( 7,Galacticus_Signal_Handler_SIGBUS )
    call Signal( 8,Galacticus_Signal_Handler_SIGFPE )
    call Signal(11,Galacticus_Signal_Handler_SIGSEGV)
    call Signal(15,Galacticus_Signal_Handler_SIGINT )
    call Signal(24,Galacticus_Signal_Handler_SIGXCPU)
    !$omp critical(gslErrorHandler)
    standardGslErrorHandler=gslSetErrorHandler(Galacticus_GSL_Error_Handler)
    !$omp end critical(gslErrorHandler)
    return
  end subroutine Galacticus_Error_Handler_Register

  subroutine Galacticus_Signal_Handler_SIGINT()
    !% Handle {\normalfont \ttfamily SIGINT} signals, by flushing all data and then aborting.
#ifndef UNCLEANEXIT
    use    :: HDF5   , only : H5Close_F
#endif
#ifdef USEMPI
    use    :: MPI    , only : MPI_COmm_Rank     , MPI_Comm_World
#endif
    !$ use :: OMP_Lib, only : OMP_Get_Thread_Num, OMP_In_Parallel
    implicit none
    integer            :: error
#ifdef USEMPI
    integer            :: mpiRank
    character(len=128) :: hostName
    logical            :: flag
#endif

    write (0,*) 'Galacticus was interrupted - will try to flush data before exiting.'
    !$ if (omp_in_parallel()) then
    !$    write (0,*) " => Error occurred in thread ",omp_get_thread_num()
    !$ else
    !$    write (0,*) " => Error occurred in master thread"
    !$ end if
    call BackTrace                (     )
    call Galacticus_Warn_Review   (     )
    call Flush                    (    0)
#ifdef UNCLEANEXIT
    call Exit(1)
#else
#ifdef USEMPI
    call MPI_Initialized(flag,error)
    if (flag) then
       call MPI_Comm_Rank(MPI_Comm_World,mpiRank,error)
       call hostnm(hostName)
       write (0,*) " => Error occurred in MPI process ",mpiRank,"; PID ",getPID(),"; host ",trim(hostName)
       write (0,'(a,i8,a)') " => Sleeping for ",errorWaitTime,"s to allow for attachment of debugger"
       call Flush(0)
       call Sleep(errorWaitTime)
    end if
#endif
    call H5Close_F(error)
    call H5Close_C(     )
    call Abort    (     )
#endif
    return
  end subroutine Galacticus_Signal_Handler_SIGINT

  subroutine Galacticus_Signal_Handler_SIGSEGV()
    !% Handle {\normalfont \ttfamily SIGSEGV} signals, by flushing all data and then aborting.
#ifndef UNCLEANEXIT
    use    :: HDF5   , only : H5Close_F
#endif
#ifdef USEMPI
    use    :: MPI    , only : MPI_Comm_Rank     , MPI_Comm_World
#endif
    !$ use :: OMP_Lib, only : OMP_Get_Thread_Num, OMP_In_Parallel
    implicit none
    integer            :: error
#ifdef USEMPI
    integer            :: mpiRank
    character(len=128) :: hostName
    logical            :: flag
#endif

    write (0,*) 'Galacticus experienced a segfault - will try to flush data before exiting.'
    !$ if (omp_in_parallel()) then
    !$    write (0,*) " => Error occurred in thread ",omp_get_thread_num()
    !$ else
    !$    write (0,*) " => Error occurred in master thread"
    !$ end if
    call BackTrace                (     )
    call Galacticus_Warn_Review   (     )
    call Flush                    (    0)
#ifdef UNCLEANEXIT
    call Exit(1)
#else
#ifdef USEMPI
    call MPI_Initialized(flag,error)
    if (flag) then
       call MPI_Comm_Rank(MPI_Comm_World,mpiRank,error)
       call hostnm(hostName)
       write (0,*) " => Error occurred in MPI process ",mpiRank,"; PID ",getPID(),"; host ",trim(hostName)
       write (0,'(a,i8,a)') " => Sleeping for ",errorWaitTime,"s to allow for attachment of debugger"
       call Flush(0)
       call Sleep(errorWaitTime)
    end if
#endif
    call H5Close_F(error)
    call H5Close_C(     )
    call Abort    (     )
#endif
    return
  end subroutine Galacticus_Signal_Handler_SIGSEGV

  subroutine Galacticus_Signal_Handler_SIGFPE()
    !% Handle {\normalfont \ttfamily SIGFPE} signals, by flushing all data and then aborting.
#ifndef UNCLEANEXIT
    use    :: HDF5   , only : H5Close_F
#endif
#ifdef USEMPI
    use    :: MPI    , only : MPI_Comm_Rank     , MPI_Comm_World
#endif
    !$ use :: OMP_Lib, only : OMP_Get_Thread_Num, OMP_In_Parallel
    implicit none
    integer            :: error
#ifdef USEMPI
    integer            :: mpiRank
    character(len=128) :: hostName
    logical            :: flag
#endif

    write (0,*) 'Galacticus experienced a floating point exception - will try to flush data before exiting.'
    !$ if (omp_in_parallel()) then
    !$    write (0,*) " => Error occurred in thread ",omp_get_thread_num()
    !$ else
    !$    write (0,*) " => Error occurred in master thread"
    !$ end if
    call BackTrace                (     )
    call Galacticus_Warn_Review   (     )
    call Flush                    (    0)
#ifdef UNCLEANEXIT
    call Exit(1)
#else
#ifdef USEMPI
    call MPI_Initialized(flag,error)
    if (flag) then
       call MPI_Comm_Rank(MPI_Comm_World,mpiRank,error)
       call hostnm(hostName)
       write (0,*) " => Error occurred in MPI process ",mpiRank,"; PID ",getPID(),"; host ",trim(hostName)
       write (0,'(a,i8,a)') " => Sleeping for ",errorWaitTime,"s to allow for attachment of debugger"
       call Flush(0)
       call Sleep(errorWaitTime)
    end if
#endif
    call H5Close_F(error)
    call H5Close_C(     )
    call Abort    (     )
#endif
    return
  end subroutine Galacticus_Signal_Handler_SIGFPE

  subroutine Galacticus_Signal_Handler_SIGBUS()
    !% Handle {\normalfont \ttfamily SIGBUS} signals, by flushing all data and then aborting.
#ifndef UNCLEANEXIT
    use    :: HDF5   , only : H5Close_F
#endif
#ifdef USEMPI
    use    :: MPI    , only : MPI_Comm_Rank     , MPI_Comm_World
#endif
    !$ use :: OMP_Lib, only : OMP_Get_Thread_Num, OMP_In_Parallel
    implicit none
    integer            :: error
#ifdef USEMPI
    integer            :: mpiRank
    character(len=128) :: hostName
    logical            :: flag
#endif

    write (0,*) 'Galacticus experienced a bus error - will try to flush data before exiting.'
    !$ if (omp_in_parallel()) then
    !$    write (0,*) " => Error occurred in thread ",omp_get_thread_num()
    !$ else
    !$    write (0,*) " => Error occurred in master thread"
    !$ end if
    call BackTrace                (     )
    call Galacticus_Warn_Review   (     )
    call Flush                    (    0)
#ifdef UNCLEANEXIT
    call Exit(1)
#else
#ifdef USEMPI
    call MPI_Initialized(flag,error)
    if (flag) then
       call MPI_Comm_Rank(MPI_Comm_World,mpiRank,error)
       call hostnm(hostName)
       write (0,*) " => Error occurred in MPI process ",mpiRank,"; PID ",getPID(),"; host ",trim(hostName)
       write (0,'(a,i8,a)') " => Sleeping for ",errorWaitTime,"s to allow for attachment of debugger"
       call Flush(0)
       call Sleep(errorWaitTime)
    end if
#endif
    call H5Close_F(error)
    call H5Close_C(     )
    call Abort    (     )
#endif
    return
  end subroutine Galacticus_Signal_Handler_SIGBUS

  subroutine Galacticus_Signal_Handler_SIGILL()
    !% Handle {\normalfont \ttfamily SIGILL} signals, by flushing all data and then aborting.
#ifndef UNCLEANEXIT
    use    :: HDF5   , only : H5Close_F
#endif
#ifdef USEMPI
    use    :: MPI    , only : MPI_Comm_Rank     , MPI_Comm_World
#endif
    !$ use :: OMP_Lib, only : OMP_Get_Thread_Num, OMP_In_Parallel
    implicit none
    integer            :: error
#ifdef USEMPI
    integer            :: mpiRank
    character(len=128) :: hostName
    logical            :: flag
#endif

    write (0,*) 'Galacticus experienced an illegal instruction - will try to flush data before exiting.'
    !$ if (omp_in_parallel()) then
    !$    write (0,*) " => Error occurred in thread ",omp_get_thread_num()
    !$ else
    !$    write (0,*) " => Error occurred in master thread"
    !$ end if
    call BackTrace                (     )
    call Galacticus_Warn_Review   (     )
    call Flush                    (    0)
#ifdef UNCLEANEXIT
    call Exit(1)
#else
#ifdef USEMPI
    call MPI_Initialized(flag,error)
    if (flag) then
       call MPI_Comm_Rank(MPI_Comm_World,mpiRank,error)
       call hostnm(hostName)
       write (0,*) " => Error occurred in MPI process ",mpiRank,"; PID ",getPID(),"; host ",trim(hostName)
       write (0,'(a,i8,a)') " => Sleeping for ",errorWaitTime,"s to allow for attachment of debugger"
       call Flush(0)
       call Sleep(errorWaitTime)
    end if
#endif
    call H5Close_F                (error)
    call H5Close_C                (     )
    call Abort                    (     )
#endif
    return
  end subroutine Galacticus_Signal_Handler_SIGILL

  subroutine Galacticus_Signal_Handler_SIGXCPU()
    !% Handle {\normalfont \ttfamily SIGXCPU} signals, by flushing all data and then aborting.
#ifndef UNCLEANEXIT
    use :: HDF5, only : H5Close_F
#endif
    implicit none
    integer :: error

    write (0,*) 'Galacticus exceeded available CPU time - will try to flush data before exiting.'
    call Flush(0)
#ifndef UNCLEANEXIT
    call H5Close_F(error)
    call H5Close_C()
#endif
    call Exit(errorStatusXCPU)
    return
  end subroutine Galacticus_Signal_Handler_SIGXCPU

  subroutine Galacticus_GSL_Error_Handler(reason,file,line,errorNumber) bind(c)
    !% Handle errors from the GSL library, by flushing all data and then aborting.
#ifndef UNCLEANEXIT
    use               :: HDF5              , only : H5Close_F
#endif
    use   , intrinsic :: ISO_C_Binding     , only : c_char
    use               :: ISO_Varying_String, only : char
#ifdef USEMPI
    use               :: MPI               , only : MPI_Comm_Rank      , MPI_Comm_World , MPI_Initialized
#endif
    !$ use            :: OMP_Lib           , only : OMP_Get_Thread_Num , OMP_In_Parallel
    use               :: String_Handling   , only : String_C_To_Fortran
    character(c_char), dimension(*) :: file       , reason
    integer  (c_int ), value        :: errorNumber, line
    integer                         :: error
#ifdef USEMPI
    integer                         :: mpiRank
    character(len=128)              :: hostName
    logical                         :: flag
#endif

    if (abortOnErrorGSL) then
       write (0,*) 'Galacticus experienced an error in the GSL library - will try to flush data before exiting.'
       write (0,*) ' => Error occurred in ',char(String_C_to_Fortran(file  )),' at line ',line
       write (0,*) ' => Reason was: '      ,char(String_C_to_Fortran(reason))
       !$ if (omp_in_parallel()) then
       !$    write (0,*) " => Error occurred in thread ",omp_get_thread_num()
       !$ else
       !$    write (0,*) " => Error occurred in master thread"
       !$ end if
    call BackTrace                (     )
    call Galacticus_Warn_Review   (     )
    call Flush                    (    0)
#ifdef UNCLEANEXIT
    call Exit(1)
#else
#ifdef USEMPI
    call MPI_Initialized(flag,error)
    if (flag) then
       call MPI_Comm_Rank(MPI_Comm_World,mpiRank,error)
       call hostnm(hostName)
       write (0,*) " => Error occurred in MPI process ",mpiRank,"; PID ",getPID(),"; host ",trim(hostName)
       write (0,'(a,i8,a)') " => Sleeping for ",errorWaitTime,"s to allow for attachment of debugger"
       call Flush(0)
       call Sleep(errorWaitTime)
    end if
#endif
    call H5Close_F(error)
    call H5Close_C(     )
    call Abort    (     )
#endif
    else
       errorStatusGSL=errorNumber
    end if
    return
  end subroutine Galacticus_GSL_Error_Handler

  subroutine Galacticus_GSL_Error_Handler_Abort_On()
    !% Record that we should abort on GSL errors.
    implicit none

    abortOnErrorGSL=.true.
    return
  end subroutine Galacticus_GSL_Error_Handler_Abort_On

  subroutine Galacticus_GSL_Error_Handler_Abort_Off()
    !% Record that we should not abort on GSL errors.
    implicit none

    abortOnErrorGSL=.false.
    return
  end subroutine Galacticus_GSL_Error_Handler_Abort_Off

  integer function Galacticus_GSL_Error_Status()
    !% Return current GSL error status.
    implicit none

    Galacticus_GSL_Error_Status=errorStatusGSL
    return
  end function Galacticus_GSL_Error_Status

  function Galacticus_Component_List(className,componentList)
    !% Construct a message describing which implementations of a component class provide required functionality.
    use :: ISO_Varying_String, only : assignment(=), operator(//)
    use :: String_Handling   , only : String_Join
    implicit none
    type     (varying_string)                                           :: Galacticus_Component_List
    character(len=*         ), intent(in   )                            :: className
    type     (varying_string), intent(in   ), dimension(:), allocatable :: componentList

    if (allocated(componentList)) then
       Galacticus_Component_List=char(10)//'Implementations of the "'   //className//'" class that provide this functionality are:'// &
            & char(10)//'   '//String_Join(componentList,char(10)//'   ')
    else
       Galacticus_Component_List=char(10)//'No implementations of the "'//className//'" class currently provide this functionality.'
    end if
    return
  end function Galacticus_Component_List

  subroutine Galacticus_Error_Wait_Set(errorWaitTimeNew)
    !% Set the time to wait after an error occurs.
    implicit none
    integer, intent(in   ) :: errorWaitTimeNew

    errorWaitTime=errorWaitTimeNew
    return
  end subroutine Galacticus_Error_Wait_Set

end module Galacticus_Error
