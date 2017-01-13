!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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
  use, intrinsic :: ISO_C_Binding
  use HDF5
  use Semaphores
  use FGSL
  use Semaphores
  use ISO_Varying_String
  implicit none
  private
  public :: Galacticus_Error_Report               , Galacticus_Error_Handler_Register    , &
       &    Galacticus_Component_List             , Galacticus_GSL_Error_Handler_Abort_On, &
       &    Galacticus_GSL_Error_Handler_Abort_Off, Galacticus_GSL_Error_Status          , &
       &    Galacticus_Warn

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
  integer, parameter, public :: errorStatusSuccess     =FGSL_Success  ! Success.
  integer, parameter, public :: errorStatusFail        =FGSL_Failure  ! Generic failure.
  integer, parameter, public :: errorStatusInputDomain =FGSL_eDom     ! Input domain error.
  integer, parameter, public :: errorStatusOutOfRange  =FGSL_eRange   ! Output range error.
  integer, parameter, public :: errorStatusDivideByZero=FGSL_eZeroDiv ! Divide by zero.
  integer, parameter, public :: errorStatusXCPU        =1025          ! CPU time limit exceeded.

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

  subroutine Galacticus_Error_Report_VarStr(unitName,message)
    !% Display an error message.
    implicit none
    character(len=*         ), intent(in   ) :: unitName
    type     (varying_string), intent(in   ) :: message

    call Galacticus_Error_Report_Char(unitName,char(message))

    return
  end subroutine Galacticus_Error_Report_VarStr

  subroutine Galacticus_Error_Report_Char(unitName,message)
    !% Display an error message (optionally reporting the unit name in which the error originated) and stop.
    !$ use OMP_Lib
    implicit none
    character(len=*), intent(in   ), optional :: message, unitName
    integer                                   :: error

    if (present(unitName)) write (0,'(a,a,a)') 'Fatal error in ',trim(unitName),'():'
    if (present(message )) write (0,'(a)'    ) trim(message)
    !$ if (omp_in_parallel()) then
    !$    write (0,*) " => Error occurred in thread ",omp_get_thread_num()
    !$ else
    !$    write (0,*) " => Error occurred in master thread"
    !$ end if
    call Galacticus_Warn_Review   (     )
    call Flush                    (    0)
#ifdef UNCLEANEXIT
    call Exit(1)
#else
    call H5Close_F                (error)
    call H5Close_C                (     )
    call Semaphore_Post_On_Error  (     )
    call Abort                    (     )
#endif
    return
  end subroutine Galacticus_Error_Report_Char

  subroutine Galacticus_Warn_VarStr(message)
    !% Display a warning message
    implicit none
    type(varying_string), intent(in   ) :: message

    call Galacticus_Warn_Char(char(message))
    return
  end subroutine Galacticus_Warn_VarStr

  subroutine Galacticus_Warn_Char(message)
    !% Display a warning message.
    use Galacticus_Display
    implicit none
    character(len=*  ), intent(in   ) :: message
    type     (warning), pointer       :: newWarning

    ! Display the message.
    call Galacticus_Display_Message(message,verbosity=verbosityWarn)
    ! Add this warning message to the list of warnings in case we need to display them on an
    ! error condition.
    !$omp critical (Galacticus_Warn)
    if (Galacticus_Verbosity_Level() < verbosityWarn) then
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
    use FGSL
    implicit none
    type(fgsl_error_handler_t) :: galacticusGslErrorHandler, standardGslErrorHandler

    call Signal( 2,Galacticus_Signal_Handler_SIGINT )
    call Signal( 8,Galacticus_Signal_Handler_SIGFPE )
    call Signal(11,Galacticus_Signal_Handler_SIGSEGV)
    call Signal(15,Galacticus_Signal_Handler_SIGINT )
    call Signal(24,Galacticus_Signal_Handler_SIGXCPU)
    galacticusGslErrorHandler=FGSL_Error_Handler_Init(Galacticus_GSL_Error_Handler)
    standardGslErrorHandler  =FGSL_Set_Error_Handler (galacticusGslErrorHandler   )
   return
  end subroutine Galacticus_Error_Handler_Register

  subroutine Galacticus_Signal_Handler_SIGINT()
    !% Handle {\normalfont \ttfamily SIGINT} signals, by flushing all data and then aborting.
    !$ use OMP_Lib
#ifdef USEMPI
    use MPI
#endif
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
       write (0,*) " => Sleeping for 86400s to allow for attachment of debugger"
       call Flush(0)
       call Sleep(86400)
    end if
#endif
    call H5Close_F                (error)
    call H5Close_C                (     )
    call Semaphore_Post_On_Error  (     )
    call Abort                    (     )
#endif
    return
  end subroutine Galacticus_Signal_Handler_SIGINT

  subroutine Galacticus_Signal_Handler_SIGSEGV()
    !% Handle {\normalfont \ttfamily SIGSEGV} signals, by flushing all data and then aborting.
    !$ use OMP_Lib
#ifdef USEMPI
    use MPI
#endif
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
       write (0,*) " => Sleeping for 86400s to allow for attachment of debugger"
       call Flush(0)
       call Sleep(86400)
    end if
#endif
    call H5Close_F                (error)
    call H5Close_C                (     )
    call Semaphore_Post_On_Error  (     )
    call Abort                    (     )
#endif
    return
  end subroutine Galacticus_Signal_Handler_SIGSEGV

  subroutine Galacticus_Signal_Handler_SIGFPE()
    !% Handle {\normalfont \ttfamily SIGFPE} signals, by flushing all data and then aborting.
    !$ use OMP_Lib
#ifdef USEMPI
    use MPI
#endif
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
       write (0,*) " => Sleeping for 86400s to allow for attachment of debugger"
       call Flush(0)
       call Sleep(86400)
    end if
#endif
    call H5Close_F                (error)
    call H5Close_C                (     )
    call Semaphore_Post_On_Error  (     )
    call Abort                    (     )
#endif
    return
  end subroutine Galacticus_Signal_Handler_SIGFPE

  subroutine Galacticus_Signal_Handler_SIGXCPU()
    !% Handle {\normalfont \ttfamily SIGXCPU} signals, by flushing all data and then aborting.
    implicit none
    integer :: error

    write (0,*) 'Galacticus exceeded available CPU time - will try to flush data before exiting.'
    call Semaphore_Post_On_Error()
    call Flush(0)
#ifndef UNCLEANEXIT
    call H5Close_F(error)
    call H5Close_C()
    call Exit(errorStatusXCPU)
#endif
    return
  end subroutine Galacticus_Signal_Handler_SIGXCPU

  subroutine Galacticus_GSL_Error_Handler(reason,file,line,errorNumber) bind(c)
    !% Handle errors from the GSL library, by flushing all data and then aborting.
    !$ use OMP_Lib
#ifdef USEMPI
    use MPI
#endif
    use FGSL
    type     (c_ptr                         ), value :: file       , reason
    integer  (kind=c_int                    ), value :: errorNumber, line
    character(kind=FGSL_Char,len=FGSL_StrMax)        :: message
    integer                                          :: error
#ifdef USEMPI
    integer                                          :: mpiRank
    character(len=128                       )        :: hostName
    logical                                          :: flag
#endif
    
    if (abortOnErrorGSL) then
       message=FGSL_StrError(errorNumber)
       write (0,*) 'Galacticus experienced an error in the GSL library - will try to flush data before exiting.'
       write (0,*) ' => Error occurred in ',trim(FGSL_Name(file  )),' at line ',line
       write (0,*) ' => Reason was: '      ,trim(FGSL_Name(reason))
       !$ if (omp_in_parallel()) then
       !$    write (0,*) " => Error occurred in thread ",omp_get_thread_num()
       !$ else
       !$    write (0,*) " => Error occurred in master thread"
       !$ end if
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
          write (0,*) " => Sleeping for 86400s to allow for attachment of debugger"
          call Flush(0)
          call Sleep(86400)
       end if
#endif
    call H5Close_F                (error)
    call H5Close_C                (     )
    call Semaphore_Post_On_Error  (     )
    call Abort                    (     )
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
    use String_Handling
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
  
end module Galacticus_Error
