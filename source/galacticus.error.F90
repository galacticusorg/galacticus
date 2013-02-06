!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements error reporting for the {\sc Galacticus} package.

module Galacticus_Error
  !% Implements error reporting for the {\sc Galacticus} package.
  use HDF5
  implicit none
  private
  public :: Galacticus_Error_Report, Galacticus_Error_Handler_Register
  
  interface Galacticus_Error_Report
     module procedure Galacticus_Error_Report_Char
     module procedure Galacticus_Error_Report_VarStr
  end interface

  ! Specify an explicit dependence on the hdf5_cFuncs.o object file.
  !: ./work/build/hdf5_cFuncs.o
  interface
     subroutine H5Close_C() bind(c,name='H5Close_C')
     end subroutine H5Close_C
  end interface
 
contains

  subroutine Galacticus_Error_Report_VarStr(unitName,message)
    !% Display an error message.
    use ISO_Varying_String
    implicit none
    character(len=*),     intent(in) :: unitName
    type(varying_string), intent(in) :: message

    call Galacticus_Error_Report_Char(unitName,char(message))
    
    return
  end subroutine Galacticus_Error_Report_VarStr

  subroutine Galacticus_Error_Report_Char(unitName,message)
    !% Display an error message (optionally reporting the unit name in which the error originated) and stop.
    !$ use OMP_Lib
    implicit none
    character(len=*), intent(in), optional :: unitName,message
    integer                                :: error

    if (present(unitName)) write (0,'(a,a,a)') 'Fatal error in ',trim(unitName),'():'
    if (present(message )) write (0,'(a)'    ) trim(message)
    !$ if (omp_in_parallel()) then
    !$    write (0,*) " => Error occurred in thread ",omp_get_thread_num()
    !$ else
    !$    write (0,*) " => Error occurred in master thread"
    !$ end if
    call Flush(0)
    call H5Close_F(error)
    call H5Close_C()
    call Abort()
    return
  end subroutine Galacticus_Error_Report_Char

  subroutine Galacticus_Error_Handler_Register()
    !% Register signal handlers.
    use FGSL
    implicit none
    type(fgsl_error_handler_t) :: galacticusGslErrorHandler,standardGslErrorHandler

    call Signal( 2,Galacticus_Signal_Handler_SIGINT )
    call Signal( 8,Galacticus_Signal_Handler_SIGFPE )
    call Signal(11,Galacticus_Signal_Handler_SIGSEGV)
    call Signal(15,Galacticus_Signal_Handler_SIGINT )
    galacticusGslErrorHandler=FGSL_Error_Handler_Init(Galacticus_GSL_Error_Handler)
    standardGslErrorHandler  =FGSL_Set_Error_Handler (galacticusGslErrorHandler   )
   return
  end subroutine Galacticus_Error_Handler_Register
  
  subroutine Galacticus_Signal_Handler_SIGINT()
    !% Handle {\tt SIGINT} signals, by flushing all data and then aborting.
    !$ use OMP_Lib
    implicit none
    integer :: error

    write (0,*) 'Galacticus was interrupted - will try to flush data before exiting.'
    !$ if (omp_in_parallel()) then
    !$    write (0,*) " => Error occurred in thread ",omp_get_thread_num()
    !$ else
    !$    write (0,*) " => Error occurred in master thread"
    !$ end if
    call Flush(0)
    call H5Close_F(error)
    call H5Close_C()
    call Abort()
    return
  end subroutine Galacticus_Signal_Handler_SIGINT

  subroutine Galacticus_Signal_Handler_SIGSEGV()
    !% Handle {\tt SIGSEGV} signals, by flushing all data and then aborting.
    !$ use OMP_Lib
    implicit none
    integer :: error

    write (0,*) 'Galacticus experienced a segfault - will try to flush data before exiting.'
    !$ if (omp_in_parallel()) then
    !$    write (0,*) " => Error occurred in thread ",omp_get_thread_num()
    !$ else
    !$    write (0,*) " => Error occurred in master thread"
    !$ end if
    call Flush(0)
    call H5Close_F(error)
    call H5Close_C()
    call Abort()
    return
  end subroutine Galacticus_Signal_Handler_SIGSEGV

  subroutine Galacticus_Signal_Handler_SIGFPE()
    !% Handle {\tt SIGFPE} signals, by flushing all data and then aborting.
    !$ use OMP_Lib
    implicit none
    integer :: error

    write (0,*) 'Galacticus experienced a floating point exception - will try to flush data before exiting.'
    !$ if (omp_in_parallel()) then
    !$    write (0,*) " => Error occurred in thread ",omp_get_thread_num()
    !$ else
    !$    write (0,*) " => Error occurred in master thread"
    !$ end if
    call Flush(0)
    call H5Close_F(error)
    call H5Close_C()
    call Abort()
    return
  end subroutine Galacticus_Signal_Handler_SIGFPE

  subroutine Galacticus_GSL_Error_Handler(reason,file,line,errorNumber) bind(c)
    !% Handle errors from the GSL library, by flushing all data and then aborting.
    !$ use OMP_Lib
    use FGSL
    use, intrinsic :: ISO_C_Binding
    type     (c_ptr                         ), value :: reason, file
    integer  (c_int                         ), value :: line, errorNumber
    character(kind=FGSL_Char,len=FGSL_StrMax)        :: message
    integer                                          :: error

    message=FGSL_StrError(errorNumber)
    write (0,*) 'Galacticus experienced an error in the GSL library - will try to flush data before exiting.'
    write (0,*) ' => Error occurred in ',trim(FGSL_Name(file  )),' at line ',line
    write (0,*) ' => Reason was: '      ,trim(FGSL_Name(reason))
    !$ if (omp_in_parallel()) then
    !$    write (0,*) " => Error occurred in thread ",omp_get_thread_num()
    !$ else
    !$    write (0,*) " => Error occurred in master thread"
    !$ end if
    call Flush(0)
    call H5Close_F(error)
    call H5Close_C()
    call Abort()
    return
  end subroutine Galacticus_GSL_Error_Handler

end module Galacticus_Error
