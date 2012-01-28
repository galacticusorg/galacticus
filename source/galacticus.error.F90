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
    implicit none
    character(len=*), intent(in), optional :: unitName,message
    integer                                :: error

    if (present(unitName)) write (0,'(a,a,a)') 'Fatal error in ',trim(unitName),'():'
    if (present(message )) write (0,'(a)'    ) trim(message)
    call Flush(0)
    call H5Close_F(error)
    call H5Close_C()
    call Abort()
    return
  end subroutine Galacticus_Error_Report_Char

  subroutine Galacticus_Error_Handler_Register()
    !% Register signal handlers.
    implicit none
    
    call Signal( 2,Galacticus_Signal_Handler_SIGINT )
    call Signal( 8,Galacticus_Signal_Handler_SIGFPE )
    call Signal(11,Galacticus_Signal_Handler_SIGSEGV)
    call Signal(15,Galacticus_Signal_Handler_SIGINT )
    return
  end subroutine Galacticus_Error_Handler_Register
  
  subroutine Galacticus_Signal_Handler_SIGINT()
    !% Handle {\tt SIGINT} signals, by flushing all data and then aborting.
    implicit none
    integer :: error

    write (0,*) 'Galacticus was interrupted - will try to flush data before exiting.'
    call Flush(0)
    call H5Close_F(error)
    call H5Close_C()
    call Abort()
    return
  end subroutine Galacticus_Signal_Handler_SIGINT

  subroutine Galacticus_Signal_Handler_SIGSEGV()
    !% Handle {\tt SIGSEGV} signals, by flushing all data and then aborting.
    implicit none
    integer :: error

    write (0,*) 'Galacticus experienced a segfault - will try to flush data before exiting.'
    call Flush(0)
    call H5Close_F(error)
    call H5Close_C()
    call Abort()
    return
  end subroutine Galacticus_Signal_Handler_SIGSEGV

  subroutine Galacticus_Signal_Handler_SIGFPE()
    !% Handle {\tt SIGFPE} signals, by flushing all data and then aborting.
    implicit none
    integer :: error

    write (0,*) 'Galacticus experienced a floating point exception - will try to flush data before exiting.'
    call Flush(0)
    call H5Close_F(error)
    call H5Close_C()
    call Abort()
    return
  end subroutine Galacticus_Signal_Handler_SIGFPE

end module Galacticus_Error
