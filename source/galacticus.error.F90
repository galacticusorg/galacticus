!! Copyright 2009, 2010, Andrew Benson <abenson@caltech.edu>
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
  private
  public :: Galacticus_Error_Report
  
  interface Galacticus_Error_Report
     module procedure Galacticus_Error_Report_Char
     module procedure Galacticus_Error_Report_VarStr
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
    
    if (present(unitName)) write (0,'(a,a,a)') 'Fatal error in ',trim(unitName),'():'
    if (present(message )) write (0,'(a)'    ) trim(message)
    call Abort()
    return
  end subroutine Galacticus_Error_Report_Char

end module Galacticus_Error
