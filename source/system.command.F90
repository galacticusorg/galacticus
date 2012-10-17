!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which executes system commands.

module System_Command
  !% Executes system commands.
  implicit none
  private
  public :: System_Command_Do
  
contains
  
  subroutine System_Command_Do(command,iStatus)
    !% Executes the system command {\tt command}, optionally returning the resulting status in {\tt iStatus}.
    use ISO_Varying_String
    use Galacticus_Error
    implicit none
    type(varying_string), intent(in)            :: command
    integer,              intent(out), optional :: iStatus
    integer                                     :: iStatusActual

    call System(char(command),iStatusActual)
    if (present(iStatus)) then
       iStatus=iStatusActual
    else
       if (iStatusActual /= 0) call Galacticus_Error_Report('System_Command_Do','failed to execute system command')
    end if
    return
  end subroutine System_Command_Do
  
end module System_Command
