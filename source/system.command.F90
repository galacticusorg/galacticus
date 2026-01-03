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
Contains a module which executes system commands.
!!}

module System_Command
  !!{
  Executes system commands.
  !!}
  implicit none
  private
  public :: System_Command_Do

  interface System_Command_Do
     module procedure System_Command_Char
     module procedure System_Command_VarStr
  end interface System_Command_Do

contains

  subroutine System_Command_VarStr(command,iStatus)
    !!{
    Executes the system command {\normalfont \ttfamily command}, optionally returning the resulting status in {\normalfont \ttfamily iStatus}.
    !!}
    use :: ISO_Varying_String, only : char, varying_string
    implicit none
    type   (varying_string), intent(in   )           :: command
    integer                , intent(  out), optional :: iStatus

    call System_Command_Char(char(command),iStatus)
    return
  end subroutine System_Command_VarStr

  subroutine System_Command_Char(command,iStatus)
    !!{
    Executes the system command {\normalfont \ttfamily command}, optionally returning the resulting status in {\normalfont \ttfamily iStatus}.
    !!}
    use :: Error, only : Error_Report
    implicit none
    character(len=*), intent(in   )           :: command
    integer         , intent(  out), optional :: iStatus
    integer                                   :: iStatusActual

    call System(command,iStatusActual)
    if (present(iStatus)) then
       iStatus=iStatusActual
    else
       if (iStatusActual /= 0) call Error_Report('failed to execute system command:'//char(10)//' --> '//command//{introspection:location})
    end if
    return
  end subroutine System_Command_Char

end module System_Command
