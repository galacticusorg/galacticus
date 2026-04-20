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
Contains a module which wraps the system \mono{which} command to allow finding of other tools.
!!}

module System_Which
  !!{
  Wraps the system \mono{which} command to allow finding of other tools.
  !!}
  implicit none
  private
  public :: which

  interface which
     module procedure whichChar
     module procedure whichVarStr
  end interface which

contains

  function whichVarStr(command,status) result(commandFull)
    !!{
    Find the path to the given \mono{command}, optionally returning status.
    !!}
    use :: ISO_Varying_String, only : char, varying_string
    implicit none
    type   (varying_string)                          :: commandFull
    type   (varying_string), intent(in   )           :: command
    integer                , intent(  out), optional :: status
    
    commandFull=which(char(command),status)
    return
  end function whichVarStr

  function whichChar(command,status) result(commandFull)
    !!{
    Find the path to the given \mono{command}, optionally returning status.
    !!}
    use :: Error             , only : Error_Report       , errorStatusFail, errorStatusSuccess
    use :: System_Command    , only : System_Command_Do
    use :: File_Utilities    , only : File_Name_Temporary, File_Remove
    use :: ISO_Varying_String, only : varying_string     , char           , trim              , assignment(=)
    implicit none
    type     (varying_string)                          :: commandFull
    character(len=*         ), intent(in   )           :: command
    integer                  , intent(  out), optional :: status
    character(len=1024      )                          :: commandFull_
    integer                                            :: status_     , commandUnit
    type     (varying_string)                          :: tempFile
    
    if (present(status)) status=errorStatusSuccess
    tempFile=File_Name_Temporary("which")
    call System_Command_Do('which '//trim(command)//' > '//char(tempFile),status_)
    if (status_ == 0) then
       open(newUnit=commandUnit,file=char(tempFile),status='unknown',form='formatted')
       read (commandUnit,'(a)') commandFull_ 
       close(commandUnit)
       commandFull=trim(commandFull_)
    else
       commandFull=""
    end if
    call File_Remove(tempFile)
    if (     present(status)                   ) status=status_
    if (.not.present(status) .and. status_ /= 0) call Error_Report("command '"//trim(command)//"' does not exist"//{introspection:location})
    return
  end function whichChar

end module System_Which
