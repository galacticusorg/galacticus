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

!% Contains a module which provides the path for \glc\ inputs and scripts.

module Galacticus_Input_Paths
  !% Provides the path for \glc\ inputs and scripts.
  use ISO_Varying_String
  implicit none
  private
  public :: Galacticus_Input_Path

  ! The input path.
  type(varying_string) :: galacticusInputPath

  ! Flag indicating if we've retrieved the path yet.
  logical              :: pathRetrieved=.false.

contains

  function Galacticus_Input_Path()
    !% Returns the path that \glc\ should use as the root for all input data reads.
    implicit none
    type(varying_string) :: Galacticus_Input_Path
    integer              :: pathLength,pathStatus

    ! Check if the path has been retrieved.
    if (.not.pathRetrieved) then
       !$omp critical (Galacticus_Path_Initialize)
       if (.not.pathRetrieved) then
          ! Get the status and path length.
          call Get_Environment_Variable("GALACTICUS_ROOT_V091",length=pathLength,status=pathStatus) 
          if (pathStatus == 0) then
             ! Path is defined, retrieve it.
             call Get_Path(pathLength)
          else
             ! No path is defined, default to current working directory.
             galacticusInputPath="./"
          end if
          ! Flag that the path has been retrieved.
          pathRetrieved=.true.
       end if
       !$omp end critical (Galacticus_Path_Initialize)
    end if

    ! Return the input path.
    Galacticus_Input_Path=trim(char(galacticusInputPath))
    return
  end function Galacticus_Input_Path
  
  subroutine Get_Path(pathLength)
    !% Retrieve the \glc\ input data path from the environment.
    implicit none
    integer,                  intent(in) :: pathLength
    character(len=pathLength+1)            :: pathName
    
    ! Get the path.
    call Get_Environment_Variable("GALACTICUS_ROOT_V091",value=pathName)
    ! Append a final "/" if necessary.
    if (pathName(pathLength:pathLength) /= "/") pathName=trim(pathName)//"/"
    ! Store the path.
    galacticusInputPath=pathName
    return
  end subroutine Get_Path

end module Galacticus_Input_Paths
