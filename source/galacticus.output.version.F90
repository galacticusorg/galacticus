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






!% Contains a module which implements writing of the version number and run time to the \glc\ output file.

module Galacticus_Version
  !% Implements writing of the version number and run time to the \glc\ output file.
  private
  public :: Galacticus_Version_Output

contains

  !# <outputFileOpenTask>
  !#  <unitName>Galacticus_Version_Output</unitName>
  !# </outputFileOpenTask>
  subroutine Galacticus_Version_Output
    !% Output version information to the main output file.
    use Galacticus_HDF5_Groups
    use ISO_Varying_String
    use HDF5
    use Galacticus_Error
    use Dates_and_Times
    implicit none
    integer(kind=HID_T)                :: versionGroupID=0,versionMajorID=0,versionMinorID=0,versionRevisionID=0,runTimeID=0
    type(varying_string)               :: groupName,groupComment
    type(varying_string), dimension(1) :: runTime

    groupName='Version'
    groupComment='Version and timestamp for this model.'
    versionGroupID=Galacticus_Output_Make_Group(groupName,groupComment)
    call Galacticus_Output_Dataset(versionGroupID,versionMajorID   ,'versionMajor'   ,'Major version number',[0])
    call Galacticus_Output_Dataset(versionGroupID,versionMinorID   ,'versionMinor'   ,'Minor version number',[9])
    call Galacticus_Output_Dataset(versionGroupID,versionRevisionID,'versionRevision','Revision number'     ,[0])
    runTime(1)=Formatted_Date_and_Time()
    call Galacticus_Output_Dataset(versionGroupID,runTimeID        ,'runTime'        ,'Time at which model was run',runTime)
    return
  end subroutine Galacticus_Version_Output
  
end module Galacticus_Version
