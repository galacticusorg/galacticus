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

!!{RST
Contains a module which provides the path for Galacticus inputs and scripts.
!!}

module Input_Paths
  !!{RST
  Provides the path for Galacticus inputs and scripts.
  !!}
  use :: ISO_Varying_String, only : varying_string
  implicit none
  private
  public :: inputPath

  !![
  <enumeration docformat="rst">
   <name>pathType</name>
   <description>
   Enumeration of various paths used by Galacticus.
   </description>
   <validator>yes</validator>
   <entry label="exec"       />
   <entry label="dataStatic" />
   <entry label="dataDynamic"/>
   <entry label="tools"      />
  </enumeration>
  !!]

  type   (varying_string), dimension(:), allocatable :: paths
  logical                                            :: pathsRetrieved=.false.

contains

  function inputPath(pathType)
    !!{RST
    Returns the path to various Galacticus resources.
    !!}
    use :: ISO_Varying_String, only : assignment(=)  , operator(//), trim, char
    implicit none
    type   (varying_string         )                :: inputPath
    type   (enumerationPathTypeType), intent(in   ) :: pathType
    integer                                         :: pathLength, pathStatus

    ! Retrieve the paths if necessary.
    if (.not.pathsRetrieved) then
       !$omp critical (Input_Path_Initialize)
       if (.not.pathsRetrieved) then
          allocate(paths(pathTypeMin:pathTypeMax))
          call Get_Environment_Variable("GALACTICUS_EXEC_PATH",length=pathLength,status=pathStatus)
          if (pathStatus == 0) then
             call pathsRetrieve(pathTypeExec,"GALACTICUS_EXEC_PATH",pathLength)
          else
             paths(pathTypeExec      %ID)="./"
          end if
          call Get_Environment_Variable("GALACTICUS_DATA_PATH",length=pathLength,status=pathStatus)
          if (pathStatus == 0) then
             call pathsRetrieve(pathTypeDataStatic,"GALACTICUS_DATA_PATH",pathLength)
          else
             paths(pathTypeDataStatic%ID)="./"
          end if
          paths(pathTypeDataDynamic%ID)=paths(pathTypeDataStatic%ID)//"dynamic/"
          paths(pathTypeDataStatic %ID)=paths(pathTypeDataStatic%ID)//"static/"
          call Get_Environment_Variable("GALACTICUS_DYNAMIC_DATA_PATH",length=pathLength,status=pathStatus)
          if (pathStatus == 0)                                                                     &
               & call pathsRetrieve(pathTypeDataDynamic,"GALACTICUS_DYNAMIC_DATA_PATH",pathLength)
          ! Pre-built external tools (CAMB, CLASS, Cloudy, FSPS, RecFast,
          ! AxionCAMB, mangle) default to living under the dynamic data path, but
          ! can be relocated independently via the GALACTICUS_TOOLS_PATH
          ! environment variable. This allows a managed (e.g. pip) install to keep
          ! immutable, pre-built tool binaries separate from the regenerable
          ! dynamic data cache, so that the cache can be purged without losing the
          ! tools (which a binary-only install cannot rebuild).
          paths(pathTypeTools%ID)=paths(pathTypeDataDynamic%ID)
          call Get_Environment_Variable("GALACTICUS_TOOLS_PATH",length=pathLength,status=pathStatus)
          if (pathStatus == 0)                                                                     &
               & call pathsRetrieve(pathTypeTools,"GALACTICUS_TOOLS_PATH",pathLength)
          pathsRetrieved=.true.
       end if
       !$omp end critical (Input_Path_Initialize)
    end if
    inputPath=trim(char(paths(pathType%ID)))
    return
  end function inputPath

  subroutine pathsRetrieve(pathType,pathName,pathLength)
    !!{RST
    Retrieve the Galacticus input data path from the environment.
    !!}
    use :: ISO_Varying_String, only : assignment(=)
    implicit none
    type(enumerationPathTypeType), intent(in   ) :: pathType
    integer                      , intent(in   ) :: pathLength
    character(len=*             ), intent(in   ) :: pathName
    character(len=pathLength+1  )                :: pathValue

    ! Get the path.
    call Get_Environment_Variable(pathName,value=pathValue)
    ! Append a final "/" if necessary.
    if (pathValue(pathLength:pathLength) /= "/") pathValue=trim(pathValue)//"/"
    ! Store the path.
    paths(pathType%ID)=trim(pathValue)
    return
  end subroutine pathsRetrieve

end module Input_Paths
