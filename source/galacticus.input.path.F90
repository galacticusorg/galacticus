!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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
Contains a module which provides the path for \glc\ inputs and scripts.
!!}

module Galacticus_Paths
  !!{
  Provides the path for \glc\ inputs and scripts.
  !!}
  use :: ISO_Varying_String, only : varying_string
  implicit none
  private
  public :: galacticusPath

  !![
  <enumeration>
   <name>pathType</name>
   <description>Enumeration of various paths used by Galacticus.</description>
   <validator>yes</validator>
   <entry label="exec"       />
   <entry label="dataStatic" />
   <entry label="dataDynamic"/>
  </enumeration>
  !!]

  type   (varying_string), dimension(:), allocatable :: paths
  logical                                            :: pathsRetrieved=.false.

contains

  function galacticusPath(pathType)
    !!{
    Returns the path to various \glc\ resources.
    !!}
    use :: ISO_Varying_String, only : assignment(=), operator(//), trim, char
    implicit none
    type   (varying_string)                :: galacticusPath
    integer                , intent(in   ) :: pathType
    integer                                :: pathLength    , pathStatus

    ! Retrieve the paths if necessary.
    if (.not.pathsRetrieved) then
       !$omp critical (Galacticus_Path_Initialize)
       if (.not.pathsRetrieved) then
          allocate(paths(pathTypeMin:pathTypeMax))
          call Get_Environment_Variable("GALACTICUS_EXEC_PATH",length=pathLength,status=pathStatus)
          if (pathStatus == 0) then
             call pathsRetrieve(pathTypeExec,"GALACTICUS_EXEC_PATH",pathLength)
          else
             paths(pathTypeExec)="./"
          end if
          call Get_Environment_Variable("GALACTICUS_DATA_PATH",length=pathLength,status=pathStatus)
          if (pathStatus == 0) then
             call pathsRetrieve(pathTypeDataStatic,"GALACTICUS_DATA_PATH",pathLength)
          else
             paths(pathTypeDataStatic)="./"
          end if
          paths(pathTypeDataDynamic)=paths(pathTypeDataStatic)//"dynamic/"
          paths(pathTypeDataStatic )=paths(pathTypeDataStatic)//"static/"
          call Get_Environment_Variable("GALACTICUS_DYNAMIC_DATA_PATH",length=pathLength,status=pathStatus)
          if (pathStatus == 0)                                                                     &
               & call pathsRetrieve(pathTypeDataDynamic,"GALACTICUS_DYNAMIC_DATA_PATH",pathLength)
          pathsRetrieved=.true.
       end if
       !$omp end critical (Galacticus_Path_Initialize)
    end if
    galacticusPath=trim(char(paths(pathType)))
    return
  end function galacticusPath

  subroutine pathsRetrieve(pathType,pathName,pathLength)
    !!{
    Retrieve the \glc\ input data path from the environment.
    !!}
    use :: ISO_Varying_String, only : assignment(=)
    implicit none
    integer                    , intent(in   ) :: pathType, pathLength
    character(len=*           ), intent(in   ) :: pathName
    character(len=pathLength+1)                :: pathValue

    ! Get the path.
    call Get_Environment_Variable(pathName,value=pathValue)
    ! Append a final "/" if necessary.
    if (pathValue(pathLength:pathLength) /= "/") pathValue=trim(pathValue)//"/"
    ! Store the path.
    paths(pathType)=pathValue
    return
  end subroutine pathsRetrieve

end module Galacticus_Paths
