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

!+    Contributions to this file made by: Andrew Benson.

!!{
Contains a module which implements calculations of gravitational accelerations.
!!}

module Galactic_Structure_Accelerations
  !!{
  Implements calculations of the gravitational accelerations.
  !!}
  use :: Kind_Numbers, only : kind_int8
  implicit none
  private
  public :: Galactic_Structure_Acceleration

  ! Module scope variables used in mapping over components.
  integer                                        :: componentTypeShared                    , massTypeShared
  double precision                , dimension(3) :: positionCartesianShared
  !$omp threadprivate(massTypeShared,componentTypeShared,positionCartesianShared)

contains

  function Galactic_Structure_Acceleration(node,positionCartesian,componentType,massType)
    !!{
    Compute the gravitational acceleration at a given position.
    !!}
    use :: Galactic_Structure_Options, only : componentTypeAll                , massTypeAll
    use :: Galacticus_Nodes          , only : optimizeForAccelerationSummation, reductionSummation, treeNode
    !![
    <include directive="accelerationTask" type="moduleUse">
    !!]
    include 'galactic_structure.acceleration.tasks.modules.inc'
    !![
    </include>
    !!]
    implicit none
    double precision                                       , dimension(3) :: Galactic_Structure_Acceleration
    type            (treeNode              ), intent(inout)               :: node
    double precision                        , intent(in   ), dimension(3) :: positionCartesian
    integer                                 , intent(in   ), optional     :: componentType                    , massType
    integer                                 , parameter                   :: accelerationSize               =3
    procedure       (Component_Acceleration), pointer                     :: componentAccelerationFunction
    double precision                                       , dimension(3) :: componentAcceleration

    ! Copy the position to module-scope.
    positionCartesianShared=positionCartesian
    ! Determine which mass type to use.
    if (present(massType)) then
       massTypeShared=massType
    else
       massTypeShared=massTypeAll
    end if
    ! Determine which component type to use.
    if (present(componentType)) then
       componentTypeShared=componentType
    else
       componentTypeShared=componentTypeAll
    end if
    ! Initialize pointer to function that supplies the acceleration for all components.
    componentAccelerationFunction => Component_Acceleration
    ! Determine which component type to use.
    if (present(componentType)) then
       componentTypeShared=componentType
    else
       componentTypeShared=componentTypeAll
    end if
    ! Determine which mass type to use.
    if (present(massType)) then
       massTypeShared     =massType
    else
       massTypeShared     =massTypeAll
    end if
    ! Compute the acceleration.
    Galactic_Structure_Acceleration=node%mapDouble1(componentAccelerationFunction,accelerationSize,reductionSummation,optimizeFor=optimizeForAccelerationSummation)
    !![
    <include directive="accelerationTask" type="functionCall" functionType="function" returnParameter="componentAcceleration">
     <functionArgs>node,positionCartesianShared,componentTypeShared,massTypeShared</functionArgs>
     <onReturn>Galactic_Structure_Acceleration=Galactic_Structure_Acceleration+componentAcceleration</onReturn>
    !!]
    include 'galactic_structure.acceleration.tasks.inc'
    !![
    </include>
    !!]
    return
  end function Galactic_Structure_Acceleration

  function Component_Acceleration(component,resultSize)
    !!{
    Function returning the acceleration in a component. Suitable for mapping over components.
    !!}
    use :: Galacticus_Nodes, only : nodeComponent
    implicit none
    integer                        , intent(in   )         :: resultSize
    class           (nodeComponent), intent(inout)         :: component
    double precision               , dimension(resultSize) :: Component_Acceleration

    Component_Acceleration=component%acceleration(positionCartesianShared,componentTypeShared,massTypeShared)
    return
  end function Component_Acceleration

end module Galactic_Structure_Accelerations












