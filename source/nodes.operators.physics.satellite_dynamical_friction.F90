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
  Implements a node operator class that applies dynamical friction to orbiting satellite halos.
  !!}

  use :: Satellite_Dynamical_Friction, only : satelliteDynamicalFrictionClass

  !![
  <nodeOperator name="nodeOperatorSatelliteDynamicalFriction">
   <description>A node operator class that applies dynamical friction to orbiting satellite halos.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorSatelliteDynamicalFriction
     !!{
     A node operator class that applies dynamical friction to orbiting satellite halos.
     !!}
     private
     class(satelliteDynamicalFrictionClass), pointer :: satelliteDynamicalFriction_ => null()
   contains
     final     ::                          satelliteDynamicalFrictionDestructor
     procedure :: differentialEvolution => satelliteDynamicalFrictionDifferentialEvolution
  end type nodeOperatorSatelliteDynamicalFriction
  
  interface nodeOperatorSatelliteDynamicalFriction
     !!{
     Constructors for the \refClass{nodeOperatorSatelliteDynamicalFriction} node operator class.
     !!}
     module procedure satelliteDynamicalFrictionConstructorParameters
     module procedure satelliteDynamicalFrictionConstructorInternal
  end interface nodeOperatorSatelliteDynamicalFriction
  
contains

  function satelliteDynamicalFrictionConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorSatelliteDynamicalFriction} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorSatelliteDynamicalFriction)                :: self
    type (inputParameters                       ), intent(inout) :: parameters
    class(satelliteDynamicalFrictionClass       ), pointer       :: satelliteDynamicalFriction_
    
    !![
    <objectBuilder class="satelliteDynamicalFriction" name="satelliteDynamicalFriction_" source="parameters"/>
    !!]
    self=nodeOperatorSatelliteDynamicalFriction(satelliteDynamicalFriction_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="satelliteDynamicalFriction_"/>
    !!]
    return
  end function satelliteDynamicalFrictionConstructorParameters

  function satelliteDynamicalFrictionConstructorInternal(satelliteDynamicalFriction_) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorSatelliteDynamicalFriction} node operator class.
    !!}
    implicit none
    type (nodeOperatorSatelliteDynamicalFriction)                        :: self
    class(satelliteDynamicalFrictionClass       ), intent(in   ), target :: satelliteDynamicalFriction_
    !![
    <constructorAssign variables="*satelliteDynamicalFriction_"/>
    !!]

    return
  end function satelliteDynamicalFrictionConstructorInternal

  subroutine satelliteDynamicalFrictionDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorSatelliteDynamicalFriction} node operator class.
    !!}
    implicit none
    type(nodeOperatorSatelliteDynamicalFriction), intent(inout) :: self

    !![
    <objectDestructor name="self%satelliteDynamicalFriction_"/>
    !!]
    return
  end subroutine satelliteDynamicalFrictionDestructor
  
  subroutine satelliteDynamicalFrictionDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Perform deceleration of a satellite due to dynamical friction.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSatellite
    implicit none
    class    (nodeOperatorSatelliteDynamicalFriction), intent(inout), target  :: self
    type     (treeNode                              ), intent(inout), target  :: node
    logical                                          , intent(inout)          :: interrupt
    procedure(interruptTask                         ), intent(inout), pointer :: functionInterrupt
    integer                                          , intent(in   )          :: propertyType
    class    (nodeComponentSatellite                )               , pointer :: satellite
    !$GLC attributes unused :: interrupt, functionInterrupt, propertyType

    if (.not.node%isSatellite()) return
    satellite => node%satellite()
    call satellite%velocityRate(self%satelliteDynamicalFriction_%acceleration(node))
    return
  end subroutine satelliteDynamicalFrictionDifferentialEvolution
  
