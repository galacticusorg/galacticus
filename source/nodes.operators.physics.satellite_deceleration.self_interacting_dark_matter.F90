!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
  Implements a node operator class that performs deceleration of satellites due to dark matter self-interactions.
  !!}

  use :: Satellite_Deceleration_SIDM, only : satelliteDecelerationSIDMClass

  !![
  <nodeOperator name="nodeOperatorSatelliteDecelerationSIDM">
   <description>A node operator class that performs deceleration of satellites due to dark matter self-interactions.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorSatelliteDecelerationSIDM
     !!{
     A node operator class that performs deceleration of satellites due to dark matter self-interactions.
     !!}
     private
     class(satelliteDecelerationSIDMClass), pointer :: satelliteDecelerationSIDM_ => null()
   contains
     final     ::                          satelliteDecelerationSIDMDestructor
     procedure :: differentialEvolution => satelliteDecelerationSIDMDifferentialEvolution
  end type nodeOperatorSatelliteDecelerationSIDM
  
  interface nodeOperatorSatelliteDecelerationSIDM
     !!{
     Constructors for the \refClass{nodeOperatorSatelliteDecelerationSIDM} node operator class.
     !!}
     module procedure satelliteDecelerationSIDMConstructorParameters
     module procedure satelliteDecelerationSIDMConstructorInternal
  end interface nodeOperatorSatelliteDecelerationSIDM
  
contains

  function satelliteDecelerationSIDMConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorSatelliteDecelerationSIDM} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorSatelliteDecelerationSIDM)                :: self
    type (inputParameters                      ), intent(inout) :: parameters
    class(satelliteDecelerationSIDMClass       ), pointer       :: satelliteDecelerationSIDM_
    
    !![
    <objectBuilder class="satelliteDecelerationSIDM" name="satelliteDecelerationSIDM_" source="parameters"/>
    !!]
    self=nodeOperatorSatelliteDecelerationSIDM(satelliteDecelerationSIDM_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="satelliteDecelerationSIDM_"/>
    !!]
    return
  end function satelliteDecelerationSIDMConstructorParameters

  function satelliteDecelerationSIDMConstructorInternal(satelliteDecelerationSIDM_) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorSatelliteDecelerationSIDM} node operator class.
    !!}
    implicit none
    type (nodeOperatorSatelliteDecelerationSIDM)                        :: self
    class(satelliteDecelerationSIDMClass       ), intent(in   ), target :: satelliteDecelerationSIDM_
    !![
    <constructorAssign variables="*satelliteDecelerationSIDM_"/>
    !!]

    return
  end function satelliteDecelerationSIDMConstructorInternal

  subroutine satelliteDecelerationSIDMDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorSatelliteDecelerationSIDM} node operator class.
    !!}
    implicit none
    type(nodeOperatorSatelliteDecelerationSIDM), intent(inout) :: self

    !![
    <objectDestructor name="self%satelliteDecelerationSIDM_"/>
    !!]
    return
  end subroutine satelliteDecelerationSIDMDestructor
  
  subroutine satelliteDecelerationSIDMDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Perform deceleration of a satellite due to dark matter self-interactions.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSatellite
    implicit none
    class    (nodeOperatorSatelliteDecelerationSIDM), intent(inout), target  :: self
    type     (treeNode                             ), intent(inout), target  :: node
    logical                                         , intent(inout)          :: interrupt
    procedure(interruptTask                        ), intent(inout), pointer :: functionInterrupt
    integer                                         , intent(in   )          :: propertyType
    class    (nodeComponentSatellite               )               , pointer :: satellite
    !$GLC attributes unused :: interrupt, functionInterrupt, propertyType

    if (.not.node%isSatellite()) return
    satellite => node%satellite()
    call satellite%velocityRate(self%satelliteDecelerationSIDM_%acceleration(node))
    return
  end subroutine satelliteDecelerationSIDMDifferentialEvolution
  
