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
  Implements a node operator class that performs evaporation of satellites due to dark matter self-interactions.
  !!}

  use :: Satellite_Evaporation_SIDM, only : satelliteEvaporationSIDMClass

  !![
  <nodeOperator name="nodeOperatorSatelliteEvaporationSIDM">
   <description>A node operator class that performs evaporation of satellites due to dark matter self-interactions.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorSatelliteEvaporationSIDM
     !!{
     A node operator class that performs evaporation of satellites due to dark matter self-interactions.
     !!}
     private
     class(satelliteEvaporationSIDMClass), pointer :: satelliteEvaporationSIDM_ => null()
   contains
     final     ::                          satelliteEvaporationSIDMDestructor
     procedure :: differentialEvolution => satelliteEvaporationSIDMDifferentialEvolution
  end type nodeOperatorSatelliteEvaporationSIDM
  
  interface nodeOperatorSatelliteEvaporationSIDM
     !!{
     Constructors for the \refClass{nodeOperatorSatelliteEvaporationSIDM} node operator class.
     !!}
     module procedure satelliteEvaporationSIDMConstructorParameters
     module procedure satelliteEvaporationSIDMConstructorInternal
  end interface nodeOperatorSatelliteEvaporationSIDM
  
contains

  function satelliteEvaporationSIDMConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorSatelliteEvaporationSIDM} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorSatelliteEvaporationSIDM)                :: self
    type (inputParameters                     ), intent(inout) :: parameters
    class(satelliteEvaporationSIDMClass       ), pointer       :: satelliteEvaporationSIDM_
    
    !![
    <objectBuilder class="satelliteEvaporationSIDM" name="satelliteEvaporationSIDM_" source="parameters"/>
    !!]
    self=nodeOperatorSatelliteEvaporationSIDM(satelliteEvaporationSIDM_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="satelliteEvaporationSIDM_"/>
    !!]
    return
  end function satelliteEvaporationSIDMConstructorParameters

  function satelliteEvaporationSIDMConstructorInternal(satelliteEvaporationSIDM_) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorSatelliteEvaporationSIDM} node operator class.
    !!}
    implicit none
    type (nodeOperatorSatelliteEvaporationSIDM)                        :: self
    class(satelliteEvaporationSIDMClass       ), intent(in   ), target :: satelliteEvaporationSIDM_
    !![
    <constructorAssign variables="*satelliteEvaporationSIDM_"/>
    !!]

    return
  end function satelliteEvaporationSIDMConstructorInternal

  subroutine satelliteEvaporationSIDMDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorSatelliteEvaporationSIDM} node operator class.
    !!}
    implicit none
    type(nodeOperatorSatelliteEvaporationSIDM), intent(inout) :: self

    !![
    <objectDestructor name="self%satelliteEvaporationSIDM_"/>
    !!]
    return
  end subroutine satelliteEvaporationSIDMDestructor
  
  subroutine satelliteEvaporationSIDMDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Perform evaporation of a satellite due to dark matter self-interactions.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSatellite
    implicit none
    class    (nodeOperatorSatelliteEvaporationSIDM), intent(inout), target  :: self
    type     (treeNode                            ), intent(inout), target  :: node
    logical                                        , intent(inout)          :: interrupt
    procedure(interruptTask                       ), intent(inout), pointer :: functionInterrupt
    integer                                        , intent(in   )          :: propertyType
    class    (nodeComponentSatellite              )               , pointer :: satellite
    !$GLC attributes unused :: interrupt, functionInterrupt, propertyType
    
    if (.not.node%isSatellite()) return
    satellite => node%satellite()
    call satellite%boundMassRate(self%satelliteEvaporationSIDM_%massLossRate(node))
    return
  end subroutine satelliteEvaporationSIDMDifferentialEvolution
  
