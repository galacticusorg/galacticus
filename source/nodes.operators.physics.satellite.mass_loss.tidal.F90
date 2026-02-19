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
  Implements a node operator class that applies tidal mass loss to orbiting satellite halos.
  !!}

  use :: Satellite_Tidal_Stripping, only : satelliteTidalStrippingClass

  !![
  <nodeOperator name="nodeOperatorSatelliteTidalMassLoss">
   <description>A node operator class that applies tidal mass loss to orbiting satellite halos.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorSatelliteTidalMassLoss
     !!{
     A node operator class that applies tidal mass loss to orbiting satellite halos.
     !!}
     private
     class  (satelliteTidalStrippingClass), pointer :: satelliteTidalStripping_ => null()
     logical                                        :: applyPreInfall
   contains
     final     ::                          satelliteTidalStrippingDestructor
     procedure :: differentialEvolution => satelliteTidalStrippingDifferentialEvolution
  end type nodeOperatorSatelliteTidalMassLoss
  
  interface nodeOperatorSatelliteTidalMassLoss
     !!{
     Constructors for the \refClass{nodeOperatorSatelliteTidalMassLoss} node operator class.
     !!}
     module procedure satelliteTidalStrippingConstructorParameters
     module procedure satelliteTidalStrippingConstructorInternal
  end interface nodeOperatorSatelliteTidalMassLoss
  
contains

  function satelliteTidalStrippingConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorSatelliteTidalMassLoss} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (nodeOperatorSatelliteTidalMassLoss)                :: self
    type   (inputParameters                   ), intent(inout) :: parameters
    class  (satelliteTidalStrippingClass      ), pointer       :: satelliteTidalStripping_
    logical                                                    :: applyPreInfall
    
    !![
     <inputParameter>
      <name>applyPreInfall</name>
      <defaultValue>.false.</defaultValue>
      <description>If true, tidal mass loss is applied pre-infall.</description>
      <source>parameters</source>
    </inputParameter>
   <objectBuilder class="satelliteTidalStripping" name="satelliteTidalStripping_" source="parameters"/>
    !!]
    self=nodeOperatorSatelliteTidalMassLoss(applyPreInfall,satelliteTidalStripping_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="satelliteTidalStripping_"/>
    !!]
    return
  end function satelliteTidalStrippingConstructorParameters

  function satelliteTidalStrippingConstructorInternal(applyPreInfall,satelliteTidalStripping_) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorSatelliteTidalMassLoss} node operator class.
    !!}
    implicit none
    type   (nodeOperatorSatelliteTidalMassLoss)                        :: self
    class  (satelliteTidalStrippingClass      ), intent(in   ), target :: satelliteTidalStripping_
    logical                                    , intent(in   )         :: applyPreInfall
    !![ 
    <constructorAssign variables="applyPreInfall, *satelliteTidalStripping_"/>
    !!]

    return
  end function satelliteTidalStrippingConstructorInternal

  subroutine satelliteTidalStrippingDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorSatelliteTidalMassLoss} node operator class.
    !!}
    implicit none
    type(nodeOperatorSatelliteTidalMassLoss), intent(inout) :: self

    !![
    <objectDestructor name="self%satelliteTidalStripping_"/>
    !!]
    return
  end subroutine satelliteTidalStrippingDestructor
  
  subroutine satelliteTidalStrippingDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Perform mass loss from a satellite due to tidal stripping.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSatellite
    implicit none
    class    (nodeOperatorSatelliteTidalMassLoss), intent(inout), target  :: self
    type     (treeNode                          ), intent(inout), target  :: node
    logical                                      , intent(inout)          :: interrupt
    procedure(interruptTask                     ), intent(inout), pointer :: functionInterrupt
    integer                                      , intent(in   )          :: propertyType
    class    (nodeComponentSatellite            )               , pointer :: satellite
    !$GLC attributes unused :: interrupt, functionInterrupt, propertyType
    
    if     (                                   &
         &          node%isOnMainBranch     () &
         &  .or.                               &
         &   (                                 &
         &     .not.self%applyPreInfall        &
         &    .and.                            &
         &     .not.node%isSatellite        () &
         &   )                                 &
         & ) return
    satellite => node%satellite()
    call satellite%boundMassRate(self%satelliteTidalStripping_%massLossRate(node))
    return
  end subroutine satelliteTidalStrippingDifferentialEvolution
  
