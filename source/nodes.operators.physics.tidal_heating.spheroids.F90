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
  Implements a node operator class that performs tidal heating in spheroids.
  !!}
  
  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass
  use :: Satellites_Tidal_Fields, only : satelliteTidalFieldClass

  !![
  <nodeOperator name="nodeOperatorTidalHeatingSpheroids">
   <description>A node operator class that performs tidal heating in spheroids.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorTidalHeatingSpheroids
     !!{
     A node operator class that performs tidal heating in spheroids.
     !!}
     private
     class(satelliteTidalFieldClass), pointer :: satelliteTidalField_ => null()
     class(darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_ => null()
   contains
     final     ::                          tidalHeatingSpheroidsDestructor
     procedure :: differentialEvolution => tidalHeatingSpheroidsDifferentialEvolution
  end type nodeOperatorTidalHeatingSpheroids

  interface nodeOperatorTidalHeatingSpheroids
     !!{
     Constructors for the \refClass{nodeOperatorTidalHeatingSpheroids} node operator class.
     !!}
     module procedure tidalHeatingSpheroidsConstructorParameters
     module procedure tidalHeatingSpheroidsConstructorInternal
  end interface nodeOperatorTidalHeatingSpheroids

contains
  
  function tidalHeatingSpheroidsConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorTidalHeatingSpheroids} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorTidalHeatingSpheroids)                :: self
    type (inputParameters                  ), intent(inout) :: parameters
    class(darkMatterHaloScaleClass         ), pointer       :: darkMatterHaloScale_
    class(satelliteTidalFieldClass         ), pointer       :: satelliteTidalField_

    !![
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    <objectBuilder class="satelliteTidalField" name="satelliteTidalField_" source="parameters"/>
    !!]
    self=nodeOperatorTidalHeatingSpheroids(darkMatterHaloScale_,satelliteTidalField_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    <objectDestructor name="satelliteTidalField_"/>
    !!]
    return
  end function tidalHeatingSpheroidsConstructorParameters

  function tidalHeatingSpheroidsConstructorInternal(darkMatterHaloScale_,satelliteTidalField_) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorTidalHeatingSpheroids} node operator class.
    !!}
    implicit none
    type (nodeOperatorTidalHeatingSpheroids)                        :: self
    class(darkMatterHaloScaleClass         ), intent(in   ), target :: darkMatterHaloScale_
    class(satelliteTidalFieldClass         ), intent(in   ), target :: satelliteTidalField_
    !![
    <constructorAssign variables="*darkMatterHaloScale_, *satelliteTidalField_"/>
    !!]

    return
  end function tidalHeatingSpheroidsConstructorInternal

  subroutine tidalHeatingSpheroidsDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorTidalHeatingSpheroids} node operator class.
    !!}
    implicit none
    type(nodeOperatorTidalHeatingSpheroids), intent(inout) :: self
    
    !![
    <objectDestructor name="self%darkMatterHaloScale_"/>
    <objectDestructor name="self%satelliteTidalField_"/>
    !!]
    return
  end subroutine tidalHeatingSpheroidsDestructor
  
  subroutine tidalHeatingSpheroidsDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Perform star formation in a spheroid.
    !!}
    use :: Galacticus_Nodes, only : propertyInactive, nodeComponentSpheroid
    implicit none
    class           (nodeOperatorTidalHeatingSpheroids), intent(inout), target  :: self
    type            (treeNode                         ), intent(inout), target  :: node
    logical                                            , intent(inout)          :: interrupt
    procedure       (interruptTask                    ), intent(inout), pointer :: functionInterrupt
    integer                                            , intent(in   )          :: propertyType
    class           (nodeComponentSpheroid            )               , pointer :: spheroid
    double precision                                                            :: tidalField       , tidalTorque
    
    ! Do nothing during inactive property solving.
    if     (propertyInactive(propertyType)      ) return
    ! Return if the node is not a satellite.
    if     (.not.node%isSatellite()             ) return
    ! Return if the spheroid is unphysical.
    spheroid => node%spheroid()
    if     (    spheroid%angularMomentum() < 0.0d0 &
         & .or. spheroid%radius         () < 0.0d0 &
         & .or. spheroid%massGas        () < 0.0d0 &
         & .or. spheroid%massStellar    () < 0.0d0 &
         & ) return
    ! Return if the spheroid has excessive angular momentum or radius.
    if     (                                                                                                                                                                      &
         &   spheroid%angularMomentum() > (spheroid%massGas()+spheroid%massStellar())*self%darkMatterHaloScale_%radiusVirial(node)*self%darkMatterHaloScale_%velocityVirial(node) &
         &  .or.                                                                                                                                                                  &
         &   spheroid%radius         () >                                             self%darkMatterHaloScale_%radiusVirial(node)                                                &
         & ) return
    tidalField =self%satelliteTidalField_%tidalTensorRadial(node)
    tidalTorque=+abs(tidalField)              &
         &      *(                           &
         &        +spheroid%massGas    ()    &
         &        +spheroid%massStellar()    &
         &       )                           &
         &      *  spheroid%radius     ()**2
    call spheroid%angularMomentumRate(tidalTorque)
    return
  end subroutine tidalHeatingSpheroidsDifferentialEvolution

