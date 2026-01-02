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
  Implements a node operator class that performs ram pressure mass loss in spheroids.
  !!}

  use :: Ram_Pressure_Stripping_Mass_Loss_Rate, only : ramPressureStrippingClass
  
  !![
  <nodeOperator name="nodeOperatorRamPressureMassLossSpheroids">
   <description>A node operator class that performs ram pressure mass loss in spheroids.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorRamPressureMassLossSpheroids
     !!{
     A node operator class that performs ram pressure mass loss in spheroids.
     !!}
     private
     class(ramPressureStrippingClass), pointer :: ramPressureStripping_ => null()
   contains
     final     ::                          ramPressureMassLossSpheroidsDestructor
     procedure :: differentialEvolution => ramPressureMassLossSpheroidsDifferentialEvolution
  end type nodeOperatorRamPressureMassLossSpheroids
  
  interface nodeOperatorRamPressureMassLossSpheroids
     !!{
     Constructors for the \refClass{nodeOperatorRamPressureMassLossSpheroids} node operator class.
     !!}
     module procedure ramPressureMassLossSpheroidsConstructorParameters
     module procedure ramPressureMassLossSpheroidsConstructorInternal
  end interface nodeOperatorRamPressureMassLossSpheroids
  
contains

  function ramPressureMassLossSpheroidsConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorRamPressureMassLossSpheroids} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorRamPressureMassLossSpheroids)                :: self
    type (inputParameters                         ), intent(inout) :: parameters
    class(ramPressureStrippingClass               ), pointer       :: ramPressureStripping_
    
    !![
    <objectBuilder class="ramPressureStripping" name="ramPressureStripping_" source="parameters"/>
    !!]
    self=nodeOperatorRamPressureMassLossSpheroids(ramPressureStripping_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="ramPressureStripping_"/>
    !!]
    return
  end function ramPressureMassLossSpheroidsConstructorParameters

  function ramPressureMassLossSpheroidsConstructorInternal(ramPressureStripping_) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorRamPressureMassLossSpheroids} node operator class.
    !!}
    implicit none
    type (nodeOperatorRamPressureMassLossSpheroids)                        :: self
    class(ramPressureStrippingClass               ), intent(in   ), target :: ramPressureStripping_
    !![
    <constructorAssign variables="*ramPressureStripping_"/>
    !!]

    return
  end function ramPressureMassLossSpheroidsConstructorInternal

  subroutine ramPressureMassLossSpheroidsDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorRamPressureMassLossSpheroids} node operator class.
    !!}
    implicit none
    type(nodeOperatorRamPressureMassLossSpheroids), intent(inout) :: self

    !![
    <objectDestructor name="self%ramPressureStripping_"/>
    !!]
    return
  end subroutine ramPressureMassLossSpheroidsDestructor
  
  subroutine ramPressureMassLossSpheroidsDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Perform ram pressure stripping-induced in a spheroid.
    !!}
    use :: Galacticus_Nodes    , only : propertyInactive, nodeComponentSpheroid, nodeComponentHotHalo
    use :: Abundances_Structure, only : operator(*)
    implicit none
    class           (nodeOperatorRamPressureMassLossSpheroids), intent(inout), target  :: self
    type            (treeNode                                ), intent(inout), target  :: node
    logical                                                   , intent(inout)          :: interrupt
    procedure       (interruptTask                           ), intent(inout), pointer :: functionInterrupt
    integer                                                   , intent(in   )          :: propertyType
    class           (nodeComponentSpheroid                   )               , pointer :: spheroid
    class           (nodeComponentHotHalo                    )               , pointer :: hotHalo
    double precision                                                                   :: massLossRate

    ! Do nothing during inactive property solving.
    if (propertyInactive(propertyType)) return
    spheroid => node%spheroid()
    if (spheroid%massGas() <= 0.0d0   ) return
    massLossRate=self%ramPressureStripping_%rateMassLoss(spheroid)
    if (massLossRate       <= 0.0d0   ) return
    hotHalo => node%hotHalo()
    call spheroid%                  massGasRate(-massLossRate                                                                       )
    call spheroid%          angularMomentumRate(-massLossRate*spheroid%angularMomentum()/(spheroid%massGas()+spheroid%massStellar()))
    call spheroid%            abundancesGasRate(-massLossRate*spheroid%abundancesGas  ()/ spheroid%massGas()                        )
    call hotHalo %           outflowingMassRate(+massLossRate                                                                       )
    call hotHalo %outflowingAngularMomentumRate(+massLossRate*spheroid%angularMomentum()/(spheroid%massGas()+spheroid%massStellar()))
    call hotHalo %     outflowingAbundancesRate(+massLossRate*spheroid%abundancesGas  ()/ spheroid%massGas()                        )
    return
  end subroutine ramPressureMassLossSpheroidsDifferentialEvolution

