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
  Implements a node operator class that performs ram pressure mass loss in disks.
  !!}

  use :: Ram_Pressure_Stripping_Mass_Loss_Rate, only : ramPressureStrippingClass
  
  !![
  <nodeOperator name="nodeOperatorRamPressureMassLossDisks">
   <description>A node operator class that performs ram pressure mass loss in disks.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorRamPressureMassLossDisks
     !!{
     A node operator class that performs ram pressure mass loss in disks.
     !!}
     private
     class(ramPressureStrippingClass), pointer :: ramPressureStripping_ => null()
   contains
     final     ::                          ramPressureMassLossDisksDestructor
     procedure :: differentialEvolution => ramPressureMassLossDisksDifferentialEvolution
  end type nodeOperatorRamPressureMassLossDisks
  
  interface nodeOperatorRamPressureMassLossDisks
     !!{
     Constructors for the \refClass{nodeOperatorRamPressureMassLossDisks} node operator class.
     !!}
     module procedure ramPressureMassLossDisksConstructorParameters
     module procedure ramPressureMassLossDisksConstructorInternal
  end interface nodeOperatorRamPressureMassLossDisks
  
contains

  function ramPressureMassLossDisksConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorRamPressureMassLossDisks} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorRamPressureMassLossDisks)                :: self
    type (inputParameters                     ), intent(inout) :: parameters
    class(ramPressureStrippingClass           ), pointer       :: ramPressureStripping_
    
    !![
    <objectBuilder class="ramPressureStripping" name="ramPressureStripping_" source="parameters"/>
    !!]
    self=nodeOperatorRamPressureMassLossDisks(ramPressureStripping_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="ramPressureStripping_"/>
    !!]
    return
  end function ramPressureMassLossDisksConstructorParameters

  function ramPressureMassLossDisksConstructorInternal(ramPressureStripping_) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorRamPressureMassLossDisks} node operator class.
    !!}
    implicit none
    type (nodeOperatorRamPressureMassLossDisks)                        :: self
    class(ramPressureStrippingClass           ), intent(in   ), target :: ramPressureStripping_
    !![
    <constructorAssign variables="*ramPressureStripping_"/>
    !!]

    return
  end function ramPressureMassLossDisksConstructorInternal

  subroutine ramPressureMassLossDisksDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorRamPressureMassLossDisks} node operator class.
    !!}
    implicit none
    type(nodeOperatorRamPressureMassLossDisks), intent(inout) :: self

    !![
    <objectDestructor name="self%ramPressureStripping_"/>
    !!]
    return
  end subroutine ramPressureMassLossDisksDestructor
  
  subroutine ramPressureMassLossDisksDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Perform ram pressure stripping-induced mass loss in a disk.
    !!}
    use :: Galacticus_Nodes    , only : propertyInactive, nodeComponentDisk, nodeComponentHotHalo
    use :: Abundances_Structure, only : operator(*)
    implicit none
    class           (nodeOperatorRamPressureMassLossDisks), intent(inout), target  :: self
    type            (treeNode                            ), intent(inout), target  :: node
    logical                                               , intent(inout)          :: interrupt
    procedure       (interruptTask                       ), intent(inout), pointer :: functionInterrupt
    integer                                               , intent(in   )          :: propertyType
    class           (nodeComponentDisk                   )               , pointer :: disk
    class           (nodeComponentHotHalo                )               , pointer :: hotHalo
    double precision                                                               :: massLossRate

    ! Do nothing during inactive property solving.
    if (propertyInactive(propertyType)        ) return
    disk => node%disk()
    if (disk%massGas() <= 0.0d0               ) return
    massLossRate=self%ramPressureStripping_%rateMassLoss(disk)
    if (massLossRate   <= 0.0d0               ) return
    hotHalo => node%hotHalo()
    call    disk%                  massGasRate(-massLossRate                                                           )
    call    disk%          angularMomentumRate(-massLossRate*disk%angularMomentum()/(disk%massGas()+disk%massStellar()))
    call    disk%            abundancesGasRate(-massLossRate*disk%abundancesGas  ()/ disk%massGas()                    )
    call hotHalo%           outflowingMassRate(+massLossRate                                                           )
    call hotHalo%outflowingAngularMomentumRate(+massLossRate*disk%angularMomentum()/(disk%massGas()+disk%massStellar()))
    call hotHalo%outflowingAbundancesRate     (+massLossRate*disk%abundancesGas  ()/ disk%massGas()                    )
    return
  end subroutine ramPressureMassLossDisksDifferentialEvolution

