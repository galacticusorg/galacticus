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
  Implements a node operator class that performs tidal mass loss in disks.
  !!}

  use :: Tidal_Stripping_Mass_Loss_Rate, only : tidalStrippingClass
  
  !![
  <nodeOperator name="nodeOperatorTidalMassLossDisks">
   <description>A node operator class that performs tidal mass loss in disks.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorTidalMassLossDisks
     !!{
     A node operator class that performs tidal mass loss in disks.
     !!}
     private
     class(tidalStrippingClass), pointer :: tidalStripping_ => null()
   contains
     final     ::                          tidalMassLossDisksDestructor
     procedure :: differentialEvolution => tidalMassLossDisksDifferentialEvolution
  end type nodeOperatorTidalMassLossDisks
  
  interface nodeOperatorTidalMassLossDisks
     !!{
     Constructors for the \refClass{nodeOperatorTidalMassLossDisks} node operator class.
     !!}
     module procedure tidalMassLossDisksConstructorParameters
     module procedure tidalMassLossDisksConstructorInternal
  end interface nodeOperatorTidalMassLossDisks
  
contains

  function tidalMassLossDisksConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorTidalMassLossDisks} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorTidalMassLossDisks)                :: self
    type (inputParameters               ), intent(inout) :: parameters
    class(tidalStrippingClass           ), pointer       :: tidalStripping_
    
    !![
    <objectBuilder class="tidalStripping" name="tidalStripping_" source="parameters"/>
    !!]
    self=nodeOperatorTidalMassLossDisks(tidalStripping_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="tidalStripping_"/>
    !!]
    return
  end function tidalMassLossDisksConstructorParameters

  function tidalMassLossDisksConstructorInternal(tidalStripping_) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorTidalMassLossDisks} node operator class.
    !!}
    implicit none
    type (nodeOperatorTidalMassLossDisks)                        :: self
    class(tidalStrippingClass           ), intent(in   ), target :: tidalStripping_
    !![
    <constructorAssign variables="*tidalStripping_"/>
    !!]

    return
  end function tidalMassLossDisksConstructorInternal

  subroutine tidalMassLossDisksDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorTidalMassLossDisks} node operator class.
    !!}
    implicit none
    type(nodeOperatorTidalMassLossDisks), intent(inout) :: self

    !![
    <objectDestructor name="self%tidalStripping_"/>
    !!]
    return
  end subroutine tidalMassLossDisksDestructor
  
  subroutine tidalMassLossDisksDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Perform star formation in a disk.
    !!}
    use :: Galacticus_Nodes              , only : propertyInactive, nodeComponentDisk  , nodeComponentHotHalo
    use :: Abundances_Structure          , only : operator(*)
    use :: Histories                     , only : operator(*)     , history
    use :: Stellar_Luminosities_Structure, only : operator(*)     , stellarLuminosities, zeroStellarLuminosities, max
    implicit none
    class           (nodeOperatorTidalMassLossDisks), intent(inout), target  :: self
    type            (treeNode                      ), intent(inout), target  :: node
    logical                                         , intent(inout)          :: interrupt
    procedure       (interruptTask                 ), intent(inout), pointer :: functionInterrupt
    integer                                         , intent(in   )          :: propertyType
    class           (nodeComponentDisk             )               , pointer :: disk
    class           (nodeComponentHotHalo          )               , pointer :: hotHalo
    type            (stellarLuminosities           ), save                   :: luminositiesTransferRate
    !$omp threadprivate(luminositiesTransferRate)
    double precision                                                         :: fractionGas             , fractionStellar, &
         &                                                                      massLossRate
    type            (history                       )                         :: historyTransferRate

    ! Do nothing during inactive property solving.
    if (propertyInactive(propertyType)            ) return
    ! Return if the disk has no mass.
    disk => node%disk()
    if (disk%massGas()+disk%massStellar() <= 0.0d0) return
    ! Return if the tidal mass loss rate is zero.
    massLossRate=self%tidalStripping_%rateMassLoss(disk)
    if (massLossRate                      <= 0.0d0) return
    ! Transfer stripped material from the disk.
    !! Gas is moved to the hot halo component.
    hotHalo         => node%hotHalo()
    fractionGas     =  min(1.0d0,max(0.0d0,disk%massGas()/(disk%massGas()+disk%massStellar())))
    fractionStellar =  1.0d0-fractionGas
    if (fractionGas     > 0.0d0 .and. disk%massGas    () > 0.0d0) then
       call        disk%                  massGasRate(-fractionGas    *massLossRate                                                             )
       call        disk%            abundancesGasRate(-fractionGas    *massLossRate*disk%abundancesGas    ()/ disk%massGas()                    )
       call     hotHalo%           outflowingMassRate(+fractionGas    *massLossRate                                                             )
       call     hotHalo%outflowingAbundancesRate     (+fractionGas    *massLossRate*disk%abundancesGas    ()/ disk%massGas()                    )
       call     hotHalo%outflowingAngularMomentumRate(+fractionGas    *massLossRate*disk%angularMomentum  ()/(disk%massGas()+disk%massStellar()))
    end if
    ! Stellar mass is simply removed.
    if (fractionStellar > 0.0d0 .and. disk%massStellar() > 0.0d0) then
       ! If luminosities are being treated as inactive properties this is an error - they appear on the right-hand side
       ! of the following ODE terms so are not inactive. (An approach similar to what is used for transfer of
       ! luminosities to the spheroid by bar instabilities could work here.)
       !! Stellar mass and metals.
       call        disk%              massStellarRate(-fractionStellar*massLossRate                                                             )
       call        disk%        abundancesStellarRate(-fractionStellar*massLossRate*disk%abundancesStellar()/                disk%massStellar() )
       !! Stellar luminosities.
       luminositiesTransferRate=max(zeroStellarLuminosities,disk%luminositiesStellar())
       call        disk%      luminositiesStellarRate(-fractionStellar*massLossRate*luminositiesTransferRate/                disk%massStellar() )
       !! Stellar properties history.
       historyTransferRate=disk%stellarPropertiesHistory()
       if (historyTransferRate%exists()) &
            & call disk%stellarPropertiesHistoryRate (-fractionStellar*massLossRate*historyTransferRate     /                disk%massStellar() )
       call historyTransferRate%destroy()
       !! Star formation history.
       historyTransferRate=disk%starFormationHistory()
       if (historyTransferRate%exists()) &
            & call disk    %starFormationHistoryRate (-fractionStellar*massLossRate*historyTransferRate     /                disk%massStellar() )
    end if
    ! Angular momentum is lost.
    call           disk%          angularMomentumRate(-                massLossRate*disk%angularMomentum  ()/(disk%massGas()+disk%massStellar()))
    return
  end subroutine tidalMassLossDisksDifferentialEvolution

