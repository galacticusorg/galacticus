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
  Implements a node operator class that performs tidal mass loss in spheroids.
  !!}

  use :: Tidal_Stripping_Mass_Loss_Rate, only : tidalStrippingClass
  
  !![
  <nodeOperator name="nodeOperatorTidalMassLossSpheroids">
   <description>A node operator class that performs tidal mass loss in spheroids.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorTidalMassLossSpheroids
     !!{
     A node operator class that performs tidal mass loss in spheroids.
     !!}
     private
     class(tidalStrippingClass), pointer :: tidalStripping_ => null()
   contains
     final     ::                          tidalMassLossSpheroidsDestructor
     procedure :: differentialEvolution => tidalMassLossSpheroidsDifferentialEvolution
  end type nodeOperatorTidalMassLossSpheroids
  
  interface nodeOperatorTidalMassLossSpheroids
     !!{
     Constructors for the \refClass{nodeOperatorTidalMassLossSpheroids} node operator class.
     !!}
     module procedure tidalMassLossSpheroidsConstructorParameters
     module procedure tidalMassLossSpheroidsConstructorInternal
  end interface nodeOperatorTidalMassLossSpheroids
  
contains

  function tidalMassLossSpheroidsConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorTidalMassLossSpheroids} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorTidalMassLossSpheroids)                :: self
    type (inputParameters                   ), intent(inout) :: parameters
    class(tidalStrippingClass               ), pointer       :: tidalStripping_
    
    !![
    <objectBuilder class="tidalStripping" name="tidalStripping_" source="parameters"/>
    !!]
    self=nodeOperatorTidalMassLossSpheroids(tidalStripping_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="tidalStripping_"/>
    !!]
    return
  end function tidalMassLossSpheroidsConstructorParameters

  function tidalMassLossSpheroidsConstructorInternal(tidalStripping_) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorTidalMassLossSpheroids} node operator class.
    !!}
    implicit none
    type (nodeOperatorTidalMassLossSpheroids)                        :: self
    class(tidalStrippingClass               ), intent(in   ), target :: tidalStripping_
    !![
    <constructorAssign variables="*tidalStripping_"/>
    !!]

    return
  end function tidalMassLossSpheroidsConstructorInternal

  subroutine tidalMassLossSpheroidsDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorTidalMassLossSpheroids} node operator class.
    !!}
    implicit none
    type(nodeOperatorTidalMassLossSpheroids), intent(inout) :: self

    !![
    <objectDestructor name="self%tidalStripping_"/>
    !!]
    return
  end subroutine tidalMassLossSpheroidsDestructor
  
  subroutine tidalMassLossSpheroidsDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Perform star formation in a spheroid.
    !!}
    use :: Galacticus_Nodes              , only : propertyInactive, nodeComponentSpheroid  , nodeComponentHotHalo
    use :: Abundances_Structure          , only : operator(*)
    use :: Histories                     , only : operator(*)     , history
    use :: Stellar_Luminosities_Structure, only : operator(*)     , stellarLuminosities, zeroStellarLuminosities, max
    implicit none
    class           (nodeOperatorTidalMassLossSpheroids), intent(inout), target  :: self
    type            (treeNode                          ), intent(inout), target  :: node
    logical                                             , intent(inout)          :: interrupt
    procedure       (interruptTask                     ), intent(inout), pointer :: functionInterrupt
    integer                                             , intent(in   )          :: propertyType
    class           (nodeComponentSpheroid             )               , pointer :: spheroid
    class           (nodeComponentHotHalo              )               , pointer :: hotHalo
    type            (stellarLuminosities               ), save                   :: luminositiesTransferRate
    !$omp threadprivate(luminositiesTransferRate)
    double precision                                                             :: fractionGas             , fractionStellar, &
         &                                                                          massLossRate
    type            (history                           )                         :: historyTransferRate

    ! Do nothing during inactive property solving.
    if (propertyInactive(propertyType)                    ) return
    ! Return if the spheroid has no mass.
    spheroid => node%spheroid()
    if (spheroid%massGas()+spheroid%massStellar() <= 0.0d0) return
    ! Return if the tidal mass loss rate is zero.
    massLossRate=self%tidalStripping_%rateMassLoss(spheroid)
    if (massLossRate                              <= 0.0d0) return
    ! Transfer stripped material from the spheroid.
    !! Gas is moved to the hot halo component.
    hotHalo         => node%hotHalo()
    fractionGas     =  min(1.0d0,max(0.0d0,spheroid%massGas()/(spheroid%massGas()+spheroid%massStellar())))
    fractionStellar =  1.0d0-fractionGas
    if (fractionGas     > 0.0d0 .and. spheroid%massGas    () > 0.0d0) then
       call        spheroid%                  massGasRate(-fractionGas    *massLossRate                                                                         )
       call        spheroid%            abundancesGasRate(-fractionGas    *massLossRate*spheroid%abundancesGas    ()/ spheroid%massGas()                        )
       call         hotHalo%           outflowingMassRate(+fractionGas    *massLossRate                                                                         )
       call         hotHalo%outflowingAbundancesRate     (+fractionGas    *massLossRate*spheroid%abundancesGas    ()/ spheroid%massGas()                        )
       call         hotHalo%outflowingAngularMomentumRate(+fractionGas    *massLossRate*spheroid%angularMomentum  ()/(spheroid%massGas()+spheroid%massStellar()))
    end if
    ! Stellar mass is simply removed.
    if (fractionStellar > 0.0d0 .and. spheroid%massStellar() > 0.0d0) then
       ! If luminosities are being treated as inactive properties this is an error - they appear on the right-hand side
       ! of the following ODE terms so are not inactive. (An approach similar to what is used for transfer of
       ! luminosities to the spheroid by bar instabilities could work here.)
       !! Stellar mass and metals.
       call        spheroid%              massStellarRate(-fractionStellar*massLossRate                                                                         )
       call        spheroid%        abundancesStellarRate(-fractionStellar*massLossRate*spheroid%abundancesStellar()/                    spheroid%massStellar() )
       !! Stellar luminosities.
       luminositiesTransferRate=max(zeroStellarLuminosities,spheroid%luminositiesStellar())
       call        spheroid%      luminositiesStellarRate(-fractionStellar*massLossRate*luminositiesTransferRate    /                    spheroid%massStellar() )
       !! Stellar properties history.
       historyTransferRate=spheroid%stellarPropertiesHistory()
       if (historyTransferRate%exists()) &
            & call spheroid%stellarPropertiesHistoryRate (-fractionStellar*massLossRate*historyTransferRate         /                    spheroid%massStellar() )
       call historyTransferRate%destroy()
       !! Star formation history.
       historyTransferRate=spheroid%starFormationHistory()
       if (historyTransferRate%exists()) &
            & call spheroid    %starFormationHistoryRate (-fractionStellar*massLossRate*historyTransferRate         /                    spheroid%massStellar() )
    end if
    ! Angular momentum is lost.
    call           spheroid%          angularMomentumRate(-                massLossRate*spheroid%angularMomentum  ()/(spheroid%massGas()+spheroid%massStellar()))
    return
  end subroutine tidalMassLossSpheroidsDifferentialEvolution

