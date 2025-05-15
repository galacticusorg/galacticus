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
  Implements a node operator class that performs stellar feedback in disks.
  !!}

  use :: Star_Formation_Rates_Disks   , only : starFormationRateDisksClass
  use :: Stellar_Population_Properties, only : stellarPopulationPropertiesClass
  use :: Stellar_Feedback_Outflows    , only : stellarFeedbackOutflowsClass
  
  !![
  <nodeOperator name="nodeOperatorStellarFeedbackDisks">
   <description>A node operator class that performs stellar feedback in disks.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorStellarFeedbackDisks
     !!{
     A node operator class that performs stellar feedback in disks.
     !!}
     private
     class(starFormationRateDisksClass     ), pointer :: starFormationRateDisks_      => null()
     class(stellarPopulationPropertiesClass), pointer :: stellarPopulationProperties_ => null()
     class(stellarFeedbackOutflowsClass    ), pointer :: stellarFeedbackOutflows_     => null()
   contains
     final     ::                          stellarFeedbackDisksDestructor
     procedure :: differentialEvolution => stellarFeedbackDisksDifferentialEvolution
  end type nodeOperatorStellarFeedbackDisks
  
  interface nodeOperatorStellarFeedbackDisks
     !!{
     Constructors for the {\normalfont \ttfamily stellarFeedbackDisks} node operator class.
     !!}
     module procedure stellarFeedbackDisksConstructorParameters
     module procedure stellarFeedbackDisksConstructorInternal
  end interface nodeOperatorStellarFeedbackDisks
  
contains

  function stellarFeedbackDisksConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily starFormation} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorStellarFeedbackDisks)                :: self
    type (inputParameters                 ), intent(inout) :: parameters
    class(starFormationRateDisksClass     ), pointer       :: starFormationRateDisks_
    class(stellarPopulationPropertiesClass), pointer       :: stellarPopulationProperties_
    class(stellarFeedbackOutflowsClass    ), pointer       :: stellarFeedbackOutflows_
    
    !![
    <objectBuilder class="starFormationRateDisks"      name="starFormationRateDisks_"      source="parameters"/>
    <objectBuilder class="stellarPopulationProperties" name="stellarPopulationProperties_" source="parameters"/>
    <objectBuilder class="stellarFeedbackOutflows"     name="stellarFeedbackOutflows_"     source="parameters"/>
    !!]
    self=nodeOperatorStellarFeedbackDisks(starFormationRateDisks_,stellarPopulationProperties_,stellarFeedbackOutflows_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="starFormationRateDisks_"     />
    <objectDestructor name="stellarPopulationProperties_"/>
    <objectDestructor name="stellarFeedbackOutflows_"    />
    !!]
    return
  end function stellarFeedbackDisksConstructorParameters

  function stellarFeedbackDisksConstructorInternal(starFormationRateDisks_,stellarPopulationProperties_,stellarFeedbackOutflows_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily stellarFeedbackDisks} node operator class.
    !!}
    implicit none
    type (nodeOperatorStellarFeedbackDisks)                        :: self
    class(starFormationRateDisksClass     ), intent(in   ), target :: starFormationRateDisks_
    class(stellarPopulationPropertiesClass), intent(in   ), target :: stellarPopulationProperties_
    class(stellarFeedbackOutflowsClass    ), intent(in   ), target :: stellarFeedbackOutflows_
    !![
    <constructorAssign variables="*starFormationRateDisks_, *stellarPopulationProperties_, *stellarFeedbackOutflows_"/>
    !!]

    return
  end function stellarFeedbackDisksConstructorInternal

  subroutine stellarFeedbackDisksDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily stellarFeedbackDisks} node operator class.
    !!}
    implicit none
    type(nodeOperatorStellarFeedbackDisks), intent(inout) :: self

    !![
    <objectDestructor name="self%starFormationRateDisks_"     />
    <objectDestructor name="self%stellarPopulationProperties_"/>
    <objectDestructor name="self%stellarFeedbackOutflows_"    />
    !!]
    return
  end subroutine stellarFeedbackDisksDestructor
  
  subroutine stellarFeedbackDisksDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Perform feedback from stars in a disk.
    !!}
    use :: Abundances_Structure          , only : abundances         , zeroAbundances
    use :: Galacticus_Nodes              , only : propertyInactive   , nodeComponentDisk, nodeComponentHotHalo
    use :: Histories                     , only : history
    use :: Stellar_Luminosities_Structure, only : stellarLuminosities
    implicit none
    class           (nodeOperatorStellarFeedbackDisks), intent(inout), target  :: self
    type            (treeNode                        ), intent(inout), target  :: node
    logical                                           , intent(inout)          :: interrupt
    procedure       (interruptTask                   ), intent(inout), pointer :: functionInterrupt
    integer                                           , intent(in   )          :: propertyType
    class           (nodeComponentDisk               )               , pointer :: disk
    class           (nodeComponentHotHalo            )               , pointer :: hotHalo
    double precision                                                           :: rateStarFormation       , massGas                   , &
         &                                                                        angularMomentum         , fractionEjected           , &
         &                                                                        rateMassStellar         , rateEnergyInput           , &
         &                                                                        rateMassFuel            , rateMassOutflowEjected    , &
         &                                                                        rateMassOutflowExpelled , rateMassOutflowTotal      , &
         &                                                                        massDisk                , rateAngularMomentumOutflow
    type            (abundances                      )                         :: abundancesGas           , rateAbundancesFuels       , &
         &                                                                        rateAbundancesStellar   , rateAbundancesOutflow
    type            (history                         )                         :: ratePropertiesStellar
    type            (stellarLuminosities             )                         :: rateLuminositiesStellar

    ! Do nothing during inactive property solving.
    if (propertyInactive(propertyType)) return
    ! Check for a realistic disk, return immediately if disk is unphysical.
    disk => node%disk()
    if     (     disk%angularMomentum() < 0.0d0 &
         &  .or. disk%radius         () < 0.0d0 &
         &  .or. disk%massGas        () < 0.0d0 &
         & ) return
    ! Get the star formation rate.
    rateStarFormation=self%starFormationRateDisks_%rate(node)   
    ! Compute abundances of star forming gas.
    massGas      =disk%massGas      ()
    abundancesGas=disk%abundancesGas()
    call abundancesGas%massToMassFraction(massGas)
    ! Find rates of change of stellar mass, gas mass, abundances and luminosities.
    ratePropertiesStellar=disk%stellarPropertiesHistory()
    call self%stellarPopulationProperties_%rates(                                                      &
         &                                                                    rateStarFormation      , &
         &                                                                    abundancesGas          , &
         &                                                                    disk                   , &
         &                                                                    node                   , &
         &                                                                    ratePropertiesStellar  , &
         &                                                                    rateMassStellar        , &
         &                                                                    rateMassFuel           , &
         &                                                                    rateEnergyInput        , &
         &                                                                    rateAbundancesFuels    , &
         &                                                                    rateAbundancesStellar  , &
         &                                                                    rateLuminositiesStellar, &
         &                                       computeRateLuminosityStellar=.false.                  &
         &                                      )
    ! Find rate of outflow of material from the disk.
    call self%stellarFeedbackOutflows_%outflowRate(disk,rateStarFormation,rateEnergyInput,rateMassOutflowEjected,rateMassOutflowExpelled)
    rateMassOutflowTotal=+rateMassOutflowEjected  &
         &               +rateMassOutflowExpelled
    if (rateMassOutflowTotal > 0.0d0) then
       ! Find the fraction of material which outflows to the hot halo.
       fractionEjected=+rateMassOutflowEjected &
            &          /rateMassOutflowTotal
       ! Get the total mass of the disk.
       massDisk=+     massGas       &
            &   +disk%massStellar()
       ! Compute the angular momentum outflow rate.
       if (massDisk > 0.0d0) then
          angularMomentum           =disk%angularMomentum()
          rateAngularMomentumOutflow=angularMomentum*(rateMassOutflowTotal/massDisk)
       else
          rateAngularMomentumOutflow=0.0d0
       end if
       if (massGas  > 0.0d0) then
          rateAbundancesOutflow=abundancesGas*rateMassOutflowTotal
       else
          rateAbundancesOutflow=zeroAbundances
       end if
       hotHalo => node%hotHalo()
       call hotHalo%           outflowingMassRate( rateMassOutflowTotal      *fractionEjected)
       call disk   %                  massGasRate(-rateMassOutflowTotal                      )
       call hotHalo%outflowingAngularMomentumRate( rateAngularMomentumOutflow*fractionEjected)
       call disk   %          angularMomentumRate(-rateAngularMomentumOutflow                )
       call hotHalo%     outflowingAbundancesRate( rateAbundancesOutflow     *fractionEjected)
       call disk   %            abundancesGasRate(-rateAbundancesOutflow                     )
    end if
    return
  end subroutine stellarFeedbackDisksDifferentialEvolution

