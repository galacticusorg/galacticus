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
  Implements a node operator class that performs stellar feedback in spheroids.
  !!}

  use :: Star_Formation_Rates_Spheroids, only : starFormationRateSpheroidsClass
  use :: Stellar_Population_Properties , only : stellarPopulationPropertiesClass
  use :: Stellar_Feedback_Outflows     , only : stellarFeedbackOutflowsClass
  
  !![
  <nodeOperator name="nodeOperatorStellarFeedbackSpheroids">
   <description>A node operator class that performs stellar feedback in spheroids.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorStellarFeedbackSpheroids
     !!{
     A node operator class that performs stellar feedback in spheroids.
     !!}
     private
     class(starFormationRateSpheroidsClass ), pointer :: starFormationRateSpheroids_  => null()
     class(stellarPopulationPropertiesClass), pointer :: stellarPopulationProperties_ => null()
     class(stellarFeedbackOutflowsClass    ), pointer :: stellarFeedbackOutflows_     => null()
   contains
     final     ::                          stellarFeedbackSpheroidsDestructor
     procedure :: differentialEvolution => stellarFeedbackSpheroidsDifferentialEvolution
  end type nodeOperatorStellarFeedbackSpheroids
  
  interface nodeOperatorStellarFeedbackSpheroids
     !!{
     Constructors for the {\normalfont \ttfamily stellarFeedbackSpheroids} node operator class.
     !!}
     module procedure stellarFeedbackSpheroidsConstructorParameters
     module procedure stellarFeedbackSpheroidsConstructorInternal
  end interface nodeOperatorStellarFeedbackSpheroids
  
contains

  function stellarFeedbackSpheroidsConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily starFormation} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorStellarFeedbackSpheroids)                :: self
    type (inputParameters                     ), intent(inout) :: parameters
    class(starFormationRateSpheroidsClass     ), pointer       :: starFormationRateSpheroids_
    class(stellarPopulationPropertiesClass    ), pointer       :: stellarPopulationProperties_
    class(stellarFeedbackOutflowsClass        ), pointer       :: stellarFeedbackOutflows_
    
    !![
    <objectBuilder class="starFormationRateSpheroids"  name="starFormationRateSpheroids_"  source="parameters"/>
    <objectBuilder class="stellarPopulationProperties" name="stellarPopulationProperties_" source="parameters"/>
    <objectBuilder class="stellarFeedbackOutflows"     name="stellarFeedbackOutflows_"     source="parameters"/>
    !!]
    self=nodeOperatorStellarFeedbackSpheroids(starFormationRateSpheroids_,stellarPopulationProperties_,stellarFeedbackOutflows_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="starFormationRateSpheroids_" />
    <objectDestructor name="stellarPopulationProperties_"/>
    <objectDestructor name="stellarFeedbackOutflows_"    />
    !!]
    return
  end function stellarFeedbackSpheroidsConstructorParameters

  function stellarFeedbackSpheroidsConstructorInternal(starFormationRateSpheroids_,stellarPopulationProperties_,stellarFeedbackOutflows_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily stellarFeedbackSpheroids} node operator class.
    !!}
    implicit none
    type (nodeOperatorStellarFeedbackSpheroids)                        :: self
    class(starFormationRateSpheroidsClass     ), intent(in   ), target :: starFormationRateSpheroids_
    class(stellarPopulationPropertiesClass    ), intent(in   ), target :: stellarPopulationProperties_
    class(stellarFeedbackOutflowsClass        ), intent(in   ), target :: stellarFeedbackOutflows_
    !![
    <constructorAssign variables="*starFormationRateSpheroids_, *stellarPopulationProperties_, *stellarFeedbackOutflows_"/>
    !!]

    return
  end function stellarFeedbackSpheroidsConstructorInternal

  subroutine stellarFeedbackSpheroidsDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily stellarFeedbackSpheroids} node operator class.
    !!}
    implicit none
    type(nodeOperatorStellarFeedbackSpheroids), intent(inout) :: self

    !![
    <objectDestructor name="self%starFormationRateSpheroids_" />
    <objectDestructor name="self%stellarPopulationProperties_"/>
    <objectDestructor name="self%stellarFeedbackOutflows_"    />
    !!]
    return
  end subroutine stellarFeedbackSpheroidsDestructor
  
  subroutine stellarFeedbackSpheroidsDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Perform feedback from stars in a spheroid.
    !!}
    use :: Abundances_Structure          , only : abundances         , zeroAbundances
    use :: Galacticus_Nodes              , only : propertyInactive   , nodeComponentSpheroid, nodeComponentHotHalo
    use :: Histories                     , only : history
    use :: Stellar_Luminosities_Structure, only : stellarLuminosities
    implicit none
    class           (nodeOperatorStellarFeedbackSpheroids), intent(inout), target  :: self
    type            (treeNode                            ), intent(inout), target  :: node
    logical                                               , intent(inout)          :: interrupt
    procedure       (interruptTask                       ), intent(inout), pointer :: functionInterrupt
    integer                                               , intent(in   )          :: propertyType
    class           (nodeComponentSpheroid               )               , pointer :: spheroid
    class           (nodeComponentHotHalo                )               , pointer :: hotHalo
    double precision                                      , parameter              :: radiusMinimum           =1.0d-12 ! ⎫
    double precision                                      , parameter              :: massMinimum             =1.0d-06 ! ⎬ Minimum absolute scales for physically plausible spheroids.
    double precision                                      , parameter              :: angularMomentumMinimum  =1.0d-20 ! ⎭
    double precision                                                               :: rateStarFormation               , massGas                   , &
         &                                                                            angularMomentum                 , fractionEjected           , &
         &                                                                            rateMassStellar                 , rateEnergyInput           , &
         &                                                                            rateMassFuel                    , rateMassOutflowEjected    , &
         &                                                                            rateMassOutflowExpelled         , rateMassOutflowTotal      , &
         &                                                                            massSpheroid                    , rateAngularMomentumOutflow
    type            (abundances                          )                         :: abundancesGas                   , rateAbundancesFuels       , &
         &                                                                            rateAbundancesStellar           , rateAbundancesOutflow
    type            (history                             )                         :: ratePropertiesStellar
    type            (stellarLuminosities                 )                         :: rateLuminositiesStellar

    ! Do nothing during inactive property solving.
    if (propertyInactive(propertyType)) return
    ! Check for a realistic spheroid, return immediately if spheroid is unphysical.
    spheroid => node%spheroid()
    if     (     spheroid%angularMomentum() < angularMomentumMinimum &
         &  .or. spheroid%radius         () <          radiusMinimum &
         &  .or. spheroid%massGas        () <            massMinimum &
         & ) return
    ! Get the star formation rate.
    rateStarFormation=self%starFormationRateSpheroids_%rate(node)   
    ! Compute abundances of star forming gas.
    massGas      =spheroid%massGas      ()
    abundancesGas=spheroid%abundancesGas()
    call abundancesGas%massToMassFraction(massGas)
    ! Find rates of change of stellar mass, gas mass, abundances and luminosities.
    ratePropertiesStellar=spheroid%stellarPropertiesHistory()
    call self%stellarPopulationProperties_%rates(                                                      &
         &                                                                    rateStarFormation      , &
         &                                                                    abundancesGas          , &
         &                                                                    spheroid               , &
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
    ! Find rate of outflow of material from the spheroid.
    call self%stellarFeedbackOutflows_%outflowRate(spheroid,rateStarFormation,rateEnergyInput,rateMassOutflowEjected,rateMassOutflowExpelled)
    rateMassOutflowTotal=+rateMassOutflowEjected  &
         &               +rateMassOutflowExpelled
    if (rateMassOutflowTotal > 0.0d0) then
       ! Find the fraction of material which outflows to the hot halo.
       fractionEjected=+rateMassOutflowEjected &
            &          /rateMassOutflowTotal
       ! Get the total mass of the spheroid.
       massSpheroid=+     massGas       &
            &       +spheroid%massStellar()
       ! Compute the angular momentum outflow rate.
       if (massSpheroid > 0.0d0) then
          angularMomentum           =spheroid%angularMomentum()
          rateAngularMomentumOutflow=angularMomentum*(rateMassOutflowTotal/massSpheroid)
       else
          rateAngularMomentumOutflow=0.0d0
       end if
       if (massGas  > 0.0d0) then
          rateAbundancesOutflow=abundancesGas*rateMassOutflowTotal
       else
          rateAbundancesOutflow=zeroAbundances
       end if
       hotHalo => node%hotHalo()
       call hotHalo %           outflowingMassRate( rateMassOutflowTotal      *fractionEjected)
       call spheroid%                  massGasRate(-rateMassOutflowTotal                      )
       call hotHalo %outflowingAngularMomentumRate( rateAngularMomentumOutflow*fractionEjected)
       call spheroid%          angularMomentumRate(-rateAngularMomentumOutflow                )
       call hotHalo %     outflowingAbundancesRate( rateAbundancesOutflow     *fractionEjected)
       call spheroid%            abundancesGasRate(-rateAbundancesOutflow                     )
    end if
    return
  end subroutine stellarFeedbackSpheroidsDifferentialEvolution

