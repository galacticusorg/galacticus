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
  Implements a node operator class that performs star formation in spheroids.
  !!}

  use :: Star_Formation_Rates_Spheroids, only : starFormationRateSpheroidsClass
  use :: Stellar_Population_Properties , only : stellarPopulationPropertiesClass
  use :: Star_Formation_Histories      , only : starFormationHistoryClass

  !![
  <nodeOperator name="nodeOperatorStarFormationSpheroids">
   <description>A node operator class that performs star formation.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorStarFormationSpheroids
     !!{
     A node operator class that shifts node indices at node promotion.
     !!}
     private
     class  (starFormationRateSpheroidsClass ), pointer :: starFormationRateSpheroids_  => null()
     class  (stellarPopulationPropertiesClass), pointer :: stellarPopulationProperties_ => null()
     class  (starFormationHistoryClass       ), pointer :: starFormationHistory_        => null()
     logical                                            :: luminositiesStellarInactive
   contains
     final     ::                                   starFormationSpheroidsDestructor
     procedure :: differentialEvolutionAnalytics => starFormationSpheroidsDifferentialEvolutionAnalytics
     procedure :: differentialEvolution          => starFormationSpheroidsDifferentialEvolution
  end type nodeOperatorStarFormationSpheroids
  
  interface nodeOperatorStarFormationSpheroids
     !!{
     Constructors for the \refClass{nodeOperatorStarFormationSpheroids} node operator class.
     !!}
     module procedure starFormationSpheroidsConstructorParameters
     module procedure starFormationSpheroidsConstructorInternal
  end interface nodeOperatorStarFormationSpheroids
  
contains

  function starFormationSpheroidsConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorStarFormationSpheroids} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (nodeOperatorStarFormationSpheroids)                :: self
    type   (inputParameters                   ), intent(inout) :: parameters
    class  (starFormationRateSpheroidsClass   ), pointer       :: starFormationRateSpheroids_
    class  (stellarPopulationPropertiesClass  ), pointer       :: stellarPopulationProperties_
    class  (starFormationHistoryClass         ), pointer       :: starFormationHistory_
    logical                                                    :: luminositiesStellarInactive
    
    !![
    <inputParameter>
      <name>luminositiesStellarInactive</name>
      <defaultValue>.false.</defaultValue>
      <source>parameters</source>
      <description>If true, stellar luminosities will be treated as inactive properties.</description>
    </inputParameter>
    <objectBuilder class="starFormationRateSpheroids"  name="starFormationRateSpheroids_"  source="parameters"/>
    <objectBuilder class="stellarPopulationProperties" name="stellarPopulationProperties_" source="parameters"/>
    <objectBuilder class="starFormationHistory"        name="starFormationHistory_"        source="parameters"/>
    !!]
    self=nodeOperatorStarFormationSpheroids(luminositiesStellarInactive,starFormationRateSpheroids_,stellarPopulationProperties_,starFormationHistory_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="starFormationRateSpheroids_" />
    <objectDestructor name="stellarPopulationProperties_"/>
    <objectDestructor name="starFormationHistory_"       />
    !!]
    return
  end function starFormationSpheroidsConstructorParameters

  function starFormationSpheroidsConstructorInternal(luminositiesStellarInactive,starFormationRateSpheroids_,stellarPopulationProperties_,starFormationHistory_) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorStarFormationSpheroids} node operator class.
    !!}
    implicit none
    type   (nodeOperatorStarFormationSpheroids)                        :: self
    class  (starFormationRateSpheroidsClass   ), intent(in   ), target :: starFormationRateSpheroids_
    class  (stellarPopulationPropertiesClass  ), intent(in   ), target :: stellarPopulationProperties_
    class(starFormationHistoryClass           ), intent(in   ), target :: starFormationHistory_
    logical                                    , intent(in   )         :: luminositiesStellarInactive
    !![
    <constructorAssign variables="luminositiesStellarInactive, *starFormationRateSpheroids_, *stellarPopulationProperties_, *starFormationHistory_"/>
    !!]

    return
  end function starFormationSpheroidsConstructorInternal

  subroutine starFormationSpheroidsDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorStarFormationSpheroids} node operator class.
    !!}
    implicit none
    type(nodeOperatorStarFormationSpheroids), intent(inout) :: self

    !![
    <objectDestructor name="self%starFormationRateSpheroids_" />
    <objectDestructor name="self%stellarPopulationProperties_"/>
    <objectDestructor name="self%starFormationHistory_"       />
    !!]
    return
  end subroutine starFormationSpheroidsDestructor

  subroutine starFormationSpheroidsDifferentialEvolutionAnalytics(self,node)
    !!{
    Initialize the mass transfer fraction.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSpheroid
    implicit none
    class(nodeOperatorStarFormationSpheroids), intent(inout) :: self
    type (treeNode                          ), intent(inout) :: node
    class(nodeComponentSpheroid             ), pointer       :: spheroid

    ! Mark the formed stellar mass as analytically-solvable (it is always zero) if we are not solving for luminosities as inactive
    ! properties.
    if (.not.self%luminositiesStellarInactive) then
       spheroid => node%spheroid()
       select type (spheroid)
       type is (nodeComponentSpheroid)
          ! Spheroid does not yet exist.
       class default
          call spheroid%massStellarFormedAnalytic()
       end select
    end if
    return
  end subroutine starFormationSpheroidsDifferentialEvolutionAnalytics

  subroutine starFormationSpheroidsDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Perform star formation in a spheroid.
    !!}
    use :: Abundances_Structure          , only : abundances
    use :: Galacticus_Nodes              , only : propertyInactive   , propertyTypeActive, propertyEvaluate, nodeComponentSpheroid
    use :: Histories                     , only : history
    use :: Stellar_Luminosities_Structure, only : stellarLuminosities
    implicit none
    class           (nodeOperatorStarFormationSpheroids), intent(inout), target  :: self
    type            (treeNode                          ), intent(inout), target  :: node
    logical                                             , intent(inout)          :: interrupt
    procedure       (interruptTask                     ), intent(inout), pointer :: functionInterrupt
    integer                                             , intent(in   )          :: propertyType
    class           (nodeComponentSpheroid             )               , pointer :: spheroid
    double precision                                    , parameter              :: radiusMinimum           =1.0d-12 ! ⎫
    double precision                                    , parameter              :: massMinimum             =1.0d-06 ! ⎬ Minimum absolute scales for physically plausible spheroids.
    double precision                                    , parameter              :: angularMomentumMinimum  =1.0d-20 ! ⎭
    double precision                                                             :: rateStarFormation               , massFuel             , &
         &                                                                          rateMassStellar                 , rateEnergyInput      , &
         &                                                                          rateMassFuel
    logical                                                                      :: luminositiesCompute
    type            (abundances                        )                         :: abundancesFuel                  , rateAbundancesFuels  , &
         &                                                                          rateAbundancesStellar
    type            (history                           )                         :: rateHistoryStarFormation        , ratePropertiesStellar
    type            (stellarLuminosities               )                         :: rateLuminositiesStellar
    
    ! Check for a realistic spheroid, return immediately if spheroid is unphysical.
    spheroid => node%spheroid()
    if     (     spheroid%angularMomentum() < angularMomentumMinimum &
         &  .or. spheroid%radius         () <          radiusMinimum &
         &  .or. spheroid%massGas        () <            massMinimum &
         &  .or.                                                     &
         &   (                                                       &
         &          propertyInactive(propertyType)                   &
         &    .and.                                                  & 
         &     .not.self%luminositiesStellarInactive                 &
         &   )                                                       &
        & ) return
    if (propertyInactive(propertyType)) then
       ! For inactive property solution make use of the "massStellarFormed" property to determine the star formation rate.
       rateStarFormation=spheroid%massStellarFormedRateGet()
    else
       ! During active property solution, integrate the star formation rate so that we will have a solution for the total mass
       ! of stars formed as a function of time. This differs from the stellar mass due to recycling, and possibly transfer of
       ! stellar mass to other components.
       rateStarFormation=self%starFormationRateSpheroids_%rate(node)   
       if (self%luminositiesStellarInactive) call spheroid%massStellarFormedRate(rateStarFormation)
    end if
    ! Compute abundances of star forming gas.
    massFuel      =spheroid%massGas      ()
    abundancesFuel=spheroid%abundancesGas()
    call abundancesFuel%massToMassFraction(massFuel)
    ! Determine if luminosities must be computed.
    luminositiesCompute=propertyEvaluate(propertyType,self%luminositiesStellarInactive)
    ! Find rates of change of stellar mass, gas mass, abundances and luminosities.
    ratePropertiesStellar=spheroid%stellarPropertiesHistory()
    call self%stellarPopulationProperties_%rates(                         &
         &                                       rateStarFormation      , &
         &                                       abundancesFuel         , &
         &                                       spheroid               , &
         &                                       node                   , &
         &                                       ratePropertiesStellar  , &
         &                                       rateMassStellar        , &
         &                                       rateMassFuel           , &
         &                                       rateEnergyInput        , &
         &                                       rateAbundancesFuels    , &
         &                                       rateAbundancesStellar  , &
         &                                       rateLuminositiesStellar, &
         &                                       luminositiesCompute      &
         &                                      )
    ! Adjust rates.
    if (propertyEvaluate(propertyTypeActive,propertyIsInactive=.false.)) then
       rateHistoryStarFormation=spheroid%starFormationHistory()
       call        rateHistoryStarFormation%reset                       (                                                              )
       call self  %starFormationHistory_   %                        rate(node,rateHistoryStarFormation,abundancesFuel,rateStarFormation)
       call        spheroid                %             massStellarRate(     rateMassStellar                                          )
       call        spheroid                %                 massGasRate(     rateMassFuel                                             )
       call        spheroid                %       abundancesStellarRate(     rateAbundancesStellar                                    )
       call        spheroid                %           abundancesGasRate(     rateAbundancesFuels                                      )
       if (ratePropertiesStellar   %exists())                                                                                            &
            & call spheroid                %stellarPropertiesHistoryRate(     ratePropertiesStellar                                    )
       if (rateHistoryStarFormation%exists())                                                                                            &
            & call spheroid                %    starFormationHistoryRate(     rateHistoryStarFormation                                 )
    end if
    if    (luminositiesCompute              )                                                                                            &
         &    call spheroid                %    luminositiesStellarRate (     rateLuminositiesStellar                                  )
    return
  end subroutine starFormationSpheroidsDifferentialEvolution

