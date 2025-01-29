!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024
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
  Implements a node operator class that performs star formation in nuclear star cluster.
  !!}

  use :: Star_Formation_Rates_NSCs     , only : starFormationRateNSCsClass
  use :: Stellar_Population_Properties , only : stellarPopulationPropertiesClass
  use :: Star_Formation_Histories      , only : starFormationHistoryClass

  !![
  <nodeOperator name="nodeOperatorStarFormationNSCs">
   <description>A node operator class that performs star formation.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorStarFormationNSCs
     !!{
     A node operator class that shifts node indices at node promotion.
     !!}
     private
     class           (starFormationRateNSCsClass      ), pointer :: starFormationRateNSCs_        => null()
     class           (stellarPopulationPropertiesClass), pointer :: stellarPopulationProperties_  => null()
     class           (starFormationHistoryClass       ), pointer :: starFormationHistory_         => null()
     logical                                                     :: luminositiesStellarInactive
   contains
     final     ::                                        starFormationNSCsDestructor
     procedure :: differentialEvolution               => starFormationNSCsDifferentialEvolution
     procedure :: differentialEvolutionStepFinalState => starFormationNSCsDifferentialEvolutionAnalytics
  end type nodeOperatorStarFormationNSCs
  
  interface nodeOperatorStarFormationNSCs
     !!{
     Constructors for the {\normalfont \ttfamily starFormationNSCs} node operator class.
     !!}
     module procedure starFormationNSCsConstructorParameters
     module procedure starFormationNSCsConstructorInternal
  end interface nodeOperatorStarFormationNSCs
  
contains

  function starFormationNSCsConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily starFormation} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (nodeOperatorStarFormationNSCs   )                :: self
    type   (inputParameters                 ), intent(inout) :: parameters
    class  (starFormationRateNSCsClass      ), pointer       :: starFormationRateNSCs_
    class  (stellarPopulationPropertiesClass), pointer       :: stellarPopulationProperties_
    class  (starFormationHistoryClass       ), pointer       :: starFormationHistory_
    logical                                                  :: luminositiesStellarInactive
    
    !![
    <inputParameter>
      <name>luminositiesStellarInactive</name>
      <defaultValue>.false.</defaultValue>
      <source>parameters</source>
      <description>If true, stellar luminosities will be treated as inactive properties.</description>
    </inputParameter>
    <objectBuilder class="starFormationRateNSCs"         name="starFormationRateNSCs_"       source="parameters"/>
    <objectBuilder class="stellarPopulationProperties"   name="stellarPopulationProperties_" source="parameters"/>
    <objectBuilder class="starFormationHistory"          name="starFormationHistory_"        source="parameters"/>
    !!]
    self=nodeOperatorStarFormationNSCs(luminositiesStellarInactive,starFormationRateNSCs_,stellarPopulationProperties_,starFormationHistory_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="starFormationRateNSCs_"         />
    <objectDestructor name="stellarPopulationProperties_"  />
    <objectDestructor name="starFormationHistory_"         />
    !!]
    return
  end function starFormationNSCsConstructorParameters

  function starFormationNSCsConstructorInternal(luminositiesStellarInactive,starFormationRateNSCs_,stellarPopulationProperties_,starFormationHistory_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily starFormationNSC} node operator class.
    !!}
    implicit none
    type   (nodeOperatorStarFormationNSCs   )                        :: self
    class  (starFormationRateNSCsClass      ), intent(in   ), target :: starFormationRateNSCs_
    class  (stellarPopulationPropertiesClass), intent(in   ), target :: stellarPopulationProperties_
    class  (starFormationHistoryClass       ), intent(in   ), target :: starFormationHistory_
    logical                                  , intent(in   )         :: luminositiesStellarInactive
    !![
    <constructorAssign variables="luminositiesStellarInactive, *starFormationRateNSCs_, *stellarPopulationProperties_, *starFormationHistory_"/>
    !!]
    return
  end function starFormationNSCsConstructorInternal

  subroutine starFormationNSCsDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily starFormationNSCs} node operator class.
    !!}
    implicit none
    type(nodeOperatorStarFormationNSCs), intent(inout) :: self

    !![
    <objectDestructor name="self%starFormationRateNSCs_"      />
    <objectDestructor name="self%stellarPopulationProperties_"/>
    <objectDestructor name="self%starFormationHistory_"       />
    !!]
    return
  end subroutine starFormationNSCsDestructor
  
  subroutine starFormationNSCsDifferentialEvolutionAnalytics(self,node)
    !!{
    Initialize the mass transfer fraction.    
    !!}
    use :: Galacticus_Nodes            , only : nodeComponentNSC
    implicit none
    class(nodeOperatorStarFormationNSCs), intent(inout) :: self
    type (treeNode                     ), intent(inout) :: node
    class(nodeComponentNSC             ), pointer       :: NSC

    ! Mark the formed stellar mass as analytically-solvable (it is always zero) if we are not solving for luminosities as inactive
    ! properties.
    if (.not.self%luminositiesStellarInactive) then
       NSC => node%NSC()
       select type (NSC)
       type is (nodeComponentNSC)
          ! NSC does not yet exist.
       class default
          call NSC%massStellarFormedAnalytic()
       end select
    end if  
    return
  end subroutine starFormationNSCsDifferentialEvolutionAnalytics

  subroutine starFormationNSCsDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Perform star formation in a nuclear star cluster.
    !!}
    use :: Abundances_Structure          , only : abundances
    use :: Galacticus_Nodes              , only : propertyInactive   , propertyTypeActive, propertyEvaluate, nodeComponentNSC
    use :: Histories                     , only : history
    use :: Stellar_Luminosities_Structure, only : stellarLuminosities
    implicit none
    class           (nodeOperatorStarFormationNSCs), intent(inout), target  :: self
    type            (treeNode                     ), intent(inout), target  :: node
    logical                                        , intent(inout)          :: interrupt
    procedure       (interruptTask                ), intent(inout), pointer :: functionInterrupt
    integer                                        , intent(in   )          :: propertyType
    class           (nodeComponentNSC             )               , pointer :: NSC
    double precision                                                        :: rateStarFormation       , massFuel                , &
         &                                                                     rateMassStellar         , rateEnergyInput         , &
         &                                                                     rateMassFuel            
    logical                                                                 :: luminositiesCompute
    type            (abundances                   )                         :: abundancesFuel          , rateAbundancesFuels     , &
         &                                                                     rateAbundancesStellar
    type            (history                      )                         :: rateHistoryStarFormation, ratePropertiesStellar
    type            (stellarLuminosities          )                         :: rateLuminositiesStellar 
    
    ! Check for a realistic nuclear star cluster, return immediately if nuclear star cluster is unphysical.
    NSC => node%NSC()
    if     (    NSC%radius () < 0.0d0 &
         &  .or.NSC%massGas() < 0.0d0 &
         & ) return
    if (propertyInactive(propertyType)) then
       ! For inactive property solution make use of the "massStellarFormed" property to determine the star formation rate.
       rateStarFormation=NSC%massStellarFormedRateGet()
    else
       ! During active property solution, integrate the star formation rate so that we will have a solution for the total mass
       ! of stars formed as a function of time. This differs from the stellar mass due to recycling, and possibly transfer of
       ! stellar mass to other components.
       rateStarFormation=self%starFormationRateNSCs_%rate(node)   
       call NSC%massStellarFormedRate(rateStarFormation)
       if (self%luminositiesStellarInactive) call NSC%massStellarFormedRate(rateStarFormation)
    end if
    ! Compute abundances of star forming gas.
    massFuel      =NSC%massGas      ()
    abundancesFuel=NSC%abundancesGas()
    call abundancesFuel%massToMassFraction(massFuel)
    ! Determine if luminosities must be computed.
    luminositiesCompute=propertyEvaluate(propertyType,self%luminositiesStellarInactive)
    ! Find rates of change of stellar mass, gas mass, abundances and luminosities.
    ratePropertiesStellar=NSC%stellarPropertiesHistory()
    call self%stellarPopulationProperties_%rates(                         &
         &                                       rateStarFormation      , &
         &                                       abundancesFuel         , &
         &                                       NSC                    , &
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
       rateHistoryStarFormation=NSC%starFormationHistory()
       call        rateHistoryStarFormation%reset                      (                                                              )
       call self  %starFormationHistory_   %                        rate(node,rateHistoryStarFormation,abundancesFuel,rateStarFormation)
       call        NSC                     %             massStellarRate(     rateMassStellar                                          )
       call        NSC                     %                 massGasRate(     rateMassFuel                                             )
       call        NSC                     %       abundancesStellarRate(     rateAbundancesStellar                                    )
       call        NSC                     %           abundancesGasRate(     rateAbundancesFuels                                      )
       if (ratePropertiesStellar   %exists())                                                                                            &
            & call NSC                     %stellarPropertiesHistoryRate(     ratePropertiesStellar                                    )
       if (rateHistoryStarFormation%exists())                                                                                            &
            & call NSC                     %    starFormationHistoryRate(     rateHistoryStarFormation                                 )
    end if
    if (luminositiesCompute                 )                                                                                            &
        & call NSC                         %luminositiesStellarRate(rateLuminositiesStellar                                            )
    return
  end subroutine starFormationNSCsDifferentialEvolution
  
