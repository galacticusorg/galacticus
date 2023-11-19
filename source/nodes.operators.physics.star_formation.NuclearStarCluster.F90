!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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

  use :: Star_Formation_Rates_NSC      , only : starFormationRateNSCClass
  use :: Stellar_Population_Properties , only : stellarPopulationPropertiesClass
  use :: Star_Formation_Histories      , only : starFormationHistoryClass

  !![
  <nodeOperator name="nodeOperatorStarFormationNSC">
   <description>A node operator class that performs star formation.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorStarFormationNSC
     !!{
     A node operator class that shifts node indices at node promotion.
     !!}
     private
     class           (starFormationRateNSCClass       ), pointer :: starFormationRateNSC_         => null()
     class           (stellarPopulationPropertiesClass), pointer :: stellarPopulationProperties_  => null()
     class           (starFormationHistoryClass       ), pointer :: starFormationHistory_         => null()
     logical                                                     :: luminositiesStellarInactive
     double precision                                            :: fractionMassRetainedInitial            , fractionMassRetainedFinal
   contains
     final     ::                                        starFormationNSCDestructor
     procedure :: differentialEvolutionPre            => starFormationNSCDifferentialEvolutionPre
     procedure :: differentialEvolution               => starFormationNSCDifferentialEvolution
     procedure :: differentialEvolutionStepFinalState => starFormationNSCDifferentialEvolutionStepFinalState
  end type nodeOperatorStarFormationNSC
  
  interface nodeOperatorStarFormationNSC
     !!{
     Constructors for the {\normalfont \ttfamily starFormationNSC} node operator class.
     !!}
     module procedure starFormationNSCConstructorParameters
     module procedure starFormationNSCConstructorInternal
  end interface nodeOperatorStarFormationNSC
  
contains

  function starFormationNSCConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily starFormation} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (nodeOperatorStarFormationNSC    )                :: self
    type   (inputParameters                 ), intent(inout) :: parameters
    class  (starFormationRateNSCClass       ), pointer       :: starFormationRateNSC_
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
    <objectBuilder class="starFormationRateNSC"          name="starFormationRateNSC_"        source="parameters"/>
    <objectBuilder class="stellarPopulationProperties"   name="stellarPopulationProperties_" source="parameters"/>
    <objectBuilder class="starFormationHistory"          name="starFormationHistory_"        source="parameters"/>
    !!]
    self=nodeOperatorStarFormationNSC(luminositiesStellarInactive,starFormationRateNSC_,stellarPopulationProperties_,starFormationHistory_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="starFormationRateNSC_"         />
    <objectDestructor name="stellarPopulationProperties_"  />
    <objectDestructor name="starFormationHistory_"         />
    !!]
    return
  end function starFormationNSCConstructorParameters

  function starFormationNSCConstructorInternal(luminositiesStellarInactive,starFormationRateNSC_,stellarPopulationProperties_,starFormationHistory_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily starFormationNSC} node operator class.
    !!}
    implicit none
    type   (nodeOperatorStarFormationNSC    )                        :: self
    class  (starFormationRateNSCClass       ), intent(in   ), target :: starFormationRateNSC_
    class  (stellarPopulationPropertiesClass), intent(in   ), target :: stellarPopulationProperties_
    class  (starFormationHistoryClass       ), intent(in   ), target :: starFormationHistory_
    logical                                  , intent(in   )         :: luminositiesStellarInactive
    !![
    <constructorAssign variables="luminositiesStellarInactive, *starFormationRateNSC_, *stellarPopulationProperties_, *starFormationHistory_"/>
    !!]

    ! Initialize values.
    self%fractionMassRetainedInitial=1.0d0
    self%fractionMassRetainedFinal  =1.0d0
    return
  end function starFormationNSCConstructorInternal

  subroutine starFormationNSCDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily starFormationNSC} node operator class.
    !!}
    implicit none
    type(nodeOperatorStarFormationNSC), intent(inout) :: self

    !![
    <objectDestructor name="self%starFormationRateNSC_"       />
    <objectDestructor name="self%stellarPopulationProperties_"/>
    <objectDestructor name="self%starFormationHistory_"       />
    !!]
    return
  end subroutine starFormationNSCDestructor
  
  subroutine starFormationNSCDifferentialEvolutionPre(self,node)
    !!{
    Initialize the mass transfer fraction.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentNSC
    implicit none
    class(nodeOperatorStarFormationNSC), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node
    class(nodeComponentNSC            ), pointer       :: NSC

    ! Initialize the mass transferred fraction to unity. The value is arbitrary as only ratios of this quantity are used, but
    ! must be non-zero.
    NSC => node%NSC()
    select type (NSC)
    type is (nodeComponentNSC)
       ! NSC does not yet exist.
    class default
       self%fractionMassRetainedInitial=1.0d0
       self%fractionMassRetainedFinal  =1.0d0
       if (NSC%fractionMassRetainedIsSettable()) call NSC%fractionMassRetainedSet(1.0d0)
    end select
    return
  end subroutine starFormationNSCDifferentialEvolutionPre

  subroutine starFormationNSCDifferentialEvolutionStepFinalState(self,node)
    !!{
    Record the final state of the Nuclear Star Cluster at the end of the timestep prior to begin evaluation of integrals for inactive
    properties.
    !!}
    use :: Galacticus_Nodes            , only : nodeComponentNSC
    implicit none
    class(nodeOperatorStarFormationNSC), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node
    class(nodeComponentNSC            ), pointer       :: NSC

    ! The retained mass fraction at the start of this step is just the fraction at the end of the previous step. Then update the
    ! retained fraction at the end of the current step.
    NSC                              =>           node%NSC                      ()
    self%fractionMassRetainedInitial =            self%fractionMassRetainedFinal
    self%fractionMassRetainedFinal   =  max(0.0d0,NSC %fractionMassRetained     ())
    return
  end subroutine starFormationNSCDifferentialEvolutionStepFinalState

  subroutine starFormationNSCDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Perform star formation in a nuclear star cluster.
    !!}
    use :: Abundances_Structure          , only : abundances
    use :: Galacticus_Nodes              , only : propertyInactive     , propertyTypeActive, propertyEvaluate, nodeComponentNSC, &
         &                                        nodeComponentSpheroid
    use :: Histories                     , only : history
    use :: Stellar_Luminosities_Structure, only : stellarLuminosities
    implicit none
    class           (nodeOperatorStarFormationNSC), intent(inout), target  :: self
    type            (treeNode                    ), intent(inout), target  :: node
    logical                                       , intent(inout)          :: interrupt
    procedure       (interruptTask               ), intent(inout), pointer :: functionInterrupt
    integer                                       , intent(in   )          :: propertyType
    class           (nodeComponentNSC            )               , pointer :: NSC
    class           (nodeComponentSpheroid       )               , pointer :: spheroid
    double precision                                                       :: rateStarFormation       , massFuel                , &
         &                                                                    rateMassStellar         , rateEnergyInput         , &
         &                                                                    rateMassFuel            , fractionMassRetained    , &
         &                                                                    fractionMassRetainedRate
    logical                                                                :: luminositiesCompute
    type            (abundances                  )                         :: abundancesFuel          , rateAbundancesFuels     , &
         &                                                                    rateAbundancesStellar
    type            (history                     )                         :: rateHistoryStarFormation, ratePropertiesStellar
    type            (stellarLuminosities         )                         :: rateLuminositiesStellar , rateLuminositiesTransfer
    
    ! Check for a realistic nuclear star cluster, return immediately if nuclear star cluster is unphysical.
    NSC => node%NSC()
    if     (     NSC%radius         () < 0.0d0 &
         &  .or. NSC%massGas        () < 0.0d0 &
         & ) return
    if (propertyInactive(propertyType)) then
       ! For inactive property solution make use of the "massStellarFormed" property to determine the star formation rate.
       rateStarFormation=NSC%massStellarFormedRateGet()
    else
       ! During active property solution, integrate the star formation rate so that we will have a solution for the total mass
       ! of stars formed as a function of time. This differs from the stellar mass due to recycling, and possibly transfer of
       ! stellar mass to other components.
       rateStarFormation=self%starFormationRateNSC_%rate(node)   
       call NSC%massStellarFormedRate(rateStarFormation)
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
       call self  %starFormationHistory_  %                        rate(node,rateHistoryStarFormation,abundancesFuel,rateStarFormation)
       call        NSC                    %             massStellarRate(     rateMassStellar                                          )
       call        NSC                    %                 massGasRate(     rateMassFuel                                             )
       call        NSC                    %       abundancesStellarRate(     rateAbundancesStellar                                    )
       call        NSC                    %           abundancesGasRate(     rateAbundancesFuels                                      )
       if (ratePropertiesStellar   %exists())                                                                                            &
            & call NSC                    %stellarPropertiesHistoryRate(     ratePropertiesStellar                                    )
       if (rateHistoryStarFormation%exists())                                                                                            &
            & call NSC                    %    starFormationHistoryRate(     rateHistoryStarFormation                                 )
    end if
    if (luminositiesCompute) then
       ! For inactive property calculations we must check if any mass (and, therefore, light) is being transferred to the
       ! spheroid component. If it is, our integrand must account for this mass transfer. The fractions of mass retained and
       ! transferred are determined from the "fractionMassRetained" property which is computed during differential evolution.
       if (propertyInactive(propertyType) .and. self%fractionMassRetainedFinal < self%fractionMassRetainedInitial) then
          spheroid => node%spheroid()
          ! Determine the fraction of mass (and light) formed at this time which will be retained in the NSC at the final time in the step.
          if      (self%fractionMassRetainedFinal   == 0.0d0                      ) then
             ! The retained fraction reached zero by the end of the step, so no mass is retained.
             fractionMassRetained    =                                                                          0.0d0
          else if (self%fractionMassRetainedFinal   >  NSC%fractionMassRetained()) then
             fractionMassRetained    =                                                                          1.0d0
          else
             ! Limit the retained fraction to unity (to avoid any rounding errors).
             fractionMassRetained    =    self%fractionMassRetainedFinal    /NSC%fractionMassRetained       ()
          end if
          ! Determine the rate at which mass (and light) that was pre-existing at the start of this timestep is being transferred.
          if      (self%fractionMassRetainedInitial == 0.0d0                      ) then
             ! The initial retained fraction was zero, so there should be no light to transfer - set a transfer rate of zero.
             fractionMassRetainedRate=                                                                        0.0d0
          else
             ! Limit the transfer rate of pre-existing light to be negative - it is not possible to transfer light *to* the
             ! NSC, so any positive value here can arise only via rounding errors.
             fractionMassRetainedRate=min(NSC%fractionMassRetainedRateGet()/self%fractionMassRetainedInitial  ,0.0d0)
          end if
          ! Find the rate of transfer of pre-existing light.
          rateLuminositiesTransfer=+NSC%luminositiesStellar() &
               &                   *fractionMassRetainedRate
          ! Evaluate the integrand for the NSC, and the corresponding one for the spheroid to account for the transfer of light.
          call    NSC  %luminositiesStellarRate(rateLuminositiesStellar*       fractionMassRetained +rateLuminositiesTransfer                            )
          call spheroid%luminositiesStellarRate(rateLuminositiesStellar*(1.0d0-fractionMassRetained)-rateLuminositiesTransfer,interrupt,functionInterrupt)
       else
          ! In this case we do not need to account for transfer of light to the spheroid because either:
          !  a) there is none, or;
          !  b) we are solving for luminosities as active properties in which case transfer to the spheroid is handled directly in the ODE.
          call     NSC%luminositiesStellarRate(rateLuminositiesStellar                                                                                  )
       end if
    end if
    return
  end subroutine starFormationNSCDifferentialEvolution
  
