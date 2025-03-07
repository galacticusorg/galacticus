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
  Implements a node operator class that performs star formation in \glspl{nsc}.
  !!}

  use :: Star_Formation_Rates_Nuclear_Star_Clusters, only : starFormationRateNuclearStarClustersClass
  use :: Stellar_Population_Properties             , only : stellarPopulationPropertiesClass
  use :: Star_Formation_Histories                  , only : starFormationHistoryClass

  !![
  <nodeOperator name="nodeOperatorStarFormationNuclearStarClusters">
   <description>A node operator class that performs star formation in \glspl{nsc}.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorStarFormationNuclearStarClusters
     !!{
     A node operator class that performs star formation in \glspl{nsc}.
     !!}
     private
     class  (starFormationRateNuclearStarClustersClass), pointer :: starFormationRateNuclearStarClusters_ => null()
     class  (stellarPopulationPropertiesClass         ), pointer :: stellarPopulationProperties_          => null()
     class  (starFormationHistoryClass                ), pointer :: starFormationHistory_                 => null()
     logical                                                     :: luminositiesStellarInactive
   contains
     final     ::                                        starFormationNuclearStarClustersDestructor
     procedure :: differentialEvolution               => starFormationNuclearStarClustersDifferentialEvolution
     procedure :: differentialEvolutionStepFinalState => starFormationNuclearStarClustersDifferentialEvolutionAnalytics
  end type nodeOperatorStarFormationNuclearStarClusters
  
  interface nodeOperatorStarFormationNuclearStarClusters
     !!{
     Constructors for the {\normalfont \ttfamily starFormationNuclearStarClusters} node operator class.
     !!}
     module procedure starFormationNuclearStarClustersConstructorParameters
     module procedure starFormationNuclearStarClustersConstructorInternal
  end interface nodeOperatorStarFormationNuclearStarClusters
  
contains

  function starFormationNuclearStarClustersConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily starFormation} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (nodeOperatorStarFormationNuclearStarClusters)                :: self
    type   (inputParameters                             ), intent(inout) :: parameters
    class  (starFormationRateNuclearStarClustersClass   ), pointer       :: starFormationRateNuclearStarClusters_
    class  (stellarPopulationPropertiesClass            ), pointer       :: stellarPopulationProperties_
    class  (starFormationHistoryClass                   ), pointer       :: starFormationHistory_
    logical                                                              :: luminositiesStellarInactive
    
    !![
    <inputParameter>
      <name>luminositiesStellarInactive</name>
      <defaultValue>.false.</defaultValue>
      <source>parameters</source>
      <description>If true, stellar luminosities will be treated as inactive properties.</description>
    </inputParameter>
    <objectBuilder class="starFormationRateNuclearStarClusters" name="starFormationRateNuclearStarClusters_" source="parameters"/>
    <objectBuilder class="stellarPopulationProperties"          name="stellarPopulationProperties_"          source="parameters"/>
    <objectBuilder class="starFormationHistory"                 name="starFormationHistory_"                 source="parameters"/>
    !!]
    self=nodeOperatorStarFormationNuclearStarClusters(luminositiesStellarInactive,starFormationRateNuclearStarClusters_,stellarPopulationProperties_,starFormationHistory_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="starFormationRateNuclearStarClusters_"/>
    <objectDestructor name="stellarPopulationProperties_"         />
    <objectDestructor name="starFormationHistory_"                />
    !!]
    return
  end function starFormationNuclearStarClustersConstructorParameters

  function starFormationNuclearStarClustersConstructorInternal(luminositiesStellarInactive,starFormationRateNuclearStarClusters_,stellarPopulationProperties_,starFormationHistory_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily starFormationNSC} node operator class.
    !!}
    implicit none
    type   (nodeOperatorStarFormationNuclearStarClusters)                        :: self
    class  (starFormationRateNuclearStarClustersClass   ), intent(in   ), target :: starFormationRateNuclearStarClusters_
    class  (stellarPopulationPropertiesClass            ), intent(in   ), target :: stellarPopulationProperties_
    class  (starFormationHistoryClass                   ), intent(in   ), target :: starFormationHistory_
    logical                                              , intent(in   )         :: luminositiesStellarInactive
    !![
    <constructorAssign variables="luminositiesStellarInactive, *starFormationRateNuclearStarClusters_, *stellarPopulationProperties_, *starFormationHistory_"/>
    !!]
    
    return
  end function starFormationNuclearStarClustersConstructorInternal

  subroutine starFormationNuclearStarClustersDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily starFormationNuclearStarClusters} node operator class.
    !!}
    implicit none
    type(nodeOperatorStarFormationNuclearStarClusters), intent(inout) :: self

    !![
    <objectDestructor name="self%starFormationRateNuclearStarClusters_"/>
    <objectDestructor name="self%stellarPopulationProperties_"         />
    <objectDestructor name="self%starFormationHistory_"                />
    !!]
    return
  end subroutine starFormationNuclearStarClustersDestructor
  
  subroutine starFormationNuclearStarClustersDifferentialEvolutionAnalytics(self,node)
    !!{
    Mark the formed stellar mass as analytically-solvable.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentNSC
    implicit none
    class(nodeOperatorStarFormationNuclearStarClusters), intent(inout) :: self
    type (treeNode                                    ), intent(inout) :: node
    class(nodeComponentNSC                            ), pointer       :: nuclearStarCluster

    ! Mark the formed stellar mass as analytically-solvable (it is always zero) if we are not solving for luminosities as inactive
    ! properties.
    if (.not.self%luminositiesStellarInactive) then
       nuclearStarCluster => node%NSC()
       select type (nuclearStarCluster)
       type is (nodeComponentNSC)
          ! Nuclear star cluster does not yet exist.
       class default
          call nuclearStarCluster%massStellarFormedAnalytic()
       end select
    end if  
    return
  end subroutine starFormationNuclearStarClustersDifferentialEvolutionAnalytics

  subroutine starFormationNuclearStarClustersDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Perform star formation in a nuclear star cluster.
    !!}
    use :: Abundances_Structure          , only : abundances
    use :: Galacticus_Nodes              , only : propertyInactive   , propertyTypeActive, propertyEvaluate, nodeComponentNSC
    use :: Histories                     , only : history
    use :: Stellar_Luminosities_Structure, only : stellarLuminosities
    implicit none
    class           (nodeOperatorStarFormationNuclearStarClusters), intent(inout), target  :: self
    type            (treeNode                                    ), intent(inout), target  :: node
    logical                                                       , intent(inout)          :: interrupt
    procedure       (interruptTask                               ), intent(inout), pointer :: functionInterrupt
    integer                                                       , intent(in   )          :: propertyType
    class           (nodeComponentNSC                            )               , pointer :: nuclearStarCluster
    double precision                                                                       :: rateStarFormation       , massFuel             , &
         &                                                                                    rateMassStellar         , rateEnergyInput      , &
         &                                                                                    rateMassFuel            
    logical                                                                                :: luminositiesCompute
    type            (abundances                                  )                         :: abundancesFuel          , rateAbundancesFuels  , &
         &                                                                                    rateAbundancesStellar
    type            (history                                     )                         :: rateHistoryStarFormation, ratePropertiesStellar
    type            (stellarLuminosities                         )                         :: rateLuminositiesStellar 
    
    ! Check for a realistic nuclear star cluster, return immediately if nuclear star cluster is unphysical.
    nuclearStarCluster => node%NSC()
    if     (    nuclearStarCluster%radius () < 0.0d0 &
         &  .or.nuclearStarCluster%massGas() < 0.0d0 &
         & ) return
    if (propertyInactive(propertyType)) then
       ! For inactive property solution make use of the "massStellarFormed" property to determine the star formation rate.
       rateStarFormation=nuclearStarCluster%massStellarFormedRateGet()
    else
       ! During active property solution, integrate the star formation rate so that we will have a solution for the total mass
       ! of stars formed as a function of time. This differs from the stellar mass due to recycling, and possibly transfer of
       ! stellar mass to other components.
       rateStarFormation=self%starFormationRateNuclearStarClusters_%rate(node)   
       call nuclearStarCluster%massStellarFormedRate(rateStarFormation)
       if (self%luminositiesStellarInactive) call nuclearStarCluster%massStellarFormedRate(rateStarFormation)
    end if
    ! Compute abundances of star forming gas.
    massFuel      =nuclearStarCluster%massGas      ()
    abundancesFuel=nuclearStarCluster%abundancesGas()
    call abundancesFuel%massToMassFraction(massFuel)
    ! Determine if luminosities must be computed.
    luminositiesCompute=propertyEvaluate(propertyType,self%luminositiesStellarInactive)
    ! Find rates of change of stellar mass, gas mass, abundances and luminosities.
    ratePropertiesStellar=nuclearStarCluster%stellarPropertiesHistory()
    call self%stellarPopulationProperties_%rates(                         &
         &                                       rateStarFormation      , &
         &                                       abundancesFuel         , &
         &                                       nuclearStarCluster     , &
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
       rateHistoryStarFormation=nuclearStarCluster%starFormationHistory()
       call        rateHistoryStarFormation%reset                       (                                                              )
       call self  %starFormationHistory_   %                        rate(node,rateHistoryStarFormation,abundancesFuel,rateStarFormation)
       call        nuclearStarCluster      %             massStellarRate(     rateMassStellar                                          )
       call        nuclearStarCluster      %                 massGasRate(     rateMassFuel                                             )
       call        nuclearStarCluster      %       abundancesStellarRate(     rateAbundancesStellar                                    )
       call        nuclearStarCluster      %           abundancesGasRate(     rateAbundancesFuels                                      )
       if (ratePropertiesStellar   %exists())                                                                                            &
            & call nuclearStarCluster      %stellarPropertiesHistoryRate(     ratePropertiesStellar                                    )
       if (rateHistoryStarFormation%exists())                                                                                            &
            & call nuclearStarCluster      %    starFormationHistoryRate(     rateHistoryStarFormation                                 )
    end if
    if (luminositiesCompute                 )                                                                                            &
        & call nuclearStarCluster          %luminositiesStellarRate     (     rateLuminositiesStellar                                  )
    return
  end subroutine starFormationNuclearStarClustersDifferentialEvolution
  
