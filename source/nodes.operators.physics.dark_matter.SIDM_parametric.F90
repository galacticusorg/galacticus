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

  !+    Contributions to this file made by: Niusha Ahvazi

  !!{RST
  Implements a node operator class that maps the CDM solution to SIDM based on the parametric model of :cite:t:`yang_parametric_2024`.
  !!}
  
  use :: Dark_Matter_Halo_Mass_Accretion_Histories, only : darkMatterHaloMassAccretionHistory       , darkMatterHaloMassAccretionHistoryClass
  use :: Cosmology_Functions                      , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters                     , only : cosmologyParametersClass
  use :: Merger_Trees_Build_Mass_Resolution       , only : mergerTreeMassResolutionClass            , mergerTreeMassResolutionFixed
  use :: Merger_Trees_Builders                    , only : mergerTreeBuilderSmoothAccretion
  use :: Dark_Matter_Profile_Scales               , only : darkMatterProfileScaleRadiusConcentration
  use :: Dark_Matter_Halo_Scales                  , only : darkMatterHaloScaleClass
  use :: Virial_Density_Contrast                  , only : virialDensityContrastClass
  use :: Dark_Matter_Profiles_Concentration       , only : darkMatterProfileConcentrationClass
  use :: SIDM_Parametric_Model                    , only : timescaleCollapse                        , velocityMaximumRateTau                  , radiusMaximumRateTau, radiusMaximumNFW, &
       &                                                   velocityMaximumNFW                       , radiusScaleNFW                          , densityScaleNFW     , densityScale    , &
       &                                                   radiusScale                              , radiusCore

  !![
  <nodeOperator name="nodeOperatorSIDMParametric" docformat="rst">
   <description>
   A node operator class that maps the CDM solution to SIDM based on the parametric model of :cite:t:`yang_parametric_2024`.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorSIDMParametric
     !!{RST
     A node operator class that maps the CDM solution to SIDM based on the parametric model of :cite:t:`yang_parametric_2024`.
     !!}
     private
     class(darkMatterParticleClass                ), pointer :: darkMatterParticle_                 => null()
     class(darkMatterHaloMassAccretionHistoryClass), pointer :: darkMatterHaloMassAccretionHistory_ => null()
     class(darkMatterProfileDMOClass              ), pointer :: darkMatterProfileDMO_               => null()
     class(cosmologyFunctionsClass                ), pointer :: cosmologyFunctions_                 => null()
     class(cosmologyParametersClass               ), pointer :: cosmologyParameters_                => null()
     class(darkMatterHaloScaleClass               ), pointer :: darkMatterHaloScale_                => null()
     class(virialDensityContrastClass             ), pointer :: virialDensityContrast_              => null()
     class(darkMatterProfileConcentrationClass    ), pointer :: darkMatterProfileConcentration_     => null()
     integer                                                 :: tauID                                        , velocityMaximumSIDMID  , &
          &                                                     radiusMaximumSIDMID                          , nodeFormationTimeSIDMID, &
          &                                                     densityScaleSIDMID                           , radiusScaleSIDMID      , &
          &                                                     radiusCoreSIDMID
     double precision                                        :: alpha                                        , C
   contains
     final     ::                                        SIDMParametricDestructor
     procedure :: nodeTreeInitialize                  => SIDMParametricNodeTreeInitialize
     procedure :: nodePromote                         => SIDMParametricNodePromote
     procedure :: differentialEvolutionScales         => SIDMParametriCalculateTauDifferentialEvolutionScale
     procedure :: differentialEvolution               => SIDMParametriCalculateTauDifferentialEvolution
     procedure :: differentialEvolutionAnalytics      => SIDMParametriDifferentialVmaxAnalytics
     procedure :: differentialEvolutionSolveAnalytics => SIDMParametriDifferentialVmaxSolveAnalytics
  end type nodeOperatorSIDMParametric
  
  interface nodeOperatorSIDMParametric
     !!{RST
     Constructors for the ``SIDMParametric`` node operator class.
     !!}
     module procedure SIDMParametricConstructorParameters
     module procedure SIDMParametricConstructorInternal
  end interface nodeOperatorSIDMParametric

contains
  
  function SIDMParametricConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the ``SIDMParametric`` node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (nodeOperatorSIDMParametric             )                :: self
    type (inputParameters                        ), intent(inout) :: parameters
    class(darkMatterParticleClass                ), pointer       :: darkMatterParticle_
    class(darkMatterHaloMassAccretionHistoryClass), pointer       :: darkMatterHaloMassAccretionHistory_
    class(darkMatterProfileDMOClass              ), pointer       :: darkMatterProfileDMO_
    class(cosmologyFunctionsClass                ), pointer       :: cosmologyFunctions_
    class(cosmologyParametersClass               ), pointer       :: cosmologyParameters_
    class(darkMatterHaloScaleClass               ), pointer       :: darkMatterHaloScale_
    class(virialDensityContrastClass             ), pointer       :: virialDensityContrast_
    class(darkMatterProfileConcentrationClass    ), pointer       :: darkMatterProfileConcentration_
    double precision                                              :: alpha                              , C

    !![
    <inputParameter docformat="rst">
      <name>alpha</name>
      <defaultValue>2.0d0</defaultValue>
      <description>
      The coefficient :math:`\alpha` of the halo mass-growth term in the gravothermal :math:`\tau` evolution, :math:`\dot\tau = 1/t_\mathrm{c} - \alpha \, (\dot{M}/M) \, \tau`. The default value of :math:`2.0` is the best-fit value found by :cite:t:`raut_extended_2026`.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>C</name>
      <defaultValue>0.75d0</defaultValue>
      <description>
      The calibration constant :math:`C` relating the gravothermal collapse timescale :math:`t_\mathrm{c}` to the relaxation time (eqn. 2.2 of :cite:t:`yang_parametric_2024`). The default value of :math:`0.75` is the best-fit value found by :cite:t:`raut_extended_2026`.
      </description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="darkMatterParticle"                 name="darkMatterParticle_"                 source="parameters"/>
    <objectBuilder class="darkMatterHaloMassAccretionHistory" name="darkMatterHaloMassAccretionHistory_" source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO"               name="darkMatterProfileDMO_"               source="parameters"/>
    <objectBuilder class="cosmologyFunctions"                 name="cosmologyFunctions_"                 source="parameters"/>
    <objectBuilder class="cosmologyParameters"                name="cosmologyParameters_"                source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"                name="darkMatterHaloScale_"                source="parameters"/>
    <objectBuilder class="virialDensityContrast"              name="virialDensityContrast_"              source="parameters"/>
    <objectBuilder class="darkMatterProfileConcentration"     name="darkMatterProfileConcentration_"     source="parameters"/>
    !!]
    self=nodeOperatorSIDMParametric(alpha,C,darkMatterParticle_,darkMatterHaloMassAccretionHistory_,darkMatterProfileDMO_,cosmologyFunctions_,cosmologyParameters_,darkMatterHaloScale_,virialDensityContrast_,darkMatterProfileConcentration_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterParticle_"                />
    <objectDestructor name="darkMatterHaloMassAccretionHistory_"/>
    <objectDestructor name="darkMatterProfileDMO_"              />
    <objectDestructor name="cosmologyFunctions_"                />
    <objectDestructor name="cosmologyParameters_"               />
    <objectDestructor name="darkMatterHaloScale_"               />
    <objectDestructor name="virialDensityContrast_"             />
    <objectDestructor name="darkMatterProfileConcentration_"    />
    !!]
    return
  end function SIDMParametricConstructorParameters

  function SIDMParametricConstructorInternal(alpha,C,darkMatterParticle_,darkMatterHaloMassAccretionHistory_,darkMatterProfileDMO_,cosmologyFunctions_,cosmologyParameters_,darkMatterHaloScale_,virialDensityContrast_,darkMatterProfileConcentration_) result(self)
    !!{RST
    Internal constructor for the ``SIDMParametric`` node operator class.
    !!}

    implicit none
    type            (nodeOperatorSIDMParametric             )                        :: self
    double precision                                         , intent(in   )         :: alpha                              , C
    class           (darkMatterParticleClass                ), intent(in   ), target :: darkMatterParticle_
    class           (darkMatterHaloMassAccretionHistoryClass), intent(in   ), target :: darkMatterHaloMassAccretionHistory_
    class           (darkMatterProfileDMOClass              ), intent(in   ), target :: darkMatterProfileDMO_
    class           (cosmologyFunctionsClass                ), intent(in   ), target :: cosmologyFunctions_
    class           (cosmologyParametersClass               ), intent(in   ), target :: cosmologyParameters_
    class           (darkMatterHaloScaleClass               ), intent(in   ), target :: darkMatterHaloScale_
    class           (virialDensityContrastClass             ), intent(in   ), target :: virialDensityContrast_
    class           (darkMatterProfileConcentrationClass    ), intent(in   ), target :: darkMatterProfileConcentration_

    !![
    <constructorAssign variables="alpha, C, *darkMatterParticle_, *darkMatterHaloMassAccretionHistory_, *darkMatterProfileDMO_, *cosmologyFunctions_, *cosmologyParameters_, *darkMatterHaloScale_, *virialDensityContrast_, *darkMatterProfileConcentration_"/>
    !!]

    !![
    <addMetaProperty component="darkMatterProfile" name="tau"                   id="self%tauID"                   isEvolvable="yes" isCreator="yes"/>
    <addMetaProperty component="darkMatterProfile" name="velocityMaximumSIDM"   id="self%velocityMaximumSIDMID"   isEvolvable="yes" isCreator="yes"/>
    <addMetaProperty component="darkMatterProfile" name="radiusMaximumSIDM"     id="self%radiusMaximumSIDMID"     isEvolvable="yes" isCreator="yes"/>
    <addMetaProperty component="basic"             name="nodeFormationTimeSIDM" id="self%nodeFormationTimeSIDMID" isEvolvable="no"  isCreator="yes"/>
    <addMetaProperty component="darkMatterProfile" name="densityScaleSIDM"      id="self%densityScaleSIDMID"      isEvolvable="yes" isCreator="yes"/>
    <addMetaProperty component="darkMatterProfile" name="radiusScaleSIDM"       id="self%radiusScaleSIDMID"       isEvolvable="yes" isCreator="yes"/>
    <addMetaProperty component="darkMatterProfile" name="radiusCoreSIDM"        id="self%radiusCoreSIDMID"        isEvolvable="yes" isCreator="yes"/>
    !!]    
    return
  end function SIDMParametricConstructorInternal

  subroutine SIDMParametricDestructor(self)
    !!{RST
    Destructor for the ``SIDMParametric`` node operator class.
    !!}
    implicit none
    type(nodeOperatorSIDMParametric), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterParticle_"                />
    <objectDestructor name="self%darkMatterHaloMassAccretionHistory_"/>
    <objectDestructor name="self%darkMatterProfileDMO_"              />
    <objectDestructor name="self%cosmologyFunctions_"                />
    <objectDestructor name="self%cosmologyParameters_"               />
    <objectDestructor name="self%darkMatterHaloScale_"               />
    <objectDestructor name="self%virialDensityContrast_"             />
    <objectDestructor name="self%darkMatterProfileConcentration_"    />
    !!]
    return
  end subroutine SIDMParametricDestructor

  subroutine SIDMParametricNodeTreeInitialize(self,node)
    !!{RST
    Initialize the SIDMParametric of all nodes in the tree.
    !!}

    use :: Galacticus_Nodes                , only : mergerTree                     , nodeComponentBasic         , nodeComponentDarkMatterProfile
    use :: Dark_Matter_Halo_Formation_Times, only : Dark_Matter_Halo_Formation_Time
    use :: Mass_Distributions              , only : massDistributionClass
    use :: Galactic_Structure_Options      , only : componentTypeDarkHalo          , componentTypeDarkMAtterOnly, massTypeDark                  , weightByMass
    implicit none
    class           (nodeOperatorSIDMParametric               ), intent(inout), target  :: self
    class           (massDistributionClass                    ),                pointer :: massDistribution_        
    type            (treeNode                                 ), intent(inout), target  :: node
    type            (treeNode                                 )               , pointer :: nodeParent                              , nodeChild             , &
         &                                                                                 nodeNew
    type            (mergerTree                               )                         :: treeNew
    class           (nodeComponentBasic                       )               , pointer :: basic                                   , basicParent           , &
         &                                                                                 basicNode                               , basicChild            , &
         &                                                                                 basicNew                                , basicBase             , &
         &                                                                                 basicNewChild
    class           (nodeComponentDarkMatterProfile           )               , pointer :: darkMatterProfile                       , darkMatterProfileChild, &
         &                                                                                 darkMatterProfileCopy
    type            (mergerTreeMassResolutionFixed            )               , pointer :: mergerTreeMassResolutionFixed_
    type            (mergerTreeBuilderSmoothAccretion         )               , pointer :: mergerTreeBuilderSmoothAccretion_
    type            (darkMatterProfileScaleRadiusConcentration)               , pointer :: darkMatterProfileScaleRadius_
    double precision                                           , parameter              :: formationMassFraction            =0.50d0
    double precision                                           , parameter              :: formationTimeFraction            =0.50d0
    ! Fractional halo mass decline per step used when extrapolating the tree below resolution.
    double precision                                           , parameter              :: massHaloDeclineFactor            =0.99d0
    double precision                                                                    :: timeFormation
    double precision                                                                    :: velocityMaximumSIDMPrevious             , timescaleCollapse_    , &
         &                                                                                 velocityMaximumSIDM                     , radiusMaximumSIDM     , &
         &                                                                                 radiusMaximumCDM                        , radiusMaximumNFW_     , &
         &                                                                                 velocityMaximumNFW_                     , radiusScaleNFW_       , &
         &                                                                                 densityScaleNFW_                        , densityScale_         , &
         &                                                                                 radiusScale_                            , radiusCore_           , &
         &                                                                                 velocityMaximumCDM                      , rateMassSpecific      , &
         &                                                                                 timeStep                                , tauPrevious           , &
         &                                                                                 massHalo                                , rateMass              , &
         &                                                                                 timeEarliest                            , massResolution        , &
         &                                                                                 tau

    call self%nodeInitialize(node)
    nodeParent => node
    basic      => node%basic()
    do while (associated(nodeParent))
       basicParent   => nodeParent%basic()
       timeFormation =  Dark_Matter_Halo_Formation_Time(                                                                              &
            &                                           node                               =nodeParent                              , &
            &                                           formationMassFraction              =formationMassFraction                   , &
            &                                           darkMatterHaloMassAccretionHistory_=self%darkMatterHaloMassAccretionHistory_  &
            &                                          )
       if (nodeParent%isPrimaryProgenitor()) then
          nodeParent => nodeParent%parent
       else
          nodeParent => null()
       end if
    end do
    call basic%floatRank0MetaPropertySet(self%nodeFormationTimeSIDMID,timeFormation)
    ! Initialize each node as an NFW profile.
    darkMatterProfile => node             %darkMatterProfile            (autoCreate=     .true.)
    massDistribution_ => self             %darkMatterProfileDMO_    %get(                node  )
    tau               =darkMatterProfile%floatRank0MetaPropertyGet(self%tauID)
    radiusMaximumCDM  =massDistribution_%radiusRotationCurveMaximum  ()
    velocityMaximumCDM=massDistribution_%velocityRotationCurveMaximum()
    radiusScaleNFW_   =radiusScaleNFW (                massDistribution_%radiusRotationCurveMaximum  ())
    densityScaleNFW_  =densityScaleNFW(radiusScaleNFW_,massDistribution_%velocityRotationCurveMaximum())
    densityScale_     =densityScale(densityScaleNFW_,tau)
    radiusScale_      =radiusScale(radiusScaleNFW_,tau)
    radiusCore_       =radiusCore(radiusScaleNFW_,tau)
    call darkMatterProfile%floatRank0MetaPropertySet(self%densityScaleSIDMID,densityScale_)
    call darkMatterProfile%floatRank0MetaPropertySet(self%radiusScaleSIDMID ,radiusScale_ )
    call darkMatterProfile%floatRank0MetaPropertySet(self%radiusCoreSIDMID  ,radiusCore_  )
    ! For leaf nodes, grow a sub-resolution merger tree to allow us to determine the formation time if necessary.
    if (.not.associated(node%firstChild)) then
       if (basic%time() > basic%floatRank0MetaPropertyGet(self%nodeFormationTimeSIDMID)) then
          ! If this node is at the tip of the branch and the time associated with it is larger than the formation time of the
          ! branch then we need to build the rest of the tree (extrapolate back in time) to start the parametric SIDM
          ! calculation from the formation time. We build the tree with smooth accretion only. 
          !! The mass resolution is set to be the current node mass times the formation mass fraction.          
          !! The earliest allowed time is set to a fraction of the formation time to ensure that the full history for which SIDM
          !! evolution is needed is captured.
          massResolution=+formationMassFraction*basic%mass                     (                            )
          timeEarliest  =+formationTimeFraction*basic%floatRank0MetaPropertyGet(self%nodeFormationTimeSIDMID)
          ! Construct objects needed to build the merger tree.
          allocate(mergerTreeMassResolutionFixed_   )
          allocate(mergerTreeBuilderSmoothAccretion_)
          allocate(darkMatterProfileScaleRadius_    )
          !![
          <referenceConstruct object="mergerTreeMassResolutionFixed_"   >
	    <constructor>
	      mergerTreeMassResolutionFixed            (                                                                             &amp;
	       &amp;                                    massResolution                     =     massResolution                      &amp;
	       &amp;                                   )
     	    </constructor>
	  </referenceConstruct>
	  <referenceConstruct object="mergerTreeBuilderSmoothAccretion_">
	    <constructor>
	      mergerTreeBuilderSmoothAccretion         (                                                                              &amp;
	       &amp;                                    massHaloDeclineFactor              =     massHaloDeclineFactor              , &amp;
	       &amp;                                    timeEarliest                       =     timeEarliest                       , &amp;
	       &amp;                                    cosmologyFunctions_                =self%cosmologyFunctions_                , &amp;
	       &amp;                                    darkMatterHaloMassAccretionHistory_=self%darkMatterHaloMassAccretionHistory_, &amp;
	       &amp;                                    mergerTreeMassResolution_          =     mergerTreeMassResolutionFixed_       &amp;
	       &amp;                                   )
     	    </constructor>
	  </referenceConstruct>
          <referenceConstruct object="darkMatterProfileScaleRadius_"    >
	    <constructor>
	      darkMatterProfileScaleRadiusConcentration(                                                                              &amp;
	       &amp;                                    correctForConcentrationDefinition  =.false.                                 , &amp;
	       &amp;                                    useMeanConcentration               =.true.                                  , &amp;
	       &amp;                                    cosmologyParameters_               =self%cosmologyParameters_               , &amp;
	       &amp;                                    cosmologyFunctions_                =self%cosmologyFunctions_                , &amp;
	       &amp;                                    darkMatterHaloScale_               =self%darkMatterHaloScale_               , &amp;
	       &amp;                                    darkMatterProfileDMO_              =self%darkMatterProfileDMO_              , &amp;
	       &amp;                                    virialDensityContrast_             =self%virialDensityContrast_             , &amp;
	       &amp;                                    darkMatterProfileConcentration_    =self%darkMatterProfileConcentration_      &amp;
	       &amp;                                   )
     	    </constructor>
	  </referenceConstruct>
          !!]
          ! Build the tree.
          treeNew%randomNumberGenerator_ => node%hostTree%randomNumberGenerator_
          allocate(treeNew%nodeBase)
          call node                             %copyNodeTo(treeNew%nodeBase)
          call mergerTreeBuilderSmoothAccretion_%build     (treeNew         )
          nodeNew   => null()
          nodeChild => treeNew%nodeBase
          do while (associated(nodeChild))
             basicChild        => nodeChild%basic            (                 )
             darkMatterProfile => nodeChild%darkMatterProfile(autoCreate=.true.)
             call darkMatterProfile%scaleSet(darkMatterProfileScaleRadius_%radius(nodeChild))
             if (basicChild%time() > timeFormation) then
                nodeNew   => nodeChild
                nodeChild => nodeChild%firstChild
             else
                nodeChild => null()
             end if
          end do
          ! Perform SIDM evolution along the subresolution tree.
          basicBase => treeNew%nodeBase%basic()
          do while (associated(nodeNew))
             basicNew                    =>  nodeNew                      %basic            (       )
             darkMatterProfile           =>  nodeNew                      %darkMatterProfile(       )
             basicNewChild               =>  nodeNew%firstChild           %basic            (       )
             darkMatterProfileChild      =>  nodeNew%firstChild           %darkMatterProfile(       )   
             massDistribution_           =>  self   %darkMatterProfileDMO_%get              (nodeNew)
             velocityMaximumSIDMPrevious =  +darkMatterProfileChild%floatRank0MetaPropertyGet   (self%velocityMaximumSIDMID) &
                  &                         +massDistribution_     %velocityRotationCurveMaximum(                          )
             ! Compute the evolution for this step.
             timeStep          =basicNew%time()-basicNewChild%time()
             timescaleCollapse_=timescaleCollapse(self%darkMatterParticle_,self%C,massDistribution_%velocityRotationCurveMaximum(),massDistribution_%radiusRotationCurveMaximum(),velocityMaximumSIDMPrevious)
             massHalo          =basicNewChild%mass()
             rateMass          =basicNewChild%accretionRate()
             rateMassSpecific  =rateMass/massHalo
             tauPrevious       =darkMatterProfileChild%floatRank0MetaPropertyGet(self%tauID)
             tau               =tauPrevious+timeStep*(1.d0/timescaleCollapse_-self%alpha*rateMassSpecific*tauPrevious) 
             ! Store the updated τ and profile structure parameters.
             call darkMatterProfile%floatRank0MetaPropertySet(self%tauID,tau)
             velocityMaximumSIDM=darkMatterProfileChild%floatRank0MetaPropertyGet(self%velocityMaximumSIDMID)+velocityMaximumRateTau(tau,massDistribution_%velocityRotationCurveMaximum())*timeStep*(1.d0/timescaleCollapse_-self%alpha*rateMassSpecific*tauPrevious)
             radiusMaximumSIDM  =darkMatterProfileChild%floatRank0MetaPropertyGet(self%radiusMaximumSIDMID  )+radiusMaximumRateTau  (tau,massDistribution_%radiusRotationCurveMaximum  ())*timeStep*(1.d0/timescaleCollapse_-self%alpha*rateMassSpecific*tauPrevious)
             call darkMatterProfile%floatRank0MetaPropertySet(self%velocityMaximumSIDMID,velocityMaximumSIDM)
             call darkMatterProfile%floatRank0MetaPropertySet(self%radiusMaximumSIDMID  ,radiusMaximumSIDM  )
             radiusMaximumCDM   =massDistribution_%radiusRotationCurveMaximum()
             radiusMaximumSIDM  =radiusMaximumSIDM+radiusMaximumCDM
             radiusMaximumNFW_  =radiusMaximumNFW  (radiusMaximumSIDM  ,tau                )
             velocityMaximumNFW_=velocityMaximumNFW(velocityMaximumSIDM,tau                )
             radiusScaleNFW_    =radiusScaleNFW    (                    radiusMaximumNFW_  )
             densityScaleNFW_   =densityScaleNFW   (radiusScaleNFW_    ,velocityMaximumNFW_)
             densityScale_      =densityScale      (densityScaleNFW_   ,tau                )
             radiusScale_       =radiusScale       (radiusScaleNFW_    ,tau                )
             radiusCore_        =radiusCore        (radiusScaleNFW_    ,tau                )
             call darkMatterProfile%floatRank0MetaPropertySet(self%densityScaleSIDMID,densityScale_)
             call darkMatterProfile%floatRank0MetaPropertySet(self%radiusScaleSIDMID ,radiusScale_ )
             call darkMatterProfile%floatRank0MetaPropertySet(self%radiusCoreSIDMID  ,radiusCore_  )
             ! Copy the values from the new tree to the tip of the original branch.
             darkMatterProfile     => node            %darkMatterProfile()
             darkMatterProfileCopy => treeNew%nodeBase%darkMatterProfile()
             call darkMatterProfile%floatRank0MetaPropertySet(self%tauID                ,darkMatterProfileCopy%floatRank0MetaPropertyGet(self%tauID                ))
             call darkMatterProfile%floatRank0MetaPropertySet(self%velocityMaximumSIDMID,darkMatterProfileCopy%floatRank0MetaPropertyGet(self%velocityMaximumSIDMID))
             call darkMatterProfile%floatRank0MetaPropertySet(self%radiusMaximumSIDMID  ,darkMatterProfileCopy%floatRank0MetaPropertyGet(self%radiusMaximumSIDMID  ))
             call darkMatterProfile%floatRank0MetaPropertySet(self%densityScaleSIDMID   ,darkMatterProfileCopy%floatRank0MetaPropertyGet(self%densityScaleSIDMID   ))
             call darkMatterProfile%floatRank0MetaPropertySet(self%radiusScaleSIDMID    ,darkMatterProfileCopy%floatRank0MetaPropertyGet(self%radiusScaleSIDMID    ))
             call darkMatterProfile%floatRank0MetaPropertySet(self%radiusCoreSIDMID     ,darkMatterProfileCopy%floatRank0MetaPropertyGet(self%radiusCoreSIDMID     ))
             ! Move to the next node unless we have reached the original leaf node.
             if (.not.nodeNew%index() == treeNew%nodeBase%index()) then
                nodeNew => nodeNew%parent
             else
                nodeNew => null()
             end if
          end do
          ! Destroy the tree and builder objects.
          call treeNew%nodeBase%destroyBranch()
          deallocate(treeNew%nodeBase)
          !![
          <objectDestructor name="mergerTreeMassResolutionFixed_"   />
          <objectDestructor name="mergerTreeBuilderSmoothAccretion_"/>
          <objectDestructor name="darkMatterProfileScaleRadius_"    />
          <objectDestructor name="massDistribution_"                />
          !!]
       end if
    end if
    return
  end subroutine SIDMParametricNodeTreeInitialize

  subroutine SIDMParametriDifferentialVmaxAnalytics(self,node)
    !!{RST
    Mark analytically-solvable properties.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentDarkMatterProfile
    implicit none
    class(nodeOperatorSIDMParametric    ), intent(inout) :: self
    type (treeNode                      ), intent(inout) :: node
    class(nodeComponentDarkMatterProfile), pointer       :: darkMatterProfile

    darkMatterProfile => node%darkMatterProfile()        
    call darkMatterProfile%floatRank0MetaPropertyAnalytic(self%densityScaleSIDMID)
    call darkMatterProfile%floatRank0MetaPropertyAnalytic(self%radiusScaleSIDMID  )
    call darkMatterProfile%floatRank0MetaPropertyAnalytic(self%radiusCoreSIDMID  )
    return
  end subroutine SIDMParametriDifferentialVmaxAnalytics

  subroutine SIDMParametriDifferentialVmaxSolveAnalytics(self,node,time)
    !!{RST
    Evolve ":term:`dark matter-only universe`" mass at a constant rate, to achieve linear interpolation in time.
    !!}
    use :: Galacticus_Nodes  , only : nodeComponentBasic   , nodeComponentDarkMatterProfile
    use :: Mass_Distributions, only : massDistributionClass
    implicit none
    class           (nodeOperatorSIDMParametric  ), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    double precision                              , intent(in   ) :: time
    class           (nodeComponentBasic          ), pointer       :: basic              , basicParent
    class         (nodeComponentDarkMatterProfile), pointer       :: darkMatterProfile
    class           (massDistributionClass       ), pointer       :: massDistribution_
    double precision                                              :: tau                , radiusScaleNFW_    , &
         &                                                           densityScaleNFW_   , densityScale_      , &
         &                                                           radiusScale_       , radiusCore_        , &
         &                                                           velocityMaximumSIDM, radiusMaximumSIDM  , &
         &                                                           radiusMaximumNFW_  , velocityMaximumNFW_

    darkMatterProfile => node             %darkMatterProfile            (    )
    massDistribution_ => self             %darkMatterProfileDMO_    %get(node)
    tau                = darkMatterProfile%floatRank0MetaPropertyGet(self%tauID                )
    velocityMaximumSIDM= darkMatterProfile%floatRank0MetaPropertyGet(self%velocityMaximumSIDMID)+massDistribution_%velocityRotationCurveMaximum()
    radiusMaximumSIDM  = darkMatterProfile%floatRank0MetaPropertyGet(self%radiusMaximumSIDMID  )+massDistribution_%radiusRotationCurveMaximum  ()
    radiusMaximumNFW_  = radiusMaximumNFW  (radiusMaximumSIDM,tau)
    velocityMaximumNFW_= velocityMaximumNFW(velocityMaximumSIDM,tau                )
    radiusScaleNFW_    = radiusScaleNFW    (radiusMaximumNFW_                      )
    densityScaleNFW_   = densityScaleNFW   (radiusScaleNFW_    ,velocityMaximumNFW_)
    densityScale_      = densityScale      (densityScaleNFW_   ,tau                )
    radiusScale_       = radiusScale       (radiusScaleNFW_    ,tau                )
    radiusCore_        = radiusCore        (radiusScaleNFW_    ,tau                )
    call darkMatterProfile%floatRank0MetaPropertySet(self%densityScaleSIDMID,densityScale_)
    call darkMatterProfile%floatRank0MetaPropertySet(self%radiusScaleSIDMID ,radiusScale_ )
    call darkMatterProfile%floatRank0MetaPropertySet(self%radiusCoreSIDMID  ,radiusCore_  )
    return
  end subroutine SIDMParametriDifferentialVmaxSolveAnalytics

  subroutine SIDMParametricNodePromote(self,node)
    !!{RST
    Ensure that ``node`` is ready for promotion to its parent.
    !!}
    use :: Error             , only : Error_Report
    use :: Galacticus_Nodes  , only : nodeComponentBasic, nodeComponentDarkMatterProfile
    use :: ISO_Varying_String, only : var_str           , varying_string                , operator(//)
    use :: String_Handling   , only : operator(//)
    implicit none
    class    (nodeOperatorSIDMParametric    ), intent(inout) :: self
    type     (treeNode                      ), intent(inout) :: node
    type     (treeNode                      ), pointer       :: nodeParent
    class    (nodeComponentBasic            ), pointer       :: basicParent      , basic
    class    (nodeComponentDarkMatterProfile), pointer       :: darkMatterProfile, darkMatterProfileParent
    type     (varying_string                )                :: message
    character(len=12                        )                :: label
    !$GLC attributes unused :: self
    
    nodeParent              => node      %parent
    basic                   => node      %basic            ()
    basicParent             => nodeParent%basic            ()
    darkMatterProfile       => node      %darkMatterProfile()    
    darkMatterProfileParent => nodeParent%darkMatterProfile()
    ! Ensure the two halos exist at the same time.
    if (basic%time() /= basicParent%time()) then
       message=var_str("node [")//node%index()//"] has not been evolved to its parent ["//nodeParent%index()//"]"//char(10)
       write (label,'(f12.6)') basic%time()
       message=message//"    node is at time: "//label//" Gyr"//char(10)
       write (label,'(f12.6)') basicParent%time()
       message=message//"  parent is at time: "//label//" Gyr"
       call Error_Report(message//{introspection:location})
    end if
    return
  end subroutine SIDMParametricNodePromote

  subroutine SIDMParametriCalculateTauDifferentialEvolutionScale(self, node)
    !!{RST
    Set the ODE solver scales for the :math:`\tau` parameter.
    !!}
    use Galacticus_Nodes, only: treeNode, nodeComponentBasic, nodeComponentDarkMatterProfile
    implicit none
    type (treeNode                      ), intent(inout) :: node
    class(nodeOperatorSIDMParametric    ), intent(inout) :: self
    class(nodeComponentDarkMatterProfile), pointer       :: darkMatterProfile

    darkMatterProfile => node%darkMatterProfile()
    call darkMatterProfile%floatRank0MetaPropertyScale(self%tauID,1.0d0)
    return
  end subroutine SIDMParametriCalculateTauDifferentialEvolutionScale

  subroutine SIDMParametriCalculateTauDifferentialEvolution(self, node,interrupt,functionInterrupt,propertyType)
    !!{RST
    Perform differential evolution of the parameters of the SIDM parametric model.
    !!}
    use :: Galacticus_Nodes                , only : treeNode                       , nodeComponentBasic         , nodeComponentDarkMatterProfile
    use :: Dark_Matter_Halo_Formation_Times, only : Dark_Matter_Halo_Formation_Time
    use :: Mass_Distributions              , only : massDistributionClass
    use :: Galactic_Structure_Options      , only : componentTypeDarkHalo          , componentTypeDarkMAtterOnly, massTypeDark                  , weightByMass
    implicit none
    type            (treeNode                      ), intent(inout), target  :: node
    logical                                         , intent(inout)          :: interrupt
    procedure       (interruptTask                 ), intent(inout), pointer :: functionInterrupt
    integer                                         , intent(in   )          :: propertyType
    class           (nodeOperatorSIDMParametric    ), intent(inout), target  :: self
    class           (massDistributionClass         )               , pointer :: massDistribution_
    class           (nodeComponentBasic            )               , pointer :: basic              , basicParent
    class           (nodeComponentDarkMatterProfile)               , pointer :: darkMatterProfile
    double precision                                                         :: timeFormation      , timescaleCollapse_, &
         &                                                                      tau                , time              , &
         &                                                                      velocityMaximumSIDM, radiusMaximumSIDM , &
         &                                                                      rateVelocity       , rateRadius        , &
         &                                                                      massHalo           , rateMass          , &
         &                                                                      rateMassSpecific   , rateTau
    !$GLC attributes unused :: interrupt, functionInterrupt, propertyType

    basic             => node %basic                    (                            )
    darkMatterProfile => node %darkMatterProfile        (                            )
    time              =  basic%time                     (                            )
    timeFormation     =  basic%floatRank0MetaPropertyGet(self%nodeFormationTimeSIDMID)
    if (time > timeFormation) then
       massDistribution_   => self%darkMatterProfileDMO_%get(node)
       velocityMaximumSIDM =  darkMatterProfile%floatRank0MetaPropertyGet(self%velocityMaximumSIDMID)+massDistribution_%velocityRotationCurveMaximum()
       timescaleCollapse_  =  timescaleCollapse(self%darkMatterParticle_,self%C,massDistribution_%velocityRotationCurveMaximum(),massDistribution_%radiusRotationCurveMaximum(),velocityMaximumSIDM)
       tau                 =  darkMatterProfile%floatRank0MetaPropertyGet(self%tauID)
       massHalo            =  basic%mass         ()
       rateMass            =  basic%accretionRate()
       rateMassSpecific    =  rateMass/massHalo
       rateTau             =  1.0d0/timescaleCollapse_-self%alpha*rateMassSpecific*tau
       call darkMatterProfile%floatRank0MetaPropertyRate(self%tauID,rateTau)
       tau         =darkMatterProfile%floatRank0MetaPropertyGet(self%tauID)
       rateVelocity=velocityMaximumRateTau(tau,massDistribution_%velocityRotationCurveMaximum())*rateTau
       rateRadius  =radiusMaximumRateTau  (tau,massDistribution_%radiusRotationCurveMaximum  ())*rateTau
       call darkMatterProfile%floatRank0MetaPropertyRate(self%velocityMaximumSIDMID,rateVelocity) 
       call darkMatterProfile%floatRank0MetaPropertyRate(self%radiusMaximumSIDMID  ,rateRadius  )
       velocityMaximumSIDM=darkMatterProfile%floatRank0MetaPropertyGet(self%velocityMaximumSIDMID)+massDistribution_%velocityRotationCurveMaximum()
       radiusMaximumSIDM  =darkMatterProfile%floatRank0MetaPropertyGet(self%radiusMaximumSIDMID  )+massDistribution_%radiusRotationCurveMaximum  ()
    else
       call darkMatterProfile%floatRank0MetaPropertyRate(self%tauID                ,0.0d0)
       call darkMatterProfile%floatRank0MetaPropertyRate(self%velocityMaximumSIDMID,0.0d0)
       call darkMatterProfile%floatRank0MetaPropertyRate(self%radiusMaximumSIDMID  ,0.0d0)       
    end if    
    return
  end subroutine SIDMParametriCalculateTauDifferentialEvolution
