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

  !!{
  Implements a node operator class that maps the CDM solution to SIDM based on the parametric model of \cite{yang_parametric_2024}.
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
  use :: SIDM_Parametric_Model                    , only : get_tc                                   , dvmaxt   , drmaxt   , Rmax_NFW, Vmax_NFW, &
       &                                                    r_s0                                     , rho_s0   , get_rho_s, get_r_s , get_r_c

  !![
  <nodeOperator name="nodeOperatorSIDMParametric">
   <description>
     A node operator class that maps the CDM solution to SIDM based on the parametric model of \cite{yang_parametric_2024}.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorSIDMParametric
     !!{
     A node operator class that maps the CDM solution to SIDM based on the parametric model of \cite{yang_parametric_2024}.
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
     integer                                                 :: tauID                                        , VmaxSIDMID             , &
          &                                                     RmaxSIDMID                                   , nodeFormationTimeSIDMID, &
          &                                                     RhosSIDMID                                   , RsSIDMID               , &
          &                                                     RcSIDMID
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
     !!{
     Constructors for the \mono{SIDMParametric} node operator class.
     !!}
     module procedure SIDMParametricConstructorParameters
     module procedure SIDMParametricConstructorInternal
  end interface nodeOperatorSIDMParametric

contains
  
  function SIDMParametricConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \mono{SIDMParametric} node operator class which takes a parameter set as input.
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
    <inputParameter>
      <name>alpha</name>
      <defaultValue>2.0d0</defaultValue>
      <description>The coefficient $\alpha$ of the halo mass-growth term in the gravothermal $\tau$ evolution, $\dot\tau = 1/t_\mathrm{c} - \alpha \, (\dot{M}/M) \, \tau$. The default value of $2.0$ is the best-fit value found by \cite{raut_extended_2026}.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>C</name>
      <defaultValue>0.75d0</defaultValue>
      <description>The calibration constant $C$ relating the gravothermal collapse timescale $t_\mathrm{c}$ to the relaxation time (eqn.~2.2 of \cite{yang_parametric_2024}). The default value of $0.75$ is the best-fit value found by \cite{raut_extended_2026}.</description>
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
    !!{
    Internal constructor for the \mono{SIDMParametric} node operator class.
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
    <addMetaProperty component="darkMatterProfile" name="VmaxSIDM"              id="self%VmaxSIDMID"              isEvolvable="yes" isCreator="yes"/>
    <addMetaProperty component="darkMatterProfile" name="RmaxSIDM"              id="self%RmaxSIDMID"              isEvolvable="yes" isCreator="yes"/>
    <addMetaProperty component="basic"             name="nodeFormationTimeSIDM" id="self%nodeFormationTimeSIDMID" isEvolvable="no"  isCreator="yes"/>
    <addMetaProperty component="darkMatterProfile" name="RhosSIDM"              id="self%RhosSIDMID"              isEvolvable="yes" isCreator="yes"/>
    <addMetaProperty component="darkMatterProfile" name="RsSIDM"                id="self%RsSIDMID"                isEvolvable="yes" isCreator="yes"/>
    <addMetaProperty component="darkMatterProfile" name="RcSIDM"                id="self%RcSIDMID"                isEvolvable="yes" isCreator="yes"/>
    !!]    
    return
  end function SIDMParametricConstructorInternal

  subroutine SIDMParametricDestructor(self)
    !!{
    Destructor for the \mono{SIDMParametric} node operator class.
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
    !!{
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
    double precision                                                                    :: VmaxSIDMPrevious                        , tc                    , &
         &                                                                                 VmaxSIDM                                , RmaxSIDM              , &
         &                                                                                 RmaxCDM                                 , RmaxNFW0              , &
         &                                                                                 VmaxNFW0                                , r_sNFW0               , &
         &                                                                                 rho_sNFW0                               , rho_s                 , &
         &                                                                                 r_s                                     , r_c                   , &
         &                                                                                 VmaxCDM                                 , Gamma                 , &
         &                                                                                 dt                                      , tau_old               , &
         &                                                                                 Mhalo                                   , dMdt                  , &
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
    tau               =  darkMatterProfile%floatRank0MetaPropertyGet    (           self%tauID )
    RmaxCDM           =  massDistribution_%radiusRotationCurveMaximum   (                      )
    VmaxCDM           =  massDistribution_%velocityRotationCurveMaximum (                      )
    r_sNFW0           =  r_s0     (        massDistribution_%radiusRotationCurveMaximum  ())
    rho_sNFW0         =  rho_s0   (r_sNFW0,massDistribution_%velocityRotationCurveMaximum())
    rho_s             =  get_rho_s(rho_sNFW0,tau)
    r_s               =  get_r_s  (r_sNFW0  ,tau)
    r_c               =  get_r_c  (r_sNFW0  ,tau)
    call darkMatterProfile%floatRank0MetaPropertySet(self%RhosSIDMID,rho_s)
    call darkMatterProfile%floatRank0MetaPropertySet(self%RsSIDMID  ,r_s  )
    call darkMatterProfile%floatRank0MetaPropertySet(self%RcSIDMID  ,r_c  )
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
             basicNew               =>  nodeNew                      %basic            (       )
             darkMatterProfile      =>  nodeNew                      %darkMatterProfile(       )
             basicNewChild          =>  nodeNew%firstChild           %basic            (       )
             darkMatterProfileChild =>  nodeNew%firstChild           %darkMatterProfile(       )   
             massDistribution_      =>  self   %darkMatterProfileDMO_%get              (nodeNew)
             VmaxSIDMPrevious       =  +darkMatterProfileChild%floatRank0MetaPropertyGet   (self%VmaxSIDMID) &
                  &                    +massDistribution_     %velocityRotationCurveMaximum(               )
             ! Compute the evolution for this step.
             dt     =basicNew%time()-basicNewChild%time()
             tc     =get_tc(self%darkMatterParticle_,self%C,massDistribution_%velocityRotationCurveMaximum(),massDistribution_%radiusRotationCurveMaximum(),VmaxSIDMPrevious)
             Mhalo  =basicNewChild%mass         ()
             dMdt   =basicNewChild%accretionRate()
             Gamma  =dMdt/Mhalo
             tau_old=darkMatterProfileChild%floatRank0MetaPropertyGet(self%tauID)
             tau    =tau_old+dt*(1.d0/tc-self%alpha*Gamma*tau_old) 
             ! Store the updated τ and profile structure parameters.
             call darkMatterProfile%floatRank0MetaPropertySet(self%tauID,tau)
             VmaxSIDM=darkMatterProfileChild%floatRank0MetaPropertyGet(self%VmaxSIDMID)+dvmaxt(tau,massDistribution_%velocityRotationCurveMaximum())*dt*(1.d0/tc-self%alpha*Gamma*tau_old) 
             RmaxSIDM=darkMatterProfileChild%floatRank0MetaPropertyGet(self%RmaxSIDMID)+drmaxt(tau,massDistribution_%radiusRotationCurveMaximum  ())*dt*(1.d0/tc-self%alpha*Gamma*tau_old)
             call darkMatterProfile%floatRank0MetaPropertySet(self%VmaxSIDMID,VmaxSIDM)
             call darkMatterProfile%floatRank0MetaPropertySet(self%RmaxSIDMID,RmaxSIDM)

             RmaxCDM =massDistribution_%radiusRotationCurveMaximum()
             RmaxSIDM=RmaxSIDM+RmaxCDM
             RmaxNFW0=Rmax_NFW(RmaxSIDM,tau)
             VmaxNFW0=Vmax_NFW(VmaxSIDM,tau)
             r_sNFW0  =r_s0(RmaxNFW0)
             rho_sNFW0=rho_s0   (r_sNFW0  ,VmaxNFW0)
             rho_s    =get_rho_s(rho_sNFW0,tau     )
             r_s      =get_r_s  (r_sNFW0  ,tau     )
             r_c      =get_r_c  (r_sNFW0  ,tau     )
             call darkMatterProfile%floatRank0MetaPropertySet(self%RhosSIDMID,rho_s)
             call darkMatterProfile%floatRank0MetaPropertySet(self%RsSIDMID  ,r_s  )
             call darkMatterProfile%floatRank0MetaPropertySet(self%RcSIDMID  ,r_c  )
             ! Copy the values from the new tree to the tip of the original branch.
             darkMatterProfile     => node            %darkMatterProfile()
             darkMatterProfileCopy => treeNew%nodeBase%darkMatterProfile()
             call darkMatterProfile%floatRank0MetaPropertySet(self%tauID     ,darkMatterProfileCopy%floatRank0MetaPropertyGet(self%tauID     ))
             call darkMatterProfile%floatRank0MetaPropertySet(self%VmaxSIDMID,darkMatterProfileCopy%floatRank0MetaPropertyGet(self%VmaxSIDMID))
             call darkMatterProfile%floatRank0MetaPropertySet(self%RmaxSIDMID,darkMatterProfileCopy%floatRank0MetaPropertyGet(self%RmaxSIDMID))
             call darkMatterProfile%floatRank0MetaPropertySet(self%RhosSIDMID,darkMatterProfileCopy%floatRank0MetaPropertyGet(self%RhosSIDMID))
             call darkMatterProfile%floatRank0MetaPropertySet(self%RsSIDMID  ,darkMatterProfileCopy%floatRank0MetaPropertyGet(self%RsSIDMID  ))
             call darkMatterProfile%floatRank0MetaPropertySet(self%RcSIDMID  ,darkMatterProfileCopy%floatRank0MetaPropertyGet(self%RcSIDMID  ))
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
    !!{
    Mark analytically-solvable properties.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentDarkMatterProfile
    implicit none
    class(nodeOperatorSIDMParametric    ), intent(inout) :: self
    type (treeNode                      ), intent(inout) :: node
    class(nodeComponentDarkMatterProfile), pointer       :: darkMatterProfile

    darkMatterProfile => node%darkMatterProfile()        
    call darkMatterProfile%floatRank0MetaPropertyAnalytic(self%RhosSIDMID)
    call darkMatterProfile%floatRank0MetaPropertyAnalytic(self%RsSIDMID  )
    call darkMatterProfile%floatRank0MetaPropertyAnalytic(self%RcSIDMID  )
    return
  end subroutine SIDMParametriDifferentialVmaxAnalytics

  subroutine SIDMParametriDifferentialVmaxSolveAnalytics(self,node,time)
    !!{
    Evolve ``\gls{dmou}'' mass at a constant rate, to achieve linear interpolation in time.
    !!}
    use :: Galacticus_Nodes  , only : nodeComponentBasic   , nodeComponentDarkMatterProfile
    use :: Mass_Distributions, only : massDistributionClass
    implicit none
    class           (nodeOperatorSIDMParametric  ), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    double precision                              , intent(in   ) :: time
    class           (nodeComponentBasic          ), pointer       :: basic            , basicParent
    class         (nodeComponentDarkMatterProfile), pointer       :: darkMatterProfile
    class           (massDistributionClass       ), pointer       :: massDistribution_
    double precision                                              :: tau              , r_sNFW0    , &
         &                                                           rho_sNFW0        , rho_s      , &
         &                                                           r_s              , r_c        , &
         &                                                           VmaxSIDM         , RmaxSIDM   , &
         &                                                           RmaxNFW0         , VmaxNFW0

    darkMatterProfile => node             %darkMatterProfile            (               )
    massDistribution_ => self             %darkMatterProfileDMO_    %get(node           )
    tau               =  darkMatterProfile%floatRank0MetaPropertyGet    (self%tauID     )   
    VmaxSIDM          =  darkMatterProfile%floatRank0MetaPropertyGet    (self%VmaxSIDMID)+massDistribution_%velocityRotationCurveMaximum()
    RmaxSIDM          =  darkMatterProfile%floatRank0MetaPropertyGet    (self%RmaxSIDMID)+massDistribution_%radiusRotationCurveMaximum  ()
    RmaxNFW0          =  Rmax_NFW (RmaxSIDM ,tau     )
    VmaxNFW0          =  Vmax_NFW (VmaxSIDM ,tau     )
    r_sNFW0           =  r_s0     (RmaxNFW0          )
    rho_sNFW0         =  rho_s0   (r_sNFW0  ,VmaxNFW0)
    rho_s             =  get_rho_s(rho_sNFW0, tau    )
    r_s               =  get_r_s  (r_sNFW0  , tau    )
    r_c               =  get_r_c  (r_sNFW0  , tau    )
    call darkMatterProfile%floatRank0MetaPropertySet(self%RhosSIDMID,rho_s)
    call darkMatterProfile%floatRank0MetaPropertySet(self%RsSIDMID  ,r_s  )
    call darkMatterProfile%floatRank0MetaPropertySet(self%RcSIDMID  ,r_c  )
    return
  end subroutine SIDMParametriDifferentialVmaxSolveAnalytics

  subroutine SIDMParametricNodePromote(self,node)
    !!{
    Ensure that \mono{node} is ready for promotion to its parent.
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
    !!{
    Set the ODE solver scales for the $\tau$ parameter.
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
    !!{
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
    class           (nodeComponentBasic            )               , pointer :: basic            , basicParent
    class           (nodeComponentDarkMatterProfile)               , pointer :: darkMatterProfile
    double precision                                                         :: timeFormation    , tc         , &
         &                                                                      tau              , time       , &
         &                                                                      VmaxSIDM         , RmaxSIDM   , &
         &                                                                      dvdt             , drdt       , &
         &                                                                      Mhalo            , dMdt       , &
         &                                                                      Gamma            , dtaudt
    !$GLC attributes unused :: interrupt, functionInterrupt, propertyType

    basic             => node %basic                    (                            )
    darkMatterProfile => node %darkMatterProfile        (                            )
    time              =  basic%time                     (                            )
    timeFormation     =  basic%floatRank0MetaPropertyGet(self%nodeFormationTimeSIDMID)
    if (time > timeFormation) then
       massDistribution_ => self%darkMatterProfileDMO_%get(node)
       VmaxSIDM          =  darkMatterProfile%floatRank0MetaPropertyGet(self%VmaxSIDMID)+massDistribution_%velocityRotationCurveMaximum()
       tc                =  get_tc(self%darkMatterParticle_,self%C,massDistribution_%velocityRotationCurveMaximum(),massDistribution_%radiusRotationCurveMaximum(),VmaxSIDM)
       tau               =  darkMatterProfile%floatRank0MetaPropertyGet(self%tauID)
       Mhalo             =  basic%mass()
       dMdt              =  basic%accretionRate()
       Gamma             =  dMdt/Mhalo
       dtaudt            =  1.0d0/tc-self%alpha*Gamma*tau
       call darkMatterProfile%floatRank0MetaPropertyRate(self%tauID,dtaudt)
       tau =darkMatterProfile%floatRank0MetaPropertyGet(self%tauID)
       dvdt=dvmaxt(tau,massDistribution_%velocityRotationCurveMaximum())*dtaudt
       drdt=drmaxt(tau,massDistribution_%radiusRotationCurveMaximum  ())*dtaudt
       call darkMatterProfile%floatRank0MetaPropertyRate(self%VmaxSIDMID,dvdt) 
       call darkMatterProfile%floatRank0MetaPropertyRate(self%RmaxSIDMID,drdt)
       VmaxSIDM=darkMatterProfile%floatRank0MetaPropertyGet(self%VmaxSIDMID)+massDistribution_%velocityRotationCurveMaximum()
       RmaxSIDM=darkMatterProfile%floatRank0MetaPropertyGet(self%RmaxSIDMID)+massDistribution_%radiusRotationCurveMaximum  ()
    else
       call darkMatterProfile%floatRank0MetaPropertyRate(self%tauID     ,0.0d0)
       call darkMatterProfile%floatRank0MetaPropertyRate(self%VmaxSIDMID,0.0d0)
       call darkMatterProfile%floatRank0MetaPropertyRate(self%RmaxSIDMID,0.0d0)       
    end if    
    return
  end subroutine SIDMParametriCalculateTauDifferentialEvolution
