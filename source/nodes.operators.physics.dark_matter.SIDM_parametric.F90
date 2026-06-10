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
    type            (treeNode                                 )               , pointer :: nodeParent                              , nodeBase              , &
         &                                                                                 nodeChild                               , nodeNew
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
    ! Fractional halo mass decline per step used when extrapolating the tree below resolution.
    double precision                                           , parameter              :: massHaloDeclineFactor            =0.99d0
    double precision                                                                    :: timeFormation
    double precision                                                                    :: VmaxSIDMPrevious                        , tc                    , &
         &                                                                                 tau                                     , dtr                   , &
         &                                                                                 VmaxSIDM                                , RmaxSIDM              , &
         &                                                                                 RmaxCDM                                 , RmaxNFW0              , &
         &                                                                                 VmaxNFW0                                , r_sNFW0               , &
         &                                                                                 rho_sNFW0                               , rho_s                 , &
         &                                                                                 r_s                                     , r_c                   , &
         &                                                                                 VmaxCDM                                 , Gamma                 , &
         &                                                                                 dt                                      , tau_old               , &
         &                                                                                 Mhalo                                   , dMdt                  , &
         &                                                                                 timeEarliest                            , massResolution


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
          ! If this node is at the tip of the branch and the time associated to it is larger than the formation time of the
          ! branch then we need to build the rest of the tree (extrapolate back in time) to start the parametric SIDM
          ! calculation from the formation time. We build the tree with smooth accretion only. 
          ! The mass resolution is set to be the current node mass times the formation mass fraction.
          massResolution=formationMassFraction*basic%mass()
          timeEarliest  =basic%floatRank0MetaPropertyGet(self%nodeFormationTimeSIDMID) - 1
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
             tc     =get_tc(self,nodeNew,massDistribution_%velocityRotationCurveMaximum(),massDistribution_%radiusRotationCurveMaximum(),VmaxSIDMPrevious)
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
    type            (treeNode                      )               , pointer :: nodeFinal        , nodeWork        , &
         &                                                                      nodeStart        , nodeParent
    class           (massDistributionClass         )               , pointer :: massDistribution_
    class           (nodeComponentBasic            )               , pointer :: basic            , basicParent
    class           (nodeComponentDarkMatterProfile)               , pointer :: darkMatterProfile
    double precision                                                         :: timeFormation    , dtr             , &
         &                                                                      tau              , time            , &
         &                                                                      timePrevious     , tc              , &
         &                                                                      RmaxNFW0         , VmaxNFW0        , &
         &                                                                      r_sNFW0          , rho_sNFW0       , &
         &                                                                      rho_s            , r_s             , &
         &                                                                      VmaxSIDM         , RmaxSIDM        , &
         &                                                                      VmaxSIDMPrevious , RmaxSIDMPrevious, &
         &                                                                      VmaxCDM          , RmaxCDM         , &
         &                                                                      VmaxCDMPrevious  , RmaxCDMPrevious , &
         &                                                                      dvdt             , drdt            , &
         &                                                                      Mhalo            , dMdt            , &
         &                                                                      Gamma            , dtaudt          , &
         &                                                                      r_c
    !$GLC attributes unused :: interrupt, functionInterrupt, propertyType

    basic             => node %basic                    (                            )
    darkMatterProfile => node %darkMatterProfile        (                            )
    time              =  basic%time                     (                            )
    timeFormation     =  basic%floatRank0MetaPropertyGet(self%nodeFormationTimeSIDMID)
    if (time > timeFormation) then
       massDistribution_ => self%darkMatterProfileDMO_%get(node)
       VmaxSIDM          =  darkMatterProfile%floatRank0MetaPropertyGet(self%VmaxSIDMID)+massDistribution_%velocityRotationCurveMaximum()
       tc                =  get_tc(self,node,massDistribution_%velocityRotationCurveMaximum(),massDistribution_%radiusRotationCurveMaximum(),VmaxSIDM)
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

  double precision function get_tc(self,node,Vmax,Rvmax, VmaxSIDM)
    !!{
    Evaluate the gravothermal evolution timescale $t_\mathrm{c}$ following eqn.~(2.2) of \cite{yang_parametric_2024}.
    !!}
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    use :: Numerical_Constants_Prefixes    , only : kilo
    use :: Error                           , only : Error_Report
    use :: Dark_Matter_Particles           , only : darkMatterParticleSelfInteractingDarkMatter
    implicit none
    class           (nodeOperatorSIDMParametric), intent(inout) :: self
    type            (treeNode                  ), intent(in   ) :: node
    double precision                            , intent(in   ) :: Vmax                  , Rvmax      , &
         &                                                         VmaxSIDM
    double precision                            , parameter     :: timescaleNormalization=150.0000d+00 ! Numerical coefficient in the gravothermal timescale (eqn. 2.2).
    double precision                            , parameter     :: radiusMaximumToScale  =  2.1626d+00 ! r_max/r_s for an NFW profile.
    double precision                            , parameter     :: velocityMaximumToScale=  1.6480d+00 ! V_max normalization coefficient for an NFW profile.
    double precision                            , parameter     :: crossSectionConversion=  2.0900d-10 ! Converts the cross section per unit mass from cm^2 g^-1 to kpc^2 M_sun^-1.
    double precision                                            :: sigmaeff              , reff       , &
         &                                                         gravitationalConstant , rhoeff
    !$GLC attributes unused :: node

    select type (darkMatterParticle_ => self%darkMatterParticle_)
    class is (darkMatterParticleSelfInteractingDarkMatter)
       sigmaeff=darkMatterParticle_%crossSectionEffective(VmaxSIDM)
    class default
       call Error_Report('unexpected class'//{introspection:location})
    end select
    ! Work in kpc length units (the internal length unit is Mpc, so apply a factor of kilo to lengths and to G).
    gravitationalConstant=gravitationalConstant_internal*kilo
    reff                 =Rvmax*kilo/radiusMaximumToScale
    rhoeff               =(Vmax/(velocityMaximumToScale*reff))**2/gravitationalConstant
    get_tc               =(timescaleNormalization/self%C)                                &
         &                *(1.0d0/(sigmaeff*crossSectionConversion*rhoeff*reff))         &
         &                *(1.0d0/sqrt(4.0d0*Pi*gravitationalConstant*rhoeff))
    return
  end function get_tc

  double precision function dvmaxt(tau, Vmaxt)
    !!{
    Return the derivative $\mathrm{d}V_\mathrm{max}/\mathrm{d}\tau$ of the SIDM maximum circular velocity with respect to the
    dimensionless gravothermal time $\tau$, scaled by the CDM maximum circular velocity \mono{Vmaxt}, using the
    polynomial fit of \cite{yang_parametric_2024}. The result is zero for $\tau>1$ (the core-collapsed regime).
    !!}
    implicit none
    double precision, intent(in   ) :: tau, Vmaxt

    if (tau > 1.0d0) then
      dvmaxt = 0.0d0
    else if (tau <= 1.0d0) then
      dvmaxt = 0.17774902d0 - 13.195824689999998d0 * tau ** 2 + 66.62092676d0 * tau ** 3 - 94.33706049999999d0 * tau ** 4 + 63.53766111d0 * tau ** 6 - 21.925108889999997d0 * tau ** 8
    end if
    dvmaxt=dvmaxt*Vmaxt
    return
  end function dvmaxt

  double precision function drmaxt(tau, Rmaxt)
    !!{
    Return the derivative $\mathrm{d}R_\mathrm{max}/\mathrm{d}\tau$ of the SIDM maximum-circular-velocity radius with respect to
    the dimensionless gravothermal time $\tau$, scaled by the CDM radius \mono{Rmaxt}, using the polynomial fit of
    \cite{yang_parametric_2024}. The result is zero for $\tau>1$ (the core-collapsed regime).
    !!}
    implicit none
    double precision, intent(in   ) :: tau, Rmaxt

    if (tau > 1.0d0) then
      drmaxt=0.0d0
    else if (tau <= 1.0d0) then
      drmaxt=0.00762288d0 - 1.43996392d0 * tau + 1.01282643d0 * tau ** 2 - 0.55015288d0 * tau ** 3
    end if
    drmaxt=drmaxt*Rmaxt
    return
  end function drmaxt

  double precision function Rmax_NFW(RmaxSIDM, tau)
    !!{
    Map an SIDM maximum-circular-velocity radius \mono{RmaxSIDM} back to the equivalent CDM (NFW) value by
    inverting the $R_\mathrm{max}(\tau)$ evolution fit of \cite{yang_parametric_2024}. \mono{tau} is clamped to
    the range $[0,1]$.
    !!}
    implicit none
    double precision, intent(in) :: RmaxSIDM, tau
    double precision             :: tau_

    tau_=min(max(tau,0.0d0),1.0d0)
    Rmax_NFW=RmaxSIDM/(1.0d0+0.007623d0*tau_-0.7200d0*tau_**2+0.3376d0*tau_**3-0.1375d0*tau_**4)
    return
  end function Rmax_NFW

  double precision function Vmax_NFW(VmaxSIDM, tau)
    !!{
    Map an SIDM maximum circular velocity \mono{VmaxSIDM} back to the equivalent CDM (NFW) value by inverting the
    $V_\mathrm{max}(\tau)$ evolution fit of \cite{yang_parametric_2024}. \mono{tau} is clamped to the range
    $[0,1]$.
    !!}
    implicit none
    double precision, intent(in) :: VmaxSIDM, tau
    double precision             :: tau_

    tau_=min(max(tau,0.0d0),1.0d0)
    Vmax_NFW = VmaxSIDM / (1 + 0.1777d0 * tau_ - 4.399d0 * tau_ ** 3 + 16.66d0 * tau_ ** 4 - 18.87d0 * tau_ ** 5 + 9.077d0 * tau_ ** 7 - 2.436d0 * tau_ ** 9)
    return
  end function Vmax_NFW

  double precision function r_s0(Rmax)
    !!{
    Return the NFW scale radius corresponding to a maximum-circular-velocity radius \mono{Rmax}, using the
    standard NFW relation $r_\mathrm{s}=R_\mathrm{max}/2.163$.
    !!}
    implicit none
    double precision, intent(in) :: Rmax

    r_s0=Rmax / 2.163d0
    return
  end function r_s0

  double precision function rho_s0(Rs, Vmax)
    !!{
    Return the NFW characteristic density corresponding to a scale radius \mono{Rs} and maximum circular velocity
    \mono{Vmax}.
    !!}
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_astronomical, only : gravitationalConstant_internal
    implicit none
    double precision, intent(in) :: Rs, Vmax

    rho_s0 = Vmax ** 2 / (0.465d0 ** 2 * 4.0d0 * Pi * gravitationalConstant_internal * Rs ** 2)
    return
  end function rho_s0

  double precision function get_rho_s(rho_s0, tau)
    !!{
    Return the SIDM parametric-profile characteristic density as a function of the dimensionless gravothermal time $\tau$,
    normalized to the initial NFW characteristic density \mono{rho\_s0}, using the fit of
    \cite{yang_parametric_2024}. \mono{tau} is clamped to the range $[0,1]$.
    !!}
    implicit none
    double precision, intent(in) :: rho_s0, tau
    double precision             :: tau_

    tau_=min(max(tau,0.0d0),1.0d0)
    get_rho_s = rho_s0 * (2.033d0 + 0.7381d0 * tau_ + 7.264d0 * tau_ ** 5 - 12.73d0 * tau_ ** 7 + 9.915d0 * tau_ ** 9 + (1.0d0 - 2.033d0) * log(tau_ + 0.001d0) / log(0.001d0))
    return
  end function get_rho_s

  double precision function get_r_s(r_s0, tau)
    !!{
    Return the SIDM parametric-profile scale radius as a function of the dimensionless gravothermal time $\tau$, normalized to the
    initial NFW scale radius \mono{r\_s0}, using the fit of \cite{yang_parametric_2024}. \mono{tau} is clamped to the range $[0,1]$.
    !!}
    implicit none
    double precision, intent(in) :: r_s0, tau
    double precision             :: tau_

    tau_=min(max(tau,0.0d0),1.0d0)
    get_r_s = r_s0 * (0.7178d0 - 0.1026d0 * tau_ + 0.2474d0 * tau_ ** 2 - 0.4079d0 * tau_ ** 3 + (1.0d0 - 0.7178d0) * log(tau_ + 0.001d0) / log(0.001d0))
    return
  end function get_r_s

  double precision function get_r_c(r_s0, tau)
    !!{
    Return the SIDM parametric-profile core radius as a function of the dimensionless gravothermal time $\tau$, normalized to the
    initial NFW scale radius \mono{r\_s0}, using the fit of \cite{yang_parametric_2024}. \mono{tau} is clamped to the range $[0,1]$.
    !!}
    implicit none
    double precision, intent(in) :: r_s0, tau
    double precision             :: tau_

    tau_=min(max(tau,0.0d0),1.0d0)
    get_r_c = r_s0 * (2.555d0 * sqrt(tau_) - 3.632d0 * tau_ + 2.131d0 * tau_ ** 2 - 1.415d0 * tau_ ** 3 + 0.4683d0 * tau_ ** 4)
    return
  end function get_r_c


 
