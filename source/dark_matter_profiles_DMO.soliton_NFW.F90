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

  !+    Contributions to this file made by: Yu Zhao

  !!{
  An implementation of fuzzy dark matter halo profiles using the soliton and NFW mass distribution.
  !!}

  use :: Dark_Matter_Halo_Scales , only : darkMatterHaloScaleClass
  use :: Dark_Matter_Particles   , only : darkMatterParticleClass
  use :: Cosmology_Parameters    , only : cosmologyParametersClass
  use :: Cosmology_Functions     , only : cosmologyFunctionsClass
  use :: Virial_Density_Contrast , only : virialDensityContrastClass
  use :: Statistics_Distributions, only : distributionFunction1DNormal
  !![
  <darkMatterProfileDMO name="darkMatterProfileDMOSolitonNFW">
   <description>
    A dark matter profile DMO class which builds \refClass{massDistributionSolitonNFW} objects to implement the \gls{fdm}
    profile. The core-halo mass relation and core radius are computed following \cite{chan_diversity_2022}, while the core
    density normalization follows \cite{schive_understanding_2014}.
   </description>
  </darkMatterProfileDMO>
  !!]
  type, extends(darkMatterProfileDMOClass) :: darkMatterProfileDMOSolitonNFW
     !!{
     A dark matter halo profile class implementing \gls{fdm} dark matter halos.
     !!}
     private
     double precision                                        :: massParticle
     class           (darkMatterHaloScaleClass    ), pointer :: darkMatterHaloScale_               => null()
     class           (cosmologyParametersClass    ), pointer :: cosmologyParameters_               => null()
     class           (cosmologyFunctionsClass     ), pointer :: cosmologyFunctions_                => null()
     class           (darkMatterParticleClass     ), pointer :: darkMatterParticle_                => null()
     class           (virialDensityContrastClass  ), pointer :: virialDensityContrast_             => null()
     type            (distributionFunction1DNormal)          :: massCoreScatter
     double precision                                        :: toleranceRelativeVelocityDispersion         , toleranceRelativeVelocityDispersionMaximum
     double precision                                        :: radiusVirialPrevious                        , radiusScalePrevious                       , &
          &                                                     radiusCorePrevious                          , radiusSolitonPrevious                     , &
          &                                                     densityCorePrevious                         , densityScalePrevious                      , &
          &                                                     massCorePrevious
     integer          (kind_int8                  )          :: lastUniqueID
     integer                                                 :: randomOffsetID                              , densityCoreID                             , &
          &                                                     radiusCoreID                                , densityCoreAccretionID
   contains
     !![
     <methods>
       <method method="computeProperties" description="Compute properties of the mass distribution."/>
       <method method="calculationReset"  description="Reset memoized calculations."                />
     </methods>
     !!]
     final     ::                      solitonNFWDestructor
     procedure :: get               => solitonNFWGet
     procedure :: autoHook          => solitonNFWAutoHook
     procedure :: calculationReset  => solitonNFWCalculationReset
     procedure :: computeProperties => solitonNFWComputeProperties
  end type darkMatterProfileDMOSolitonNFW

  interface darkMatterProfileDMOSolitonNFW
     !!{
     Constructors for the {\normalfont \ttfamily solitonNFW} dark matter halo profile class.
     !!}
     module procedure solitonNFWConstructorParameters
     module procedure solitonNFWConstructorInternal
  end interface darkMatterProfileDMOSolitonNFW

  ! Sub-module scope variable used in root finding.
  double precision :: radiusCore_, radiusScale_, densityScale_, densityCore_
  !$omp threadprivate(radiusCore_, radiusScale_, densityScale_, densityCore_)

contains

  function solitonNFWConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{darkMatterProfileDMOSolitonNFW} dark matter halo profile class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (darkMatterProfileDMOSolitonNFW)                :: self
    type            (inputParameters               ), intent(inout) :: parameters
    class           (darkMatterHaloScaleClass      ), pointer       :: darkMatterHaloScale_
    class           (darkMatterParticleClass       ), pointer       :: darkMatterParticle_
    class           (cosmologyFunctionsClass       ), pointer       :: cosmologyFunctions_
    class           (cosmologyParametersClass      ), pointer       :: cosmologyParameters_
    class           (virialDensityContrastClass    ), pointer       :: virialDensityContrast_
    double precision                                                :: toleranceRelativeVelocityDispersion, toleranceRelativeVelocityDispersionMaximum

    !![
    <inputParameter>
      <name>toleranceRelativeVelocityDispersion</name>
      <defaultValue>1.0d-6</defaultValue>
      <source>parameters</source>
      <description>The relative tolerance to use in numerical solutions for the velocity dispersion.</description>
    </inputParameter>
    <inputParameter>
      <name>toleranceRelativeVelocityDispersionMaximum</name>
      <defaultValue>1.0d-3</defaultValue>
      <source>parameters</source>
      <description>The maximum relative tolerance to use in numerical solutions for the velocity dispersion.</description>
    </inputParameter>
    <objectBuilder class="darkMatterHaloScale"   name="darkMatterHaloScale_"   source="parameters"/>
    <objectBuilder class="darkMatterParticle"    name="darkMatterParticle_"    source="parameters"/>
    <objectBuilder class="cosmologyFunctions"    name="cosmologyFunctions_"    source="parameters"/>
    <objectBuilder class="cosmologyParameters"   name="cosmologyParameters_"   source="parameters"/>
    <objectBuilder class="virialDensityContrast" name="virialDensityContrast_" source="parameters"/>
    !!]
    self = darkMatterProfileDMOSolitonNFW(darkMatterHaloScale_,darkMatterParticle_,cosmologyFunctions_,cosmologyParameters_,virialDensityContrast_,toleranceRelativeVelocityDispersion,toleranceRelativeVelocityDispersionMaximum)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"  />
    <objectDestructor name="darkMatterParticle_"   />
    <objectDestructor name="cosmologyFunctions_"   />
    <objectDestructor name="cosmologyParameters_"  />
    <objectDestructor name="virialDensityContrast_"/>
    !!]
    return
  end function solitonNFWConstructorParameters

  function solitonNFWConstructorInternal(darkMatterHaloScale_,darkMatterParticle_,cosmologyFunctions_,cosmologyParameters_,virialDensityContrast_,toleranceRelativeVelocityDispersion,toleranceRelativeVelocityDispersionMaximum) result(self)
    !!{
    Generic constructor for the \refClass{darkMatterProfileDMOSolitonNFW} dark matter halo profile class.
    !!}
    use :: Error                         , only : Component_List                   , Error_Report
    use :: Dark_Matter_Particles         , only : darkMatterParticleFuzzyDarkMatter
    use :: Numerical_Constants_Prefixes  , only : kilo
    use :: Galacticus_Nodes              , only : defaultDarkMatterProfileComponent
    implicit none
    type            (darkMatterProfileDMOSolitonNFW)                     :: self
    class           (darkMatterHaloScaleClass      ), intent(in), target :: darkMatterHaloScale_
    class           (darkMatterParticleClass       ), intent(in), target :: darkMatterParticle_
    class           (cosmologyFunctionsClass       ), intent(in), target :: cosmologyFunctions_
    class           (cosmologyParametersClass      ), intent(in), target :: cosmologyParameters_
    class           (virialDensityContrastClass    ), intent(in), target :: virialDensityContrast_
    double precision                                , intent(in)         :: toleranceRelativeVelocityDispersion, toleranceRelativeVelocityDispersionMaximum
    !![
    <constructorAssign variables="*darkMatterHaloScale_, *darkMatterParticle_, *cosmologyFunctions_, *cosmologyParameters_, *virialDensityContrast_, toleranceRelativeVelocityDispersion, toleranceRelativeVelocityDispersionMaximum"/>
    <addMetaProperty component="darkMatterProfile" name="randomOffset" id="self%randomOffsetID" isEvolvable="no" isCreator="yes"/>
    <addMetaProperty component="basic"             name="densityCore"  id="self%densityCoreID"  isEvolvable="no" isCreator="yes"/>
    <addMetaProperty component="basic"             name="radiusCore"   id="self%radiusCoreID"   isEvolvable="no" isCreator="yes"/>
    <addMetaProperty component="basic"             name="densityCoreAccretion" id="self%densityCoreAccretionID" isEvolvable="yes" isCreator="no"/>
    !!]

    self%lastUniqueID=-huge(1_kind_int8)
    self%massCoreScatter = distributionFunction1DNormal(mean=0.0d0,variance = log10(1.5d0)**2) ! 50% log-normal scatter from Eq.(15) of Chan et al. (2022; MNRAS; 551; 943; https://ui.adsabs.harvard.edu/abs/2022MNRAS.511..943C).

    select type (darkMatterParticle__ => self%darkMatterParticle_)
    class is (darkMatterParticleFuzzyDarkMatter)
       self%massParticle=+darkMatterParticle__%mass()*kilo
    class default
       call Error_Report('expected a `darkMatterParticleFuzzyDarkMatter` dark matter particle object'//{introspection:location})
    end select
    if (.not.defaultDarkMatterProfileComponent%scaleIsGettable())                                                             &
        & call Error_Report                                                                                                   &
        &      (                                                                                                              &
        &       'solitonNFW dark matter profile requires a dark matter profile component with a gettable "scale" property.'// &
        &       Component_List(                                                                                               &
        &                      'darkMatterProfile'                                                                         ,  &
        &                      defaultDarkMatterProfileComponent%scaleAttributeMatch(requireGettable=.true.)                  &
        &                     )                                                                                            // &
        &      {introspection:location}                                                                                       &
        &      )
    return
  end function solitonNFWConstructorInternal

  subroutine solitonNFWAutoHook(self)
    !!{
    Attach to the calculation reset event.
    !!}
    use :: Events_Hooks, only : calculationResetEvent, openMPThreadBindingAllLevels
    implicit none
    class(darkMatterProfileDMOSolitonNFW), intent(inout) :: self

    call calculationResetEvent%attach(self,solitonNFWCalculationReset,openMPThreadBindingAllLevels,label='darkMatterProfileDMOSolitonNFW')
    return
  end subroutine solitonNFWAutoHook

  subroutine solitonNFWDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily solitonNFW} dark matter halo profile class.
    !!}
    use :: Events_Hooks, only : calculationResetEvent
    implicit none
    type(darkMatterProfileDMOSolitonNFW), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_"  />
    <objectDestructor name="self%darkMatterParticle_"   />
    <objectDestructor name="self%cosmologyFunctions_"   />
    <objectDestructor name="self%cosmologyParameters_"  />
    <objectDestructor name="self%virialDensityContrast_"/>
    !!]
    if (calculationResetEvent%isAttached(self,solitonNFWCalculationReset)) call calculationResetEvent%detach(self,solitonNFWCalculationReset)
    return
  end subroutine solitonNFWDestructor

  subroutine solitonNFWCalculationReset(self,node,uniqueID)
    !!{
    Reset the dark matter profile calculation.
    !!}
    use :: Kind_Numbers, only : kind_int8
    implicit none
    class  (darkMatterProfileDMOSolitonNFW), intent(inout) :: self
    type   (treeNode                      ), intent(inout) :: node
    integer(kind_int8                     ), intent(in   ) :: uniqueID
    !$GLC attributes unused :: node

    self%lastUniqueID            =uniqueID
    self%radiusVirialPrevious    =-huge(0.0d0)
    self%radiusScalePrevious     =-huge(0.0d0)
    self%radiusCorePrevious      =-huge(0.0d0)
    self%radiusSolitonPrevious   =-huge(0.0d0)
    self%densityCorePrevious     =-huge(0.0d0)
    self%densityScalePrevious    =-huge(0.0d0)
    self%massCorePrevious        =-huge(0.0d0)
    return
  end subroutine solitonNFWCalculationReset

  function solitonNFWGet(self,node,weightBy,weightIndex) result(massDistribution_)
    !!{
    Return the soliton plus NFW fuzzy dark matter mass distribution for the given {\normalfont \ttfamily node}.
    !!}
    use :: Galactic_Structure_Options, only : componentTypeDarkHalo     , massTypeDark                    , weightByMass
    use :: Mass_Distributions        , only : massDistributionSolitonNFW, kinematicsDistributionSolitonNFW, kinematicsDistributionClass
    implicit none
    class           (darkMatterProfileDMOSolitonNFW), intent(inout)           :: self
    class           (massDistributionClass         ), pointer                 :: massDistribution_
    class           (kinematicsDistributionClass   ), pointer                 :: kinematicsDistribution_
    type            (treeNode                      ), intent(inout)           :: node
    type            (enumerationWeightByType       ), intent(in   ), optional :: weightBy
    integer                                         , intent(in   ), optional :: weightIndex
    type            (enumerationWeightByType       )                          :: weightBy_
    double precision                                                          :: radiusCore             , radiusScale , &
         &                                                                       radiusSoliton          , radiusVirial, &
         &                                                                       densityScale           , densityCore , &
         &                                                                       massCore
    !![
    <optionalArgument name="weightBy" defaultsTo="weightByMass" />
    !!]

    ! Return a null distribution if weighting is not by mass.
    massDistribution_ => null()
    if (weightBy_ /= weightByMass) return
    ! Compute properties of the distribution.
    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node,node%uniqueID())
    if (self%radiusCorePrevious < 0.0d0) then
       call self%computeProperties(node,radiusVirial,radiusScale,radiusCore,radiusSoliton,densityCore,densityScale,massCore)
       self%radiusVirialPrevious =radiusVirial
       self%radiusScalePrevious  =radiusScale
       self%radiusCorePrevious   =radiusCore
       self%radiusSolitonPrevious=radiusSoliton
       self%densityCorePrevious  =densityCore
       self%densityScalePrevious =densityScale
       self%massCorePrevious     =massCore
    end if
    radiusVirial =self%radiusVirialPrevious
    radiusScale  =self%radiusScalePrevious
    radiusCore   =self%radiusCorePrevious
    radiusSoliton=self%radiusSolitonPrevious
    densityCore  =self%densityCorePrevious
    densityScale =self%densityScalePrevious
    massCore     =self%massCorePrevious
    
    ! Construct the distribution.
    allocate(massDistributionSolitonNFW       :: massDistribution_      )
    allocate(kinematicsDistributionSolitonNFW :: kinematicsDistribution_)
    select type(massDistribution_)
    type is (massDistributionSolitonNFW)
       !![
       <referenceConstruct object="massDistribution_">
	 <constructor>
           massDistributionSolitonNFW(                                                &amp;
            &amp;                     radiusScale            = radiusScale          , &amp;
            &amp;                     radiusCore             = radiusCore           , &amp;
            &amp;                     radiusSoliton          = radiusSoliton        , &amp;
            &amp;                     densitySolitonCentral  = densityCore          , &amp;
            &amp;                     densityNormalizationNFW= densityScale         , &amp;
            &amp;                     radiusVirial           = radiusVirial         , &amp;
            &amp;                     componentType          = componentTypeDarkHalo, &amp;
            &amp;                     massType               = massTypeDark           &amp;
            &amp;                    )
	 </constructor>
       </referenceConstruct>
       !!]
    end select
    select type (kinematicsDistribution_)
    type is (kinematicsDistributionSolitonNFW)
       !![
       <referenceConstruct object="kinematicsDistribution_">
	 <constructor>
           kinematicsDistributionSolitonNFW(                                                                                            &amp;
            &amp;                           toleranceRelativeVelocityDispersion       =self%toleranceRelativeVelocityDispersion       , &amp;
            &amp;                           toleranceRelativeVelocityDispersionMaximum=self%toleranceRelativeVelocityDispersionMaximum  &amp;
	    &amp;                          )
	 </constructor>
       </referenceConstruct>
       !!]
    end select
    call massDistribution_%setKinematicsDistribution(kinematicsDistribution_)
    !![
    <objectDestructor name="kinematicsDistribution_"/>
    !!]
    return
  end function solitonNFWGet

  subroutine solitonNFWComputeProperties(self,node,radiusVirial,radiusScale,radiusCore,radiusSoliton,densityCore,densityScale,massCore)
    use :: Galacticus_Nodes                , only : treeNode           , nodeComponentBasic       , nodeComponentDarkMatterProfile
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Units       , only : electronVolt
    use :: Numerical_Constants_Astronomical, only : megaParsec
    use :: Numerical_Constants_Physical    , only : speedLight         , plancksConstant
    use :: Numerical_Constants_Prefixes    , only : kilo
    use :: Cosmology_Parameters            , only : hubbleUnitsStandard, hubbleUnitsLittleH
    use :: Root_Finder                     , only : rootFinder         , rangeExpandMultiplicative, rangeExpandSignExpectPositive , rangeExpandSignExpectNegative
    implicit none
    class           (darkMatterProfileDMOSolitonNFW), intent(inout) :: self
    type            (treeNode                      ), intent(inout) :: node
    double precision                                , intent(  out) :: radiusVirial                        , radiusScale                , &
         &                                                             radiusCore                          , radiusSoliton              , &
         &                                                             densityCore                         , densityScale               , &
         &                                                             massCore
    class           (nodeComponentBasic            ), pointer       :: basic
    class           (nodeComponentDarkMatterProfile), pointer       :: darkMatterProfile
    type            (rootFinder                    ), save          :: finder
    logical                                         , save          :: finderInitialized =.false.
    !$omp threadprivate(finder, finderInitialized)
    double precision                                , parameter     :: toleranceAbsolute = 0.0d0           , toleranceRelative  =1.0d-3
    double precision                                , parameter     :: plancksConstantBar=+plancksConstant                                & ! ℏ in units of eV s.
         &                                                                                /2.0d0                                          &
         &                                                                                /Pi                                             &
         &                                                                                /electronVolt
    double precision                                                :: massHalo                            , expansionFactor            , &
         &                                                             redshift                            , concentration              , &
         &                                                             hubbleConstant                      , hubbleConstantLittle       , &
         &                                                             OmegaMatter                         , densityMatter              , &
         &                                                             zeta_0                              , zeta_z                     , &
         &                                                             randomOffset                        , massCoreNormal
    double precision                                , parameter     :: alpha             =0.515            , beta               =8.0d6  , &
         &                                                             gamma             =10.0d0**(-5.73d0)                                 ! Best-fitting parameters from Chan et al. (2022; MNRAS; 551; 943; https://ui.adsabs.harvard.edu/abs/2022MNRAS.511..943C).
    integer                                                         :: status                              , sampleCount                , &
         &                                                             maxSamples = 50

    ! Get required components.
    basic             => node%basic            ()
    darkMatterProfile => node%darkMatterProfile()
    ! Extract basic properties of the node.
    expansionFactor=+self             %cosmologyFunctions_% expansionFactor            (basic%time           ())
    redshift       =+self             %cosmologyFunctions_ %redshiftFromExpansionFactor(      expansionFactor  )
    massHalo       =+basic                                 %mass                       (                       )
    radiusScale    =+darkMatterProfile                     %scale                      (                       )
    radiusVirial   =+self             %darkMatterHaloScale_%radiusVirial               (node                   )
    concentration  =+                                       radiusVirial                                         &
         &          /                                       radiusScale
    densityScale   =+1.0d0                                    &
         &          /4.0d0                                    &
         &          /Pi                                       &
         &          *massHalo                                 &
         &          /radiusScale**3                           &
         &          /(                                        &
         &            +              log(1.0d0+concentration) &
         &            -concentration/   (1.0d0+concentration) &
         &          )
    ! Extract cosmological parameters for later use.
    hubbleConstant      =+self%cosmologyParameters_%HubbleConstant (hubbleUnitsStandard)
    hubbleConstantLittle=+self%cosmologyParameters_%HubbleConstant (hubbleUnitsLittleH )
    OmegaMatter         =+self%cosmologyParameters_%OmegaMatter    (                   )
    densityMatter       =+self%cosmologyParameters_%densityCritical(                   ) &
         &                      *                   OmegaMatter
    ! Compute the core mass.
    zeta_0             =+self% virialDensityContrast_%densityContrast(massHalo,expansionFactor=1.0d0          )
    zeta_z             =+self% virialDensityContrast_%densityContrast(massHalo,expansionFactor=expansionFactor)
    massCoreNormal     =+basic%floatRank0MetaPropertyGet(self%densityCoreAccretionID)

    radiusScale_       =radiusScale
    densityScale_      =densityScale
    ! Solve for the soliton radius.
    if (.not.finderInitialized) then
       finder=rootFinder(                                        &
            &            rootFunction     =radiusTransitionRoot, &
            &            toleranceAbsolute=toleranceAbsolute   , &
            &            toleranceRelative=toleranceRelative     &
            &           )
       finderInitialized=.true.
    end if

    do sampleCount=1,maxSamples
       if (sampleCount == 1) then
           randomOffset     = darkMatterProfile%floatRank0MetaPropertyGet(self%randomOffsetID)
           if (randomOffset == 0.0d0) then
               randomOffset = self%massCoreScatter%sample(randomNumberGenerator_=node%hostTree%randomNumberGenerator_)
           end if
       else
           randomOffset     = self%massCoreScatter%sample(randomNumberGenerator_=node%hostTree%randomNumberGenerator_)
       end if
       massCore             = massCoreNormal*10.0d0**randomOffset
       ! Compute the core radius.
       radiusCore         =+5.5d6                           & ! Equation (14) of Chan et al. (2022; MNRAS; 551; 943; https://ui.adsabs.harvard.edu/abs/2022MNRAS.511..943C).
            &              /(self%massParticle/1.0d-23)**2  &
            &              /expansionFactor                 &
            &              /massCore
       ! Compute the core density normalization.
       densityCore       =+massCore                         & ! Equation (3) of Schive et al. (2014; PRL; 113; 1302; https://ui.adsabs.harvard.edu/abs/2014PhRvL.113z1302S).
            &             /0.413d0                          &
            &             /(radiusCore               **3)   &
            &             /Pi
       radiusCore_        =radiusCore
       densityCore_       =densityCore
       call finder%rangeExpand(                                                              &
         &                     rangeExpandUpward            =2.0d0                        , &
         &                     rangeExpandDownward          =0.5d0                        , &
         &                     rangeDownwardLimit           =1.0d0*radiusCore             , &
         &                     rangeUpwardLimit             =1.0d1*radiusCore             , &
         &                     rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive, &
         &                     rangeExpandUpwardSignExpect  =rangeExpandSignExpectNegative, &
         &                     rangeExpandType              =rangeExpandMultiplicative      &
         &                    )
       radiusSoliton=finder%find(rootGuess=3.0d0*radiusCore,status=status)
       if (status == errorStatusSuccess) then
           call darkMatterProfile%floatRank0MetaPropertySet(self%randomOffsetID,randomOffset)
           exit
       end if
    end do
    call basic%floatRank0MetaPropertySet(self%densityCoreID,densityCore)
    call basic%floatRank0MetaPropertySet(self%radiusCoreID ,radiusCore )
    return
  end subroutine solitonNFWComputeProperties

  double precision function radiusTransitionRoot(radius) result(f)
    !!{
    Root function used in seeking the transition radius in fuzzy dark matter profiles.
    !!}
    implicit none
    double precision, intent(in   ) :: radius

    f      =+                radius/radiusScale_         &
         &  *(1.0d0+         radius/radiusScale_)**2     &
         &  /(1.0d0+0.091d0*(radius/radiusCore_ )**2)**8 &
         &  -densityScale_                               &
         &  /densityCore_
    return
  end function radiusTransitionRoot
