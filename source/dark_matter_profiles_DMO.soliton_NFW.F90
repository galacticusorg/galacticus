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

  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass
  use :: Dark_Matter_Particles  , only : darkMatterParticleClass
  use :: Cosmology_Parameters   , only : cosmologyParametersClass
  use :: Cosmology_Functions    , only : cosmologyFunctionsClass
  use :: Virial_Density_Contrast, only : virialDensityContrastClass
  !![
  <darkMatterProfileDMO name="darkMatterProfileDMOSolitonNFW">
   <description>
    A dark matter profile DMO class which builds \refClass{massDistributionSolitonNFW} objects to implement the \gls{fdm}
    profile. The calculation of soliton properties follows \cite{schive_understanding_2014}.
   </description>
  </darkMatterProfileDMO>
  !!]
  type, extends(darkMatterProfileDMOClass) :: darkMatterProfileDMOSolitonNFW
     !!{
     A dark matter halo profile class implementing \gls{fdm} dark matter halos.
     !!}
     private
     double precision                                      :: massParticle
     class           (darkMatterHaloScaleClass  ), pointer :: darkMatterHaloScale_   => null()
     class           (cosmologyParametersClass  ), pointer :: cosmologyParameters_   => null()
     class           (cosmologyFunctionsClass   ), pointer :: cosmologyFunctions_    => null()
     class           (darkMatterParticleClass   ), pointer :: darkMatterParticle_    => null()
     class           (virialDensityContrastClass), pointer :: virialDensityContrast_ => null()
   contains
     !![
     <methods>
       <method method="computeProperties" description="Compute properties of the mass distribution."/>
     </methods>
     !!]
     final     ::                      solitonNFWDestructor
     procedure :: get               => solitonNFWGet
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
    type (darkMatterProfileDMOSolitonNFW)                :: self
    type (inputParameters               ), intent(inout) :: parameters
    class(darkMatterHaloScaleClass      ), pointer       :: darkMatterHaloScale_
    class(darkMatterParticleClass       ), pointer       :: darkMatterParticle_
    class(cosmologyFunctionsClass       ), pointer       :: cosmologyFunctions_
    class(cosmologyParametersClass      ), pointer       :: cosmologyParameters_
    class(virialDensityContrastClass    ), pointer       :: virialDensityContrast_

    !![
    <objectBuilder class="darkMatterHaloScale"   name="darkMatterHaloScale_"   source="parameters"/>
    <objectBuilder class="darkMatterParticle"    name="darkMatterParticle_"    source="parameters"/>
    <objectBuilder class="cosmologyFunctions"    name="cosmologyFunctions_"    source="parameters"/>
    <objectBuilder class="cosmologyParameters"   name="cosmologyParameters_"   source="parameters"/>
    <objectBuilder class="virialDensityContrast" name="virialDensityContrast_" source="parameters"/>
    !!]
    self = darkMatterProfileDMOSolitonNFW(darkMatterHaloScale_,darkMatterParticle_,cosmologyFunctions_,cosmologyParameters_,virialDensityContrast_)
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

  function solitonNFWConstructorInternal(darkMatterHaloScale_,darkMatterParticle_,cosmologyFunctions_,cosmologyParameters_,virialDensityContrast_) result(self)
    !!{
    Generic constructor for the \refClass{darkMatterProfileDMOSolitonNFW} dark matter halo profile class.
    !!}
    use :: Error                         , only : Component_List                   , Error_Report
    use :: Dark_Matter_Particles         , only : darkMatterParticleFuzzyDarkMatter
    use :: Numerical_Constants_Prefixes  , only : kilo
    use :: Galacticus_Nodes              , only : defaultDarkMatterProfileComponent
    implicit none
    type (darkMatterProfileDMOSolitonNFW)                     :: self
    class(darkMatterHaloScaleClass      ), intent(in), target :: darkMatterHaloScale_
    class(darkMatterParticleClass       ), intent(in), target :: darkMatterParticle_
    class(cosmologyFunctionsClass       ), intent(in), target :: cosmologyFunctions_
    class(cosmologyParametersClass      ), intent(in), target :: cosmologyParameters_
    class(virialDensityContrastClass    ), intent(in), target :: virialDensityContrast_
    !![
    <constructorAssign variables="*darkMatterHaloScale_, *darkMatterParticle_, *cosmologyFunctions_, *cosmologyParameters_, *virialDensityContrast_"/>
    !!]

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

  subroutine solitonNFWDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily solitonNFW} dark matter halo profile class.
    !!}
    implicit none
    type(darkMatterProfileDMOSolitonNFW), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_"  />
    <objectDestructor name="self%darkMatterParticle_"   />
    <objectDestructor name="self%cosmologyFunctions_"   />
    <objectDestructor name="self%cosmologyParameters_"  />
    <objectDestructor name="self%virialDensityContrast_"/>
    !!]
    return
  end subroutine solitonNFWDestructor

  function solitonNFWGet(self,node,weightBy,weightIndex) result(massDistribution_)
    !!{
    Return the soliton plus NFW fuzzy dark matter mass distribution for the given {\normalfont \ttfamily node}.
    !!}
    use :: Galactic_Structure_Options, only : componentTypeDarkHalo     , massTypeDark                                , weightByMass
    use :: Mass_Distributions        , only : massDistributionSolitonNFW, kinematicsDistributionCollisionlessTabulated
    implicit none
    class           (darkMatterProfileDMOSolitonNFW              ), intent(inout)           :: self
    class           (massDistributionClass                       ), pointer                 :: massDistribution_
    type            (kinematicsDistributionCollisionlessTabulated), pointer                 :: kinematicsDistribution_
    type            (treeNode                                    ), intent(inout)           :: node
    type            (enumerationWeightByType                     ), intent(in   ), optional :: weightBy
    integer                                                       , intent(in   ), optional :: weightIndex
    type            (enumerationWeightByType                     )                          :: weightBy_
    double precision                                                                        :: radiusCore             , radiusScale , &
         &                                                                                     radiusSoliton          , radiusVirial, &
         &                                                                                     densityScale           , densityCore  , &
         &                                                                                     massHaloMinimum0       , massCore
    !![
    <optionalArgument name="weightBy" defaultsTo="weightByMass" />
    !!]
    
    ! Return a null distribution if weighting is not by mass.
    massDistribution_ => null()
    if (weightBy_ /= weightByMass) return
    ! Compute properties of the distribution.
    call self%computeProperties(node,radiusVirial,radiusScale,radiusCore,radiusSoliton,densityCore,densityScale,massHaloMinimum0,massCore)
    ! Construct the distribution.
    allocate(massDistributionSolitonNFW :: massDistribution_)
    select type(massDistribution_)
    type is (massDistributionSolitonNFW)
       !![
       <referenceConstruct object="massDistribution_">
	 <constructor>
           massDistributionSolitonNFW(                           &amp;
           &amp; radiusScale            = radiusScale          , &amp;
           &amp; radiusCore             = radiusCore           , &amp;
           &amp; radiusSoliton          = radiusSoliton        , &amp;
           &amp; densitySolitonCentral  = densityCore          , &amp;
           &amp; densityNormalizationNFW= densityScale         , &amp;
           &amp; radiusVirial           = radiusVirial         , &amp;
           &amp; componentType          = componentTypeDarkHalo, &amp;
           &amp; massType               = massTypeDark           &amp;
           &amp;                     )
	 </constructor>
       </referenceConstruct>
       !!]
    end select
    allocate(kinematicsDistribution_)
    !![
    <referenceConstruct object="kinematicsDistribution_">
      <constructor>
        kinematicsDistributionCollisionlessTabulated()
      </constructor>
    </referenceConstruct>
    !!]
    call massDistribution_%setKinematicsDistribution(kinematicsDistribution_)
    !![
    <objectDestructor name="kinematicsDistribution_"/>
    !!]
    return
  end function solitonNFWGet

  subroutine solitonNFWComputeProperties(self,node,radiusVirial,radiusScale,radiusCore,radiusSoliton,densityCore,densityScale,massHaloMinimum0,massCore)
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
    double precision                                , intent(  out) :: radiusVirial                       , radiusScale             , &
         &                                                             radiusCore                         , radiusSoliton           , &
         &                                                             densityCore                        , densityScale            , &
         &                                                             massHaloMinimum0                   , massCore
    class           (nodeComponentBasic            ), pointer       :: basic
    class           (nodeComponentDarkMatterProfile), pointer       :: darkMatterProfile
    type            (rootFinder                    ), save          :: finder
    logical                                         , save          :: finderInitialized =.false.
    !$omp threadprivate(finder, finderInitialized)
    double precision                                , parameter     :: toleranceAbsolute = 0.0d0          , toleranceRelative   =1.0d-3
    double precision                                , parameter     :: plancksConstantBar=+plancksConstant                               & ! ℏ in units of eV s.
         &                                                                                /2.0d0                                         &
         &                                                                                /Pi                                            &
         &                                                                                /electronVolt
    double precision                                                :: massHalo                           , expansionFactor            , &
         &                                                             redshift                           , concentration              , &
         &                                                             hubbleConstant                     , hubbleConstantLittle       , &
         &                                                             OmegaMatter                        , densityMatter              , &
         &                                                             zeta_0                             , zeta_z
    
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
    zeta_0             =+self%virialDensityContrast_%densityContrast(massHalo,expansionFactor=1.0d0          )
    zeta_z             =+self%virialDensityContrast_%densityContrast(massHalo,expansionFactor=expansionFactor)
    massHaloMinimum0   =+ 32.0d0               & ! Text after equation (6) of Schive et al. (2014; PRL; 113; 1302; https://ui.adsabs.harvard.edu/abs/2014PhRvL.113z1302S).
         &              *Pi                    &
         &              /375.0d0      **0.25d0 &
         &              *zeta_0       **0.25d0 &
         &              /OmegaMatter  **0.75d0 &
         &              *densityMatter         &
         &              *speedLight   **3      &
         &              /megaParsec   **3      &
         &              /(                     &
         &                +kilo                &
         &                /megaParsec          &
         &                *hubbleConstant      &
         &                *self%massParticle   &
         &                /plancksConstantBar  &
         &               )**1.5d0
    massCore           =+0.25d0                                     & ! Equation (6) of Schive et al. (2014; PRL; 113; 1302; https://ui.adsabs.harvard.edu/abs/2014PhRvL.113z1302S).
         &              *massHaloMinimum0                           &
         &              /sqrt(expansionFactor)                      &
         &              *(zeta_z  /zeta_0          )**(1.0d0/6.0d0) &
         &              *(massHalo/massHaloMinimum0)**(1.0d0/3.0d0)
    ! Compute the core radius.
    radiusCore         =+1.6d-3                           & ! Equation (7) of Schive et al. (2014; PRL; 113; 1302; https://ui.adsabs.harvard.edu/abs/2014PhRvL.113z1302S).
         &              /(self%massParticle/1.0d-22)      &
         &              *sqrt(expansionFactor)            &
         &              /(zeta_z  /zeta_0)**(1.0d0/6.0d0) &
         &              /(massHalo/1.0d9 )**(1.0d0/3.0d0)
    ! Compute the core density normalization. Equation (3) of Schive et al. (2014; PRL; 113; 1302; https://ui.adsabs.harvard.edu/abs/2014PhRvL.113z1302S).
    densityCore       =+1.9d6                          &
         &             /(self%massParticle/1.0d-23)**2 &
         &             /radiusCore                 **4 &
         &             /expansionFactor
    ! Solve for the soliton radius.
    if (.not.finderInitialized) then
       finder=rootFinder(                                        &
            &            rootFunction     =radiusTransitionRoot, &
            &            toleranceAbsolute=toleranceAbsolute   , &
            &            toleranceRelative=toleranceRelative     &
            &           )
       call finder%rangeExpand(                                                              &
            &                  rangeExpandUpward            = 2.0d0             , &
            &                  rangeExpandDownward          = 0.5d0                  , &
            &                  rangeDownwardLimit           = 1.0d0*radiusCore             , &
            &                  rangeUpwardLimit             = 1.0d1*radiusCore             , &
            &                  rangeExpandDownwardSignExpect= rangeExpandSignExpectPositive, &
            &                  rangeExpandUpwardSignExpect  = rangeExpandSignExpectNegative, &
            &                  rangeExpandType              = rangeExpandMultiplicative      &
            &                 )
       finderInitialized=.true.
    end if
    radiusCore_  =radiusCore
    radiusScale_ =radiusScale
    densityScale_=densityScale
    densityCore_ =densityCore
    radiusSoliton=finder%find(rootGuess=2.0d0*radiusCore)    
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
