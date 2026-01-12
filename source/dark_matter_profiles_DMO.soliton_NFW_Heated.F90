!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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
  use :: Mass_Distributions      , only : enumerationNonAnalyticSolversType
  use :: Cosmology_Parameters    , only : cosmologyParametersClass
  use :: Cosmology_Functions     , only : cosmologyFunctionsClass
  use :: Virial_Density_Contrast , only : virialDensityContrastClass
  use :: Statistics_Distributions, only : distributionFunction1DNormal
  !![
  <darkMatterProfileDMO name="darkMatterProfileDMOSolitonNFWHeated">
   <description>
    A dark matter profile DMO class which builds \refClass{massDistributionSolitonNFWHeated} objects to implement the \gls{fdm}
    profile. The inner region follows the soliton solution, while the outer region transitions to a heated NFW envelope.
    The core-halo mass relation and core radius are computed following \cite{chan_diversity_2022}, while the core
    density normalization follows \cite{schive_understanding_2014}.
   </description>
  </darkMatterProfileDMO>
  !!]
  type, extends(darkMatterProfileDMOClass) :: darkMatterProfileDMOSolitonNFWHeated
     !!{
     A dark matter halo profile class implementing \gls{fdm} dark matter halos.
     !!}
     private
     double precision                                             :: massParticle
     class           (darkMatterHaloScaleClass         ), pointer :: darkMatterHaloScale_               => null()
     class           (darkMatterProfileDMOClass        ), pointer :: darkMatterProfileDMO_              => null()
     class           (darkMatterProfileHeatingClass    ), pointer :: darkMatterProfileHeating_          => null()
     class           (massDistributionClass            ), pointer :: massDistributionHeated_            => null()
     class           (cosmologyParametersClass         ), pointer :: cosmologyParameters_               => null()
     class           (cosmologyFunctionsClass          ), pointer :: cosmologyFunctions_                => null()
     class           (darkMatterParticleClass          ), pointer :: darkMatterParticle_                => null()
     class           (virialDensityContrastClass       ), pointer :: virialDensityContrast_             => null()
     type            (distributionFunction1DNormal     )          :: massCoreScatter
     type            (enumerationNonAnalyticSolversType)          :: nonAnalyticSolver
     logical                                                      :: tolerateVelocityMaximumFailure     , tolerateEnclosedMassIntegrationFailure    , &
          &                                                          toleratePotentialIntegrationFailure, velocityDispersionApproximate
     double precision                                             :: toleranceRelativeVelocityDispersion, toleranceRelativeVelocityDispersionMaximum, &
          &                                                           fractionRadiusFinalSmall          , toleranceRelativePotential
     double precision                                             :: radiusVirialPrevious               , radiusScalePrevious                       , &
          &                                                          radiusCorePrevious                 , radiusSolitonPrevious                     , &
          &                                                          densityCorePrevious                , densityScalePrevious                      , &
          &                                                          massCorePrevious                   , scatterFractional   
     integer          (kind_int8                      )           :: lastUniqueID
     integer                                                      :: randomOffsetID                     , densityCoreID                             , &
          &                                                          radiusCoreID                       , radiusSolitonID                           , &
          &                                                          massCoreNormalID                   , massCoreID                                , &
          &                                                          zetaID
   contains    
     !![
     <methods>
       <method method="computeProperties" description="Compute properties of the mass distribution."/>
       <method method="calculationReset"  description="Reset memoized calculations."                />
     </methods>
     !!]
     final     ::                      solitonNFWHeatedDestructor
     procedure :: get               => solitonNFWHeatedGet
     procedure :: autoHook          => solitonNFWHeatedAutoHook
     procedure :: calculationReset  => solitonNFWHeatedCalculationReset
     procedure :: computeProperties => solitonNFWHeatedComputeProperties
  end type darkMatterProfileDMOSolitonNFWHeated

  interface darkMatterProfileDMOSolitonNFWHeated
     !!{
     Constructors for the {\normalfont \ttfamily solitonNFWHeated} dark matter halo profile class.
     !!}
     module procedure solitonNFWHeatedConstructorParameters
     module procedure solitonNFWHeatedConstructorInternal
  end interface darkMatterProfileDMOSolitonNFWHeated

  ! Sub-module scope variable used in root finding.
  class           (darkMatterProfileDMOSolitonNFWHeated), pointer  :: self_                    => null()
  double precision                                                 :: radiusCore_, densityCore_
  !$omp threadprivate(self_, radiusCore_, densityCore_)

contains

  function solitonNFWHeatedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{darkMatterProfileDMOSolitonNFWHeated} dark matter halo profile class which takes a parameter set as input.
    !!}
    use :: Mass_Distributions, only : enumerationNonAnalyticSolversEncode
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (darkMatterProfileDMOSolitonNFWHeated)                :: self
    type            (inputParameters                     ), intent(inout) :: parameters
    class           (darkMatterProfileHeatingClass       ), pointer       :: darkMatterProfileHeating_
    class           (darkMatterHaloScaleClass            ), pointer       :: darkMatterHaloScale_
    class           (darkMatterParticleClass             ), pointer       :: darkMatterParticle_
    class           (cosmologyFunctionsClass             ), pointer       :: cosmologyFunctions_
    class           (cosmologyParametersClass            ), pointer       :: cosmologyParameters_
    class           (virialDensityContrastClass          ), pointer       :: virialDensityContrast_
    type            (varying_string                      )                :: nonAnalyticSolver
    logical                                                               :: tolerateVelocityMaximumFailure     , tolerateEnclosedMassIntegrationFailure    , &
         &                                                                   toleratePotentialIntegrationFailure, velocityDispersionApproximate
    double precision                                                      :: toleranceRelativeVelocityDispersion, toleranceRelativeVelocityDispersionMaximum, &
         &                                                                   fractionRadiusFinalSmall           , toleranceRelativePotential                , &
         &                                                                   scatterFractional

    !![
    <inputParameter>
      <name>nonAnalyticSolver</name>
      <defaultValue>var_str('fallThrough')</defaultValue>
      <source>parameters</source>
      <description>Selects how solutions are computed when no analytic solution is available. If set to ``{\normalfont \ttfamily fallThrough}'' then the solution ignoring heating is used, while if set to ``{\normalfont \ttfamily numerical}'' then numerical solvers are used to find solutions.</description>
    </inputParameter>
    <inputParameter>
      <name>tolerateEnclosedMassIntegrationFailure</name>
      <defaultValue>.false.</defaultValue>
      <source>parameters</source>
      <description>If {\normalfont \ttfamily true}, tolerate failures to find the mass enclosed as a function of radius.</description>
    </inputParameter>
    <inputParameter>
      <name>tolerateVelocityMaximumFailure</name>
      <defaultValue>.false.</defaultValue>
      <source>parameters</source>
      <description>If {\normalfont \ttfamily true}, tolerate failures to find the radius of the maximum circular velocity.</description>
    </inputParameter>
    <inputParameter>
      <name>toleratePotentialIntegrationFailure</name>
      <defaultValue>.false.</defaultValue>
      <source>parameters</source>
      <description>If {\normalfont \ttfamily true}, tolerate failures to compute the potential.</description>
    </inputParameter>
    <inputParameter>
      <name>velocityDispersionApproximate</name>
      <defaultValue>.true.</defaultValue>
      <source>parameters</source>
      <description>If {\normalfont \ttfamily true}, radial velocity dispersion is computed using an approximate method in which we assume that $\sigma_\mathrm{r}^2(r) \rightarrow \sigma_\mathrm{r}^2(r) - (2/3) \epsilon(r)$, where $\epsilon(r)$ is the specific heating energy. If {\normalfont \ttfamily false} then radial velocity dispersion is computed by numerically solving the Jeans equation.</description>
    </inputParameter>
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
    <inputParameter>
      <name>fractionRadiusFinalSmall</name>
      <defaultValue>1.0d-3</defaultValue>
      <source>parameters</source>
      <description>The initial radius is limited to be no smaller than this fraction of the final radius. This can help avoid problems in profiles that are extremely close to being disrupted.</description>
    </inputParameter>
    <inputParameter>
      <name>toleranceRelativePotential</name>
      <defaultValue>1.0d-3</defaultValue>
      <source>parameters</source>
      <description>The maximum allowed relative tolerance to use in numerical solutions for the gravitational potential in dark-matter-only density profiles before aborting.</description>
    </inputParameter>
    <inputParameter>
      <name>scatterFractional</name>
      <defaultValue>0.5d0</defaultValue>
      <source>parameters</source>
      <description>The fractional scatter in the solitonic core-halo mass relation (default corresponds to a 50\% fractional scatter).</description>
    </inputParameter>
    <objectBuilder class="darkMatterHaloScale"      name="darkMatterHaloScale_"      source="parameters"/>
    <objectBuilder class="darkMatterParticle"       name="darkMatterParticle_"       source="parameters"/>
    <objectBuilder class="darkMatterProfileHeating" name="darkMatterProfileHeating_" source="parameters"/>
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    <objectBuilder class="cosmologyParameters"      name="cosmologyParameters_"      source="parameters"/>
    <objectBuilder class="virialDensityContrast"    name="virialDensityContrast_"    source="parameters"/>
    !!]
    self = darkMatterProfileDMOSolitonNFWHeated(enumerationNonAnalyticSolversEncode(char(nonAnalyticSolver),includesPrefix=.false.),darkMatterHaloScale_,darkMatterParticle_,darkMatterProfileHeating_,cosmologyFunctions_,cosmologyParameters_,virialDensityContrast_,toleranceRelativeVelocityDispersion,toleranceRelativeVelocityDispersionMaximum,tolerateEnclosedMassIntegrationFailure,tolerateVelocityMaximumFailure,toleratePotentialIntegrationFailure,velocityDispersionApproximate,fractionRadiusFinalSmall,toleranceRelativePotential,scatterFractional)
    !![
    <inputParametersValidate source="parameters"               />
    <objectDestructor        name  ="darkMatterHaloScale_"     />
    <objectDestructor        name  ="darkMatterParticle_"      />
    <objectDestructor        name  ="darkMatterProfileHeating_"/>
    <objectDestructor        name  ="cosmologyFunctions_"      />
    <objectDestructor        name  ="cosmologyParameters_"     />
    <objectDestructor        name  ="virialDensityContrast_"   />
    !!]
    return
  end function solitonNFWHeatedConstructorParameters

  function solitonNFWHeatedConstructorInternal(nonAnalyticSolver,darkMatterHaloScale_,darkMatterParticle_,darkMatterProfileHeating_,cosmologyFunctions_,cosmologyParameters_,virialDensityContrast_,toleranceRelativeVelocityDispersion,toleranceRelativeVelocityDispersionMaximum,tolerateEnclosedMassIntegrationFailure,tolerateVelocityMaximumFailure,toleratePotentialIntegrationFailure,velocityDispersionApproximate,fractionRadiusFinalSmall,toleranceRelativePotential,scatterFractional) result(self)
    !!{
    Generic constructor for the \refClass{darkMatterProfileDMOSolitonNFWHeated} dark matter halo profile class.
    !!}
    use :: Mass_Distributions          , only : enumerationNonAnalyticSolversIsValid
    use :: Error                       , only : Component_List                      , Error_Report
    use :: Dark_Matter_Particles       , only : darkMatterParticleFuzzyDarkMatter
    use :: Numerical_Constants_Prefixes, only : kilo
    use :: Galacticus_Nodes            , only : defaultDarkMatterProfileComponent
    implicit none
    type            (darkMatterProfileDMOSolitonNFWHeated)                     :: self
    class           (darkMatterProfileHeatingClass       ), intent(in), target :: darkMatterProfileHeating_
    class           (darkMatterHaloScaleClass            ), intent(in), target :: darkMatterHaloScale_
    class           (darkMatterParticleClass             ), intent(in), target :: darkMatterParticle_
    class           (cosmologyFunctionsClass             ), intent(in), target :: cosmologyFunctions_
    class           (cosmologyParametersClass            ), intent(in), target :: cosmologyParameters_
    class           (virialDensityContrastClass          ), intent(in), target :: virialDensityContrast_
    type            (enumerationNonAnalyticSolversType   ), intent(in)         :: nonAnalyticSolver
    logical                                               , intent(in)         :: tolerateVelocityMaximumFailure        , toleratePotentialIntegrationFailure , &
         &                                                                        tolerateEnclosedMassIntegrationFailure, velocityDispersionApproximate          
    double precision                                      , intent(in)         :: toleranceRelativeVelocityDispersion, toleranceRelativeVelocityDispersionMaximum , & 
         &                                                                        fractionRadiusFinalSmall              , toleranceRelativePotential          , &
         &                                                                        scatterFractional
    !![
    <constructorAssign variables="nonAnalyticSolver,*darkMatterHaloScale_,*darkMatterParticle_,*darkMatterProfileHeating_,*cosmologyFunctions_,*cosmologyParameters_,*virialDensityContrast_,toleranceRelativeVelocityDispersion,toleranceRelativeVelocityDispersionMaximum,tolerateEnclosedMassIntegrationFailure,tolerateVelocityMaximumFailure,toleratePotentialIntegrationFailure,velocityDispersionApproximate,fractionRadiusFinalSmall,toleranceRelativePotential,scatterFractional"/>
    <addMetaProperty component="darkMatterProfile" name="solitonRandomOffset"   id="self%randomOffsetID"   isEvolvable="no"  isCreator="yes"/>
    <addMetaProperty component="darkMatterProfile" name="solitonDensityCore"    id="self%densityCoreID"    isEvolvable="no"  isCreator="yes"/>
    <addMetaProperty component="darkMatterProfile" name="solitonRadiusCore"     id="self%radiusCoreID"     isEvolvable="no"  isCreator="yes"/>
    <addMetaProperty component="darkMatterProfile" name="solitonRadiusSoliton"  id="self%radiusSolitonID"  isEvolvable="no"  isCreator="yes"/>
    <addMetaProperty component="darkMatterProfile" name="solitonMassCoreNormal" id="self%massCoreNormalID" isEvolvable="yes" isCreator="no"/>
    <addMetaProperty component="darkMatterProfile" name="solitonMassCore"       id="self%massCoreID"       isEvolvable="no"  isCreator="yes"/>
    <addMetaProperty component="darkMatterProfile" name="solitonZeta"           id="self%zetaID"           isEvolvable="no"  isCreator="yes"/>
    !!]

    self%lastUniqueID=-huge(1_kind_int8)
    self%massCoreScatter = distributionFunction1DNormal(mean=0.0d0,variance=log10(1.0d0+scatterFractional)**2)
    select type (darkMatterParticle_ => self%darkMatterParticle_)
    class is (darkMatterParticleFuzzyDarkMatter)
       self%massParticle=+darkMatterParticle_%mass()*kilo
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
    ! Validate.
    if (.not.enumerationNonAnalyticSolversIsValid(nonAnalyticSolver)) call Error_Report('invalid non-analytic solver type'//{introspection:location})
    
    ! Construct an NFW profile.
    allocate(darkMatterProfileDMONFW :: self%darkMatterProfileDMO_)
    select type (darkMatterProfileDMO_ => self%darkMatterProfileDMO_)
    type is (darkMatterProfileDMONFW)
       !![
       <referenceConstruct isResult="yes" owner="self" nameAssociated="darkMatterProfileDMO_"  object="darkMatterProfileDMO_">
	 <constructor>
	   darkMatterProfileDMONFW(                                                           &amp;
	    &amp;                  velocityDispersionUseSeriesExpansion=.true.              , &amp;
	    &amp;                  darkMatterHaloScale_                =darkMatterHaloScale_  &amp;
	    &amp;                 )
	 </constructor>
       </referenceConstruct>
       !!]
    end select
    return
  end function solitonNFWHeatedConstructorInternal

  subroutine solitonNFWHeatedAutoHook(self)
    !!{
    Attach to the calculation reset event.
    !!}
    use :: Events_Hooks, only : calculationResetEvent, openMPThreadBindingAllLevels
    implicit none
    class(darkMatterProfileDMOSolitonNFWHeated), intent(inout) :: self

    call calculationResetEvent%attach(self,solitonNFWHeatedCalculationReset,openMPThreadBindingAllLevels,label='darkMatterProfileDMOSolitonNFWHeated')
    return
  end subroutine solitonNFWHeatedAutoHook

  subroutine solitonNFWHeatedDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily solitonNFWHeated} dark matter halo profile class.
    !!}
    use :: Events_Hooks, only : calculationResetEvent
    implicit none
    type(darkMatterProfileDMOSolitonNFWHeated), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_"     />
    <objectDestructor name="self%darkMatterParticle_"      />
    <objectDestructor name="self%darkMatterProfileDMO_"    />
    <objectDestructor name="self%darkMatterProfileHeating_"/>
    <objectDestructor name="self%massDistributionHeated_"  />
    <objectDestructor name="self%cosmologyFunctions_"      />
    <objectDestructor name="self%cosmologyParameters_"     />
    <objectDestructor name="self%virialDensityContrast_"   />
    !!]
    if (calculationResetEvent%isAttached(self,solitonNFWHeatedCalculationReset)) call calculationResetEvent%detach(self,solitonNFWHeatedCalculationReset)
    return
  end subroutine solitonNFWHeatedDestructor

  subroutine solitonNFWHeatedCalculationReset(self,node,uniqueID)
    !!{
    Reset the dark matter profile calculation.
    !!}
    use :: Kind_Numbers, only : kind_int8
    implicit none
    class  (darkMatterProfileDMOSolitonNFWHeated), intent(inout) :: self
    type   (treeNode                            ), intent(inout) :: node
    integer(kind_int8                           ), intent(in   ) :: uniqueID
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
  end subroutine solitonNFWHeatedCalculationReset

  function solitonNFWHeatedGet(self,node,weightBy,weightIndex) result(massDistribution_)
    !!{
    Return the soliton plus NFW fuzzy dark matter mass distribution for the given {\normalfont \ttfamily node}.
    !!}
    use :: Error                     , only : Error_Report
    use :: Galacticus_Nodes          , only : nodeComponentBasic
    use :: Galactic_Structure_Options, only : componentTypeDarkHalo           , massTypeDark                          , weightByMass
    use :: Mass_Distributions        , only : massDistributionSolitonNFWHeated, kinematicsDistributionSolitonNFWHeated, kinematicsDistributionNFW   , &
         &                                    kinematicsDistributionClass     , massDistributionSpherical             , massDistributionNFW         , &
         &                                    kinematicsDistributionHeated    , massDistributionSphericalHeated       , massDistributionHeatingClass
    implicit none
    class           (darkMatterProfileDMOSolitonNFWHeated), intent(inout)           :: self
    class           (massDistributionClass               ), pointer                 :: massDistribution_
    class           (kinematicsDistributionClass         ), pointer                 :: kinematicsDistribution_, kinematicsDistributionNFW_
    class           (massDistributionClass               ), pointer                 :: massDistributionNFW_
    class           (massDistributionHeatingClass        ), pointer                 :: massDistributionHeating_
    type            (treeNode                            ), intent(inout)           :: node
    class           (nodeComponentBasic                  ), pointer                 :: basic
    type            (enumerationWeightByType             ), intent(in   ), optional :: weightBy
    integer                                               , intent(in   ), optional :: weightIndex
    type            (enumerationWeightByType             )                          :: weightBy_
    double precision                                                                :: radiusCore             , radiusScale , &
         &                                                                             radiusSoliton          , radiusVirial, &
         &                                                                             densityScale           , densityCore , &
         &                                                                             massCore
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
    if (radiusSoliton <= 0.0d0 .or. densityCore <= 0.0d0) then
       ! No soliton - build a simple NFW mass distribution.
       allocate(massDistributionNFW       :: massDistributionNFW_      )
       select type(massDistributionNFW_)
       type is (massDistributionNFW)
          basic => node%basic()
          !![
	  <referenceConstruct object="massDistributionNFW_">
	    <constructor>
              massDistributionNFW(                                             &amp;
               &amp;              mass         =basic%mass                 (), &amp;
               &amp;              virialRadius =      radiusVirial           , &amp;
               &amp;              scaleLength  =      radiusScale            , &amp;
               &amp;              componentType=      componentTypeDarkHalo  , &amp;
               &amp;              massType     =      massTypeDark             &amp;
               &amp;             )
	    </constructor>
	  </referenceConstruct>
          !!]
       end select
       allocate(kinematicsDistributionNFW :: kinematicsDistributionNFW_)
       select type(kinematicsDistributionNFW_)
       type is (kinematicsDistributionNFW)
          !![
	  <referenceConstruct object="kinematicsDistributionNFW_">
	    <constructor>
              kinematicsDistributionNFW(                                                                 &amp;
	      &amp;                    useSeriesApproximation=.true. &amp;
	      &amp;                   )
	    </constructor>
	  </referenceConstruct>
          !!]
       end select
       call massDistributionNFW_%setKinematicsDistribution(kinematicsDistributionNFW_)
       !![
       <objectDestructor name="kinematicsDistributionNFW_"/>
       !!]

       allocate(massDistributionSphericalHeated :: massDistribution_)
       allocate(kinematicsDistributionHeated :: kinematicsDistribution_)
       select type(massDistribution_)
       type is (massDistributionSphericalHeated)
          massDistributionHeating_  => self%darkMatterProfileHeating_%get(node                     )
          select type (massDistributionNFW_)
          class is (massDistributionSpherical)
             !![
             <referenceConstruct object="massDistribution_">
                <constructor>
                 massDistributionSphericalHeated(                                                                                    &amp;
               &amp;                          nonAnalyticSolver                     =self%nonAnalyticSolver                     , &amp;
               &amp;                          tolerateVelocityMaximumFailure        =self%tolerateVelocityMaximumFailure        , &amp;
               &amp;                          tolerateEnclosedMassIntegrationFailure=self%tolerateEnclosedMassIntegrationFailure, &amp;
               &amp;                          toleratePotentialIntegrationFailure   =self%toleratePotentialIntegrationFailure   , &amp;
               &amp;                          fractionRadiusFinalSmall              =self%fractionRadiusFinalSmall              , &amp;
               &amp;                          toleranceRelativePotential            =self%toleranceRelativePotential            , &amp;
               &amp;                          massDistribution_                     =     massDistributionNFW_                  , &amp;
               &amp;                          massDistributionHeating_              =     massDistributionHeating_              , &amp;
               &amp;                          componentType                         =     componentTypeDarkHalo                 , &amp;
               &amp;                          massType                              =     massTypeDark                            &amp;
               &amp;                         )
                </constructor>
             </referenceConstruct>
             !!]
          class default
             call Error_Report('expected a spherical mass distribution'//{introspection:location})
          end select
       end select
       select type (kinematicsDistribution_)
       type is (kinematicsDistributionHeated)
          !![
          <referenceConstruct object="kinematicsDistribution_">
            <constructor>
               kinematicsDistributionHeated(                                                                                            &amp;
         &amp;                       nonAnalyticSolver                         =self%nonAnalyticSolver                         , &amp;
         &amp;                       velocityDispersionApproximate             =self%velocityDispersionApproximate             , &amp;
         &amp;                       toleranceRelativeVelocityDispersion       =self%toleranceRelativeVelocityDispersion       , &amp;
         &amp;                       toleranceRelativeVelocityDispersionMaximum=self%toleranceRelativeVelocityDispersionMaximum  &amp;
	 &amp;                      )
            </constructor>
          </referenceConstruct>
          !!]
       end select
       call massDistribution_%setKinematicsDistribution(kinematicsDistribution_)
       !![
       <objectDestructor name="kinematicsDistribution_"/>
       !!]
    else
       ! Build a soliton-NFW Heated mass distribution.
       allocate(massDistributionSolitonNFWHeated       :: massDistribution_      )
       allocate(kinematicsDistributionSolitonNFWHeated :: kinematicsDistribution_)
       select type(massDistribution_)
       type is (massDistributionSolitonNFWHeated)
          select type (massDistributionHeated_ => self%massDistributionHeated_)
          class is (massDistributionSpherical)
          !![
          <referenceConstruct object="massDistribution_">
            <constructor>
              massDistributionSolitonNFWHeated(                                                       &amp;
               &amp;                           radiusCore             =radiusCore                   , &amp;
               &amp;                           radiusSoliton          =radiusSoliton                , &amp;
               &amp;                           densitySolitonCentral  =densityCore                  , &amp;
               &amp;                           massDistributionHeated_=massDistributionHeated_      , &amp;
               &amp;                           componentType          =componentTypeDarkHalo        , &amp;
               &amp;                           massType               =massTypeDark                   &amp;
               &amp;                          )
            </constructor>
          </referenceConstruct>
          !!]
          class default
             call Error_Report('a spherically-symmetric mass distribution is required'//{introspection:location})
          end select
       end select
       select type (kinematicsDistribution_)
       type is (kinematicsDistributionSolitonNFWHeated)
          !![
          <referenceConstruct object="kinematicsDistribution_">
            <constructor>
              kinematicsDistributionSolitonNFWHeated(                                                                                            &amp;
            &amp;                                 toleranceRelativeVelocityDispersion       =self%toleranceRelativeVelocityDispersion       , &amp;
            &amp;                                 toleranceRelativeVelocityDispersionMaximum=self%toleranceRelativeVelocityDispersionMaximum  &amp;
            &amp;                                )
            </constructor>
          </referenceConstruct>
          !!]
       end select
       call massDistribution_%setKinematicsDistribution(kinematicsDistribution_)
       !![
       <objectDestructor name="kinematicsDistribution_"/>
       !!]
    end if
    return
  end function solitonNFWHeatedGet

  subroutine solitonNFWHeatedComputeProperties(self,node,radiusVirial,radiusScale,radiusCore,radiusSoliton,densityCore,densityScale,massCore,weightBy,weightIndex)
    use :: Galacticus_Nodes          , only : treeNode                       , nodeComponentBasic         , nodeComponentDarkMatterProfile
    use :: Galactic_Structure_Options, only : componentTypeDarkHalo          , massTypeDark               , weightByMass
    use :: Numerical_Constants_Math  , only : Pi
    use :: Root_Finder               , only : rootFinder                     , rangeExpandMultiplicative  , rangeExpandSignExpectPositive , rangeExpandSignExpectNegative
    use :: Mass_Distributions        , only : massDistributionSphericalHeated, kinematicsDistributionHeated, massDistributionSpherical    , massDistributionHeatingClass
    implicit none
    class           (darkMatterProfileDMOSolitonNFWHeated), intent(inout), target   :: self
    type            (treeNode                            ), intent(inout)           :: node
    type            (enumerationWeightByType             ), intent(in   ), optional :: weightBy
    integer                                               , intent(in   ), optional :: weightIndex
    type            (enumerationWeightByType             )                          :: weightBy_
    double precision                                      , intent(  out)           :: radiusVirial              , radiusScale               , &
         &                                                                             radiusCore                , radiusSoliton             , &
         &                                                                             densityCore               , densityScale              , &
         &                                                                             massCore
    class           (nodeComponentBasic                  ), pointer                 :: basic
    class           (nodeComponentDarkMatterProfile      ), pointer                 :: darkMatterProfile
    class           (massDistributionClass               ), pointer                 :: massDistributionDecorated
    class           (massDistributionHeatingClass        ), pointer                 :: massDistributionHeating_
    type            (rootFinder                          ), save                    :: finder
    logical                                               , save                    :: finderInitialized =.false.
    !$omp threadprivate(finder, finderInitialized)
    double precision                                      , parameter               :: toleranceAbsolute = 0.0d0 , toleranceRelative = 1.0d-3
    double precision                                                                :: massHalo                  , expansionFactor           , &
         &                                                                             redshift                  , concentration             , &
         &                                                                             randomOffset              , massCoreNormal            , &
         &                                                                             zeta_0                    , zeta_z
    integer                                                                         :: sampleCountMaximum=50
    integer                                                                         :: status                    , sampleCount
         
    !![
    <optionalArgument name="weightBy" defaultsTo="weightByMass" />
    !!]

    ! Get required components.
    basic             => node%basic            ()
    darkMatterProfile => node%darkMatterProfile()
    call darkMatterProfile%floatRank0MetaPropertySet(self%radiusSolitonID,-1.0d0)
    call darkMatterProfile%floatRank0MetaPropertySet(self%massCoreID     ,-1.0d0)
    call darkMatterProfile%floatRank0MetaPropertySet(self%radiusCoreID   ,-1.0d0)
    call darkMatterProfile%floatRank0MetaPropertySet(self%densityCoreID  ,-1.0d0)
    ! Extract basic properties of the node.
    expansionFactor=+self             %cosmologyFunctions_% expansionFactor            (basic%time           ())
    redshift       =+self             %cosmologyFunctions_ %redshiftFromExpansionFactor(      expansionFactor  )
    massHalo       =+basic                                 %mass                       (                       )
    zeta_0         =+self%virialDensityContrast_%densityContrast(massHalo,expansionFactor=1.0d0                )
    zeta_z         =+self%virialDensityContrast_%densityContrast(massHalo,expansionFactor=expansionFactor      )
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

    ! If weighting is not by mass, return a null profile.
    if (weightBy_ /= weightByMass) return
    ! Create the mass distribution.
    allocate(massDistributionSphericalHeated :: self%massDistributionHeated_)
    select type(massDistributionHeated_ => self%massDistributionHeated_)
    type is (massDistributionSphericalHeated)
       massDistributionDecorated => self%darkMatterProfileDMO_    %get(node,weightBy,weightIndex)
       massDistributionHeating_  => self%darkMatterProfileHeating_%get(node                     )
       select type (massDistributionDecorated)
       class is (massDistributionSpherical)
          !![
	  <referenceConstruct isResult="yes" owner="self" nameAssociated="massDistributionHeated_" object="massDistributionHeated_">
	    <constructor>
              massDistributionSphericalHeated(                                                                                    &amp;
               &amp;                          nonAnalyticSolver                     =self%nonAnalyticSolver                     , &amp;
	       &amp;                          tolerateVelocityMaximumFailure        =self%tolerateVelocityMaximumFailure        , &amp;
	       &amp;                          tolerateEnclosedMassIntegrationFailure=self%tolerateEnclosedMassIntegrationFailure, &amp;
	       &amp;                          toleratePotentialIntegrationFailure   =self%toleratePotentialIntegrationFailure   , &amp;
	       &amp;                          fractionRadiusFinalSmall              =self%fractionRadiusFinalSmall              , &amp;
	       &amp;                          toleranceRelativePotential            =self%toleranceRelativePotential            , &amp;
               &amp;                          massDistribution_                     =     massDistributionDecorated             , &amp;
               &amp;                          massDistributionHeating_              =     massDistributionHeating_              , &amp;
               &amp;                          componentType                         =     componentTypeDarkHalo                 , &amp;
               &amp;                          massType                              =     massTypeDark                            &amp;
               &amp;                         )
	    </constructor>
	  </referenceConstruct>
	  <objectDestructor name="massDistributionDecorated"/>
	  <objectDestructor name="massDistributionHeating_" />
          !!]
       class default
          call Error_Report('expected a spherical mass distribution'//{introspection:location})
       end select
    end select

    ! Compute the core mass.
    massCoreNormal =+darkMatterProfile%floatRank0MetaPropertyGet(self%massCoreNormalID)

    ! Solve for the soliton radius.
    self_ => self
    if (.not.finderInitialized) then
       finder=rootFinder(                                        &
            &            rootFunction     =radiusTransitionRoot, &
            &            toleranceAbsolute=toleranceAbsolute   , &
            &            toleranceRelative=toleranceRelative     &
            &           )
       finderInitialized=.true.
    end if
    do sampleCount=1,sampleCountMaximum
       ! Find the random offset in the core mass.
       if (sampleCount == 1) then
          randomOffset       =darkMatterProfile%floatRank0MetaPropertyGet(self%randomOffsetID)
          if (randomOffset == 0.0d0) &
               & randomOffset=self%massCoreScatter%sample(randomNumberGenerator_=node%hostTree%randomNumberGenerator_)
       else
          randomOffset       =self%massCoreScatter%sample(randomNumberGenerator_=node%hostTree%randomNumberGenerator_)
       end if
       ! Find the core mass, including the random offset.
       massCore           =+massCoreNormal       &
            &              *10.0d0**randomOffset
       ! Compute the core radius.
       radiusCore         =+5.5d6                           & ! Equation (14) of Chan et al. (2022; MNRAS; 551; 943; https://ui.adsabs.harvard.edu/abs/2022MNRAS.511..943C).
            &              /(self%massParticle/1.0d-23)**2  &
            &              /expansionFactor                 &
            &              /massCore
       ! Compute the core density normalization.
       densityCore       =+massCore                         & ! Equation (3) of Schive et al. (2014; PRL; 113; 1302; https://ui.adsabs.harvard.edu/abs/2014PhRvL.113z1302S).
            &             /0.413d0                          &
            &             /radiusCore                  **3  &
            &             /Pi
       radiusCore_        =radiusCore
       densityCore_       =densityCore
       call finder%rangeExpand(                                                             &
            &                  rangeExpandUpward            =2.0d0                        , &
            &                  rangeExpandDownward          =0.5d0                        , &
            &                  rangeDownwardLimit           =1.0d0*radiusCore             , &
            &                  rangeUpwardLimit             =1.0d1*radiusCore             , &
            &                  rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive, &
            &                  rangeExpandUpwardSignExpect  =rangeExpandSignExpectNegative, &
            &                  rangeExpandType              =rangeExpandMultiplicative      &
            &                 )
       radiusSoliton=finder%find(rootGuess=3.0d0*radiusCore,status=status)
       if (status == errorStatusSuccess) then
           call darkMatterProfile%floatRank0MetaPropertySet(self%randomOffsetID ,randomOffset )
           call darkMatterProfile%floatRank0MetaPropertySet(self%radiusSolitonID,radiusSoliton)
           call darkMatterProfile%floatRank0MetaPropertySet(self%massCoreID     ,massCore     )
           exit
       end if
    end do
    call darkMatterProfile%floatRank0MetaPropertySet(self%densityCoreID,densityCore       )
    call darkMatterProfile%floatRank0MetaPropertySet(self%radiusCoreID ,radiusCore        )
    call darkMatterProfile%floatRank0MetaPropertySet(self%zetaID       ,zeta_z     /zeta_0)
    return
  end subroutine solitonNFWHeatedComputeProperties

  double precision function radiusTransitionRoot(radius) result(f)
    !!{
    Root function used in seeking the transition radius in fuzzy dark matter profiles.
    !!}
    use :: Coordinates          , only          : coordinateCartesian, assignment(=)
    implicit none
    double precision            , intent(in)    :: radius
    type   (coordinateCartesian)                :: coordinates

    coordinates=[radius,0.0d0,0.0d0]

    f      =+                densityCore_                      &
         &  /(1.0d0+0.091d0*(radius/radiusCore_ )**2)**8       &
         &  -self_%massDistributionHeated_%density(coordinates)
    return
  end function radiusTransitionRoot
