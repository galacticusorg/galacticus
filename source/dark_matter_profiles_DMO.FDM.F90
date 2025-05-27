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
  An implementation of fuzzy dark matter halo profiles using the Soliton and NFW mass distribution.
  !!}

  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass
  use :: Dark_Matter_Particles  , only : darkMatterParticleClass
  use :: Cosmology_Parameters   , only : cosmologyParametersClass
  use :: Cosmology_Functions    , only : cosmologyFunctionsClass

  !![
  <darkMatterProfileDMO name="darkMatterProfileDMOSolitonNFW">
   <description>
    A dark matter profile DMO class which builds massDistributionSolitonNFW objects to implement the fuzzy dark matter (FDM) profile.
   </description>
  </darkMatterProfileDMO>
  !!]
  type, extends(darkMatterProfileDMOClass) :: darkMatterProfileDMOSolitonNFW
     !!{
     A dark matter halo profile class implementing FDM dark matter halos.
     !!}
     private
     double precision                                    :: m_psi, redshift
     double precision                                    :: r_c, r_s, r_sol, rho_s, rho_0
     class  (darkMatterHaloScaleClass)         , pointer :: darkMatterHaloScale_ => null()
     class  (cosmologyParametersClass)         , pointer :: cosmologyParameters_ => null()
     class  (cosmologyFunctionsClass)          , pointer :: cosmologyFunctions_  => null()
     class  (darkMatterParticleClass)          , pointer :: darkMatterParticle_  => null()
     logical                                             :: paramsUseInput
   contains
     final     ::               solitonNFWDestructor
     procedure :: get        => solitonNFWGet
     procedure :: paraminput => solitonNFWInputParameters
  end type darkMatterProfileDMOSolitonNFW

  interface darkMatterProfileDMOSolitonNFW
     !!{
     Constructors for the {\normalfont \ttfamily solitonNFW} dark matter halo profile class.
     !!}
     module procedure SolitonNFWConstructorParameters
     module procedure SolitonNFWConstructorInternal
  end interface darkMatterProfileDMOSolitonNFW

  class   (darkMatterProfileDMOSolitonNFW      ), pointer :: self_  => null()
  !$omp threadprivate(self_)

contains

  function SolitonNFWConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily solitonNFW} dark matter halo profile class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (darkMatterProfileDMOSolitonNFW   )                :: self
    type   (inputParameters                  ), intent(inout) :: parameters
    class  (darkMatterHaloScaleClass         ), pointer       :: darkMatterHaloScale_
    class  (darkMatterParticleClass          ), pointer       :: darkMatterParticle_
    class  (cosmologyFunctionsClass          ), pointer       :: cosmologyFunctions_
    class  (cosmologyParametersClass         ), pointer       :: cosmologyParameters_
    logical                                                   :: paramsUseInput

    !![
    <inputParameter>
      <name>paramsUseInput</name>
      <defaultValue>.true.</defaultValue>
      <source>parameters</source>
      <description>If {\normalfont \ttfamily true}, parameters of the FDM halo are using the input parameters.</description>
    </inputParameter>
    <objectBuilder class="darkMatterHaloScale"  name="darkMatterHaloScale_" source="parameters"/>
    <objectBuilder class="darkMatterParticle"   name="darkMatterParticle_"  source="parameters"/>
    <objectBuilder class="cosmologyFunctions"   name="cosmologyFunctions_"  source="parameters"/>
    <objectBuilder class="cosmologyParameters"  name="cosmologyParameters_" source="parameters"/>
    !!]
    self = darkMatterProfileDMOSolitonNFW(paramsUseInput, darkMatterHaloScale_, darkMatterParticle_, cosmologyFunctions_,cosmologyParameters_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    <objectDestructor name="darkMatterParticle_" />
    <objectDestructor name="cosmologyFunctions_" />
    <objectDestructor name="cosmologyParameters_"/>
    !!]
    return
  end function SolitonNFWConstructorParameters

  function SolitonNFWConstructorInternal(paramsUseInput, darkMatterHaloScale_, darkMatterParticle_, cosmologyFunctions_,cosmologyParameters_) result(self)
    !!{
    Generic constructor for the {\normalfont \ttfamily solitonNFW} dark matter halo profile class.
    !!}
    use :: Error                         , only : Component_List                   , Error_Report
    use :: Dark_Matter_Particles         , only : darkMatterParticleFuzzyDarkMatter
    use :: Numerical_Constants_Prefixes  , only : kilo
    use :: Galacticus_Nodes              , only : defaultDarkMatterProfileComponent
    implicit none
    type   (darkMatterProfileDMOSolitonNFW   )                           :: self
    class  (darkMatterHaloScaleClass         ), intent(in), target       :: darkMatterHaloScale_
    class  (darkMatterParticleClass          ), intent(in), target       :: darkMatterParticle_
    class  (cosmologyFunctionsClass          ), intent(in), target       :: cosmologyFunctions_
    class  (cosmologyParametersClass         ), intent(in), target       :: cosmologyParameters_
    logical                                   , intent(in)               :: paramsUseInput
    !![
    <constructorAssign variables="paramsUseInput,*darkMatterHaloScale_,*darkMatterParticle_,*cosmologyFunctions_, *cosmologyParameters_"/>
    !!]

    select type (darkMatterParticle__ => self%darkMatterParticle_)
    class is (darkMatterParticleFuzzyDarkMatter)
       self%m_psi=+darkMatterParticle__%mass()*kilo
    end select

    if (.not.defaultDarkMatterProfileComponent%scaleIsGettable())                                                        &
        & call Error_Report                                                                                             &
        &      (                                                                                                        &
        &       'SolitonNFW dark matter profile requires a dark matter profile component with a gettable "scale" property.'// &
        &       Component_List(                                                                                         &
        &                      'darkMatterProfile'                                                                   ,  &
        &                      defaultDarkMatterProfileComponent%scaleAttributeMatch(requireGettable=.true.)            &
        &                     )                                                                                      // &
        &      {introspection:location}                                                                                 &
        &      )
    return
  end function SolitonNFWConstructorInternal

  subroutine solitonNFWDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily solitonNFW} dark matter halo profile class.
    !!}
    implicit none
    type(darkMatterProfileDMOSolitonNFW), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_" />
    <objectDestructor name="self%darkMatterParticle_" />
    <objectDestructor name="self%cosmologyFunctions_" />
    <objectDestructor name="self%cosmologyParameters_"/>
    !!]
    return
  end subroutine solitonNFWDestructor

  double precision function fdmRadiusTransitionRoot(x) result(f)
    !!{
    Root function used in seeking the transition radius in fuzzy dark matter profiles.
    !!}
    implicit none
    double precision                       , intent(in   ) :: x

    f = (1.0d0 + 0.091d0 * (x/self_%r_c)**2.0d0)**(-8.0d0) * (x/self_%r_s) * (1.0d0 + x/self_%r_s)**2.0d0 - (self_%rho_s/self_%rho_0)

  end function fdmRadiusTransitionRoot

  subroutine computeDMproperties(self, M_h, r_c, r_sol, rho_0)
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Units       , only : electronVolt
    use :: Numerical_Constants_Astronomical, only : parsec
    use :: Numerical_Constants_Physical    , only : speedLight, plancksConstant
    use :: Cosmology_Parameters            , only : hubbleUnitsStandard, hubbleUnitsLittleH
    use :: Root_Finder                     , only : rootFinder, rangeExpandMultiplicative, rangeExpandSignExpectPositive, rangeExpandSignExpectNegative
    implicit none
    class  (darkMatterProfileDMOSolitonNFW), intent(inout), target :: self
    type   (rootFinder                    ), save                  :: finder
    logical                                , save                  :: finderInitialized=.false.
    double precision                       , intent(in)            :: M_h
    double precision                       , intent(out)           :: r_c, r_sol, rho_0
    double precision                                               :: z, a, H0, h, Omega_m0, Omega_m, rho_m0, rho_m, zeta_0, zeta_z, M_min0, M_c
    double precision                                               :: x_c, hbar
    double precision                       , parameter             :: toleranceAbsolute=0.0d0, toleranceRelative=1.0d-3
    
    self_          => self

    hbar     = plancksConstant/2.0d0/Pi/electronVolt
    z        = self%redshift
    H0       = self%cosmologyParameters_%HubbleConstant(hubbleUnitsStandard)
    h        = self%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)
    Omega_m0 = self%cosmologyParameters_%OmegaMatter()
    Omega_m  = self%cosmologyParameters_%OmegaMatter() * (1.0d0 + z)**3
    rho_m0   = self%cosmologyParameters_%densityCritical() * Omega_m0
    rho_m    = rho_m0 * (1.0d0 + z)**3
    a        = 1.0d0 / (1.0d0 + z)

    zeta_0   = (18.0d0*Pi**2 + 82.0d0*(Omega_m0-1.0d0) - 39.0d0* (Omega_m0 - 1.0d0)**2)/Omega_m0
    zeta_z   = (18.0d0*Pi**2 + 82.0d0*(Omega_m -1.0d0) - 39.0d0* (Omega_m - 1.0d0)**2)/Omega_m

    M_min0 = 375.0d0**(-0.25d0) * 32.0d0 * Pi * zeta_0**0.25d0 * rho_m0 * h**2 * speedLight**3 * (H0 * self%m_psi / hbar)**(-1.5d0) * Omega_m0**(-0.75d0)
    M_c    = 0.25d0 * a**(-0.5d0) * (zeta_z / zeta_0)**(1.0d0/6.0d0) * (M_h / M_min0)**(1.0d0/3.0d0) * M_min0

    r_c      = (1.6d-3/(self%m_psi / 1.0d-22))*(a**0.5d0)*((zeta_z/zeta_0)**(-1.0d0/6.0d0))*((M_h/1.0d9)**(-1.0d0/3.0d0))
    
    x_c = r_c * (3.0d0 * H0**2 * Omega_m0 / (8.0d0 * Pi))**(-0.25d0) * (self%m_psi / hbar)**(-0.5d0)
    rho_0 = 1.9d9 * 1.0d9 * (self%m_psi / 1.0d-23)**(-2.0d0) * x_c**(-4.0d0) / a

    self%r_c   = r_c
    self%rho_0 = rho_0

    if (.not.finderInitialized) then
       finder=rootFinder(                                             &
            &            rootFunction     =fdmRadiusTransitionRoot  , &
            &            toleranceAbsolute=toleranceAbsolute        , &
            &            toleranceRelative=toleranceRelative          &
            &           )
       call finder%rangeExpand(                                                              &
            &                  rangeExpandUpward            = 5.0d0*r_c                    , &
            &                  rangeExpandDownward          = r_c                          , &
            &                  rangeExpandDownwardSignExpect= rangeExpandSignExpectPositive, &
            &                  rangeExpandUpwardSignExpect  = rangeExpandSignExpectNegative, &
            &                  rangeExpandType              = rangeExpandMultiplicative    , &
            &                  rangeDownwardLimit           = 1.0d0 * r_c                  , &
            &                  rangeUpwardLimit             = 1.0d+1 * r_c                   &
            &                 )
       finderInitialized=.true.
    end if

    r_sol = finder                  %find                  (rootGuess=2.0d0*r_c)
    self%r_sol = r_sol

  end subroutine computeDMproperties

  function solitonNFWGet(self,node,weightBy,weightIndex) result(massDistribution_)
    !!{
    Return the Soliton and NFW fuzzy dark matter mass distribution for the given {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes          , only : nodeComponentBasic        , nodeComponentDarkMatterProfile
    use :: Galactic_Structure_Options, only : componentTypeDarkHalo     , massTypeDark                  , weightByMass
    use :: Mass_Distributions        , only : massDistributionSolitonNFW, kinematicsDistributionCollisionlessTabulated
    use :: Numerical_Constants_Math  , only : Pi
    implicit none
    class  (darkMatterProfileDMOSolitonNFW              ), intent(inout)           :: self
    class  (massDistributionClass                       ), pointer                 :: massDistribution_
    type   (kinematicsDistributionCollisionlessTabulated), pointer                 :: kinematicsDistribution_
    type   (treeNode                                    ), intent(inout)           :: node
    type   (enumerationWeightByType                     ), intent(in   ), optional :: weightBy
    integer                                              , intent(in   ), optional :: weightIndex
    class  (nodeComponentBasic                          ), pointer                 :: basic
    class  (nodeComponentDarkMatterProfile              ), pointer                 :: darkMatterProfile
    type   (enumerationWeightByType                     )                          :: weightBy_

    double precision :: r_c, r_s, r_sol, r_vir, rho_s, rho_0, concentration
    double precision :: M_h
    basic             => node %basic()
    darkMatterProfile => node %darkMatterProfile()

    self%redshift = self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(basic%time()))
    M_h           =+basic%mass()
    r_vir         =+self%darkMatterHaloScale_%radiusVirial(node)

    if  (self%paramsUseInput) then
            concentration =+r_vir/self%r_s
    else
            r_s           =+darkMatterProfile%scale()
            concentration =+r_vir/r_s
            rho_s         =+3.0d0*M_h/4.0d0/Pi/(+log(1.0d0+concentration)-concentration/(1.0d0+concentration))/r_s**3
            self%r_s      = r_s
            self%rho_s    = rho_s
            call computeDMproperties(self, M_h, r_c, r_sol, rho_0)
    end if

    if (present(weightBy)) then
      weightBy_ = weightBy
    else
      weightBy_ = weightByMass
    end if

    massDistribution_ => null()
    if (weightBy_ /= weightByMass) return

    allocate(massDistributionSolitonNFW :: massDistribution_)
    select type(massDistribution_)
    type is (massDistributionSolitonNFW)
      basic             => node%basic            ()
      darkMatterProfile => node%darkMatterProfile()
      !![
      <referenceConstruct object="massDistribution_">
   <constructor>
          massDistributionSolitonNFW(                                 &amp;
          &amp; radiusScale           = self%r_s,                     &amp;
          &amp; radiusCore            = self%r_c,                     &amp;
          &amp; radiusSol             = self%r_sol,                   &amp;
          &amp; solitonCentralDensity = self%rho_0,                   &amp;
          &amp; densityNormalization  = self%rho_s,                   &amp;
          &amp; radiusVirial          = r_vir,                        &amp;
          &amp; concentration         = concentration,                &amp;
          &amp; componentType         = componentTypeDarkHalo,        &amp;
          &amp; massType              = massTypeDark                  &amp;
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

  subroutine solitonNFWInputParameters(self, r_c, r_s, r_sol, rho_0, rho_s)
    implicit none
    class  (darkMatterProfileDMOSolitonNFW), intent(inout) :: self
    double precision                       , intent(in   ) :: r_c, r_s, r_sol, rho_0, rho_s
    self%r_c   = r_c
    self%r_s   = r_s
    self%r_sol = r_sol
    self%rho_0 = rho_0
    self%rho_s = rho_s
    return
  end subroutine solitonNFWInputParameters

