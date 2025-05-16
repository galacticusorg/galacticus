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
  An implementation of dark matter halo profile concentrations using the
  \cite{diemer_universal_2014} algorithm.
  !!}

  use :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass, criticalOverdensityClass
  use :: Cosmology_Functions       , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters      , only : cosmologyParametersClass
  use :: Dark_Matter_Profiles_DMO  , only : darkMatterProfileDMONFW
  use :: Power_Spectra             , only : powerSpectrumClass
  use :: Virial_Density_Contrast   , only : virialDensityContrastFixed

  !![
  <darkMatterProfileConcentration name="darkMatterProfileConcentrationDiemerKravtsov2014">
   <description>
    A dark matter profile concentration class in which the concentration is computed using a fitting function from
    \cite{diemer_universal_2014}:
    \begin{equation}
    c = {c_\mathrm{min} \over 2} \left[ \left({\nu\over\nu_\mathrm{min}}\right)^{-\alpha} +
    \left({\nu\over\nu_\mathrm{min}}\right)^{\beta} \right],
    \end{equation}
    where $c_\mathrm{min}=\phi_0+\phi_1 n$, $\nu_\mathrm{min}=\eta_0+\eta_1 n$, $n$ is the logarithmic slope of the linear
    power spectrum at wavenumber $k = \kappa 2 \pi / R$, $R$ is the comoving Lagrangian radius of the halo, $R=[3 M / 4 \pi
    \rho_\mathrm{M}(z=0)]^{1/3}$, and $\nu=\delta_\mathrm{crit}(t)/\sigma(M)$ is the peak height parameter. The numerical
    parameters $(\kappa,\phi_0,\phi_1,\eta_0,\eta_1,\alpha,\beta)$ are set by the parameters {\normalfont \ttfamily [kappa]},
    {\normalfont \ttfamily [phi0]}, {\normalfont \ttfamily [phi1]}, {\normalfont \ttfamily [eta0]}, {\normalfont \ttfamily
    [eta1]}, {\normalfont \ttfamily [alpha]}, {\normalfont \ttfamily [beta]}, respectively, and default to the values given in
    Table 3 of \cite{diemer_universal_2014} for the median relation, namely $(0.69,6.58,1.37,6.82,1.42,1.12,1.69)$.
   </description>
   <deepCopy>
    <functionClass variables="virialDensityContrastDefinition_, darkMatterProfileDMODefinition_"/>
   </deepCopy>
   <stateStorable>
    <functionClass variables="virialDensityContrastDefinition_, darkMatterProfileDMODefinition_"/>
   </stateStorable>
  </darkMatterProfileConcentration>
  !!]
  type, extends(darkMatterProfileConcentrationClass) :: darkMatterProfileConcentrationDiemerKravtsov2014
     !!{
     A dark matter halo profile concentration class implementing the algorithm of
     \cite{diemer_universal_2014}.
     !!}
     private
     class           (cosmologyFunctionsClass      ), pointer :: cosmologyFunctions_              => null()
     class           (cosmologyParametersClass     ), pointer :: cosmologyParameters_             => null()
     class           (criticalOverdensityClass     ), pointer :: criticalOverdensity_             => null()
     class           (cosmologicalMassVarianceClass), pointer :: cosmologicalMassVariance_        => null()
     class           (powerSpectrumClass           ), pointer :: powerSpectrum_                   => null()
     type            (virialDensityContrastFixed   ), pointer :: virialDensityContrastDefinition_ => null()
     type            (darkMatterProfileDMONFW      ), pointer :: darkMatterProfileDMODefinition_  => null()
     double precision                                         :: kappa                                     , scatter     , &
          &                                                      phi0                                      , phi1        , &
          &                                                      eta0                                      , eta1        , &
          &                                                      alpha                                     , beta        , &
          &                                                      timePrevious                              , massPrevious, &
          &                                                      concentrationMeanPrevious
   contains
     final     ::                                   diemerKravtsov2014Destructor
     procedure :: concentration                  => diemerKravtsov2014Concentration
     procedure :: concentrationMean              => diemerKravtsov2014ConcentrationMean
     procedure :: densityContrastDefinition      => diemerKravtsov2014DensityContrastDefinition
     procedure :: darkMatterProfileDMODefinition => diemerKravtsov2014DarkMatterProfileDefinition
  end type darkMatterProfileConcentrationDiemerKravtsov2014

  interface darkMatterProfileConcentrationDiemerKravtsov2014
     !!{
     Constructors for the \refClass{darkMatterProfileConcentrationDiemerKravtsov2014} dark matter halo profile concentration class.
     !!}
     module procedure diemerKravtsov2014ConstructorParameters
     module procedure diemerKravtsov2014ConstructorInternal
  end interface darkMatterProfileConcentrationDiemerKravtsov2014

contains

  function diemerKravtsov2014ConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily diemerKravtsov2014} dark matter halo
    profile concentration class.
    !!}
    implicit none
    type            (darkMatterProfileConcentrationDiemerKravtsov2014)                :: self
    type            (inputParameters                                 ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass                         ), pointer       :: cosmologyFunctions_
    class           (cosmologyParametersClass                        ), pointer       :: cosmologyParameters_
    class           (criticalOverdensityClass                        ), pointer       :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass                   ), pointer       :: cosmologicalMassVariance_
    class           (powerSpectrumClass                              ), pointer       :: powerSpectrum_
    double precision                                                                  :: kappa                    , phi0   , &
         &                                                                               phi1                     , eta0   , &
         &                                                                               eta1                     , alpha  , &
         &                                                                               beta                     , scatter

    ! Check and read parameters.
    !![
    <inputParameter>
      <name>kappa</name>
      <source>parameters</source>
      <variable>kappa</variable>
      <defaultValue>0.69d0</defaultValue>
      <description>The parameter $\kappa$ appearing in the halo concentration algorithm of \cite{diemer_universal_2014}.</description>
    </inputParameter>
    <inputParameter>
      <name>phi0</name>
      <source>parameters</source>
      <variable>phi0</variable>
      <defaultValue>6.58d0</defaultValue>
      <description>The parameter $\phi_0$ appearing in the halo concentration algorithm of \cite{diemer_universal_2014}.</description>
    </inputParameter>
    <inputParameter>
      <name>phi1</name>
      <source>parameters</source>
      <variable>phi1</variable>
      <defaultValue>1.37d0</defaultValue>
      <description>The parameter $\phi_1$ appearing in the halo concentration algorithm of \cite{diemer_universal_2014}.</description>
    </inputParameter>
    <inputParameter>
      <name>eta0</name>
      <source>parameters</source>
      <variable>eta0</variable>
      <defaultValue>6.82d0</defaultValue>
      <description>The parameter $\eta_0$ appearing in the halo concentration algorithm of \cite{diemer_universal_2014}.</description>
    </inputParameter>
    <inputParameter>
      <name>eta1</name>
      <source>parameters</source>
      <variable>eta1</variable>
      <defaultValue>1.42d0</defaultValue>
      <description>The parameter $\eta_1$ appearing in the halo concentration algorithm of \cite{diemer_universal_2014}.</description>
    </inputParameter>
    <inputParameter>
      <name>alpha</name>
      <source>parameters</source>
      <variable>alpha</variable>
      <defaultValue>1.12d0</defaultValue>
      <description>The parameter $\alpha$ appearing in the halo concentration algorithm of \cite{diemer_universal_2014}.</description>
    </inputParameter>
    <inputParameter>
      <name>beta</name>
      <source>parameters</source>
      <variable>beta</variable>
      <defaultValue>1.69d0</defaultValue>
      <description>The parameter $\beta$ appearing in the halo concentration algorithm of \cite{diemer_universal_2014}.</description>
    </inputParameter>
    <inputParameter>
      <name>scatter</name>
      <source>parameters</source>
      <variable>scatter</variable>
      <defaultValue>0.0d0</defaultValue>
      <description>The scatter (in dex) to assume in the halo concentration algorithm of \cite{diemer_universal_2014}.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    <objectBuilder class="cosmologyParameters"      name="cosmologyParameters_"      source="parameters"/>
    <objectBuilder class="criticalOverdensity"      name="criticalOverdensity_"      source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    <objectBuilder class="powerSpectrum"            name="powerSpectrum_"            source="parameters"/>
    !!]
     self=darkMatterProfileConcentrationDiemerKravtsov2014(kappa,phi0,phi1,eta0,eta1,alpha,beta,scatter,cosmologyFunctions_,cosmologyParameters_,criticalOverdensity_,cosmologicalMassVariance_,powerSpectrum_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"      />
    <objectDestructor name="cosmologyParameters_"     />
    <objectDestructor name="criticalOverdensity_"     />
    <objectDestructor name="cosmologicalMassVariance_"/>
    <objectDestructor name="powerSpectrum_"           />
    !!]
    return
  end function diemerKravtsov2014ConstructorParameters

  function diemerKravtsov2014ConstructorInternal(kappa,phi0,phi1,eta0,eta1,alpha,beta,scatter,cosmologyFunctions_,cosmologyParameters_,criticalOverdensity_,cosmologicalMassVariance_,powerSpectrum_) result(self)
    !!{
    Constructor for the \refClass{darkMatterProfileConcentrationDiemerKravtsov2014} dark matter halo profile
    concentration class.
    !!}
    use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleVirialDensityContrastDefinition
    use :: Virial_Density_Contrast, only : fixedDensityTypeCritical
    implicit none
    type            (darkMatterProfileConcentrationDiemerKravtsov2014  )                         :: self
    double precision                                                    , intent(in   )          :: kappa                          , phi0   , &
         &                                                                                          phi1                           , eta0   , &
         &                                                                                          eta1                           , alpha  , &
         &                                                                                          beta                           , scatter
    class           (cosmologyFunctionsClass                           ), intent(in   ), target  :: cosmologyFunctions_
    class           (cosmologyParametersClass                          ), intent(in   ), target  :: cosmologyParameters_
    class           (criticalOverdensityClass                          ), intent(in   ), target  :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass                     ), intent(in   ), target  :: cosmologicalMassVariance_
    class           (powerSpectrumClass                                ), intent(in   ), target  :: powerSpectrum_
    type            (darkMatterHaloScaleVirialDensityContrastDefinition)               , pointer :: darkMatterHaloScaleDefinition_
    !![
    <constructorAssign variables="kappa, phi0, phi1, eta0, eta1, alpha, beta, scatter, *cosmologyFunctions_, *cosmologyParameters_, *criticalOverdensity_, *cosmologicalMassVariance_, *powerSpectrum_"/>
    !!]

    self%timePrevious             =-1.0d0
    self%massPrevious             =-1.0d0
    self%concentrationMeanPrevious=-1.0d0
    allocate(     darkMatterHaloScaleDefinition_  )
    allocate(self%virialDensityContrastDefinition_)
    allocate(self%darkMatterProfileDMODefinition_ )
    !![
    <referenceConstruct owner="self" isResult="yes" object="virialDensityContrastDefinition_">
     <constructor>
      virialDensityContrastFixed                        (                                                                            &amp;
       &amp;                                             densityContrastValue                =200.0d0                              , &amp;
       &amp;                                             densityType                         =fixedDensityTypeCritical             , &amp;
       &amp;                                             turnAroundOverVirialRadius          =2.0d0                                , &amp;
       &amp;                                             cosmologyParameters_                =self%cosmologyParameters_            , &amp;
       &amp;                                             cosmologyFunctions_                 =self%cosmologyFunctions_               &amp;
       &amp;                                            )
     </constructor>
    </referenceConstruct>
    <referenceConstruct                             object="darkMatterHaloScaleDefinition_"  >
     <constructor>
      darkMatterHaloScaleVirialDensityContrastDefinition(                                                                            &amp;
       &amp;                                             cosmologyParameters_                =self%cosmologyParameters_            , &amp;
       &amp;                                             cosmologyFunctions_                 =self%cosmologyFunctions_             , &amp;
       &amp;                                             virialDensityContrast_              =self%virialDensityContrastDefinition_  &amp;
       &amp;                                            )
     </constructor>
    </referenceConstruct>
    <referenceConstruct owner="self" isResult="yes" object="darkMatterProfileDMODefinition_" >
     <constructor>
      darkMatterProfileDMONFW                           (                                                                            &amp;
       &amp;                                             velocityDispersionUseSeriesExpansion=.true.                               , &amp;
       &amp;                                             darkMatterHaloScale_                =darkMatterHaloScaleDefinition_         &amp;
       &amp;                                            )
     </constructor>
    </referenceConstruct>
    <objectDestructor                               name  ="darkMatterHaloScaleDefinition_" />
    !!]
    return
  end function diemerKravtsov2014ConstructorInternal

  subroutine diemerKravtsov2014Destructor(self)
    !!{
    Destructor for the \refClass{darkMatterProfileConcentrationDiemerKravtsov2014} dark matter halo profile concentration class.
    !!}
    implicit none
    type(darkMatterProfileConcentrationDiemerKravtsov2014), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"             />
    <objectDestructor name="self%cosmologyParameters_"            />
    <objectDestructor name="self%criticalOverdensity_"            />
    <objectDestructor name="self%cosmologicalMassVariance_"       />
    <objectDestructor name="self%powerSpectrum_"                  />
    <objectDestructor name="self%virialDensityContrastDefinition_"/>
    <objectDestructor name="self%darkMatterProfileDMODefinition_" />
    !!]
    return
  end subroutine diemerKravtsov2014Destructor

  double precision function diemerKravtsov2014Concentration(self,node)
    !!{
    Return the concentration of the dark matter halo profile of {\normalfont \ttfamily node}
    using the \cite{diemer_universal_2014} algorithm.
    !!}
    implicit none
    class(darkMatterProfileConcentrationDiemerKravtsov2014), intent(inout), target :: self
    type (treeNode                                        ), intent(inout), target :: node

    ! Get the mean concentration.
    diemerKravtsov2014Concentration=self%concentrationMean(node)
    ! Add scatter if necessary.
    if (self%scatter > 0.0d0)                                                     &
         &  diemerKravtsov2014Concentration                                       &
         & =diemerKravtsov2014Concentration                                       &
         & *10.0d0**(                                                             &
         &           +self%scatter                                                &
         &           *node%hostTree%randomNumberGenerator_%standardNormalSample() &
         &          )
    return
  end function diemerKravtsov2014Concentration

  double precision function diemerKravtsov2014ConcentrationMean(self,node)
    !!{
    Return the mean concentration of the dark matter halo profile of {\normalfont \ttfamily node}
    using the \cite{diemer_universal_2014} algorithm.
    !!}
    use :: Galacticus_Nodes        , only : nodeComponentBasic, treeNode
    use :: Math_Exponentiation     , only : cubeRoot
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (darkMatterProfileConcentrationDiemerKravtsov2014), intent(inout)          :: self
    type            (treeNode                                        ), intent(inout), target  :: node
    class           (nodeComponentBasic                              )               , pointer :: basic
    double precision                                                                           :: radiusHaloLagrangian, peakHeight        , &
         &                                                                                        wavenumber          , powerSpectrumSlope, &
         &                                                                                        concentrationMinimum, peakHeightMinimum

    basic => node%basic()
    if     (                                   &
         &   basic%mass() /= self%massPrevious &
         &  .or.                               &
         &   basic%time() /= self%timePrevious &
         & ) then
        radiusHaloLagrangian         =+cubeRoot(                                                                         &
            &                                  +3.0d0                                                                    &
            &                                  *basic%mass()                                                             &
            &                                  /4.0d0                                                                    &
            &                                  /Pi                                                                       &
            &                                  /self%cosmologyParameters_%densityCritical()                              &
            &                                  /self%cosmologyParameters_%OmegaMatter    ()                              &
            &                                 )
       peakHeight                    =+self%criticalOverdensity_     %value       (time=basic%time(),mass=basic%mass())  &
            &                         /self%cosmologicalMassVariance_%rootVariance(time=basic%time(),mass=basic%mass())
       wavenumber                    =+self%kappa                                                                        &
            &                         *2.0d0                                                                             &
            &                         *Pi                                                                                &
            &                         /radiusHaloLagrangian
       powerSpectrumSlope            =+self%powerSpectrum_%powerLogarithmicDerivative(wavenumber,basic%time())
       concentrationMinimum          =+max(                                                                              &
            &                              +self%phi0                                                                    &
            &                              +self%phi1                                                                    &
            &                              *powerSpectrumSlope                                                         , &
            &                              +2.0d0                                                                        &
            &                             )
       peakHeightMinimum             =+max(                                                                              &
            &                              +self%eta0                                                                    &
            &                              +self%eta1                                                                    &
            &                              *powerSpectrumSlope                                                         , &
            &                              +1.0d0                                                                        &
            &                             )
       self%concentrationMeanPrevious=+0.5d0                                                                             &
            &                         *concentrationMinimum                                                              &
            &                         *(                                                                                 &
            &                           +(peakHeight/peakHeightMinimum)**(-self%alpha)                                   &
            &                           +(peakHeight/peakHeightMinimum)**(+self%beta )                                   &
            &                          )
       self%massPrevious             = basic%mass()
       self%timePrevious             = basic%time()
    end if
    diemerKravtsov2014ConcentrationMean=self%concentrationMeanPrevious
    return
  end function diemerKravtsov2014ConcentrationMean

  function diemerKravtsov2014DensityContrastDefinition(self)
    !!{
    Return a virial density contrast object defining that used in the definition of concentration in the
    \cite{diemer_universal_2014} algorithm.
    !!}
    implicit none
    class(virialDensityContrastClass                      ), pointer       :: diemerKravtsov2014DensityContrastDefinition
    class(darkMatterProfileConcentrationDiemerKravtsov2014), intent(inout) :: self

    diemerKravtsov2014DensityContrastDefinition => self%virialDensityContrastDefinition_
    return
  end function diemerKravtsov2014DensityContrastDefinition

  function diemerKravtsov2014DarkMatterProfileDefinition(self)
    !!{
    Return a dark matter density profile object defining that used in the definition of concentration in the
    \cite{diemer_universal_2014} algorithm.
    !!}
    implicit none
    class(darkMatterProfileDMOClass                       ), pointer       :: diemerKravtsov2014DarkMatterProfileDefinition
    class(darkMatterProfileConcentrationDiemerKravtsov2014), intent(inout) :: self

    diemerKravtsov2014DarkMatterProfileDefinition => self%darkMatterProfileDMODefinition_
    return
  end function diemerKravtsov2014DarkMatterProfileDefinition
