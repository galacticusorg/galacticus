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

  !+    Contributions to this file made by: Xiaolong Du
  
  !!{
  An implementation of dark matter halo profile concentrations using the
  \cite{diemer_accurate_2019} algorithm.
  !!}

  use :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass, criticalOverdensityClass
  use :: Cosmology_Functions       , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters      , only : cosmologyParametersClass
  use :: Dark_Matter_Profiles_DMO  , only : darkMatterProfileDMONFW
  use :: Virial_Density_Contrast   , only : virialDensityContrastFixed
  use :: Linear_Growth             , only : linearGrowthClass
  use :: Root_Finder               , only : rootFinder

  !![
  <darkMatterProfileConcentration name="darkMatterProfileConcentrationDiemerJoyce2019">
   <description>Dark matter halo concentrations are computed using the algorithm of \cite{diemer_accurate_2019}.</description>
   <deepCopy>
    <functionClass variables="virialDensityContrastDefinition_, darkMatterProfileDMODefinition_"/>
   </deepCopy>
   <stateStorable>
    <functionClass variables="virialDensityContrastDefinition_, darkMatterProfileDMODefinition_"/>
   </stateStorable>
  </darkMatterProfileConcentration>
  !!]
  type, extends(darkMatterProfileConcentrationClass) :: darkMatterProfileConcentrationDiemerJoyce2019
     !!{
     A dark matter halo profile concentration class implementing the algorithm of
     \cite{diemer_accurate_2019}.
     !!}
     private
     class           (cosmologyFunctionsClass      ), pointer     :: cosmologyFunctions_              => null()
     class           (cosmologyParametersClass     ), pointer     :: cosmologyParameters_             => null()
     class           (criticalOverdensityClass     ), pointer     :: criticalOverdensity_             => null()
     class           (cosmologicalMassVarianceClass), pointer     :: cosmologicalMassVariance_        => null()
     class           (linearGrowthClass            ), pointer     :: linearGrowth_                    => null()
     type            (virialDensityContrastFixed   ), pointer     :: virialDensityContrastDefinition_ => null()
     type            (darkMatterProfileDMONFW      ), pointer     :: darkMatterProfileDMODefinition_  => null()
     type            (rootFinder                   )              :: finder
     double precision                                             :: kappa                                     , scatter      , &
          &                                                          a0                                        , a1           , &
          &                                                          b0                                        , b1           , &
          &                                                          cAlpha                                    , timePrevious , &
          &                                                          concentrationMeanPrevious                 , massPrevious , &
          &                                                          GTildePrevious
     logical                                                      :: truncateConcentration                     , includeUpturn, &
          &                                                          truncateUpturn
   contains
     final     ::                                   diemerJoyce2019Destructor
     procedure :: concentration                  => diemerJoyce2019Concentration
     procedure :: concentrationMean              => diemerJoyce2019ConcentrationMean
     procedure :: densityContrastDefinition      => diemerJoyce2019DensityContrastDefinition
     procedure :: darkMatterProfileDMODefinition => diemerJoyce2019DarkMatterProfileDefinition
  end type darkMatterProfileConcentrationDiemerJoyce2019

  interface darkMatterProfileConcentrationDiemerJoyce2019
     !!{
     Constructors for the \refClass{darkMatterProfileConcentrationDiemerJoyce2019} dark matter halo profile concentration class.
     !!}
     module procedure diemerJoyce2019ConstructorParameters
     module procedure diemerJoyce2019ConstructorInternal
  end interface darkMatterProfileConcentrationDiemerJoyce2019

  ! Module global variables used in root finding.
  double precision :: GRoot_, powerSpectrumSlope
  !$omp threadprivate(GRoot_, powerSpectrumSlope)

contains

  function diemerJoyce2019ConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily diemerJoyce2019} dark matter halo
    profile concentration class.
    !!}
    implicit none
    type            (darkMatterProfileConcentrationDiemerJoyce2019   )                :: self
    type            (inputParameters                                 ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass                         ), pointer       :: cosmologyFunctions_
    class           (cosmologyParametersClass                        ), pointer       :: cosmologyParameters_
    class           (criticalOverdensityClass                        ), pointer       :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass                   ), pointer       :: cosmologicalMassVariance_
    class           (linearGrowthClass                               ), pointer       :: linearGrowth_
    double precision                                                                  :: kappa                    , a0           , &
         &                                                                               a1                       , b0           , &
         &                                                                               b1                       , cAlpha       , &
         &                                                                               scatter
    logical                                                                           :: truncateConcentration    , includeUpturn, &
         &                                                                               truncateUpturn

    ! Check and read parameters.
    !![
    <inputParameter>
      <name>kappa</name>
      <source>parameters</source>
      <defaultValue>0.41d0</defaultValue>
      <description>The parameter $\kappa$ appearing in the halo concentration algorithm of \cite{diemer_accurate_2019}.</description>
    </inputParameter>
    <inputParameter>
      <name>a0</name>
      <source>parameters</source>
      <defaultValue>2.45d0</defaultValue>
      <description>The parameter $a_0$ appearing in the halo concentration algorithm of \cite{diemer_accurate_2019}.</description>
    </inputParameter>
    <inputParameter>
      <name>a1</name>
      <source>parameters</source>
      <defaultValue>1.82d0</defaultValue>
      <description>The parameter $a_1$ appearing in the halo concentration algorithm of \cite{diemer_accurate_2019}.</description>
    </inputParameter>
    <inputParameter>
      <name>b0</name>
      <source>parameters</source>
      <defaultValue>3.20d0</defaultValue>
      <description>The parameter $b_0$ appearing in the halo concentration algorithm of \cite{diemer_accurate_2019}.</description>
    </inputParameter>
    <inputParameter>
      <name>b1</name>
      <source>parameters</source>
      <defaultValue>2.30d0</defaultValue>
      <description>The parameter $b_1$ appearing in the halo concentration algorithm of \cite{diemer_accurate_2019}.</description>
    </inputParameter>
    <inputParameter>
      <name>cAlpha</name>
      <source>parameters</source>
      <defaultValue>0.21d0</defaultValue>
      <description>The parameter $c_{\alpha}$ appearing in the halo concentration algorithm of \cite{diemer_accurate_2019}.</description>
    </inputParameter>
    <inputParameter>
      <name>scatter</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>The scatter (in dex) to assume in the halo concentration algorithm of \cite{diemer_accurate_2019}.</description>
    </inputParameter>
    <inputParameter>
      <name>truncateConcentration</name>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
      <description>If false, solutions to equation~(30) of \cite{diemer_accurate_2019} requiring $x&lt;1$ will cause a fatal error. If true, such cases are simply truncated to $x=1$.</description>
    </inputParameter>
    <inputParameter>
      <name>includeUpturn</name>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
      <description>If true, the term modeling the upturn in the $c(M)$ relation at high masses (i.e. $[1+\nu^2/B(n)]$ in equation~(28) of \citealt{diemer_accurate_2019}) is included. Otherwise this term is set equal to $1$ so that no upturn occurs.</description>
    </inputParameter>
    <inputParameter>
      <name>truncateUpturn</name>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
      <description>If true, the term modeling the upturn in the $c(M)$ relation at high masses (i.e. $[1+\nu^2/B(n)]$ in equation~(28) of \citealt{diemer_accurate_2019}) will be truncated at $\nu=\sqrt{B(n)}$ where the right-hand side of equation~(28) reaches the minimum, i.e. for any $\nu>\sqrt{B(n)}$, the right-hand side of equation~(28) will be truncated to $2 A/\sqrt{B}$.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    <objectBuilder class="cosmologyParameters"      name="cosmologyParameters_"      source="parameters"/>
    <objectBuilder class="criticalOverdensity"      name="criticalOverdensity_"      source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    <objectBuilder class="linearGrowth"             name="linearGrowth_"             source="parameters"/>
    !!]
    self=darkMatterProfileConcentrationDiemerJoyce2019(kappa,a0,a1,b0,b1,cAlpha,scatter,truncateConcentration,includeUpturn,truncateUpturn,cosmologyFunctions_,cosmologyParameters_,criticalOverdensity_,cosmologicalMassVariance_,linearGrowth_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"      />
    <objectDestructor name="cosmologyParameters_"     />
    <objectDestructor name="criticalOverdensity_"     />
    <objectDestructor name="cosmologicalMassVariance_"/>
    <objectDestructor name="linearGrowth_"            />
    !!]
    return
  end function diemerJoyce2019ConstructorParameters

  function diemerJoyce2019ConstructorInternal(kappa,a0,a1,b0,b1,cAlpha,scatter,truncateConcentration,includeUpturn,truncateUpturn,cosmologyFunctions_,cosmologyParameters_,criticalOverdensity_,cosmologicalMassVariance_,linearGrowth_) result(self)
    !!{
    Constructor for the \refClass{darkMatterProfileConcentrationDiemerJoyce2019} dark matter halo profile
    concentration class.
    !!}
    use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleVirialDensityContrastDefinition
    use :: Virial_Density_Contrast, only : fixedDensityTypeCritical
    use :: Root_Finder            , only : rangeExpandMultiplicative, rangeExpandSignExpectNegative, rangeExpandSignExpectPositive
    implicit none
    type            (darkMatterProfileConcentrationDiemerJoyce2019     )                         :: self
    double precision                                                    , intent(in   )          :: kappa                          , a0           , &
         &                                                                                          a1                             , b0           , &
         &                                                                                          b1                             , cAlpha       , &
         &                                                                                          scatter
    logical                                                             , intent(in   )          :: truncateConcentration          , includeUpturn, &
         &                                                                                          truncateUpturn
    class           (cosmologyFunctionsClass                           ), intent(in   ), target  :: cosmologyFunctions_
    class           (cosmologyParametersClass                          ), intent(in   ), target  :: cosmologyParameters_
    class           (criticalOverdensityClass                          ), intent(in   ), target  :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass                     ), intent(in   ), target  :: cosmologicalMassVariance_
    class           (linearGrowthClass                                 ), intent(in   ), target  :: linearGrowth_
    type            (darkMatterHaloScaleVirialDensityContrastDefinition)               , pointer :: darkMatterHaloScaleDefinition_
    !![
    <constructorAssign variables="kappa, a0, a1, b0, b1, cAlpha, scatter, truncateConcentration, includeUpturn, truncateUpturn, *cosmologyFunctions_, *cosmologyParameters_, *criticalOverdensity_, *cosmologicalMassVariance_, *linearGrowth_"/>
    !!]

    self%timePrevious             =-1.0d0
    self%massPrevious             =-1.0d0
    self%concentrationMeanPrevious=-1.0d0
    self%GTildePrevious           =-1.0d0
    self%finder                   =rootFinder(                                                             &
         &                                    rootFunction                 =GRoot                        , &
         &                                    rangeExpandDownward          =0.5d0                        , &
         &                                    rangeExpandUpward            =2.0d0                        , &
         &                                    rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive, &
         &                                    rangeExpandUpwardSignExpect  =rangeExpandSignExpectNegative, &
         &                                    rangeExpandType              =rangeExpandMultiplicative    , &
         &                                    rangeDownwardLimit           =1.0d0                        , &
         &                                    toleranceAbsolute            =0.0d0                        , &
         &                                    toleranceRelative            =1.0d-3                         &
         &                                   )
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
  end function diemerJoyce2019ConstructorInternal

  subroutine diemerJoyce2019Destructor(self)
    !!{
    Destructor for the \refClass{darkMatterProfileConcentrationDiemerJoyce2019} dark matter halo profile concentration class.
    !!}
    implicit none
    type(darkMatterProfileConcentrationDiemerJoyce2019), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"             />
    <objectDestructor name="self%cosmologyParameters_"            />
    <objectDestructor name="self%criticalOverdensity_"            />
    <objectDestructor name="self%cosmologicalMassVariance_"       />
    <objectDestructor name="self%linearGrowth_"                   />
    <objectDestructor name="self%virialDensityContrastDefinition_"/>
    <objectDestructor name="self%darkMatterProfileDMODefinition_" />
    !!]
    return
  end subroutine diemerJoyce2019Destructor

  double precision function diemerJoyce2019Concentration(self,node)
    !!{
    Return the concentration of the dark matter halo profile of {\normalfont \ttfamily node}
    using the \cite{diemer_accurate_2019} algorithm.
    !!}
    implicit none
    class(darkMatterProfileConcentrationDiemerJoyce2019), intent(inout), target :: self
    type (treeNode                                     ), intent(inout), target :: node

    ! Get the mean concentration.
    diemerJoyce2019Concentration=self%concentrationMean(node)
    ! Add scatter if necessary.
    if (self%scatter > 0.0d0)                                                     &
         &  diemerJoyce2019Concentration                                          &
         & =diemerJoyce2019Concentration                                          &
         & *10.0d0**(                                                             &
         &           +self%scatter                                                &
         &           *node%hostTree%randomNumberGenerator_%standardNormalSample() &
         &          )
    return
  end function diemerJoyce2019Concentration

  double precision function diemerJoyce2019ConcentrationMean(self,node)
    !!{
    Return the mean concentration of the dark matter halo profile of {\normalfont \ttfamily node}
    using the \cite{diemer_accurate_2019} algorithm.
    !!}
    use :: Galacticus_Nodes        , only : nodeComponentBasic, treeNode
    use :: Math_Exponentiation     , only : cubeRoot
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (darkMatterProfileConcentrationDiemerJoyce2019), intent(inout)          :: self
    type            (treeNode                                     ), intent(inout), target  :: node
    class           (nodeComponentBasic                           )               , pointer :: basic
    double precision                                                                        :: peakHeight    , massHalo  , &
         &                                                                                     alphaEffective, A         , &
         &                                                                                     B             , C         , &
         &                                                                                     GTilde        , termUpturn

    basic => node%basic()
    if     (                                   &
         &   basic%mass() /= self%massPrevious &
         &  .or.                               &
         &   basic%time() /= self%timePrevious &
         & ) then
       peakHeight        =+self%criticalOverdensity_     %value                               (time=basic%time(),mass=basic%mass(),node=node) &
            &             /self%cosmologicalMassVariance_%rootVariance                        (time=basic%time(),mass=basic%mass()          )
       massHalo          =+self%kappa**3                                                                                                      &
            &             *basic%mass()
       powerSpectrumSlope=-6.0d0                                                                                                              &
            &             *self%cosmologicalMassVariance_%rootVarianceLogarithmicGradient     (time=basic%time(),mass=massHalo              ) &
            &             -3.0d0
       alphaEffective    =+self%linearGrowth_            %logarithmicDerivativeExpansionFactor(time=basic%time())
       A                 =+self%a0                                                                                                            &
            &             *(                                                                                                                  &
            &              +1.0d0                                                                                                             &
            &              +self%a1*(powerSpectrumSlope+3.0d0)                                                                                &
            &              )
       B                 =+self%b0                                                                                                            &
            &             *(                                                                                                                  &
            &              +1.0d0                                                                                                             &
            &              +self%b1*(powerSpectrumSlope+3.0d0)                                                                                &
            &              )
       C                 =+1.0d0                                                                                                              &
            &             -self%cAlpha*(1.0d0-alphaEffective)
       if (self%includeUpturn) then
          termUpturn     =+1.0d0+peakHeight**2/B
       else
          termUpturn     =+1.0d0
       end if
       GRoot_            =+A                                                                                                                  &
            &             /peakHeight                                                                                                         &
            &             *termUpturn
       if     ( self%includeUpturn            &
            &  .and.                          &
            &   self%truncateUpturn           &
            &  .and.                          &
            &   peakHeight          > sqrt(B) &
            & ) then
          GRoot_         =+2.0d0*A/sqrt(B)
       end if
       ! Initial guess of the concentration.
       if (self%GTildePrevious < 0.0d0) self%GTildePrevious=5.0d0
       ! Solve for the concentration.
       if (self%truncateConcentration .and. GRoot(1.0d0) <= 0.0d0) then
          GTilde=1.0d0
       else
          GTilde=self%finder%find(rootGuess=self%GTildePrevious)
       end if
       self%GTildePrevious              =   GTilde
       self%concentrationMeanPrevious   = C*GTilde
       self%massPrevious                = basic%mass()
       self%timePrevious                = basic%time()
    end if
    diemerJoyce2019ConcentrationMean=self%concentrationMeanPrevious
    return
  end function diemerJoyce2019ConcentrationMean

  function diemerJoyce2019DensityContrastDefinition(self)
    !!{
    Return a virial density contrast object defining that used in the definition of concentration in the
    \cite{diemer_accurate_2019} algorithm.
    !!}
    implicit none
    class(virialDensityContrastClass                   ), pointer       :: diemerJoyce2019DensityContrastDefinition
    class(darkMatterProfileConcentrationDiemerJoyce2019), intent(inout) :: self

    diemerJoyce2019DensityContrastDefinition => self%virialDensityContrastDefinition_
    return
  end function diemerJoyce2019DensityContrastDefinition

  function diemerJoyce2019DarkMatterProfileDefinition(self)
    !!{
    Return a dark matter density profile object defining that used in the definition of concentration in the
    \cite{diemer_accurate_2019} algorithm.
    !!}
    implicit none
    class(darkMatterProfileDMOClass                    ), pointer       :: diemerJoyce2019DarkMatterProfileDefinition
    class(darkMatterProfileConcentrationDiemerJoyce2019), intent(inout) :: self

    diemerJoyce2019DarkMatterProfileDefinition => self%darkMatterProfileDMODefinition_
    return
  end function diemerJoyce2019DarkMatterProfileDefinition

  double precision function GRoot(concentration)
     !!{
     Root function used to find the concentration.
     !!}
     implicit none
     double precision, intent(in   ) :: concentration

     GRoot  =+GRoot_                                   &
          &  -concentration                            &
          &  /(                                        &
          &    +              log(1.0d0+concentration) &
          &    -concentration/   (1.0d0+concentration) &
          &   )**((5.0d0+powerSpectrumSlope)/6.0d0)
     return
  end function GRoot
