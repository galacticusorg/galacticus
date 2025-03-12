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
  Implements a transfer function class for fuzzy dark matter using the fitting function of
  \cite{passaglia_accurate_2022}.
  !!}

  use :: Dark_Matter_Particles, only : darkMatterParticleClass
  use :: Cosmology_Functions  , only : cosmologyFunctionsClass

  !![
  <transferFunction name="transferFunctionFuzzyDMPassaglia2022">
   <description>A transfer function class for fuzzy dark matter using the fitting function of \cite{passaglia_accurate_2022}.</description>
  </transferFunction>
  !!]
  type, extends(transferFunctionClass) :: transferFunctionFuzzyDMPassaglia2022
     !!{
     A transfer function class for fuzzy dark matter using the fitting function of \cite{passaglia_accurate_2022}.
     !!}
     private
     double precision                                   :: m22                          , time, &
          &                                                A                            , B   , &
          &                                                wavenumberJeans
     logical                                            :: solvingForMode
     class           (darkMatterParticleClass), pointer :: darkMatterParticle_ => null()
     class           (cosmologyFunctionsClass), pointer :: cosmologyFunctions_ => null()
     class           (transferFunctionClass  ), pointer :: transferFunctionCDM => null()
   contains
     final     ::                          fuzzyDMPassaglia2022Destructor
     procedure :: value                 => fuzzyDMPassaglia2022Value
     procedure :: logarithmicDerivative => fuzzyDMPassaglia2022LogarithmicDerivative
     procedure :: halfModeMass          => fuzzyDMPassaglia2022HalfModeMass
     procedure :: quarterModeMass       => fuzzyDMPassaglia2022QuarterModeMass
     procedure :: fractionModeMass      => fuzzyDMPassaglia2022FractionModeMass
     procedure :: epochTime             => fuzzyDMPassaglia2022EpochTime    
  end type transferFunctionFuzzyDMPassaglia2022
   
  interface transferFunctionFuzzyDMPassaglia2022
     !!{
     Constructors for the {\normalfont \ttfamily fuzzyDMPassaglia2022} transfer function class.
     !!}
     module procedure fuzzyDMPassaglia2022ConstructorParameters
     module procedure fuzzyDMPassaglia2022ConstructorInternal
  end interface transferFunctionFuzzyDMPassaglia2022

  ! Submodule-scope variables used in root finding.
  class           (transferFunctionFuzzyDMPassaglia2022), pointer :: self_
  double precision                                                :: modeTarget
  !$omp threadprivate(self_,modeTarget)

contains

  function fuzzyDMPassaglia2022ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily fuzzyDMPassaglia2022} transfer function class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (transferFunctionFuzzyDMPassaglia2022)                :: self
    type (inputParameters                     ), intent(inout) :: parameters
    class(transferFunctionClass               ), pointer       :: transferFunctionCDM
    class(darkMatterParticleClass             ), pointer       :: darkMatterParticle_
    class(cosmologyFunctionsClass             ), pointer       :: cosmologyFunctions_
    class(cosmologyParametersClass            ), pointer       :: cosmologyParameters_
  
    ! Validate parameters.
    if (.not.parameters%isPresent('transferFunction')) call Error_Report("an explicit 'transferFunction' must be given"//{introspection:location})
    ! Read parameters.

    !![
    <objectBuilder class="transferFunction"    name="transferFunctionCDM"  source="parameters"/>
    <objectBuilder class="darkMatterParticle"  name="darkMatterParticle_"  source="parameters"/>
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    !!]
    self=transferFunctionFuzzyDMPassaglia2022(transferFunctionCDM,darkMatterParticle_,cosmologyFunctions_,cosmologyParameters_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="transferFunctionCDM" />
    <objectDestructor name="darkMatterParticle_" />
    <objectDestructor name="cosmologyFunctions_" />
    !!]
    return
  end function fuzzyDMPassaglia2022ConstructorParameters

  function fuzzyDMPassaglia2022ConstructorInternal(transferFunctionCDM,darkMatterParticle_,cosmologyFunctions_,cosmologyParameters_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily fuzzyDMPassaglia2022} transfer function class.
    !!}
    use :: Cosmology_Functions_Parameters, only : requestTypeExpansionFactor
    use :: Dark_Matter_Particles         , only : darkMatterParticleFuzzyDarkMatter
    use :: Error                         , only : Error_Report
    use :: Numerical_Constants_Prefixes  , only : kilo
    implicit none
    type (transferFunctionFuzzyDMPassaglia2022)                        :: self
    class(transferFunctionClass               ), target, intent(in   ) :: transferFunctionCDM
    class(darkMatterParticleClass             ), target, intent(in   ) :: darkMatterParticle_
    class(cosmologyFunctionsClass             ), target, intent(in   ) :: cosmologyFunctions_
    class(cosmologyParametersClass            ), target, intent(in   ) :: cosmologyParameters_
    !![
    <constructorAssign variables="*transferFunctionCDM, *darkMatterParticle_, *cosmologyFunctions_, *cosmologyParameters_"/>
    !!]

    select type (darkMatterParticle__ => self%darkMatterParticle_)
    class is (darkMatterParticleFuzzyDarkMatter)
       if (darkMatterParticle__%densityFraction() == 1.0d0) then
          self%m22=+darkMatterParticle__%mass() &
               &   *kilo                        &
               &   /1.0d-22
       else
          call Error_Report('transfer function is not implemented for a mixed CDM and fuzzy dark matter model'//{introspection:location})
       end if
    class default
       call Error_Report('transfer function expects a fuzzy dark matter particle'//{introspection:location})
    end select
    ! Transfer function is defined at the present day.
    self%time              =self%cosmologyFunctions_%cosmicTime(expansionFactor=1.0d0)
    ! Compute parameters of the fitting function.
    self%A                 =+2.22d0*     self%m22**(1.0d0/25.0d0-1.0d0/1000.0d0*log(self%m22))  ! Equation 46 of Passaglia & Hu (2022, PRD, 105, 123529; https://ui.adsabs.harvard.edu/abs/2022PhRvD.105l3529P).
    self%B                 =+0.16d0/     self%m22**(1.0d0/20.0d0                             )  ! Equation 46 of Passaglia & Hu (2022, PRD, 105, 123529; https://ui.adsabs.harvard.edu/abs/2022PhRvD.105l3529P).
    self%wavenumberJeans   =+9.00d0*sqrt(self%m22                                             ) ! Equation 45 of Passaglia & Hu (2022, PRD, 105, 123529; https://ui.adsabs.harvard.edu/abs/2022PhRvD.105l3529P).
    ! Initialize mode solving state.
    self%solvingForMode    =.false.
    return
  end function fuzzyDMPassaglia2022ConstructorInternal

  subroutine fuzzyDMPassaglia2022Destructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily fuzzyDMPassaglia2022} transfer function class.
    !!}
    implicit none
    type(transferFunctionFuzzyDMPassaglia2022), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    <objectDestructor name="self%darkMatterParticle_" />
    <objectDestructor name="self%transferFunctionCDM" />
    <objectDestructor name="self%cosmologyFunctions_" />
    !!]
    return
  end subroutine fuzzyDMPassaglia2022Destructor

  double precision function fuzzyDMPassaglia2022Value(self,wavenumber)
    !!{
    Return the transfer function at the given wavenumber.
    !!}
    implicit none
    class           (transferFunctionFuzzyDMPassaglia2022), intent(inout) :: self
    double precision                                      , intent(in   ) :: wavenumber
    double precision                                      , parameter     :: n         =2.5d0
    double precision                                                      :: x
    
    x                        =+self%A               &
         &                    *     wavenumber      &
         &                    /self%wavenumberJeans
    if (self%solvingForMode) then
       fuzzyDMPassaglia2022Value=+1.0d0
    else
       fuzzyDMPassaglia2022Value=+self%transferFunctionCDM%value(wavenumber)
    end if
    fuzzyDMPassaglia2022Value=+fuzzyDMPassaglia2022Value &
         &                    *sin(                      &
         &                            x**       n        &
         &                        )                      &
         &                    /       x**       n        &
         &                    /   (                      &
         &                         +1.0d0                &
         &                         +self%B               &
         &                         *     x**(6.0d0-n)    &
         &                        )
    return
  end function fuzzyDMPassaglia2022Value

  double precision function fuzzyDMPassaglia2022LogarithmicDerivative(self,wavenumber)
    !!{
    Return the logarithmic derivative of the transfer function at the given wavenumber.
    !!}
    implicit none
    class           (transferFunctionFuzzyDMPassaglia2022), intent(inout) :: self
    double precision                                      , intent(in   ) :: wavenumber
    double precision                                      , parameter     :: n         =2.5d0
    double precision                                                      :: x

    x                                        =+self%A               &
         &                                    *     wavenumber      &
         &                                    /self%wavenumberJeans
    fuzzyDMPassaglia2022LogarithmicDerivative=+self%transferFunctionCDM%logarithmicDerivative(wavenumber)&
    &+x**6*(-6.0d0+n)*self%B/(x**6*self%B+x**n)-n+x**n*n/tan(x**n)
    return
  end function fuzzyDMPassaglia2022LogarithmicDerivative

  double precision function fuzzyDMPassaglia2022HalfModeMass(self,status)
    !!{
    Compute the mass corresponding to the wavenumber at which the transfer function is suppressed by a factor of two relative
    to a \gls{cdm} transfer function.
    !!}
    implicit none
    class  (transferFunctionfuzzyDMPassaglia2022), intent(inout), target   :: self
    integer                                      , intent(  out), optional :: status

    fuzzyDMPassaglia2022HalfModeMass=self%fractionModeMass(0.50d0,status)
    return
  end function fuzzyDMPassaglia2022HalfModeMass

  double precision function fuzzyDMPassaglia2022QuarterModeMass(self,status)
    !!{
    Compute the mass corresponding to the wavenumber at which the transfer function is suppressed by a factor of four relative
    to a \gls{cdm} transfer function.
    !!}
    implicit none
    class  (transferFunctionfuzzyDMPassaglia2022), intent(inout), target   :: self
    integer                                      , intent(  out), optional :: status

    fuzzyDMPassaglia2022QuarterModeMass=self%fractionModeMass(0.25d0,status)
    return
  end function fuzzyDMPassaglia2022QuarterModeMass

  double precision function fuzzyDMPassaglia2022FractionModeMass(self,fraction,status)
    !!{
    Compute the mass corresponding to the wavenumber at which the transfer function is suppressed by a factor of two relative
    to a \gls{cdm} transfer function.
    !!}
    use :: Error                   , only : errorStatusSuccess
    use :: Numerical_Constants_Math, only : Pi
    use :: Root_Finder             , only : rangeExpandMultiplicative, rangeExpandSignExpectNegative, rangeExpandSignExpectPositive, rootFinder
    implicit none
    class           (transferFunctionfuzzyDMPassaglia2022), intent(inout), target   :: self
    double precision                                      , intent(in   )           :: fraction
    integer                                               , intent(  out), optional :: status
    double precision                                                                :: matterDensity, wavenumber
    type            (rootFinder                          )                          :: finder

    ! There is no analytic solution for the fraction-mode mass so we resort to numerical root finding. This is complicated by the
    ! fact that the transfer function oscillates. Our approach is to start at a wavenumber much smaller than the Jeans wavenumber,
    ! and slowly increase the wavenumber until the root is bracketed.
    finder                               =   rootFinder(                                                             &
         &                                              rootFunction                 =modeSolver                   , &
         &                                              toleranceRelative            =1.000d-3                     , &
         &                                              rangeExpandUpward            =1.001d+0                     , &
         &                                              rangeExpandDownward          =0.500d+0                     , &
         &                                              rangeExpandType              =rangeExpandMultiplicative    , &
         &                                              rangeExpandUpwardSignExpect  =rangeExpandSignExpectNegative, &
         &                                              rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive  &
         &                                             )
    self_                                =>  self 
    modeTarget                           =  +fraction
    self%solvingForMode                  =   .true.
    wavenumber                           =   finder%find(rootGuess=1.0d-2*self%wavenumberJeans)
    self%solvingForMode                  =   .false.
    matterDensity                        =  +self%cosmologyParameters_%OmegaMatter    () &
         &                                  *self%cosmologyParameters_%densityCritical()
    fuzzyDMPassaglia2022FractionModeMass =  +4.0d0         &
         &                                  *Pi            &
         &                                  /3.0d0         &
         &                                  *matterDensity &
         &                                  *(             &
         &                                    +Pi          &
         &                                    /wavenumber  &
         &                                   )**3
    if (present(status)) status=errorStatusSuccess
    return
  end function fuzzyDMPassaglia2022FractionModeMass
  
  double precision function modeSolver(wavenumber)
    !!{
    Function used in solving for half- and quarter-mode masses in the {\normalfont \ttfamily fuzzyDMPassaglia2022} transfer function.
    !!}
    implicit none
    double precision, intent(in   ) :: wavenumber
    
    modeSolver=+self_%value     (wavenumber) &
         &     -      modeTarget
    return
  end function modeSolver

  double precision function fuzzyDMPassaglia2022EpochTime(self)
    !!{
    Return the cosmic time at the epoch at which this transfer function is defined.
    !!}
    implicit none
    class(transferFunctionFuzzyDMPassaglia2022), intent(inout) :: self

    fuzzyDMPassaglia2022EpochTime=self%time
    return
  end function fuzzyDMPassaglia2022EpochTime
