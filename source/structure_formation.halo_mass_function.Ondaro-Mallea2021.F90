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
Implements the dark matter halo mass function class of \cite{ondaro-mallea_non-universality_2022} for non-universal
primordial power spectra and structure growth rates.
!!}
  use :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass
  use :: Linear_Growth             , only : linearGrowthClass
  use :: Root_Finder               , only : rootFinder
  
  !![
  <haloMassFunction name="haloMassFunctionOndaroMallea2021">
    <description>
      The dark matter halo mass function class of \cite{ondaro-mallea_non-universality_2022} for non-universal
      primordial power spectra and structure growth rates. The mass function is given by
      \begin{equation}
      n(M) = n^\prime(M) f_2(n_\mathrm{eff}) f_3(\alpha_\mathrm{eff}),
      \end{equation}
      where $n^\prime(M)$ is some other mass function,
      \begin{equation}
      f_2(n_\mathrm{eff})=n_0 n_\mathrm{eff}^2 + n_1 n_\mathrm{eff} + n_0,
      \end{equation}
      and
      \begin{equation}
      f_3(\alpha_\mathrm{eff})=a_0 \alpha_\mathrm{eff}^2 + a_1.
      \end{equation}
      Here
      \begin{equation}
      n_\mathrm{eff} = -3 -2 \frac{\mathrm{d} \log \sigma(R)}{\mathrm{d} \log R} = -3 -6 \frac{\mathrm{d} \log \sigma(M)}{\mathrm{d} \log M},
      \end{equation}
      where $M$ is halo mass, and $\sigma(M)$ is the fractional root-variance in the linear theory cosmological density field on that scale, and
      \begin{equation}
      \alpha_\mathrm{eff}(a) = \left. \frac{\mathrm{d} \log D}{\mathrm{d} \log a}\right|_{a=a_\mathrm{ev}},
      \end{equation}
      where $D(a)$ is the linear growth factor, $a$ is the expansion factor, and $D(a_\mathrm{ev})=\gamma D(a)$ with $\gamma=4/5$.
    </description>
  </haloMassFunction>
  !!]
  type, extends(haloMassFunctionClass) :: haloMassFunctionOndaroMallea2021
     !!{
     The halo mass function class of \cite{ondaro-mallea_non-universality_2022} for non-universal
     primordial power spectra and structure growth rates.
     !!}
     private
     double precision                               , dimension(0:2) :: coefficientsN
     double precision                               , dimension(0:1) :: coefficientsA
     class           (cosmologicalMassVarianceClass), pointer        :: cosmologicalMassVariance_ => null()
     class           (linearGrowthClass            ), pointer        :: linearGrowth_             => null()
     class           (haloMassFunctionClass        ), pointer        :: haloMassFunction_         => null()
     type            (rootFinder                   )                 :: finder
   contains
     final     ::                 ondaroMallea2021Destructor
     procedure :: differential => ondaroMallea2021Differential
  end type haloMassFunctionOndaroMallea2021

  interface haloMassFunctionOndaroMallea2021
     !!{
     Constructors for the \refClass{haloMassFunctionOndaroMallea2021} halo mass function class.
     !!}
     module procedure ondaroMallea2021ConstructorParameters
     module procedure ondaroMallea2021ConstructorInternal
  end interface haloMassFunctionOndaroMallea2021

  ! Sub-module-scope objects used in root-finding.
  class           (haloMassFunctionOndaroMallea2021), pointer :: self_
  double precision                                            :: linearGrowthTarget
  !$omp threadprivate(self_,linearGrowthTarget)
  
contains

  function ondaroMallea2021ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{haloMassFunctionOndaroMallea2021} halo mass function class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (haloMassFunctionOndaroMallea2021)                 :: self
    type            (inputParameters                 ), intent(inout)  :: parameters
    class           (haloMassFunctionClass           ), pointer        :: haloMassFunction_
    class           (cosmologicalMassVarianceClass   ), pointer        :: cosmologicalMassVariance_
    class           (cosmologyParametersClass        ), pointer        :: cosmologyParameters_
    class           (linearGrowthClass               ), pointer        :: linearGrowth_
    double precision                                  , dimension(0:2) :: coefficientsN
    double precision                                  , dimension(0:1) :: coefficientsA

    !![
    <objectBuilder class="cosmologyParameters"      name="cosmologyParameters_"      source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    <objectBuilder class="linearGrowth"             name="linearGrowth_"             source="parameters"/>
    <objectBuilder class="haloMassFunction"         name="haloMassFunction_"         source="parameters"/>
    <inputParameter>
      <name>coefficientsN</name>
      <source>parameters</source>
      <defaultValue>[-0.1178d0,-0.3389d0,0.3022d0]</defaultValue>
      <defaultSource>\cite[][Table~3, row 4]{ondaro-mallea_non-universality_2022}</defaultSource>
      <description>The coefficients, $n_{0\ldots2}$, appearing in equation~(7) of the \cite{ondaro-mallea_non-universality_2022} halo mass function model.</description>
    </inputParameter>
    <inputParameter>
      <name>coefficientsA</name>
      <source>parameters</source>
      <defaultValue>[-1.0785d0,2.9700d0]</defaultValue>
      <defaultSource>\cite[][Table~3, row 4]{ondaro-mallea_non-universality_2022}</defaultSource>
      <description>The coefficients, $a_{0\ldots1}$, appearing in equation~(8) of the \cite{ondaro-mallea_non-universality_2022} halo mass function model.</description>
    </inputParameter>
    !!]
    self=haloMassFunctionOndaroMallea2021(coefficientsN,coefficientsA,cosmologyParameters_,cosmologicalMassVariance_,linearGrowth_,haloMassFunction_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"     />
    <objectDestructor name="cosmologicalMassVariance_"/>
    <objectDestructor name="linearGrowth_"            />
    <objectDestructor name="haloMassFunction_"        />
    !!]
    return
  end function ondaroMallea2021ConstructorParameters

  function ondaroMallea2021ConstructorInternal(coefficientsN,coefficientsA,cosmologyParameters_,cosmologicalMassVariance_,linearGrowth_,haloMassFunction_) result(self)
    !!{
    Internal constructor for the \refClass{haloMassFunctionOndaroMallea2021} halo mass function class.
    !!}
    use :: Root_Finder, only : rangeExpandMultiplicative, rangeExpandSignExpectPositive, rangeExpandSignExpectNegative
    implicit none
    type            (haloMassFunctionOndaroMallea2021)                                :: self
    class           (cosmologyParametersClass        ), target        , intent(in   ) :: cosmologyParameters_
    class           (cosmologicalMassVarianceClass   ), target        , intent(in   ) :: cosmologicalMassVariance_
    class           (linearGrowthClass               ), target        , intent(in   ) :: linearGrowth_
    class           (haloMassFunctionClass           ), target        , intent(in   ) :: haloMassFunction_
    double precision                                  , dimension(0:2), intent(in   ) :: coefficientsN
    double precision                                  , dimension(0:1), intent(in   ) :: coefficientsA
    !![
    <constructorAssign variables="coefficientsN, coefficientsA, *cosmologyParameters_, *cosmologicalMassVariance_, *linearGrowth_, *haloMassFunction_"/>
    !!]

    self%finder=rootFinder(                                                             &
         &                 rootFunction                 =linearGrowthFactorRoot       , &
         &                 toleranceRelative            =1.000d-3                     , &
         &                 rangeExpandUpward            =1.001d+0                     , &
         &                 rangeExpandDownward          =0.500d+0                     , &
         &                 rangeExpandType              =rangeExpandMultiplicative    , &
         &                 rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
         &                 rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative  &
         &                )
    return
  end function ondaroMallea2021ConstructorInternal

  subroutine ondaroMallea2021Destructor(self)
    !!{
    Destructor for the \refClass{haloMassFunctionOndaroMallea2021} halo mass function class.
    !!}
    implicit none
    type(haloMassFunctionOndaroMallea2021), intent(inout) :: self
    
    !![
    <objectDestructor name="self%cosmologyParameters_"     />
    <objectDestructor name="self%cosmologicalMassVariance_"/>
    <objectDestructor name="self%haloMassFunction_"        />
    <objectDestructor name="self%linearGrowth_"            />
    !!]
    return
  end subroutine ondaroMallea2021Destructor

  double precision function ondaroMallea2021Differential(self,time,mass,node)
    !!{
    Return the differential halo mass function at the given time and mass.
    !!}
    implicit none
    class           (haloMassFunctionOndaroMallea2021), intent(inout), target   :: self
    double precision                                  , intent(in   )           :: time                       , mass
    type            (treeNode                        ), intent(inout), optional :: node
    double precision                                                            :: slopeEffectivePowerSpectrum, rateGrowthLinear, &
         &                                                                         factorPowerSpectrum        , factorGrowth    , &
         &                                                                         timeGrowth

    ondaroMallea2021Differential =  +self%haloMassFunction_        %differential                        (time,mass,node)
    self_                        =>  self
    linearGrowthTarget           =  +self%linearGrowth_            %value                               (time          )                                &
         &                          *4.0d0                                                                                                              &
         &                          /5.0d0
    timeGrowth                   =  +self%finder                   %find                                (rootGuess=time)
    rateGrowthLinear             =  +self%linearGrowth_            %logarithmicDerivativeExpansionFactor(timeGrowth    )
    slopeEffectivePowerSpectrum  =  -3.0d0                                                                                                              &
         &                          -6.0d0                                                                                                              &
         &                          *self%cosmologicalMassVariance_%rootVarianceLogarithmicGradient     (mass,time     )
    factorPowerSpectrum          =  +self                          %coefficientsN                       (0             )*slopeEffectivePowerSpectrum**2 &
    &                               +self                          %coefficientsN                       (1             )*slopeEffectivePowerSpectrum    &
    &                               +self                          %coefficientsN                       (2             )
    factorGrowth                 =  +self                          %coefficientsA                       (0             )*rateGrowthLinear               &
         &                          +self                          %coefficientsA                       (1             )
    ondaroMallea2021Differential =  +ondaroMallea2021Differential                                                                                       &
            &                       *factorPowerSpectrum                                                                                                &
            &                       *factorGrowth     
    return
  end function ondaroMallea2021Differential

  double precision function linearGrowthFactorRoot(time)
    !!{
    Root function used in finding the epoch at which to evaluate the growth factor for the \cite{ondaro-mallea_non-universality_2022} dark matter
    halo mass function.
    !!}
    implicit none
    double precision, intent(in   ) :: time

    linearGrowthFactorRoot=+self_%linearGrowth_     %value(time) &
         &                 -      linearGrowthTarget
    return
  end function linearGrowthFactorRoot
  
