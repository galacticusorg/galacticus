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
  Implements a hybrid top-hat/sharp $k$-space power spectrum window function class.
  !!}
  
  use :: Cosmology_Parameters, only : cosmologyParametersClass
  
  !![
  <powerSpectrumWindowFunction name="powerSpectrumWindowFunctionTopHatSharpKHybrid">
   <description>
    A hybrid top-hat/sharp $k$-space window function for filtering of power spectra. This class implements a convolution of a
    top-hat window function and sharp $k$-space window function in $k$-space:
    \begin{equation}
     W(k) = W_\mathrm{th}(k) W_\mathrm{s}(k),
    \end{equation}
    where
    \begin{equation}
     W(k) = {3 (\sin(x)-x \cos(x)) \over x^3},
    \end{equation}
    where $x = k R_\mathrm{th}$, and
    \begin{equation}
     W_\mathrm{s}(k) = \left\{ \begin{array}{ll} 1 &amp; \hbox{if } k &lt; k_\mathrm{s} \\ 0 &amp; \hbox{if } k &gt; k_\mathrm{s}, \end{array} \right.
    \end{equation}
    where $k\mathrm{s} = \alpha / R_\mathrm{s}$ if {\normalfont \ttfamily
    [normalization]} is assigned a numerical value. Alternatively, if {\normalfont
    \ttfamily [normalization]}$=${\normalfont \ttfamily natural} then the value of
    $\alpha$ is chosen such that $k_\mathrm{s} = (6 \Pi^2 \bar{\rho}/M)^{1/3}$ if $R_\mathrm{s}=3M/4\pi\bar{\rho}$.
    The radii, $R_\mathrm{th}$ and $R_\mathrm{s}$, are chosen such that:
    \begin{eqnarray}
    R_\mathrm{th}^2 + R_\mathrm{s}^2 &amp;=&amp; (3M/4\pi\bar{\rho})^{2/3} \\
    R_\mathrm{s} &amp;=&amp; \beta R_\mathrm{th},
    \end{eqnarray}
    where $\beta=${\normalfont \ttfamily [radiiRatio]}.
   </description>
  </powerSpectrumWindowFunction>
  !!]
  type, extends(powerSpectrumWindowFunctionClass) :: powerSpectrumWindowFunctionTopHatSharpKHybrid
     !!{
     A hybrid top-hat/sharp $k$-space power spectrum window function class.
     !!}
     private
     class           (cosmologyParametersClass), pointer :: cosmologyParameters_ => null()
     double precision                                    :: cutOffNormalization           , radiiRatio
     type            (varying_string          )          :: normalization
   contains
     !![
     <methods>
       <method description="Set the radii of the components of the window function." method="radii" />
     </methods>
     !!]
     final     ::                      topHatSharpKHybridDestructor
     procedure :: value             => topHatSharpKHybridValue
     procedure :: wavenumberMaximum => topHatSharpKHybridWavenumberMaximum
     procedure :: radii             => topHatSharpKHybridRadii
  end type powerSpectrumWindowFunctionTopHatSharpKHybrid

  interface powerSpectrumWindowFunctionTopHatSharpKHybrid
     !!{
     Constructors for the \refClass{powerSpectrumWindowFunctionTopHatSharpKHybrid} power spectrum window function class.
     !!}
     module procedure topHatSharpKHybridConstructorParameters
     module procedure topHatSharpKHybridConstructorInternal
  end interface powerSpectrumWindowFunctionTopHatSharpKHybrid

contains

  function topHatSharpKHybridConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{powerSpectrumWindowFunctionTopHatSharpKHybrid} power spectrum window function class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (powerSpectrumWindowFunctionTopHatSharpKHybrid)                :: self
    type            (inputParameters                              ), intent(inout) :: parameters
    class           (cosmologyParametersClass                     ), pointer       :: cosmologyParameters_
    type            (varying_string                               )                :: normalization
    character       (len=32                                       )                :: normalizationChar
    double precision                                                               :: normalizationValue                     , radiiRatio

    ! Check parameters.
    !![
    <inputParameter>
      <name>normalization</name>
      <source>parameters</source>
      <variable>normalization</variable>
      <defaultValue>var_str('natural')</defaultValue>
      <description>
        The parameter $a$ in the relation $k_\mathrm{s} = a/r_\mathrm{s}$, where $k_\mathrm{s}$ is the cut-off wavenumber for
        the sharp $k$-space window function and $r_\mathrm{s}$ is the radius of a sphere (in real-space) enclosing the
        requested smoothing mass. Alternatively, a value of {\normalfont \ttfamily natural} will be supplied in which case the normalization
        is chosen such that, in real-space, $W(r=0)=1$. This results in a contained mass of $M=6 \pi^2 \bar{\rho} k_\mathrm{s}^{-3}$.
      </description>
    </inputParameter>
    <inputParameter>
      <name>radiiRatio</name>
      <source>parameters</source>
      <defaultValue>1.0d0</defaultValue>
      <description>
        The parameter $\beta$ in the relation $r_\mathrm{s}=\beta r_\mathrm{th}$ between $k$-space sharp and top-hat window
        function radii in the hybrid window function used for computing the variance in the power spectrum.
      </description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    !!]
    if (normalization == 'natural') then
       normalizationValue=0.0d0
    else
       normalizationChar=char(normalization)
       read (normalizationChar,*) normalizationValue
    end if
    self=powerSpectrumWindowFunctionTopHatSharpKHybrid(cosmologyParameters_,normalizationValue,radiiRatio)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    !!]
    return
  end function topHatSharpKHybridConstructorParameters

  function topHatSharpKHybridConstructorInternal(cosmologyParameters_,normalization,radiiRatio) result(self)
    !!{
    Internal constructor for the \refClass{powerSpectrumWindowFunctionTopHatSharpKHybrid} power spectrum window function class.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    type            (powerSpectrumWindowFunctionTopHatSharpKHybrid)                        :: self
    class           (cosmologyParametersClass                     ), target, intent(in   ) :: cosmologyParameters_
    double precision                                                       , intent(in   ) :: normalization                        , radiiRatio
    character       (len=18                                       )                        :: normalizationText
    !![
    <constructorAssign variables="radiiRatio, *cosmologyParameters_"/>
    !!]

    ! Compute normalization.
    if (normalization <= 0.0d0) then
       ! Compute the "natural" normalization.
       self%cutOffNormalization=+(                &
            &                     +9.0d0          &
            &                     *Pi             &
            &                     /2.0d0          &
            &                    )**(1.0d0/3.0d0)
       self%normalization='natural'
    else
       ! Use provided normalization.
       self%cutOffNormalization=+normalization
       write (normalizationText,'(e17.10)') normalization
       self%normalization=trim(normalizationText)
    end if
    return
  end function topHatSharpKHybridConstructorInternal

  subroutine topHatSharpKHybridDestructor(self)
    !!{
    Destructor for the \refClass{powerSpectrumWindowFunctionTopHatSharpKHybrid} power spectrum window function class.
    !!}
    implicit none
    type(powerSpectrumWindowFunctionTopHatSharpKHybrid), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    !!]
    return
  end subroutine topHatSharpKHybridDestructor

  double precision function topHatSharpKHybridValue(self,wavenumber,smoothingMass)
    !!{
    Computes a window function for calculations of the variance in the power spectrum. Specifically, uses a convolution of
    top-hat real-space and sharp $k$-space window functions. The top-hat radius is $r_\mathrm{th}$, while the $k$-space
    cut-off wavenumber is $k_\mathrm{s}=a/r_\mathrm{s}$, where $a=${\normalfont \ttfamily [normalization]}. The two radii are
    chosen such that $r_\mathrm{th}^2 + r_\mathrm{s}^2 = (3 M / 4 \pi \bar{rho})^{1/3}$ and $r_\mathrm{s}=\beta r_{\mathrm
    th}$ where $\beta=${\normalfont \ttfamily [pradiiRatio]}.
    !!}
    implicit none
    class           (powerSpectrumWindowFunctionTopHatSharpKHybrid), intent(inout) :: self
    double precision                                               , intent(in   ) :: smoothingMass              , wavenumber
    double precision                                               , parameter     :: xSeriesMaximum      =1.0d-3
    double precision                                                               :: radiusKSpaceSharp          , radiusTopHat             , &
         &                                                                            topHatWindowFunction       , kSpaceSharpWindowFunction, &
         &                                                                            wavenumberCutOff           , x                        , &
         &                                                                            xSquared

    ! Find radii for both filters and cut-off wavenumber for the sharp k-space filter.
    call self%radii(smoothingMass,radiusTopHat,radiusKSpaceSharp)
    wavenumberCutOff=self%wavenumberMaximum(smoothingMass)
    ! Compute the top-hat window function.
    x               =+wavenumber   &
         &           *radiusTopHat
    if      (x <= 0.0d0) then
       topHatWindowFunction=0.0d0
    else if (x <= xSeriesMaximum) then
       ! Use a series expansion of the window function for small x.
       xSquared            =+x**2
       topHatWindowFunction=+1.0d0                        &
            &               +xSquared*(  -1.0d0/   10.0d0 &
            &               +xSquared* ( +1.0d0/  280.0d0 &
            &               +xSquared*  (-1.0d0/15120.0d0 &
            &                           )                 &
            &                          )                  &
            &                         )
    else
       ! For larger x, use the full expression.
       topHatWindowFunction=3.0d0*(sin(x)-x*cos(x))/(x**3)
    end if
    ! Compute k-space sharp window function.
    if      (wavenumber <=            0.0d0) then
       kSpaceSharpWindowFunction=0.0d0
    else if (wavenumber <= wavenumberCutOff) then
       kSpaceSharpWindowFunction=1.0d0
    else
       kSpaceSharpWindowFunction=0.0d0
    end if
    ! Compute the convolution (which is just the multiplication in k-space).
    topHatSharpKHybridValue=+kSpaceSharpWindowFunction &
         &                  *     topHatWindowFunction
    return
  end function topHatSharpKHybridValue

  double precision function topHatSharpKHybridWavenumberMaximum(self,smoothingMass)
    !!{
    Computes the maximum wavenumber at which the window function for calculations of the variance in the power spectrum is
    non-zero. Specifically, uses a convolution of top-hat real-space and sharp $k$-space window functions. The top-hat radius
    is $r_\mathrm{th}$, while the $k$-space cut-off wavenumber is $k_\mathrm{s}=a/r_\mathrm{s}$, where $a=${\normalfont \ttfamily
    [normalization]}. The two radii are chosen such that $r_\mathrm{th}^2 + r_\mathrm{s}^2 = (3 M / 4 \pi \bar{rho})^{1/3}$
    and $r_\mathrm{s}=\beta r_\mathrm{th}$ where $\beta=${\normalfont \ttfamily [radiiRatio]}.
    !!}
    implicit none
    class           (powerSpectrumWindowFunctionTopHatSharpKHybrid), intent(inout) :: self
    double precision                                               , intent(in   ) :: smoothingMass
    double precision                                                               :: radiusTopHat , radiusKSpaceSharp

    call self%radii(smoothingMass,radiusTopHat,radiusKSpaceSharp)
    topHatSharpKHybridWavenumberMaximum=+self%cutOffNormalization &
         &                              /radiusKSpaceSharp
    return
  end function topHatSharpKHybridWavenumberMaximum

  subroutine topHatSharpKHybridRadii(self,smoothingMass,radiusTopHat,radiusKSpaceSharp)
    !!{
    Computes the radii of the top-hat and sharp $k$-space filters. Specifically, uses a convolution of top-hat real-space and
    sharp $k$-space window functions. The top-hat radius is $r_\mathrm{th}$, while the $k$-space cut-off wavenumber is
    $k_\mathrm{s}=a/r_\mathrm{s}$, where $a=${\normalfont \ttfamily [normalization]}. The two radii are chosen such that $r_\mathrm{th}^2 +
    r_\mathrm{s}^2 = (3 M / 4 \pi \bar{rho})^{1/3}$ and $r_\mathrm{s}=\beta r_\mathrm{th}$ where $\beta=${\normalfont \ttfamily [radiiRatio]}.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (powerSpectrumWindowFunctionTopHatSharpKHybrid), intent(inout) :: self
    double precision                                               , intent(in   ) :: smoothingMass
    double precision                                               , intent(  out) :: radiusTopHat , radiusKSpaceSharp
    double precision                                                               :: radiusTotal

    ! Find the radius enclosing this mass.
    radiusTotal      =+(                                             &
         &              +3.0d0                                       &
         &              /4.0d0                                       &
         &              /Pi                                          &
         &              *smoothingMass                               &
         &              /self%cosmologyParameters_%OmegaMatter    () &
         &              /self%cosmologyParameters_%densityCritical() &
         &              )**(1.0d0/3.0d0)
    ! Find the top-hat and sharp k-space radii, and the k-space wavenumber.
    radiusTopHat     =+radiusTotal                                   &
         &            /sqrt(                                         &
         &                  +1.0d0                                   &
         &                  +self%radiiRatio**2                      &
         &                 )
    radiusKSpaceSharp=+self%radiiRatio                               &
         &            *radiusTopHat
    return
  end subroutine topHatSharpKHybridRadii

