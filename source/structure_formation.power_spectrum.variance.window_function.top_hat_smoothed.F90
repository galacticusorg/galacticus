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

  !+ Contributions to this file made by: Ivan Esteban

  !!{
  Implements a top-hat power spectrum window function class, convoluted with a Gaussian.
  !!}

  use :: Cosmology_Parameters, only : cosmologyParametersClass
  
  !![
  <powerSpectrumWindowFunction name="powerSpectrumWindowFunctionTopHatSmoothed">
    <description>
      A top-hat in real space window function for filtering of power spectra, smoothed with a Gaussian. The window function is given by:
      \begin{equation}
      W(k) = {3 (\sin(x)-x \cos(x)) \over x^3} \times \exp{-k^2\sigma^2 \over 2},
      \end{equation}
      where $x = k R$ and $R=(3M/4\pi\bar{\rho})^{1/3}$ for a smoothing scale $M$ and mean matter density $\bar{\rho}$.
      $\sigma$ is the width of the smoothing Gaussian in real space. This exponentially cuts off the window function at $k \gg 1/\sigma$.
        </description>
  </powerSpectrumWindowFunction>
  !!]
  type, extends(powerSpectrumWindowFunctionClass) :: powerSpectrumWindowFunctionTopHatSmoothed
     !!{
     A top-hat power spectrum window function class, smoothed with a Gaussian.
     !!}
     private
     class(cosmologyParametersClass), pointer :: cosmologyParameters_ => null()
     double precision                         :: sigma
    contains
     final     ::                      topHatSmoothedDestructor
     procedure :: value             => topHatSmoothedValue
     procedure :: wavenumberMaximum => topHatSmoothedWavenumberMaximum
  end type powerSpectrumWindowFunctionTopHatSmoothed

  interface powerSpectrumWindowFunctionTopHatSmoothed
     !!{
     Constructors for the \refClass{powerSpectrumWindowFunctionTopHatSmoothed} power spectrum window function class.
     !!}
     module procedure topHatSmoothedConstructorParameters
     module procedure topHatSmoothedConstructorInternal
  end interface powerSpectrumWindowFunctionTopHatSmoothed

contains

  function topHatSmoothedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{powerSpectrumWindowFunctionTopHatSmoothed} power spectrum window function class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (powerSpectrumWindowFunctionTopHatSmoothed)                :: self
    type            (inputParameters                          ), intent(inout) :: parameters
    class           (cosmologyParametersClass                 ), pointer       :: cosmologyParameters_
    double precision                                                           :: sigma
    
    !![
    <inputParameter>
      <name>sigma</name>
      <source>parameters</source>
      <defaultValue>3.0d0</defaultValue>
      <defaultSource>Corresponds roughly to the smallest scale probed by Lyman-$\alpha$ data.</defaultSource>
      <description>The parameter ``$\sigma$'' which defines the width of the smoothing Gaussian.</description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    !!]
    self=powerSpectrumWindowFunctionTopHatSmoothed(cosmologyParameters_,sigma)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    !!]
    return
  end function topHatSmoothedConstructorParameters

  function topHatSmoothedConstructorInternal(cosmologyParameters_,sigma) result(self)
    !!{
    Internal constructor for the \refClass{powerSpectrumWindowFunctionTopHatSmoothed} power spectrum window function class.
    !!}
    implicit none
    type            (powerSpectrumWindowFunctionTopHatSmoothed)                        :: self
    class           (cosmologyParametersClass                 ), target, intent(in   ) :: cosmologyParameters_
    double precision                                                   , intent(in   ) :: sigma
    !![
    <constructorAssign variables="sigma, *cosmologyParameters_"/>
    !!]
    
    return
  end function topHatSmoothedConstructorInternal

  subroutine topHatSmoothedDestructor(self)
    !!{
    Destructor for the \refClass{powerSpectrumWindowFunctionTopHatSmoothed} power spectrum window function class.
    !!}
    implicit none
    type(powerSpectrumWindowFunctionTopHatSmoothed), intent(inout) :: self
    
    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    !!]
    return
  end subroutine topHatSmoothedDestructor

  double precision function topHatSmoothedValue(self,wavenumber,smoothingMass)
    !!{
    Top hat in real space window function Fourier transformed into $k$-space used in computing the variance of the power
    spectrum. Everything is convolved with a Gaussian of real-space width $\sigma$.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (powerSpectrumWindowFunctionTopHatSmoothed), intent(inout) :: self
    double precision                                           , intent(in   ) :: smoothingMass        , wavenumber
    double precision                                           , parameter     :: xSeriesMaximum=1.0d-3
    double precision                                                           :: topHatRadius         , x         , &
         &                                                                        xSquared

    topHatRadius=+(                                             &
         &         +3.0d0                                       &
         &         /4.0d0                                       &
         &         /Pi                                          &
         &         *smoothingMass                               &
         &         /self%cosmologyParameters_%OmegaMatter    () &
         &         /self%cosmologyParameters_%densityCritical() &
         &        )**(1.0d0/3.0d0)
    x           =+wavenumber                                    &
         &       *topHatRadius
    if      (x <= 0.0d0         ) then
       topHatSmoothedValue=+0.0d0
    else if (x <= xSeriesMaximum) then
       ! Use a series expansion of the window function for small x.
       xSquared           =+x**2
       topHatSmoothedValue=+1.0d0                        &
            &              +xSquared*(  -1.0d0/   10.0d0 &
            &              +xSquared* ( +1.0d0/  280.0d0 &
            &              +xSquared*  (-1.0d0/15120.0d0 &
            &                          )                 &
            &                         )                  &
            &                        )
    else
       ! For larger x, use the full expression.
       topHatSmoothedValue=+3.0d0         &
            &              *(             &
            &                +     sin(x) &
            &                -x   *cos(x) &
            &               )             &
            &              /  x**3
    end if
    topHatSmoothedValue=+topHatSmoothedValue &
         &              *exp(                &
         &                   -wavenumber**2  &
         &                   *self%sigma**2  &
         &                   /2.0d0          &
         &                  )
    return
  end function topHatSmoothedValue

  double precision function topHatSmoothedWavenumberMaximum(self,smoothingMass)
    !!{
    Maximum wavenumber for a top hat in real space window function convoluted with a Gaussian Fourier transformed into $k$-space
    used in computing the variance of the power spectrum. It is set to $k=3.5/\sigma$, with $\sigma$ being the real-space width of
    the smoothing Gaussian. Here, the Gaussian has dropped to around $2\times 10^{-3}$.    
    !!}
    implicit none
    class           (powerSpectrumWindowFunctionTopHatSmoothed), intent(inout) :: self
    double precision                                           , intent(in   ) :: smoothingMass
    !$GLC attributes unused :: self, smoothingMass

    topHatSmoothedWavenumberMaximum=+3.5d0      &
         &                          /self%sigma
    return
  end function topHatSmoothedWavenumberMaximum
