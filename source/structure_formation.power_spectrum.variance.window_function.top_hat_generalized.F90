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
  Implements a generalized top-hat power spectrum window function class \citep{brown_towards_2022}.
  !!}

  use :: Cosmology_Parameters, only : cosmologyParametersClass

  !![
  <powerSpectrumWindowFunction name="powerSpectrumWindowFunctionTopHatGeneralized">
   <description>
    A generalized top-hat in real space window function for filtering of power spectra. The window function is given by
    \citep{brown_towards_2022}:
    \begin{equation}
     W(k) = {3 (\sin(\mu_\mathrm{g} x)-\mu_\mathrm{g} x \cos(\mu_\mathrm{g} x)) \over \mu_\mathrm{g} x^3},
    \end{equation}
    where $x = k R$ and $R=(3M/4\pi\bar{\rho})^{1/3}$ for a smoothing scale $M$ and mean matter density $\bar{\rho}$, and
    $\mu_\mathrm{g}$ is a parameter.
   </description>
  </powerSpectrumWindowFunction>
  !!]
  type, extends(powerSpectrumWindowFunctionClass) :: powerSpectrumWindowFunctionTopHatGeneralized
     !!{
     A top-hat power spectrum window function class.
     !!}
     private
     class           (cosmologyParametersClass), pointer :: cosmologyParameters_ => null()
     double precision                                    :: mu
    contains
     final     ::                      topHatGeneralizedDestructor
     procedure :: value             => topHatGeneralizedValue
     procedure :: wavenumberMaximum => topHatGeneralizedWavenumberMaximum
  end type powerSpectrumWindowFunctionTopHatGeneralized

  interface powerSpectrumWindowFunctionTopHatGeneralized
     !!{
     Constructors for the \refClass{powerSpectrumWindowFunctionTopHatGeneralized} power spectrum window function class.
     !!}
     module procedure topHatGeneralizedConstructorParameters
     module procedure topHatGeneralizedConstructorInternal
  end interface powerSpectrumWindowFunctionTopHatGeneralized

contains

  function topHatGeneralizedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{powerSpectrumWindowFunctionTopHatGeneralized} power spectrum window function class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (powerSpectrumWindowFunctionTopHatGeneralized)                :: self
    type            (inputParameters                             ), intent(inout) :: parameters
    class           (cosmologyParametersClass                    ), pointer       :: cosmologyParameters_
    double precision                                                              :: mu
    
    !![
    <inputParameter>
      <name>mu</name>
      <source>parameters</source>
      <description>The parameter $\mu_\mathrm{g}$ appearing in the generalized top-hat window function of \cite{brown_towards_2022}.</description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    !!]
    self=powerSpectrumWindowFunctionTopHatGeneralized(mu,cosmologyParameters_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    !!]
    return
  end function topHatGeneralizedConstructorParameters

  function topHatGeneralizedConstructorInternal(mu,cosmologyParameters_) result(self)
    !!{
    Internal constructor for the \refClass{powerSpectrumWindowFunctionTopHatGeneralized} power spectrum window function class.
    !!}
    implicit none
    type            (powerSpectrumWindowFunctionTopHatGeneralized)                        :: self
    class           (cosmologyParametersClass                    ), target, intent(in   ) :: cosmologyParameters_
    double precision                                                      , intent(in   ) :: mu
    !![
    <constructorAssign variables="mu, *cosmologyParameters_"/>
    !!]

    return
  end function topHatGeneralizedConstructorInternal

  subroutine topHatGeneralizedDestructor(self)
    !!{
    Destructor for the \refClass{powerSpectrumWindowFunctionTopHatGeneralized} power spectrum window function class.
    !!}
    implicit none
    type(powerSpectrumWindowFunctionTopHatGeneralized), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    !!]
    return
  end subroutine topHatGeneralizedDestructor

  double precision function topHatGeneralizedValue(self,wavenumber,smoothingMass,time)
    !!{
    Generalized top hat in real space window function Fourier transformed into $k$-space used in computing the variance of the
    power spectrum.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (powerSpectrumWindowFunctionTopHatGeneralized), intent(inout) :: self
    double precision                                              , intent(in   ) :: smoothingMass        , wavenumber, &
         &                                                                           time
    double precision                                              , parameter     :: xSeriesMaximum=1.0d-3
    double precision                                                              :: topHatRadius         , x         , &
         &                                                                           xSquared
    !$GLC attributes unused :: time

    topHatRadius=+(                                             &
         &         +3.0d0                                       &
         &         /4.0d0                                       &
         &         /Pi                                          &
         &         *smoothingMass                               &
         &         /self%cosmologyParameters_%OmegaMatter    () &
         &         /self%cosmologyParameters_%densityCritical() &
         &        )**(1.0d0/3.0d0)
    x           =+wavenumber                                    &
         &       *topHatRadius                                  &
         &       *self%mu
    if      (x <= 0.0d0         ) then
       topHatGeneralizedValue=+0.0d0
    else if (x <= xSeriesMaximum) then
       ! Use a series expansion of the window function for small x.
       xSquared              =+x**2
       topHatGeneralizedValue=+1.0d0                        &
            &                 +xSquared*(  -1.0d0/   10.0d0 &
            &                 +xSquared* ( +1.0d0/  280.0d0 &
            &                 +xSquared*  (-1.0d0/15120.0d0 &
            &                             )                 &
            &                            )                  &
            &                           )
    else
       ! For larger x, use the full expression.
       topHatGeneralizedValue=3.0d0*(sin(x)-x*cos(x))/(x**3)
    end if
    return
  end function topHatGeneralizedValue

  double precision function topHatGeneralizedWavenumberMaximum(self,smoothingMass)
    !!{
    Maximum wavenumber for a top hat in real space window function Fourier transformed into $k$-space used in computing the
    variance of the power spectrum.
    !!}
    implicit none
    class           (powerSpectrumWindowFunctionTopHatGeneralized), intent(inout) :: self
    double precision                                              , intent(in   ) :: smoothingMass
    double precision                                              , parameter     :: wavenumberLarge=huge(1.0d0) ! Effective infinity.
    !$GLC attributes unused :: self, smoothingMass

    topHatGeneralizedWavenumberMaximum=wavenumberLarge
    return
  end function topHatGeneralizedWavenumberMaximum
