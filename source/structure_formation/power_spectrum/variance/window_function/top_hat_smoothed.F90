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

!+    Contributions to this file made by: Ivan Esteban, Claude.

  !!{RST
  Implements a top-hat power spectrum window function class, convolved with a Gaussian.
  !!}

  use :: Cosmology_Parameters, only : cosmologyParametersClass
  
  !![
  <powerSpectrumWindowFunction name="powerSpectrumWindowFunctionTopHatSmoothed" docformat="rst">
    <description>
    A top-hat in real space window function for filtering of power spectra, smoothed with a Gaussian. The window function is given by:

    .. math::

       W(k) = {3 (\sin(x)-x \cos(x)) \over x^3} \times \exp{-k^2\sigma^2 \over 2},

    where :math:`x = k R` and :math:`R=(3M/4\pi\bar{\rho})^{1/3}` for a smoothing scale :math:`M` and mean matter density :math:`\bar{\rho}`. :math:`\sigma` is the width of the smoothing Gaussian in real space. This exponentially cuts off the window function at :math:`k \gg 1/\sigma`.
    </description>
  </powerSpectrumWindowFunction>
  !!]
  type, extends(powerSpectrumWindowFunctionClass) :: powerSpectrumWindowFunctionTopHatSmoothed
     !!{RST
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
     !!{RST
     Constructors for the :galacticus-class:`powerSpectrumWindowFunctionTopHatSmoothed` power spectrum window function class.
     !!}
     module procedure topHatSmoothedConstructorParameters
     module procedure topHatSmoothedConstructorInternal
  end interface powerSpectrumWindowFunctionTopHatSmoothed

contains

  function topHatSmoothedConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`powerSpectrumWindowFunctionTopHatSmoothed` power spectrum window function class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (powerSpectrumWindowFunctionTopHatSmoothed)                :: self
    type            (inputParameters                          ), intent(inout) :: parameters
    class           (cosmologyParametersClass                 ), pointer       :: cosmologyParameters_
    double precision                                                           :: sigma
    
    !![
    <inputParameter docformat="rst">
      <name>sigma</name>
      <source>parameters</source>
      <defaultValue>3.0d0</defaultValue>
      <defaultSource>
      Corresponds roughly to the smallest scale probed by Lyman-:math:`\alpha` data.
      </defaultSource>
      <description>
      The parameter ":math:`\sigma`" which defines the width of the smoothing Gaussian.
      </description>
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
    !!{RST
    Internal constructor for the :galacticus-class:`powerSpectrumWindowFunctionTopHatSmoothed` power spectrum window function class.
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
    !!{RST
    Destructor for the :galacticus-class:`powerSpectrumWindowFunctionTopHatSmoothed` power spectrum window function class.
    !!}
    implicit none
    type(powerSpectrumWindowFunctionTopHatSmoothed), intent(inout) :: self
    
    !![
    <objectDestructor name="self%cosmologyParameters_"/>
    !!]
    return
  end subroutine topHatSmoothedDestructor

  double precision function topHatSmoothedValue(self,wavenumber,smoothingMass,time)
    !!{RST
    Top hat in real space window function Fourier transformed into :math:`k`-space used in computing the variance of the power spectrum. Everything is convolved with a Gaussian of real-space width :math:`\sigma`.
    !!}
    use :: Power_Spectrum_Window_Function_Utilities, only : Window_Function_Radius_Lagrangian, Window_Function_Top_Hat
    implicit none
    class           (powerSpectrumWindowFunctionTopHatSmoothed), intent(inout) :: self
    double precision                                           , intent(in   ) :: smoothingMass, wavenumber, &
         &                                                                        time
    double precision                                                           :: topHatRadius
    !$GLC attributes unused :: time

    topHatRadius       =Window_Function_Radius_Lagrangian(smoothingMass          ,self%cosmologyParameters_)
    topHatSmoothedValue=Window_Function_Top_Hat          (wavenumber*topHatRadius                          )
    topHatSmoothedValue=+topHatSmoothedValue &
         &              *exp(                &
         &                   -wavenumber**2  &
         &                   *self%sigma**2  &
         &                   /2.0d0          &
         &                  )
    return
  end function topHatSmoothedValue

  double precision function topHatSmoothedWavenumberMaximum(self,smoothingMass)
    !!{RST
    Maximum wavenumber for a top hat in real space window function convoluted with a Gaussian Fourier transformed into :math:`k`-space used in computing the variance of the power spectrum. It is set to :math:`k=3.5/\sigma`, with :math:`\sigma` being the real-space width of the smoothing Gaussian. Here, the Gaussian has dropped to around :math:`2\times 10^{-3}`.
    !!}
    implicit none
    class           (powerSpectrumWindowFunctionTopHatSmoothed), intent(inout) :: self
    double precision                                           , intent(in   ) :: smoothingMass
    !$GLC attributes unused :: self, smoothingMass

    topHatSmoothedWavenumberMaximum=+3.5d0      &
         &                          /self%sigma
    return
  end function topHatSmoothedWavenumberMaximum
