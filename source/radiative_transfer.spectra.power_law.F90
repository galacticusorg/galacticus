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

  !![
  <radiativeTransferSpectrum name="radiativeTransferSpectrumPowerLaw">
   <description>
    A photon spectrum class for power law spectra of the form
    \begin{equation}
     L_\lambda (\lambda) = \begin{array}{cc} A (\lambda/\lambda_\mathrm{min})^\alpha &amp; \hbox{if } \lambda_\mathrm{min} \le \lambda &lt; \lambda_\mathrm{max} \\ 0 &amp; \hbox{otherwise,} \end{array}
    \end{equation}
    where $A=${\normalfont \ttfamily normalization}, $\alpha=${\normalfont \ttfamily exponent}, $\lambda_\mathrm{min}=${\normalfont
    \ttfamily wavelengthMinimum}, and $\lambda_\mathrm{max}=${\normalfont \ttfamily wavelengthMaximum}.
   </description>
  </radiativeTransferSpectrum>
  !!]
  type, extends(radiativeTransferSpectrumClass) :: radiativeTransferSpectrumPowerLaw
     !!{
     Implementation of a power-law spectrum for radiative transfer calculations.
     !!}
     private
     double precision :: wavelengthMinimum, wavelengthMaximum, &
          &              exponent         , normalization
   contains
     procedure :: luminosity => powerLawLuminosity
     procedure :: spectrum   => powerLawSpectrum
  end type radiativeTransferSpectrumPowerLaw
  
  interface radiativeTransferSpectrumPowerLaw
     !!{
     Constructors for the \refClass{radiativeTransferSpectrumPowerLaw} radiative transfer spectrum class.
     !!}
     module procedure powerLawConstructorParameters
     module procedure powerLawConstructorInternal
  end interface radiativeTransferSpectrumPowerLaw
  
contains

  function powerLawConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{radiativeTransferSpectrumPowerLaw} radiative transfer spectrum class which takes a parameter set as
    input.
    !!}
    use :: Input_Parameters, only : inputParameters, inputParameter
    implicit none
    type            (radiativeTransferSpectrumPowerLaw)                :: self
    type            (inputParameters                  ), intent(inout) :: parameters
    double precision                                                   :: wavelengthMinimum, wavelengthMaximum, &
         &                                                                exponent         , normalization

    !![
    <inputParameter>
      <name>wavelengthMinimum</name>
      <description>The minimum wavelength (in units of \AA) for the power-law spectrum.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>wavelengthMaximum</name>
      <description>The maximum wavelength (in units of \AA) for the power-law spectrum.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>exponent</name>
      <description>The exponent of the power-law spectrum.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>normalization</name>
      <description>The normalization (in units of $L_\odot / \AA$) of the power-law spectrum.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=radiativeTransferSpectrumPowerLaw(wavelengthMinimum,wavelengthMaximum,exponent,normalization)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function powerLawConstructorParameters

  function powerLawConstructorInternal(wavelengthMinimum,wavelengthMaximum,exponent,normalization) result(self)
    !!{
    Internal constructor for the \refClass{radiativeTransferSpectrumPowerLaw} radiative transfer photon packet class.
    !!}
    implicit none
    type            (radiativeTransferSpectrumPowerLaw)                :: self
    double precision                                   , intent(in   ) :: wavelengthMinimum, wavelengthMaximum, &
         &                                                                exponent         , normalization
    !![
    <constructorAssign variables="wavelengthMinimum, wavelengthMaximum, exponent, normalization"/>
    !!]

    return
  end function powerLawConstructorInternal

  double precision function powerLawLuminosity(self,wavelengthMinimum,wavelengthMaximum)
    !!{
    Compute the luminosity in the given wavelength range for a power-law spectrum.
    !!}
    use :: Numerical_Comparison, only : Values_Agree
    implicit none
    class           (radiativeTransferSpectrumPowerLaw), intent(inout) :: self
    double precision                                   , intent(in   ) :: wavelengthMinimum, wavelengthMaximum

    if (Values_Agree(self%exponent,1.0d0,absTol=1.0d-3)) then
       powerLawLuminosity=+self%normalization     &
            &             *self%wavelengthMinimum &
            &             *log(                   &
            &                  +wavelengthMaximum &
            &                  /wavelengthMinimum &
            &             )
    else
       powerLawLuminosity=+self%normalization                                                  &
            &             *self%wavelengthMinimum                                              &
            &             /                                              (self%exponent-1.0d0) &
            &             *(                                                                   &
            &               +(wavelengthMaximum/self%wavelengthMinimum)**(self%exponent-1.0d0) &
            &               -(wavelengthMinimum/self%wavelengthMinimum)**(self%exponent-1.0d0) &
            &             )
    end if
    return    
  end function powerLawLuminosity
  
  double precision function powerLawSpectrum(self,wavelength)
    !!{
    Return the spectrum of the power-law.
    !!}
    implicit none
    class           (radiativeTransferSpectrumPowerLaw), intent(inout) :: self
    double precision                                   , intent(in   ) :: wavelength

    if     (                                      &
         &   wavelength >= self%wavelengthMinimum &
         &  .and.                                 &
         &   wavelength <  self%wavelengthMaximum &
         & ) then
       powerLawSpectrum=+   self%normalization     &
            &           *(                         &
            &             +      wavelength        &
            &             / self%wavelengthMinimum &
            &            )**self%exponent
    else
       powerLawSpectrum=+0.0d0
    end if
    return
  end function powerLawSpectrum
