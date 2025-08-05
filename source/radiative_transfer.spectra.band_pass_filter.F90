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
  <radiativeTransferSpectrum name="radiativeTransferSpectrumBandPassFilter">
   <description>
    A photon spectrum class which simply truncates some other photon spectrum class outside of some band.
   </description>
  </radiativeTransferSpectrum>
  !!]
  type, extends(radiativeTransferSpectrumClass) :: radiativeTransferSpectrumBandPassFilter
     !!{
     Implementation of a spectrum band pass filter for radiative transfer calculations.
     !!}
     private
     class           (radiativeTransferSpectrumClass), pointer :: radiativeTransferSpectrum_ => null()
     double precision                                          :: wavelengthMinimum                   , wavelengthMaximum
   contains
     final     ::               bandPassFilterDestructor
     procedure :: luminosity => bandPassFilterLuminosity
     procedure :: spectrum   => bandPassFilterSpectrum
  end type radiativeTransferSpectrumBandPassFilter
  
  interface radiativeTransferSpectrumBandPassFilter
     !!{
     Constructors for the \refClass{radiativeTransferSpectrumBandPassFilter} radiative transfer spectrum class.
     !!}
     module procedure bandPassFilterConstructorParameters
     module procedure bandPassFilterConstructorInternal
  end interface radiativeTransferSpectrumBandPassFilter
  
contains

  function bandPassFilterConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{radiativeTransferSpectrumBandPassFilter} radiative transfer spectrum class which takes a parameter set as
    input.
    !!}
    use :: Input_Parameters, only : inputParameters, inputParameter
    implicit none
    type            (radiativeTransferSpectrumBandPassFilter)                :: self
    type            (inputParameters                        ), intent(inout) :: parameters
    class           (radiativeTransferSpectrumClass         ), pointer       :: radiativeTransferSpectrum_
    double precision                                                         :: wavelengthMinimum         , wavelengthMaximum

    !![
    <inputParameter>
      <name>wavelengthMinimum</name>
      <defaultValue>0.0d0</defaultValue>
      <description>The minimum wavelength (in units of \AA) to pass the spectrum.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>wavelengthMaximum</name>
      <defaultValue>huge(0.0d0)</defaultValue>
      <description>The maximum wavelength (in units of \AA) to pass the spectrum.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="radiativeTransferSpectrum" name="radiativeTransferSpectrum_" source="parameters"/>
    !!]
    self=radiativeTransferSpectrumBandPassFilter(wavelengthMinimum,wavelengthMaximum,radiativeTransferSpectrum_)
    !![
    <inputParametersValidate source="parameters"/>
     <objectDestructor name="radiativeTransferSpectrum_"/>
     !!]
   return
  end function bandPassFilterConstructorParameters

  function bandPassFilterConstructorInternal(wavelengthMinimum,wavelengthMaximum,radiativeTransferSpectrum_) result(self)
    !!{
    Internal constructor for the \refClass{radiativeTransferSpectrumBandPassFilter} radiative transfer photon packet class.
    !!}
    implicit none
    type            (radiativeTransferSpectrumBandPassFilter)                        :: self
    double precision                                         , intent(in   )         :: wavelengthMinimum         , wavelengthMaximum
    class           (radiativeTransferSpectrumClass         ), intent(in   ), target :: radiativeTransferSpectrum_
    !![
    <constructorAssign variables="wavelengthMinimum, wavelengthMaximum, *radiativeTransferSpectrum_"/>
    !!]

    return
  end function bandPassFilterConstructorInternal

  subroutine bandPassFilterDestructor(self)
    !!{
    Destructor for the \refClass{radiativeTransferSpectrumBandPassFilter} radiative transfer photon packet class.
    !!}
    implicit none
    type(radiativeTransferSpectrumBandPassFilter), intent(inout) :: self

    !![
    <objectDestructor name="self%radiativeTransferSpectrum_"/>
    !!]
    return
  end subroutine bandPassFilterDestructor

  double precision function bandPassFilterLuminosity(self,wavelengthMinimum,wavelengthMaximum)
    !!{
    Compute the luminosity in the given wavelength range for a power-law spectrum.
    !!}
    use :: Numerical_Comparison, only : Values_Agree
    implicit none
    class           (radiativeTransferSpectrumBandPassFilter), intent(inout) :: self
    double precision                                         , intent(in   ) :: wavelengthMinimum , wavelengthMaximum
    double precision                                                         :: wavelengthMinimum_, wavelengthMaximum_

    wavelengthMinimum_=max(wavelengthMinimum,self%wavelengthMinimum)
    wavelengthMaximum_=min(wavelengthMaximum,self%wavelengthMaximum)
    if (wavelengthMaximum_ > wavelengthMinimum_) then
       bandPassFilterLuminosity=self%radiativeTransferSpectrum_%luminosity(wavelengthMinimum,wavelengthMaximum)
    else
       bandPassFilterLuminosity=0.0d0
    end if
    return    
  end function bandPassFilterLuminosity
  
  double precision function bandPassFilterSpectrum(self,wavelength)
    !!{
    Return the spectrum of the power-law.
    !!}
    implicit none
    class           (radiativeTransferSpectrumBandPassFilter), intent(inout) :: self
    double precision                                         , intent(in   ) :: wavelength

    if     (                                      &
         &   wavelength >= self%wavelengthMinimum &
         &  .and.                                 &
         &   wavelength <  self%wavelengthMaximum &
         & ) then
       bandPassFilterSpectrum=self%radiativeTransferSpectrum_%spectrum(wavelength)
    else
       bandPassFilterSpectrum=0.0d0
    end if
    return
  end function bandPassFilterSpectrum
