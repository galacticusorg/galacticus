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

  !!{
  Implements calculations of zero dust attenuation of stellar spectra.
  !!}

  !![
  <stellarSpectraDustAttenuation name="stellarSpectraDustAttenuationZero">
   <description>Returns a zero dust attenuation.</description>
  </stellarSpectraDustAttenuation>
  !!]
  type, extends(stellarSpectraDustAttenuationClass) :: stellarSpectraDustAttenuationZero
     !!{
     A class implementing zero dust attenuation of stellar spectra.
     !!}
     private
   contains
     procedure :: attenuation => zeroAttenuation
  end type stellarSpectraDustAttenuationZero

  interface stellarSpectraDustAttenuationZero
     !!{
     Constructors for the \refClass{stellarSpectraDustAttenuationZero} stellar spectra dust attenuation class.
     !!}
     module procedure zeroConstructorParameters
  end interface stellarSpectraDustAttenuationZero

contains

  function zeroConstructorParameters(parameters) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily zero} stellar spectra dust attenuation class.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(stellarSpectraDustAttenuationZero)                :: self
    type(inputParameters                  ), intent(inout) :: parameters

    self=stellarSpectraDustAttenuationZero()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function zeroConstructorParameters

  double precision function zeroAttenuation(self,wavelength,age,vBandAttenuation)
    !!{
    Return a zero attenuation.
    !!}
    implicit none
    class           (stellarSpectraDustAttenuationZero), intent(inout) :: self
    double precision                                   , intent(in   ) :: wavelength      , age, &
         &                                                                vBandAttenuation
    !$GLC attributes unused :: self, wavelength, age, vBandAttenuation

    zeroAttenuation=0.0d0
    return
  end function zeroAttenuation
