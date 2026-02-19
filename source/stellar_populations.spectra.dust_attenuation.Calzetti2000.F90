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
  Implements calculations of attenuation of stellar spectra using the model of \cite{calzetti_dust_2000}.
  !!}

  !![
  <stellarSpectraDustAttenuation name="stellarSpectraDustAttenuationCalzetti2000">
   <description>Returns the dust attenuation of stellar spectra according to the model of \cite{calzetti_dust_2000}.</description>
  </stellarSpectraDustAttenuation>
  !!]
  type, extends(stellarSpectraDustAttenuationClass) :: stellarSpectraDustAttenuationCalzetti2000
     !!{
     A class implementing calculations of attenuation of stellar spectra using the model of \cite{calzetti_dust_2000}.
     !!}
     private
   contains
     procedure :: attenuation => calzetti2000Attenuation
  end type stellarSpectraDustAttenuationCalzetti2000

  interface stellarSpectraDustAttenuationCalzetti2000
     !!{
     Constructors for the \refClass{stellarSpectraDustAttenuationCalzetti2000} stellar spectra dust attenuation class.
     !!}
     module procedure calzetti2000ConstructorParameters
  end interface stellarSpectraDustAttenuationCalzetti2000

contains

  function calzetti2000ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{stellarSpectraDustAttenuationCalzetti2000} stellar spectra dust attenuation class which takes a parameter set
    as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(stellarSpectraDustAttenuationCalzetti2000)                :: self
    type(inputParameters                          ), intent(inout) :: parameters

    self=stellarSpectraDustAttenuationCalzetti2000()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function calzetti2000ConstructorParameters

  double precision function calzetti2000Attenuation(self,wavelength,age,vBandAttenuation)
    !!{
    Return attenuation of stellar spectra according to the model of \cite{calzetti_dust_2000}.
    !!}
    use :: Numerical_Constants_Units, only : micronsToAngstroms
    implicit none
    class           (stellarSpectraDustAttenuationCalzetti2000), intent(inout) :: self
    double precision                                           , intent(in   ) :: wavelength              , age  , &
         &                                                                        vBandAttenuation
    double precision                                           , parameter     :: Rv               =4.05d0           ! Eqn. (5) of Calzetti et al.
    double precision                                                           :: wavelengthMicrons       , kappa
    !$GLC attributes unused :: self, age

    ! Eqn. (4) of Calzetti et al.
    wavelengthMicrons=wavelength/micronsToAngstroms
    if      (wavelengthMicrons > 0.12d0 .and. wavelengthMicrons <= 0.63d0) then
       kappa=2.659d0*(-2.156d0+1.509d0/wavelengthMicrons-0.198d0/wavelengthMicrons**2+0.011d0/wavelengthMicrons**3)+Rv
    else if (wavelengthMicrons > 0.63d0 .and. wavelengthMicrons <= 2.20d0) then
       kappa=2.659d0*(-1.857d0+1.040d0/wavelengthMicrons                                                          )+Rv
    else
       kappa=0.0d0
    end if
    ! Compute attenuation.
    calzetti2000Attenuation=vBandAttenuation*kappa/Rv
    return
  end function calzetti2000Attenuation
