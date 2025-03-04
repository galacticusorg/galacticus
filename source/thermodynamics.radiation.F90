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

module Thermodynamics_Radiation
  !!{
  Implements calculations of thermal radiation.
  !!}
  implicit none
  private
  public :: Blackbody_Emission, Blackbody_Radiance

  ! Define labels that specify if the computed radiance is per unit wavelength or per unit frequency.
  integer, parameter, public :: radianceTypeWavelength=0
  integer, parameter, public :: radianceTypeFrequency =1

contains

  double precision function Blackbody_Emission(wavelength,temperature,radianceType)
    !!{
    Compute the Planck blackbody spectral radiance (defined per unit wavelength, in units of J s$^{-1}$ m$^{-2}$ sr$^{-1}$
    \AA$^{-1}$) or J s$^{-1}$ m$^{-2}$ sr$^{-1}$ Hz$^{-1}$ depending on the optional {\normalfont \ttfamily radianceType} argument). Input {\normalfont \ttfamily
    wavelength} is in Angstroms, input temperature is in Kelvin.
    !!}
    use :: Numerical_Constants_Physical, only : boltzmannsConstant, plancksConstant, speedLight
    use :: Numerical_Constants_Units   , only : metersToAngstroms
    implicit none
    double precision, intent(in   )           :: temperature                       , wavelength
    integer         , intent(in   ), optional :: radianceType
    double precision, parameter               :: exponentialArgumentMaximum=100.0d0
    double precision                          :: exponentialArgument
    integer                                   :: radianceTypeActual

    ! Assume zero emission by default.
    Blackbody_Emission=0.0d0

    ! Return zero emission for non-positive temperatures.
    if (temperature <= 0.0d0) return

    ! Determine what type of radiance to use.
    if (present(radianceType)) then
       radianceTypeActual=radianceType
    else
       radianceTypeActual=radianceTypeWavelength
    end if

    ! Compute the argument appearing in the exponential term.
    exponentialArgument=plancksConstant*speedLight*metersToAngstroms/wavelength/boltzmannsConstant/temperature
    ! If the exponential argument is not too large, then compute the spectrum, otherwise return zero.
    if (exponentialArgument < exponentialArgumentMaximum) then
       select case (radianceTypeActual)
       case (radianceTypeWavelength)
          Blackbody_Emission=2.0d0*plancksConstant*speedLight**2*metersToAngstroms**4/wavelength**5/(exp(exponentialArgument)-1.0d0)
       case (radianceTypeFrequency )
          Blackbody_Emission=2.0d0*plancksConstant*speedLight   *metersToAngstroms**3/wavelength**3/(exp(exponentialArgument)-1.0d0)
       end select
    else
       Blackbody_Emission=0.0d0
    end if
    return
  end function Blackbody_Emission

  double precision function Blackbody_Radiance(temperature)
    !!{
    Compute the Planck blackbody spectral radiance integrated over all frequencies (in units of J s$^{-1}$ m$^{-2}$
    sr$^{-1}$). Input temperature is in Kelvin.
    !!}
    use :: Numerical_Constants_Math    , only : Pi
    use :: Numerical_Constants_Physical, only : stefanBoltzmannConstant
    implicit none
    double precision, intent(in   ) :: temperature

    Blackbody_Radiance=+stefanBoltzmannConstant    &
         &             *temperature            **4 &
         &             /Pi
    return
  end function Blackbody_Radiance

end module Thermodynamics_Radiation
