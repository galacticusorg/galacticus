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
Contains a module which implements thermodynamic properties of ideal gases, including the adiabatic sound speed and the
Jeans length for a self-gravitating gas cloud.
!!}

module Ideal_Gases_Thermodynamics
  !!{
  Implements thermodynamic properties of ideal gases. Provides the adiabatic sound speed (assuming $\gamma = 5/3$) and the
  Jeans length, defined as the sound crossing time times the free-fall time, for gas of a given temperature, density, and
  mean atomic mass.
  !!}
  implicit none
  private
  public :: Ideal_Gas_Sound_Speed, Ideal_Gas_Jeans_Length

contains

  double precision function Ideal_Gas_Jeans_Length(temperature,density)
    !!{
    Return the Jeans length (in Mpc) for gas of given temperature (in Kelvin) and density (in internal units). The Jeans
    length is computed as the ratio of the adiabatic sound speed to the square root of the gravitational constant times
    the density, $\lambda_\mathrm{J} = c_\mathrm{s} / \sqrt{G \rho}$, which is the geometric mean of the sound crossing
    time and the free-fall time scaled to a length. The sound speed is evaluated assuming a monatomic ideal gas
    ($\gamma = 5/3$) with primordial mean atomic mass.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    double precision, intent(in   ) :: density, temperature

    Ideal_gas_Jeans_Length=Ideal_Gas_Sound_Speed(temperature)/sqrt(gravitationalConstant_internal)/sqrt(density)
    return
  end function Ideal_Gas_Jeans_Length

  double precision function Ideal_Gas_Sound_Speed(temperature,meanAtomicMass)
    !!{
    Return the sound speed (in km/s) for an ideal gas of given \mono{temperature} and (optionally) \mono{meanAtomicMass}.
    !!}
    use :: Numerical_Constants_Astronomical, only : meanAtomicMassPrimordial
    use :: Numerical_Constants_Atomic      , only : atomicMassUnit
    use :: Numerical_Constants_Physical    , only : boltzmannsConstant
    use :: Numerical_Constants_Prefixes    , only : kilo
    implicit none
    double precision, intent(in   )           :: temperature
    double precision, intent(in   ), optional :: meanAtomicMass
    double precision                          :: meanAtomicMassActual

    ! Determine what mean atomic mass to use.
    if (present(meanAtomicMass)) then
       meanAtomicMassActual=meanAtomicMass
    else
       meanAtomicMassActual=meanAtomicMassPrimordial
    end if

    ! Compute the sound speed.
    Ideal_Gas_Sound_Speed=sqrt(5.0d0*boltzmannsConstant*temperature/3.0d0/meanAtomicMassActual/atomicMassUnit)/kilo
    return
  end function Ideal_Gas_Sound_Speed

end module Ideal_Gases_Thermodynamics
