!! Copyright 2009, 2010, Andrew Benson <abenson@caltech.edu>
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






!% Contains a module of useful physical constants.

module Numerical_Constants_Physical
  !% Contains various useful physical constants.
  use FGSL
  use Numerical_Constants_Prefixes
  public

  ! Speed of light (m/s).
  double precision, parameter :: speedLight=FGSL_CONST_MKSA_SPEED_OF_LIGHT
  
  ! Newton's gravitational constant (in Galacticus' M_Solar, Mpc, km/s unit system).
  double precision, parameter :: gravitationalConstantGalacticus=FGSL_CONST_MKSA_GRAVITATIONAL_CONSTANT&
       &*FGSL_CONST_MKSA_SOLAR_MASS/(kilo**2)/FGSL_CONST_MKSA_PARSEC/mega

  ! Newton's gravitational constant (in SI units).
  double precision, parameter :: gravitationalConstant=FGSL_CONST_MKSA_GRAVITATIONAL_CONSTANT

  ! Stefan-Boltzmann constant (in units of J/s/M^2/K^4).
  double precision, parameter :: stefanBoltzmannConstant=FGSL_CONST_MKSA_STEFAN_BOLTZMANN_CONSTANT

  ! Radiation constant (in units of J/m^3/K^4).
  double precision, parameter :: radiationConstant=4.0d0*FGSL_CONST_MKSA_STEFAN_BOLTZMANN_CONSTANT/FGSL_CONST_MKSA_SPEED_OF_LIGHT

  ! Boltzmann's constant (in units of J/K).
  double precision, parameter :: boltzmannsConstant=FGSL_CONST_MKSA_BOLTZMANN

  ! Thomson cross section (in units of m^2).
  double precision, parameter :: thomsonCrossSection=FGSL_CONST_MKSA_THOMSON_CROSS_SECTION

  ! Electron mass (in units of kg).
  double precision, parameter :: electronMass=FGSL_CONST_MKSA_MASS_ELECTRON

  ! Planck's constant (in units of J s).
  double precision, parameter :: plancksConstant=FGSL_CONST_MKSA_PLANCKS_CONSTANT_H

end module Numerical_Constants_Physical
