!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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

!% Contains a module of useful physical constants.

module Numerical_Constants_Physical
  !% Contains various useful physical constants.
  use :: Numerical_Constants_Math, only : Pi
  implicit none
  public

  ! Speed of light (m/s).
  !# <gslConstant variable="speedLight" gslSymbol="GSL_CONST_MKSA_SPEED_OF_LIGHT" gslHeader="gsl_const_mksa"/>

  ! Newton's gravitational constant (in SI units).
  !# <gslConstant variable="gravitationalConstant" gslSymbol="GSL_CONST_MKSA_GRAVITATIONAL_CONSTANT" gslHeader="gsl_const_mksa"/>

  ! Stefan-Boltzmann constant (in units of J/s/M^2/K^4).
  !# <gslConstant variable="stefanBoltzmannConstant" gslSymbol="GSL_CONST_MKSA_STEFAN_BOLTZMANN_CONSTANT" gslHeader="gsl_const_mksa"/>

  ! Radiation constant (in units of J/m^3/K^4).
  double precision, parameter :: radiationConstant              =4.0d0*stefanBoltzmannConstant/speedLight

  ! Boltzmann's constant (in units of J/K).
  !# <gslConstant variable="boltzmannsConstant" gslSymbol="GSL_CONST_MKSA_BOLTZMANN" gslHeader="gsl_const_mksa"/>

  ! Thomson cross section (in units of m^2).
  !# <gslConstant variable="thomsonCrossSection" gslSymbol="GSL_CONST_MKSA_THOMSON_CROSS_SECTION" gslHeader="gsl_const_mksa"/>

  ! Electron mass (in units of kg).
  !# <gslConstant variable="electronMass" gslSymbol="GSL_CONST_MKSA_MASS_ELECTRON" gslHeader="gsl_const_mksa"/>

  ! Planck's constant (in units of J s).
  !# <gslConstant variable="plancksConstant" gslSymbol="GSL_CONST_MKSA_PLANCKS_CONSTANT_H" gslHeader="gsl_const_mksa"/>

  ! Electron Charge (in units of C).
  !# <gslConstant variable="electronCharge" gslSymbol="GSL_CONST_MKSA_ELECTRON_CHARGE" gslHeader="gsl_const_mksa"/>

  ! Permitivity of free space (in SI units).
  !# <gslConstant variable="eps0" gslSymbol="GSL_CONST_MKSA_VACUUM_PERMITTIVITY" gslHeader="gsl_const_mksa"/>

  ! Fine structure constant (unitless)
  !# <gslConstant variable="fineStructure" gslSymbol="GSL_CONST_NUM_FINE_STRUCTURE" gslHeader="gsl_const_num"/>

  ! Classical electron radius (m)
  double precision, parameter :: electronRadius                 = 1.0d0 / (4.0d0 * Pi * eps0) * electronCharge**2 / (electronMass * speedLight**2)

end module Numerical_Constants_Physical
