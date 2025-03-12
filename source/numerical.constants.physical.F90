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
Contains a module of useful physical constants.
!!}

module Numerical_Constants_Physical
  !!{
  Contains various useful physical constants.
  !!}
  use :: Numerical_Constants_Math, only : Pi
  implicit none
  public

  !![
  <constant variable="speedLight" gslSymbol="GSL_CONST_MKSA_SPEED_OF_LIGHT" gslHeader="gsl_const_mksa" symbol="\mathrm{c}" units="m/s" unitsInSI="1.0" description="The speed of light in vacuum." externalDescription="https://en.wikipedia.org/wiki/Speed_of_light" reference="Gnu Scientific Library" referenceURL="https://www.gnu.org/software/gsl/doc/html/const.html#c.GSL_CONST_MKSA_SPEED_OF_LIGHT" group="physical"/>
  <constant variable="gravitationalConstant" gslSymbol="GSL_CONST_MKSA_GRAVITATIONAL_CONSTANT" gslHeader="gsl_const_mksa" symbol="\mathrm{G}" units="N m$^2$ kg$^{-2}$" unitsInSI="1.0" description="The gravitational constant." reference="Gnu Scientific Library" referenceURL="https://www.gnu.org/software/gsl/doc/html/const.html#c.GSL_CONST_MKSA_GRAVITATIONAL_CONSTANT" group="physical"/>
  <constant variable="stefanBoltzmannConstant" gslSymbol="GSL_CONST_MKSA_STEFAN_BOLTZMANN_CONSTANT" gslHeader="gsl_const_mksa" symbol="\sigma" units="J/s/m$^2$/K$^4$" unitsInSI="1.0" description="The Stefan-Boltzmann constant." reference="Gnu Scientific Library" referenceURL="https://www.gnu.org/software/gsl/doc/html/const.html#c.GSL_CONST_MKSA_STEFAN_BOLTZMANN_CONSTANT" group="physical"/>
  <constant variable="boltzmannsConstant" gslSymbol="GSL_CONST_MKSA_BOLTZMANN" gslHeader="gsl_const_mksa" symbol="\mathrm{k}" units="J/K" unitsInSI="1.0" description="Boltzmann's constant." reference="Gnu Scientific Library" referenceURL="https://www.gnu.org/software/gsl/doc/html/const.html#c.GSL_CONST_MKSA_BOLTZMANN" group="physical"/>
  <constant variable="thomsonCrossSection" gslSymbol="GSL_CONST_MKSA_THOMSON_CROSS_SECTION" gslHeader="gsl_const_mksa" symbol="\sigma_\mathrm{T}" units="m/s" unitsInSI="1.0" description="The Thompson cross-section." reference="Gnu Scientific Library" referenceURL="https://www.gnu.org/software/gsl/doc/html/const.html#c.GSL_CONST_MKSA_THOMSON_CROSS_SECTION" group="physical"/>
  <constant variable="electronMass" gslSymbol="GSL_CONST_MKSA_MASS_ELECTRON" gslHeader="gsl_const_mksa" symbol="\mathrm{m}_\mathrm{e}" units="kg" unitsInSI="1.0" description="The mass of an electron." reference="Gnu Scientific Library" referenceURL="https://www.gnu.org/software/gsl/doc/html/const.html#c.GSL_CONST_MKSA_MASS_ELECTRON" group="physical"/>
  <constant variable="plancksConstant" gslSymbol="GSL_CONST_MKSA_PLANCKS_CONSTANT_H" gslHeader="gsl_const_mksa" symbol="\mathrm{h}" units="J s" unitsInSI="1.0" description="Planck's constant." reference="Gnu Scientific Library" referenceURL="https://www.gnu.org/software/gsl/doc/html/const.html#c.GSL_CONST_MKSA_PLANCKS_CONSTANT_H" group="physical"/>
  <constant variable="electronCharge" gslSymbol="GSL_CONST_MKSA_ELECTRON_CHARGE" gslHeader="gsl_const_mksa" symbol="\mathrm{e}" units="C" unitsInSI="1.0" description="The charge of the electron." reference="Gnu Scientific Library" referenceURL="https://www.gnu.org/software/gsl/doc/html/const.html#c.GSL_CONST_MKSA_ELECTRON_CHARGE" group="physical"/>
  <constant variable="fineStructure" gslSymbol="GSL_CONST_NUM_FINE_STRUCTURE" gslHeader="gsl_const_num" symbol="\alpha" units="dimensionless" unitsInSI="1.0" description="The electromagnetic fine structure constant." reference="Gnu Scientific Library" referenceURL="https://www.gnu.org/software/gsl/doc/html/const.html#c.GSL_CONST_NUM_FINE_STRUCTURE" group="physical"/>
  <constant variable="permittivityFreeSpace" gslSymbol="GSL_CONST_MKSA_VACUUM_PERMITTIVITY" gslHeader="gsl_const_mksa" symbol="\epsilon_0" units="m/s" unitsInSI="1.0" description="The permittivity of free space." reference="Gnu Scientific Library" referenceURL="https://www.gnu.org/software/gsl/doc/html/const.html#c.GSL_CONST_MKSA_VACUUM_PERMITTIVITY" group="physical"/>
  <constant variable="radiationConstant" value="4.0d0*stefanBoltzmannConstant/speedLight" symbol="a" units="J/m$^3$/K$^4$" unitsInSI="1.0" description="The radiation density constant." externalDescription="https://en.wikipedia.org/wiki/Stefan%E2%80%93Boltzmann_law#Energy_density" reference="Definition." group="physical"/>
  <constant variable="electronRadius" value="1.0d0/(4.0d0*Pi*permittivityFreeSpace)*electronCharge**2/(electronMass*speedLight**2)" symbol="r_\mathrm{e}" units="m" unitsInSI="1.0" description="The classical electron radius." externalDescription="https://en.wikipedia.org/wiki/Classical_electron_radius" reference="Definition." group="physical"/>
  !!]

end module Numerical_Constants_Physical
