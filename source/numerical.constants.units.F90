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
Contains a module of useful unit conversions.
!!}

module Numerical_Constants_Units
  !!{
  Contains various useful unit conversions.
  !!}
  use :: Numerical_Constants_Math, only : Pi
  implicit none
  public

  !![
  <constant variable="metersToAngstroms" value="1.0d10" symbol="\text{\AA}/\mathrm{m}" description="Number of Angstroms per meter." units="\AA/m" unitsInSI="1.0" reference="Defined." group="units"/>
  <constant variable="micronsToAngstroms" value="1.0d4" symbol="\text{\AA}/\mu\mathrm{m}" description="Number of Angstroms per micron." units="\AA/$\mu$m" unitsInSI="1.0" reference="Defined." group="units"/>
  <constant variable="rydberg" gslSymbol="GSL_CONST_MKSA_RYDBERG" gslHeader="gsl_const_mksa" symbol="\mathrm{Ry}" description="Rydberg---the ionization energy of hydrogen in its ground state." units="J" unitsInSI="1.0" reference="Gnu Scientific Library" referenceURL="https://www.gnu.org/software/gsl/doc/html/const.html#c.GSL_CONST_MKSA_RYDBERG" group="units:atomic"/>
  <constant variable="electronVolt" gslSymbol="GSL_CONST_MKSA_ELECTRON_VOLT" gslHeader="gsl_const_mksa" symbol="\mathrm{eV}" description="Electron-volt---the energy of an electron accelerated through a 1 volt electric field." units="J" unitsInSI="1.0" reference="Gnu Scientific Library" referenceURL="https://www.gnu.org/software/gsl/doc/html/const.html#c.GSL_CONST_MKSA_ELECTRON_VOLT" group="units"/>
  <constant variable="ergs" value="1.0d-7" symbol="\mathrm{ergs}" description="Unit of energy in the CGS system." externalDescription="https://en.wikipedia.org/wiki/Erg" units="J" unitsInSI="1.0" reference="\cite{jackson_classical_1999}" group="units"/>
  <constant variable="barn" value="1.0d-28" symbol="\mathrm{b}" description="Unit of area used for nuclear cross-sections." externalDescription="https://en.wikipedia.org/wiki/Barn_(unit)" units="m$^2$" unitsInSI="1.0" reference="IUPAC Gold Book" referenceURL="https://goldbook.iupac.org/terms/view/B00598" group="units"/>
  <constant variable="degree" value="Pi/180.0d0" symbol="^\circ" description="Degree---unit of angle." externalDescription="https://en.wikipedia.org/wiki/Degree_(angle)" units="rad" unitsInSI="1.0" reference="Derived" group="units"/>
  <constant variable="arcsecond" value="Pi/180.0d0/60.0d0/60.0d0" symbol="^{\prime\prime}" description="Arcsecond---unit of angle." externalDescription="https://en.wikipedia.org/wiki/Minute_and_second_of_arc" units="rad" unitsInSI="1.0" reference="Derived" group="units"/>
  !!]

end module Numerical_Constants_Units
