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
Contains a module of useful atomic constants.
!!}

module Numerical_Constants_Atomic
  !!{
  Contains various useful atomic constants.
  !!}
  use :: Numerical_Constants_Physical, only : electronMass     , plancksConstant, speedLight
  use :: Numerical_Constants_Units   , only : metersToAngstroms, rydberg
  implicit none
  public

  ! Atomic masses.
  !![
  <constant variable="atomicMassUnit" gslSymbol="GSL_CONST_MKSA_UNIFIED_ATOMIC_MASS" gslHeader="gsl_const_mksa" symbol="\mathrm{u}" units="kg" unitsInSI="1.0" description="The unified atomic mass unit."  reference="Gnu Scientific Library" referenceURL="https://www.gnu.org/software/gsl/doc/html/const.html#c.GSL_CONST_MKSA_UNIFIED_ATOMIC_MASS" group="units:atomic"/>
  <constant variable="atomicMassHydrogen" value="1.0078250322d0" symbol="A_{^1\mathrm{H}}" units="amu" unitsInSI="1.66053906892e-27" description="Atomic mass of the $^1$H isotope of hydrogen." reference="Commission on Isotopic Abundances and Atomic Weights" referenceURL="https://www.ciaaw.org/hydrogen.htm" group="atomic"/>
  <constant variable="atomicMassHelium" value="4.0026032545d0" symbol="A_{^4\mathrm{He}}" units="amu" unitsInSI="1.66053906892e-27" description="Atomic mass of the $^4$He isotope of helium." reference="Commission on Isotopic Abundances and Atomic Weights" referenceURL="https://www.ciaaw.org/helium.htm" group="atomic"/>
  <constant variable="atomicMassLithium7" value="7.01600344d0" symbol="A_{^7\mathrm{Li}}" units="amu" unitsInSI="1.66053906892e-27" description="Atomic mass of the $^7$Li isotope of lithium." reference="Commission on Isotopic Abundances and Atomic Weights" referenceURL="https://www.ciaaw.org/lithium.htm" group="atomic"/>
  <constant variable="massHydrogenAtom" value="atomicMassHydrogen*atomicMassUnit" symbol="m_\mathrm{H}" units="kg" unitsInSI="1.0" description="Mass of the $^1$H isotope of hydrogen." reference="Derived" group="atomic"/>
  <constant variable="massHeliumAtom" value="atomicMassHelium*atomicMassUnit" symbol="m_\mathrm{He}" units="kg" unitsInSI="1.0" description="Mass of the $^4$He isotope of helium." reference="Derived" group="atomic"/>
  !!]

  ! Hydrogen Lyman series limit wavelength including correction for finite mass of the atom.
  !![
  <constant variable="lymanSeriesLimitWavelengthHydrogen_atomic" value="+plancksConstant*speedLight/rydberg*metersToAngstroms/(+1.0d0-electronMass/massHydrogenAtom)" symbol="R_\mathrm{H}" units="$\AA$" unitsInSI="1.0e-10" description="Hydrogen Lyman series limit wavelength including correction for finite mass of the atom." externalDescription="https://en.wikipedia.org/wiki/Rydberg_constant#Rydberg_constant" reference="Derived." group="atomic"/>
  !!]
  
end module Numerical_Constants_Atomic
