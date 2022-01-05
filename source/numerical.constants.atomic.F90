!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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
  use :: Numerical_Constants_Units   , only : angstromsPerMeter, rydbergs
  implicit none
  public

  ! Atomic mass unit (in kg).
  !![
  <gslConstant variable="atomicMassUnit" gslSymbol="GSL_CONST_MKSA_UNIFIED_ATOMIC_MASS" gslHeader="gsl_const_mksa"/>
  !!]

  ! Atomic masses.
  double precision, parameter :: atomicMassHydrogen                =1.007825d0
  double precision, parameter :: atomicMassHelium                  =4.002602d0
  double precision, parameter :: atomicMassLithium7                =7.016004d0

  ! Mass of hydrogen and helium atom (in kg).
  double precision, parameter :: massHydrogenAtom                  =atomicMassHydrogen*atomicMassUnit
  double precision, parameter :: massHeliumAtom                    =atomicMassHelium  *atomicMassUnit

  ! Hydrogen Lyman series limit wavelength including correction for finite mass of the atom.
  double precision, parameter :: lymanSeriesLimitWavelengthHydrogen=+plancksConstant    &
       &                                                            *speedLight         &
       &                                                            /rydbergs           &
       &                                                            *angstromsPerMeter  &
       &                                                            /(                  &
       &                                                              +1.0d0            &
       &                                                              -electronMass     &
       &                                                              /massHydrogenAtom &
       &                                                             )
end module Numerical_Constants_Atomic
