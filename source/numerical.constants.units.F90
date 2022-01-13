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
Contains a module of useful unit conversions.
!!}

module Numerical_Constants_Units
  !!{
  Contains various useful unit conversions.
  !!}
  use :: Numerical_Constants_Math, only : Pi
  implicit none
  public

  ! Ergs in Joules.
  double precision, parameter :: ergs              =1.0d-7

  ! Rydberg in Joules.
  !![
  <gslConstant variable="rydbergs" gslSymbol="GSL_CONST_MKSA_RYDBERG" gslHeader="gsl_const_mksa"/>
  !!]

  ! Angstroms in microns.
  double precision, parameter :: angstromsPerMicron=1.0d4

  ! Angstroms in meters.
  double precision, parameter :: angstromsPerMeter =1.0d10

  ! Electron volt (in units of Joules).
  !![
  <gslConstant variable="electronVolt" gslSymbol="GSL_CONST_MKSA_ELECTRON_VOLT" gslHeader="gsl_const_mksa"/>
  !!]

  ! Barn (cross section unit, in units of m²).
  double precision, parameter :: barn              =1.0d-28

  ! Degree (in units of radians).
  double precision, parameter :: degree            =Pi/180.0d0

  ! Arcsecond (in units of radians).
  double precision, parameter :: arcsecond         =Pi/180.0d0/60.0d0/60.0d0

end module Numerical_Constants_Units
