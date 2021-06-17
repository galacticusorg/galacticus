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

!!{
Contains a module of useful mathematical constants.
!!}

module Numerical_Constants_Math
  !!{
  Contains various useful mathematical constants.
  !!}
  use :: Kind_Numbers , only : kind_quad
  implicit none

  ! e.
  !![
  <gslConstant variable="e" gslSymbol="M_E" gslHeader="gsl_math"/>
  !!]

  ! Pi.
  !![
  <gslConstant variable="Pi" gslSymbol="M_PI" gslHeader="gsl_math"/>
  !!]
  real(kind=kind_quad), public, parameter :: PiQuadPrecision=3.141592653589793238462643383279502884197_kind_quad

  ! ! Natural logarithm of 10.
  !![
  <gslConstant variable="ln10" gslSymbol="M_LN10" gslHeader="gsl_math"/>
  !!]

  ! Natural logarithm of 2.
  !![
  <gslConstant variable="ln2" gslSymbol="M_LN2" gslHeader="gsl_math"/>
  !!]

  ! Euler's constant.
  !![
  <gslConstant variable="eulersConstant" gslSymbol="M_EULER" gslHeader="gsl_math"/>
  !!]

  ! Riemann zeta-function values.
  !! ζ(3) - https://oeis.org/A002117
  double precision, public, parameter :: riemannZeta3=1.20205690315959428539973816151144999076d0
  
end module Numerical_Constants_Math
