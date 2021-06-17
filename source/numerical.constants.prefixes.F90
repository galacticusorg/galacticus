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
Contains a module of useful numerical prefixes.
!!}

module Numerical_Constants_Prefixes
  !!{
  Contains useful numerical prefixes.
  !!}
  implicit none
  public

  double precision, parameter :: hella=1.0d27
  double precision, parameter :: yotta=1.0d24
  double precision, parameter :: zetta=1.0d21
  double precision, parameter :: exa  =1.0d18
  double precision, parameter :: peta =1.0d15
  double precision, parameter :: tera =1.0d12
  double precision, parameter :: giga =1.0d9
  double precision, parameter :: mega =1.0d6
  double precision, parameter :: kilo =1.0d3
  double precision, parameter :: hecto=1.0d2
  double precision, parameter :: deca =1.0d1
  double precision, parameter :: deci =1.0d-1
  double precision, parameter :: centi=1.0d-2
  double precision, parameter :: milli=1.0d-3
  double precision, parameter :: micro=1.0d-6
  double precision, parameter :: nano =1.0d-9
  double precision, parameter :: pico =1.0d-12
  double precision, parameter :: femto=1.0d-15
  double precision, parameter :: atto =1.0d-18
  double precision, parameter :: zepto=1.0d-21
  double precision, parameter :: yocto=1.0d-24

end module Numerical_Constants_Prefixes
