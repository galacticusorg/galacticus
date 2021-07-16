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
Contains a module which implements dilogarithms.
!!}

! Add dependency on GSL library.
!; gsl

module Dilogarithms
  !!{
  Implements dilogarithms.
  !!}
  use, intrinsic :: ISO_C_Binding, only : c_double
  implicit none
  private
  public :: Dilogarithm

  interface
     function gsl_sf_dilog(x) bind(c,name='gsl_sf_dilog')
       !!{
       Template for the GSL dilogarithm function.
       !!}
       import c_double
       real(c_double)        :: gsl_sf_dilog
       real(c_double), value :: x
     end function gsl_sf_dilog
  end interface
  
contains

  double precision function Dilogarithm(x)
    !!{
    Evaluate the $\hbox{Si}(x)\equiv\int_0^x \d t \sin(t)/t$ sine integral.
    !!}
    implicit none
    double precision, intent(in   ) :: x

    Dilogarithm=GSL_SF_Dilog(x)
    return
  end function Dilogarithm

end module Dilogarithms
