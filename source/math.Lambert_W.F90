!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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
Contains a module which implements Lambert W functions.
!!}

! Add dependency on GSL library.
!; gsl

module Lambert_Ws
  !!{
  Implements Lambert W functions.
  !!}
  use, intrinsic :: ISO_C_Binding, only : c_double
  implicit none
  private
  public :: Lambert_W0, Lambert_Wm1

  interface
     function gsl_sf_lambert_W0(x) bind(c,name='gsl_sf_lambert_W0')
       !!{
       Template for the GSL Lambert W (primary branch) C function.
       !!}
       import
       real(c_double)        :: gsl_sf_lambert_W0
       real(c_double), value :: x
     end function gsl_sf_lambert_W0
  end interface

  interface
     function gsl_sf_lambert_Wm1(x) bind(c,name='gsl_sf_lambert_Wm1')
       !!{
       Template for the GSL Lambert W (secondary branch) C function.
       !!}
       import
       real(c_double)        :: gsl_sf_lambert_Wm1
       real(c_double), value :: x
     end function gsl_sf_lambert_Wm1
  end interface

contains

  double precision function Lambert_W0(x)
    !!{
    Evaluate the (primary branch of the) Lambert W function, $W(x)$.
    !!}
    implicit none
    double precision, intent(in   ) :: x

    Lambert_W0=GSL_SF_Lambert_W0(x)
    return
  end function Lambert_W0

  double precision function Lambert_Wm1(x)
    !!{
    Evaluate the (secondary branch of the) Lambert W function, $W(x)$.
    !!}
    implicit none
    double precision, intent(in   ) :: x

    Lambert_Wm1=GSL_SF_Lambert_Wm1(x)
    return
  end function Lambert_Wm1

end module Lambert_Ws
