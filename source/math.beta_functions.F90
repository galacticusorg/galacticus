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
Contains a module which implements beta functions.
!!}

! Add dependency on GSL library.
!; gsl

module Beta_Functions
  !!{
  Implements beta functions.
  !!}
  use, intrinsic :: ISO_C_Binding, only : c_double
  implicit none
  private
  public :: Beta_Function, Beta_Function_Incomplete_Normalized

  interface
     function gsl_sf_beta(a,b) bind(c,name='gsl_sf_beta')
       !!{
       Template for the GSL beta C function.
       !!}
       import
       real(c_double)        :: gsl_sf_beta
       real(c_double), value :: a          , b
     end function gsl_sf_beta
  end interface

  interface
     function gsl_sf_beta_inc(a,b,x) bind(c,name='gsl_sf_beta_inc')
       !!{
       Template for the GSL incomplete beta C function.
       !!}
       import
       real(c_double)        :: gsl_sf_beta_inc
       real(c_double), value :: a              , b, &
            &                   x
     end function gsl_sf_beta_inc
  end interface

contains

  double precision function Beta_Function(a,b)
    !!{
    Evaluate the beta function, $B(a,b)$.
    !!}
    implicit none
    double precision, intent(in   ) :: a, b

    Beta_Function=GSL_SF_Beta(a,b)
    return
  end function Beta_Function

  double precision function Beta_Function_Incomplete_Normalized(a,b,x)
    !!{
    Evaluate the normalized incomplete beta function, $B_x(a,b)/B(a,b)$.
    !!}
    implicit none
    double precision, intent(in   ) :: a, b, &
         &                             x

    Beta_Function_Incomplete_Normalized=GSL_SF_Beta_Inc(a,b,x)
    return
  end function Beta_Function_Incomplete_Normalized

end module Beta_Functions
