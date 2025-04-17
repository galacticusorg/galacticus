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
Contains a module which implements beta functions.
!!}

! Add dependency on GSL library.
!; gsl

module Beta_Functions
  !!{
  Implements beta functions.
  !!}
  use, intrinsic :: ISO_C_Binding, only : c_double     , c_int
  use            :: Interface_GSL, only : gsl_sf_result
  implicit none
  private
  public :: Beta_Function, Beta_Function_Incomplete_Normalized

  interface
     function gsl_sf_beta_e(a,b,r) bind(c,name='gsl_sf_beta_e')
       !!{
       Template for the GSL beta C function.
       !!}
       import
       integer(c_int        )        :: gsl_sf_beta_e
       real   (c_double     ), value :: a            , b
       type   (gsl_sf_result)        :: r
     end function gsl_sf_beta_e
  end interface

  interface
     function gsl_sf_beta_inc_e(a,b,x,r) bind(c,name='gsl_sf_beta_inc_e')
       !!{
       Template for the GSL incomplete beta C function.
       !!}
       import
       integer(c_int        )        :: gsl_sf_beta_inc_e
       real   (c_double     ), value :: a                , b, &
            &                           x
       type   (gsl_sf_result)        :: r
     end function gsl_sf_beta_inc_e
  end interface

contains

  double precision function Beta_Function(a,b) result(r)
    !!{
    Evaluate the beta function, $B(a,b)$.
    !!}
    use :: Error             , only : Error_Report
    use :: Interface_GSL     , only : GSL_Success , gslErrorDecode
    use :: ISO_Varying_String, only : operator(//)
    implicit none
    double precision               , intent(in   ) :: a     , b
    integer         (c_int        )                :: status
    type            (gsl_sf_result)                :: r_

    status=GSL_SF_Beta_E(a,b,r_    )
    r     =                  r_%val
    if (status /= GSL_Success) call Error_Report('beta function evaluation failed: '//gslErrorDecode(status)//{introspection:location})
    return
  end function Beta_Function

  double precision function Beta_Function_Incomplete_Normalized(a,b,x) result(r)
    !!{
    Evaluate the normalized incomplete beta function, $B_x(a,b)/B(a,b)$.
    !!}
    use :: Error             , only : Error_Report, GSL_Error_Handler_Abort_Off, GSL_Error_Handler_Abort_On
    use :: Interface_GSL     , only : GSL_Success , gslErrorDecode
    use :: ISO_Varying_String, only : operator(//)
    implicit none
    double precision               , intent(in   ) :: a     , b, &
         &                                            x
    integer         (c_int        )                :: status
    type            (gsl_sf_result)                :: r_

    status=GSL_SF_Beta_Inc_E(a,b,x,r_    )
    r     =                        r_%val
    if (status /= GSL_Success) call Error_Report('incomplete beta function evaluation failed: '//gslErrorDecode(status)//{introspection:location})
    return
  end function Beta_Function_Incomplete_Normalized

end module Beta_Functions
