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
Contains a module which implements binomial coefficients.
!!}

! Add dependency on GSL library.
!; gsl

module Binomial_Coefficients
  !!{
  Implements binomial coefficients.
  !!}
  use, intrinsic :: ISO_C_Binding, only : c_double, c_int
  implicit none
  private
  public :: Binomial_Coefficient

  interface
     function gsl_sf_choose(n,m) bind(c,name='gsl_sf_choose')
       !!{
       Template for the GSL choose C function.
       !!}
       import
       real   (c_double)        :: gsl_sf_choose
       integer(c_int   ), value :: n            , m
     end function gsl_sf_choose
  end interface

contains

  double precision function Binomial_Coefficient(n,m)
    !!{
    Evaluate the binomial coefficient, $\left({n \over m}\right)$.
    !!}
    implicit none
    integer, intent(in   ) :: n,m
    integer :: n_, m_, prefactor
    
    if (n > 0 .and. m > 0) then
       n_=n
       m_=m
       prefactor=1
    else if (m >= 0) then
       n_=-n+m-1
       m_=m
       prefactor=(-1)**m
    else if (m <= n) then
       n_=-m-1
       m_=n-m
       prefactor=(-1)**(n-m)
    else
       n_=1
       m_=1
       prefactor=1
    end if
    Binomial_Coefficient=dble(prefactor)*GSL_SF_Choose(n_,m_)
    return
  end function Binomial_Coefficient

end module Binomial_Coefficients
