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
Contains a module which performs numerical differentiation.
!!}

! Add dependency on GSL library.
!; gsl

module Numerical_Differentiation
  !!{
  Implements numerical differentiation.
  !!}
  use, intrinsic :: ISO_C_Binding, only : c_double   , c_int             , c_ptr              , c_null_ptr
  use            :: Interface_GSL, only : gslFunction, gslFunctionDestroy, gslFunctionTemplate
  implicit none
  private
  public :: differentiator

  type differentiator
     !!{
     Type which computes numerical derivatives.
     !!}
     private
     type(c_ptr) :: f=c_null_ptr
   contains
     !![
     <methods>
       <method description="Returns the derivative of the function at argument {\normalfont \ttfamily x}." method="derivative" />
     </methods>
     !!]
     final     ::               differentiatorDestructor
     procedure :: derivative => differentiatorDerivative
  end type differentiator

  interface differentiator
     !!{
     Constructors for the numerical derivatives class.
     !!}
     module procedure differentiatorConstructorInternal
  end interface differentiator

  interface
     function gsl_deriv_central(f,x,h,result_,abserr_) bind(c)
       import c_int, c_ptr, c_double
       integer(c_int   )                       :: gsl_deriv_central
       type   (c_ptr   ), intent(in   ), value :: f
       real   (c_double), intent(in   ), value :: x                , h
       real   (c_double), intent(  out)        :: result_          , abserr_
     end function gsl_deriv_central
  end interface

contains

  function differentiatorConstructorInternal(f) result(self)
    !!{
    Constructor for the numerical derivative class. Must be passed the function {\normalfont \ttfamily f} for which derivatives
    will be computed.
    !!}
    implicit none
    type     (differentiator     ) :: self
    procedure(gslFunctionTemplate) :: f

    self%f=gslFunction(f)
    return
  end function differentiatorConstructorInternal

  subroutine differentiatorDestructor(self)
    !!{
    Destructor for the numerical derivative class.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_associated
    implicit none
    type(differentiator), intent(inout) :: self

    if (c_associated(self%f)) call gslFunctionDestroy(self%f)
    return
  end subroutine differentiatorDestructor

  double precision function differentiatorDerivative(self,x,h,errorAbsolute)
    !!{
    Compute a numerical derivative of the function at {\normalfont \ttfamily x}. The initial stepsize is {\normalfont \ttfamily
    h}. If present, the absolute error estimate is returned in {\normalfont \ttfamily errorAbsolute}.
    !!}
    use :: Error, only : Error_Report, errorStatusSuccess
    implicit none
    class           (differentiator), intent(in   )           :: self
    double precision                , intent(in   )           :: x             , h
    double precision                , intent(  out), optional :: errorAbsolute
    double precision                                          :: errorAbsolute_
    integer         (c_int         )                          :: status

    status=gsl_deriv_central(self%f,x,h,differentiatorDerivative,errorAbsolute_)
    if (status /= errorStatusSuccess) call Error_Report('failed to compute numerical derivative'//{introspection:location})
    if (present(errorAbsolute)) errorAbsolute=errorAbsolute_
    return
  end function differentiatorDerivative

end module Numerical_Differentiation
