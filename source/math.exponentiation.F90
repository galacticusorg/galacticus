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
Contains a module which provides fast exponentiation utilizing tables.
!!}

module Math_Exponentiation
  !!{
  Provides a fast exponentiation class which utilizes tables to do rapid exponentiation in a limited range of argument for fixed
  exponent, along with other exponentiation functions.
  !!}
  use :: Tables, only : table1DLinearLinear
  implicit none
  private
  public :: cubeRoot

  type, public :: fastExponentiator
     !!{
     A fast exponentiation class which utilizes tables to do rapid exponentiation in a limited range of argument for fixed
     exponent. The tabulation used is linear in both the argument and the result to avoid expensive $\log()$ and $\exp()$
     operations. Outside of the valid range the standard exponentiation can be used, or an abort can occur.
     !!}
     double precision                      :: rangeMinimum     , rangeMaximum, &
          &                                   density
     double precision                      :: exponent
     type            (table1DLinearLinear) :: solution
     logical                               :: abortOutsideRange
   contains
     !![
     <methods>
       <method description="Evaluate {\normalfont \ttfamily x}$^y$ using table look-up." method="exponentiate" />
     </methods>
     !!]
     procedure :: exponentiate => fastExponentiatorExponentiate
  end type fastExponentiator

  interface fastExponentiator
     module procedure fastExponentiatorConstructor
  end interface fastExponentiator

contains

  function fastExponentiatorConstructor(rangeMinimum,rangeMaximum,exponent,density,abortOutsideRange) result (self)
    !!{
    Constructor for the fast exponentiator class.
    !!}
    implicit none
    type            (fastExponentiator)                :: self
    double precision                   , intent(in   ) :: rangeMinimum     , rangeMaximum, &
         &                                                exponent         , density
    logical                            , intent(in   ) :: abortOutsideRange
    integer                                            :: pointCount
    !![
    <constructorAssign variables="rangeMinimum, rangeMaximum, exponent, density, abortOutsideRange"/>
    !!]

    ! Construct the table.
    pointCount=int((rangeMaximum-rangeMinimum)*density)+1
    call self%solution%create  (rangeMinimum,rangeMaximum,pointCount)
    call self%solution%populate(self%solution%xs()**exponent)
    return
  end function fastExponentiatorConstructor

  double precision function fastExponentiatorExponentiate(self,x)
    !!{
    Evaluate the result of an exponentiation operation.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (fastExponentiator), intent(inout) :: self
    double precision                   , intent(in   ) :: x

    if     (                       &
         &   x < self%rangeMinimum &
         &  .or.                   &
         &   x > self%rangeMaximum &
         & ) then
       ! Argument is outside of range - either abort or compute using intrinsic exponentiation function.
       if (self%abortOutsideRange) then
          fastExponentiatorExponentiate=0.0d0
          call Error_Report('argument is outside range'//{introspection:location})
       else
          fastExponentiatorExponentiate=x**self%exponent
       end if
    else
       ! Argument is within range - interpolate in table to find the solution.
       fastExponentiatorExponentiate=self%solution%interpolate(x)
    end if
    return
  end function fastExponentiatorExponentiate

  pure double precision function cubeRoot(x)
    !!{
    Utilize the fast {\normalfont \ttfamily cbrt()} function from the standard C library for computing cube roots.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_double
    implicit none
    double precision, intent(in   ) :: x

    interface
       !!{
       Interface to the standard C library cube root function.
       !!}
       pure real(kind=c_double) function cbrt(x) bind(c)
         import
         real(kind=c_double), value :: x
       end function cbrt
    end interface

    cubeRoot=cbrt(x)
    return
  end function cubeRoot

end module Math_Exponentiation
