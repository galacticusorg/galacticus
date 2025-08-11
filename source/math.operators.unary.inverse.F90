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
  Implementation of an inverse unary operator.
  !!}

  !![
  <operatorUnary name="operatorUnaryInverse">
   <description>An inverse unary operator.</description>
  </operatorUnary>
  !!]
  type, extends(operatorUnaryClass) :: operatorUnaryInverse
     !!{
     Implementation of an inverse unary operator.
     !!}
     private
   contains
     procedure :: operate   => inverseOperate
     procedure :: unoperate => inverseUnoperate
     procedure :: jacobian  => inverseJacobian
  end type operatorUnaryInverse

  interface operatorUnaryInverse
     !!{
     Constructors for the \refClass{operatorUnaryInverse} 1D distribution function class.
     !!}
     module procedure inverseConstructorParameters
  end interface operatorUnaryInverse

contains

  function inverseConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{operatorUnaryInverse} 1D distribution function class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(operatorUnaryInverse)                :: self
    type(inputParameters     ), intent(inout) :: parameters

    self=operatorUnaryInverse()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function inverseConstructorParameters

  double precision function inverseOperate(self,x)
    !!{
    Apply an inverse operation.
    !!}
    implicit none
    class           (operatorUnaryInverse), intent(inout) :: self
    double precision                      , intent(in   ) :: x
    !$GLC attributes unused :: self

    inverseOperate=1.0d0/x
    return
  end function inverseOperate

  double precision function inverseUnoperate(self,f)
    !!{
    Unapply an inverse operation.
    !!}
    implicit none
    class           (operatorUnaryInverse), intent(inout) :: self
    double precision                       , intent(in   ) :: f
    !$GLC attributes unused :: self

    inverseUnoperate=1.0d0/f
    return
  end function inverseUnoperate

  double precision function inverseJacobian(self,x)
    !!{
    Comput the Jacobian of the inverse operation.
    !!}
    implicit none
    class           (operatorUnaryInverse), intent(inout) :: self
    double precision                      , intent(in   ) :: x
    !$GLC attributes unused :: self

    inverseJacobian=-1.0d0/x**2
    return
  end function inverseJacobian
