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
  Implementation of a logarithm unary operator.
  !!}

  !![
  <operatorUnary name="operatorUnaryLogarithm">
   <description>A logarithm unary operator.</description>
  </operatorUnary>
  !!]
  type, extends(operatorUnaryClass) :: operatorUnaryLogarithm
     !!{
     Implementation of a logarithm unary operator.
     !!}
     private
   contains
     procedure :: operate   => logarithmOperate
     procedure :: unoperate => logarithmUnoperate
  end type operatorUnaryLogarithm

  interface operatorUnaryLogarithm
     !!{
     Constructors for the \refClass{operatorUnaryLogarithm} 1D distribution function class.
     !!}
     module procedure logarithmConstructorParameters
  end interface operatorUnaryLogarithm

  ! Largest argument to exp().
  double precision, parameter :: expArgumentHuge=log(huge(0.0d0))

contains

  function logarithmConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{operatorUnaryLogarithm} 1D distribution function class which builds the object from a
    parameter set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(operatorUnaryLogarithm)                :: self
    type(inputParameters       ), intent(inout) :: parameters

    self=operatorUnaryLogarithm()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function logarithmConstructorParameters

  double precision function logarithmOperate(self,x)
    !!{
    Apply an logarithm operation.
    !!}
    implicit none
    class           (operatorUnaryLogarithm), intent(inout) :: self
    double precision                        , intent(in   ) :: x
    !$GLC attributes unused :: self

    logarithmOperate=log(x)
    return
  end function logarithmOperate

  double precision function logarithmUnoperate(self,f)
    !!{
    Unapply an logarithm operation.
    !!}
    implicit none
    class           (operatorUnaryLogarithm), intent(inout) :: self
    double precision                        , intent(in   ) :: f
    !$GLC attributes unused :: self

    if (f < expArgumentHuge) then
       logarithmUnoperate=exp (f)
    else
       logarithmUnoperate=huge(f)
    end if
    return
  end function logarithmUnoperate
