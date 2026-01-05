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
  Implementation of an identity unary operator.
  !!}

  !![
  <operatorUnary name="operatorUnaryIdentity">
   <description>An identity unary operator.</description>
  </operatorUnary>
  !!]
  type, extends(operatorUnaryClass) :: operatorUnaryIdentity
     !!{
     Implementation of an identity unary operator.
     !!}
     private
   contains
     procedure :: operate   => identityOperate
     procedure :: unoperate => identityUnoperate
     procedure :: jacobian  => identityJacobian
  end type operatorUnaryIdentity

  interface operatorUnaryIdentity
     !!{
     Constructors for the \refClass{operatorUnaryIdentity} 1D distribution function class.
     !!}
     module procedure identityConstructorParameters
  end interface operatorUnaryIdentity

contains

  function identityConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{operatorUnaryIdentity} 1D distribution function class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(operatorUnaryIdentity)                :: self
    type(inputParameters      ), intent(inout) :: parameters

    self=operatorUnaryIdentity()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function identityConstructorParameters

  double precision function identityOperate(self,x)
    !!{
    Apply an identity operation.
    !!}
    implicit none
    class           (operatorUnaryIdentity), intent(inout) :: self
    double precision                       , intent(in   ) :: x
    !$GLC attributes unused :: self

    identityOperate=x
    return
  end function identityOperate

  double precision function identityUnoperate(self,f)
    !!{
    Unapply an identity operation.
    !!}
    implicit none
    class           (operatorUnaryIdentity), intent(inout) :: self
    double precision                       , intent(in   ) :: f
    !$GLC attributes unused :: self

    identityUnoperate=f
    return
  end function identityUnoperate

  double precision function identityJacobian(self,x)
    !!{
    Compute the Jacobian of an identity operation.
    !!}
    implicit none
    class           (operatorUnaryIdentity), intent(inout) :: self
    double precision                       , intent(in   ) :: x
    !$GLC attributes unused :: self, x

    identityJacobian=1.0d0
    return
  end function identityJacobian
