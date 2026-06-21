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

  !!{RST
  Implementation of a logarithm unary operator.
  !!}

  !![
  <operatorUnary name="operatorUnaryLogarithm" docformat="rst">
   <description>
   A unary operator implementing the natural logarithm :math:`f(x) = \ln(x)`; applying this operator transforms a positive scalar to log-space, and its inverse is the exponential :math:`f^{-1}(f) = \exp(f)`.
   </description>
  </operatorUnary>
  !!]
  type, extends(operatorUnaryClass) :: operatorUnaryLogarithm
     !!{RST
     Implementation of a logarithm unary operator.
     !!}
     private
   contains
     procedure :: operate   => logarithmOperate
     procedure :: unoperate => logarithmUnoperate
     procedure :: jacobian  => logarithmJacobian
  end type operatorUnaryLogarithm

  interface operatorUnaryLogarithm
     !!{RST
     Constructors for the :galacticus-class:`operatorUnaryLogarithm` 1D distribution function class.
     !!}
     module procedure logarithmConstructorParameters
  end interface operatorUnaryLogarithm

  ! Largest argument to exp().
  double precision, parameter :: expArgumentHuge=log(huge(0.0d0))

contains

  function logarithmConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`operatorUnaryLogarithm` 1D distribution function class which builds the object from a parameter set.
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
    !!{RST
    Apply a logarithm operation.
    !!}
    implicit none
    class           (operatorUnaryLogarithm), intent(inout) :: self
    double precision                        , intent(in   ) :: x
    !$GLC attributes unused :: self

    logarithmOperate=log(x)
    return
  end function logarithmOperate

  double precision function logarithmUnoperate(self,f)
    !!{RST
    Unapply a logarithm operation.
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

  double precision function logarithmJacobian(self,x)
    !!{RST
    Compute the Jacobian of a logarithm operation.
    !!}
    implicit none
    class           (operatorUnaryLogarithm), intent(inout) :: self
    double precision                        , intent(in   ) :: x
    !$GLC attributes unused :: self

    logarithmJacobian=1.0d0/x
    return
  end function logarithmJacobian
