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
Contains a module of functions for root finding unit tests.
!!}

module Test_Multidimensional_Minimizer_Functions
  !!{
  Contains functions for root finding unit tests.
  !!}
  implicit none
  private
  public :: minimizerFunction_, minimizeFunctionDerivative_, minimizeFunctionBoth_

  ! Parameters of the function to be minimized.
  double precision, dimension(0:4) :: paraboloidParameter=[1.0d0, 2.0d0, 10.0d0, 20.0d0, 30.0d0]
  
contains

  double precision function minimizerFunction_(x)
    !!{
    Evaluate the value of a multidimensional function for minimization.
    !!}
    implicit none
    double precision, intent(in   ), dimension(:) :: x

    minimizerFunction_=+paraboloidParameter(2)*(x(1)-paraboloidParameter(0))**2 &
         &             +paraboloidParameter(3)*(x(2)-paraboloidParameter(1))**2 &
         &             +paraboloidParameter(4)
    return
  end function minimizerFunction_

  function minimizeFunctionDerivative_(x) result(gradient)
    !!{
    Evaluate the gradient of a multidimensional function for minimization.
    !!}
    implicit none
    double precision, intent(in   ), dimension(     : ) :: x
    double precision               , dimension(size(x)) :: gradient

    gradient(1)=+2.0d0*paraboloidParameter(2)*(x(1)-paraboloidParameter(0))
    gradient(2)=+2.0d0*paraboloidParameter(3)*(x(2)-paraboloidParameter(1))
    return
  end function minimizeFunctionDerivative_

  subroutine minimizeFunctionBoth_(x,functionValue,gradient)
    !!{
    Evaluate the value and gradient of a multidimensional function for minimization.
    !!}
    implicit none
    double precision, intent(in   ), dimension(     : ) :: x
    double precision, intent(  out)                     :: functionValue
    double precision, intent(  out), dimension(size(x)) :: gradient

    functionValue   =+      paraboloidParameter(2)*(x(1)-paraboloidParameter(0))**2 &
         &           +      paraboloidParameter(3)*(x(2)-paraboloidParameter(1))**2 &
         &           +      paraboloidParameter(4)
    gradient     (1)=+2.0d0*paraboloidParameter(2)*(x(1)-paraboloidParameter(0))
    gradient     (2)=+2.0d0*paraboloidParameter(3)*(x(2)-paraboloidParameter(1))
    return
  end subroutine minimizeFunctionBoth_

end module Test_Multidimensional_Minimizer_Functions
