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

!+    Contributions to this file made by: Andrew Benson, Claude.
!!{RST
Contains a module which provides integration of tabulated functions.
!!}

module Numerical_Integration_Tabulated
  !!{RST
  Provides integration of functions known only at a set of tabulated points, by the trapezoidal and Simpson rules.

  These are for functions which have already been tabulated, or which are cheap enough that a fixed rule on a fine grid
  is the simplest approach. For functions which can be evaluated anywhere, the adaptive integrators of
  ``Numerical_Integration`` and ``Numerical_Integration2`` should usually be preferred.
  !!}
  implicit none
  private
  public :: integrateTrapezoidal, integrateSimpson

contains

  double precision function integrateTrapezoidal(x,y) result(integral)
    !!{RST
    Return the integral of ``y`` over ``x`` by the trapezoidal rule. The abscissae ``x`` need not be uniformly spaced,
    and the integral is negative if they decrease.
    !!}
    use :: Error, only : Error_Report
    implicit none
    double precision, intent(in   ), dimension(:) :: x, y

    if (size(x) /= size(y)) call Error_Report('`x` and `y` must be the same size'//{introspection:location})
    integral=sum(0.5d0*(y(2:size(y))+y(1:size(y)-1))*(x(2:size(x))-x(1:size(x)-1)))
    return
  end function integrateTrapezoidal

  double precision function integrateSimpson(x,y) result(integral)
    !!{RST
    Return the integral of ``y`` over ``x`` by the composite Simpson rule. The abscissae ``x`` must be uniformly spaced,
    and there must be an odd number of them, at least three.
    !!}
    use :: Error, only : Error_Report
    implicit none
    double precision, intent(in   ), dimension(:) :: x                         , y
    ! Tolerance, relative to the mean spacing, on the departure of the abscissae from uniform spacing.
    double precision, parameter                   :: toleranceUniformity=1.0d-6
    integer                                       :: n

    n=size(x)
    if (size(y) /= n                                                                                  ) &
         & call Error_Report('`x` and `y` must be the same size'                                    //{introspection:location})
    if (n < 3 .or. mod(n,2) == 0                                                                      ) &
         & call Error_Report('Simpson''s rule requires an odd number of points, at least three'     //{introspection:location})
    if (any(abs(x(2:n)-x(1:n-1)-(x(n)-x(1))/dble(n-1)) > toleranceUniformity*abs(x(n)-x(1))/dble(n-1))) &
         & call Error_Report('Simpson''s rule requires uniformly spaced abscissae'                  //{introspection:location})
    integral=+(x(2)-x(1))/3.0d0       &
         &   *(                       &
         &     +      y(1)            &
         &     +4.0d0*sum(y(2:n-1:2)) &
         &     +2.0d0*sum(y(3:n-2:2)) &
         &     +      y(n)            &
         &    )
    return
  end function integrateSimpson

end module Numerical_Integration_Tabulated
