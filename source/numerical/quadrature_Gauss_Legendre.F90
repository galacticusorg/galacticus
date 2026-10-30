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
Contains a module which provides Gauss-Legendre quadrature rules.
!!}

module Numerical_Quadrature_Gauss_Legendre
  !!{RST
  Provides Gauss-Legendre quadrature rules.
  !!}
  implicit none
  private
  public :: gaussLegendreRule

contains

  subroutine gaussLegendreRule(order,abscissae,weights)
    !!{RST
    Return the abscissae and weights of the Gauss-Legendre rule of the given ``order`` on the interval
    :math:`[0,1]`.

    The nodes are the roots of the Legendre polynomial of that order, found by Newton iteration from the standard
    Chebyshev-like starting guess, with the polynomial and its derivative evaluated by the usual recurrence. Both are
    then mapped from :math:`[-1,1]` onto :math:`[0,1]`, and the weights scaled by the half-width of the interval so
    that they sum to unity---which is what makes the result an average rather than an integral.

    It is used for averages over orientation, both of dust attenuation and of the orientation of a component relative
    to its galaxy, which is why it lives in a module of its own.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    integer                       , intent(in   )               :: order
    double precision, allocatable , intent(inout), dimension(:) :: abscissae                     , weights
    double precision              , parameter                   :: toleranceRelative     =1.0d-15
    integer                       , parameter                   :: countIterationsMaximum=100
    double precision                                            :: root                          , rootPrevious , &
         &                                                         legendre                      , legendreLower, &
         &                                                         legendreLowerLower            , derivative
    integer                                                     :: i                             , j            , &
         &                                                         countIterations

    if (allocated(abscissae)) deallocate(abscissae)
    if (allocated(weights  )) deallocate(weights  )
    allocate(abscissae(order))
    allocate(weights  (order))
    do i=1,order
       ! Initial guess for the i'th root of the Legendre polynomial of this order.
       root=cos(Pi*(dble(i)-0.25d0)/(dble(order)+0.5d0))
       countIterations=0
       do
          ! Evaluate the Legendre polynomial and its derivative at the current estimate by recurrence.
          legendre     =1.0d0
          legendreLower=0.0d0
          do j=1,order
             legendreLowerLower=legendreLower
             legendreLower     =legendre
             legendre          =(dble(2*j-1)*root*legendreLower-dble(j-1)*legendreLowerLower)/dble(j)
          end do
          derivative  =dble(order)*(root*legendre-legendreLower)/(root**2-1.0d0)
          rootPrevious=root
          root        =rootPrevious-legendre/derivative
          countIterations=countIterations+1
          if (abs(root-rootPrevious) <= toleranceRelative*abs(root) .or. countIterations >= countIterationsMaximum) exit
       end do
       ! Map from [-1,1] onto [0,1]. The weights are halved along with the interval, so that they sum to unity and
       ! the quadrature returns a mean.
       abscissae(i)=0.5d0*(1.0d0-root)
       weights  (i)=1.0d0/((1.0d0-root**2)*derivative**2)
    end do
    return
  end subroutine gaussLegendreRule

end module Numerical_Quadrature_Gauss_Legendre
