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
Contains a module which provides analytic test integrands for the multi-dimensional integration unit tests.
!!}

module Test_Integration_ND_Functions
  !!{RST
  Provides analytic test integrands for the multi-dimensional integration unit tests.
  !!}
  implicit none
  private
  public :: integrandPolynomial5 , integrandPolynomial7, &
       &    integrandSinCos      , integrandXSquaredCos, &
       &    integrandRootXYSquare, integrandYOverRootX , &
       &    integrandSeparable3D , integrandSphere     , &
       &    radiusSphere         , integrandProduct4D

  ! Radius of the sphere used by the discontinuous integrand.
  double precision :: radiusSphere=1.0d0

contains

  double precision function integrandPolynomial5(x)
    !!{RST
    A polynomial of total degree five, :math:`f = 1 + 2 x_1 + 3 x_1^2 x_2 + x_1^2 x_2^2 x_3 + 4 x_2^5`. The degree-5 rule
    embedded in the degree-7 cubature integrates this exactly, so both rules agree and the error estimate must vanish.
    !!}
    implicit none
    double precision, intent(in   ), dimension(:) :: x

    integrandPolynomial5=+1.0d0                       &
         &               +2.0d0*x(1)                  &
         &               +3.0d0*x(1)**2*x(2)          &
         &               +      x(1)**2*x(2)**2*x(3)  &
         &               +4.0d0*        x(2)**5
    return
  end function integrandPolynomial5

  double precision function integrandPolynomial7(x)
    !!{RST
    A polynomial of total degree seven, :math:`f = x_1^7 + x_1^4 x_2^2 x_3 + x_1^3 x_2^2 x_3^2 + x_2^6`. The degree-7 cubature
    integrates this exactly, while the embedded degree-5 rule does not.
    !!}
    implicit none
    double precision, intent(in   ), dimension(:) :: x

    integrandPolynomial7=+x(1)**7                     &
         &               +x(1)**4*x(2)**2*x(3)        &
         &               +x(1)**3*x(2)**2*x(3)**2     &
         &               +        x(2)**6
    return
  end function integrandPolynomial7

  double precision function integrandSinCos(x)
    !!{RST
    The integrand :math:`f = \sin(x_1^2) \cos(x_2)`, as used by the two-dimensional integration unit tests.
    !!}
    implicit none
    double precision, intent(in   ), dimension(:) :: x

    integrandSinCos=sin(x(1)**2)*cos(x(2))
    return
  end function integrandSinCos

  double precision function integrandXSquaredCos(x)
    !!{RST
    The integrand :math:`f = x_1^2 \cos(x_2)`, as used by the two-dimensional integration unit tests.
    !!}
    implicit none
    double precision, intent(in   ), dimension(:) :: x

    integrandXSquaredCos=x(1)**2*cos(x(2))
    return
  end function integrandXSquaredCos

  double precision function integrandRootXYSquare(x)
    !!{RST
    The integrand :math:`f = \sqrt{x_1} x_2^2`, as used by the two-dimensional integration unit tests. Its derivative is singular
    at :math:`x_1=0`.
    !!}
    implicit none
    double precision, intent(in   ), dimension(:) :: x

    integrandRootXYSquare=sqrt(x(1))*x(2)**2
    return
  end function integrandRootXYSquare

  double precision function integrandYOverRootX(x)
    !!{RST
    The integrand :math:`f = x_2/\sqrt{x_1}`, as used by the two-dimensional integration unit tests. This is itself singular at
    :math:`x_1=0`, though integrably so. The cubature rule never places a point on the boundary of a region, so the singularity
    is never evaluated.
    !!}
    implicit none
    double precision, intent(in   ), dimension(:) :: x

    integrandYOverRootX=x(2)/sqrt(x(1))
    return
  end function integrandYOverRootX

  double precision function integrandSeparable3D(x)
    !!{RST
    A smooth, separable three-dimensional integrand, :math:`f = \exp(x_1) \sin(x_2) \cos(x_3)`.
    !!}
    implicit none
    double precision, intent(in   ), dimension(:) :: x

    integrandSeparable3D=exp(x(1))*sin(x(2))*cos(x(3))
    return
  end function integrandSeparable3D

  double precision function integrandSphere(x)
    !!{RST
    A discontinuous integrand: unity inside a sphere of radius  radiusSphere centered on the origin, and zero outside. Its
    integral is the volume of that part of the sphere lying within the region of integration.
    !!}
    implicit none
    double precision, intent(in   ), dimension(:) :: x

    if (sum(x**2) <= radiusSphere**2) then
       integrandSphere=1.0d0
    else
       integrandSphere=0.0d0
    end if
    return
  end function integrandSphere

  double precision function integrandProduct4D(x)
    !!{RST
    A smooth, separable four-dimensional integrand, :math:`f = x_1 x_2^2 x_3^3 x_4^4`.
    !!}
    implicit none
    double precision, intent(in   ), dimension(:) :: x

    integrandProduct4D=x(1)*x(2)**2*x(3)**3*x(4)**4
    return
  end function integrandProduct4D

end module Test_Integration_ND_Functions
