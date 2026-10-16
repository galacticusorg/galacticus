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
Contains a module which provides analytic test integrands for the computational domain volume integrator unit tests.
!!}

module Test_Computational_Domain_Volume_Integrators_Functions
  !!{RST
  Provides analytic test integrands for the computational domain volume integrator unit tests. Each integrand counts the number of
  times it is evaluated, allowing the cost of the integration to be reported alongside its accuracy.
  !!}
  use :: Coordinates, only : coordinate
  implicit none
  private
  public :: countEvaluations   , countEvaluationsReset, &
       &    integrandUnity     , integrandCartesian   , &
       &    integrandSpherical , integrandCylindrical , &
       &    integrandCloud     , radiusCloud

  ! Count of the number of integrand evaluations performed.
  integer          :: countEvaluations=0

  ! Radius of the spherical "cloud" used by the discontinuous integrand.
  double precision :: radiusCloud     =1.0d0

contains

  subroutine countEvaluationsReset()
    !!{RST
    Reset the count of integrand evaluations.
    !!}
    implicit none

    countEvaluations=0
    return
  end subroutine countEvaluationsReset

  double precision function integrandUnity(coordinates)
    !!{RST
    A unit integrand---its integral over any cell equals that cell's volume.
    !!}
    implicit none
    class(coordinate), intent(in   ) :: coordinates
    !$GLC attributes unused :: coordinates

    countEvaluations=countEvaluations+1
    integrandUnity  =1.0d0
    return
  end function integrandUnity

  double precision function integrandCartesian(coordinates)
    !!{RST
    A smooth integrand, :math:`f(x,y,z) = x^2+y^2+z^2`, for use in Cartesian cells.
    !!}
    use :: Coordinates, only : coordinateCartesian
    use :: Error      , only : Error_Report
    implicit none
    class(coordinate), intent(in   ) :: coordinates

    countEvaluations=countEvaluations+1
    select type (coordinates)
    type is (coordinateCartesian)
       integrandCartesian=+coordinates%x()**2 &
            &             +coordinates%y()**2 &
            &             +coordinates%z()**2
    class default
       integrandCartesian=0.0d0
       call Error_Report('coordinates must be Cartesian'//{introspection:location})
    end select
    return
  end function integrandCartesian

  double precision function integrandSpherical(coordinates)
    !!{RST
    A smooth integrand, :math:`f(r,\theta,\phi) = \cos^2 \theta`, for use in spherical cells.
    !!}
    use :: Coordinates, only : coordinateSpherical
    use :: Error      , only : Error_Report
    implicit none
    class(coordinate), intent(in   ) :: coordinates

    countEvaluations=countEvaluations+1
    select type (coordinates)
    type is (coordinateSpherical)
       integrandSpherical=cos(coordinates%theta())**2
    class default
       integrandSpherical=0.0d0
       call Error_Report('coordinates must be spherical'//{introspection:location})
    end select
    return
  end function integrandSpherical

  double precision function integrandCylindrical(coordinates)
    !!{RST
    A smooth integrand, :math:`f(r,\phi,z) = r z`, for use in cylindrical cells.
    !!}
    use :: Coordinates, only : coordinateCylindrical
    use :: Error      , only : Error_Report
    implicit none
    class(coordinate), intent(in   ) :: coordinates

    countEvaluations=countEvaluations+1
    select type (coordinates)
    type is (coordinateCylindrical)
       integrandCylindrical=+coordinates%r() &
            &               *coordinates%z()
    class default
       integrandCylindrical=0.0d0
       call Error_Report('coordinates must be cylindrical'//{introspection:location})
    end select
    return
  end function integrandCylindrical

  double precision function integrandCloud(coordinates)
    !!{RST
    A discontinuous integrand: unity inside a sphere of radius  radiusCloud centered on the origin, and zero outside. This mimics
    the constant density cloud mass distributions used in radiative transfer models, for which cells straddling the surface of the
    cloud have a discontinuous integrand.
    !!}
    implicit none
    class(coordinate), intent(in   ) :: coordinates

    countEvaluations=countEvaluations+1
    if (coordinates%rSphericalSquared() <= radiusCloud**2) then
       integrandCloud=1.0d0
    else
       integrandCloud=0.0d0
    end if
    return
  end function integrandCloud

end module Test_Computational_Domain_Volume_Integrators_Functions
