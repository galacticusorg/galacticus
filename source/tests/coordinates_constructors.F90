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
Contains a program to test the explicit constructors of the ``coordinate`` class.
!!}

program Test_Coordinates_Constructors
  !!{RST
  Tests of the explicit ``coordinateCartesian(x,y,z)``, ``coordinateSpherical(r,theta,phi)`` and
  ``coordinateCylindrical(r,phi,z)`` constructors.
  !!}
  use :: Coordinates             , only : coordinateCartesian, coordinateCylindrical , coordinateSpherical , assignment(=)
  use :: Display                 , only : displayVerbositySet, verbosityLevelStandard
  use :: Numerical_Constants_Math, only : Pi
  use :: Unit_Tests              , only : Assert             , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  type            (coordinateCartesian  )               :: coordinatesCartesian
  type            (coordinateSpherical  )               :: coordinatesSpherical
  type            (coordinateCylindrical)               :: coordinatesCylindrical
  double precision                       , dimension(3) :: cartesian

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Coordinate constructors")

  ! Cartesian constructor: components are stored directly, and the accessors return them.
  coordinatesCartesian=coordinateCartesian(1.0d0,2.0d0,3.0d0)
  call Assert('cartesian(x,y,z) components',[coordinatesCartesian%x(),coordinatesCartesian%y(),coordinatesCartesian%z()],[1.0d0,2.0d0,3.0d0],absTol=1.0d-12)

  ! Null Cartesian constructor returns the origin.
  coordinatesCartesian=coordinateCartesian()
  call Assert('cartesian() is the origin',[coordinatesCartesian%x(),coordinatesCartesian%y(),coordinatesCartesian%z()],[0.0d0,0.0d0,0.0d0],absTol=1.0d-12)

  ! Spherical constructor: components are stored directly, and conversion to Cartesian is correct. The point
  ! (r,θ,φ)=(1,π/2,0) is the Cartesian point (1,0,0).
  coordinatesSpherical=coordinateSpherical(1.0d0,0.5d0*Pi,0.0d0)
  call Assert('spherical(r,theta,phi) components',[coordinatesSpherical%r(),coordinatesSpherical%theta(),coordinatesSpherical%phi()],[1.0d0,0.5d0*Pi,0.0d0],absTol=1.0d-12)
  cartesian=coordinatesSpherical%toCartesian()
  call Assert('spherical -> cartesian conversion',cartesian,[1.0d0,0.0d0,0.0d0],absTol=1.0d-6)

  ! Cylindrical constructor: components are stored directly, and conversion to Cartesian is correct. The point
  ! (r,φ,z)=(1,0,2) is the Cartesian point (1,0,2).
  coordinatesCylindrical=coordinateCylindrical(1.0d0,0.0d0,2.0d0)
  call Assert('cylindrical(r,phi,z) components',[coordinatesCylindrical%r(),coordinatesCylindrical%phi(),coordinatesCylindrical%z()],[1.0d0,0.0d0,2.0d0],absTol=1.0d-12)
  cartesian=coordinatesCylindrical%toCartesian()
  call Assert('cylindrical -> cartesian conversion',cartesian,[1.0d0,0.0d0,2.0d0],absTol=1.0d-6)

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Test_Coordinates_Constructors
