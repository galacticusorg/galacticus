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
Contains a program to test the ``coordinate`` class.
!!}

program Test_Coordinates
  !!{RST
  Tests of the ``coordinate`` class, in particular the automatic conversion between coordinate systems performed on
  assignment.
  !!}
  use :: Coordinates             , only : coordinateCartesian                 , coordinateCylindrical , coordinateSpherical , assignment(=)
  use :: Display                 , only : displayVerbositySet                 , verbosityLevelStandard
  use :: Numerical_Constants_Math, only : Pi
  use :: Unit_Tests              , only : Assert                              , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  double precision           , dimension(3,8) :: cartesian  =reshape([0.0d0,0.0d0,0.0d0,1.0d0,1.0d0,1.0d0,1.0d0,0.0d0,0.0d0,0.0d0,1.0d0,0.0d0,0.0d0,0.0d0,1.0d0,-1.0d0,0.0d0,0.0d0,0.0d0,-1.0d0,0.0d0,0.0d0,0.0d0,-1.0d0],shape(cartesian))
  double precision           , dimension(3,8) :: spherical  =reshape([0.0d0,0.0d0,0.0d0,sqrt(3.0d0),0.955316618d0,Pi/4.0d0,1.0d0,Pi/2.0d0,0.0d0,1.0d0,Pi/2.0d0,Pi/2.0d0,1.0d0,0.0d0,0.0d0,1.0d0,Pi/2.0d0,Pi,1.0d0,Pi/2.0d0,-0.5d0*Pi,1.0d0,Pi,0.0d0],shape(spherical))
  double precision           , dimension(3,8) :: cylindrical=reshape([0.0d0,0.0d0,0.0d0,sqrt(2.0d0),Pi/4.0d0,1.0d0,1.0d0,0.0d0,0.0d0,1.0d0,Pi/2.0d0,0.0d0,0.0d0,0.0d0,1.0d0,1.0d0,Pi,0.0d0,1.0d0,-0.5d0*Pi,0.0d0,0.0d0,0.0d0,-1.0d0],shape(cylindrical))
  type            (coordinateCartesian  )     :: coordinatesCartesian
  type            (coordinateSpherical  )     :: coordinatesSpherical
  type            (coordinateCylindrical)     :: coordinatesCylindrical
  integer                                     :: i
  character       (len=43)                    :: groupName

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Coordinate system conversions")

  do i=1,size(cartesian,dim=2)
     write (groupName,'(a,3(f4.1,a))') 'convert point [',cartesian(1,i),',',cartesian(2,i),',',cartesian(3,i),'] (cartesian):'
     call Unit_Tests_Begin_Group(groupName)
     ! Cartesian to spherical.
     coordinatesCartesian=cartesian  (:,i)
     coordinatesSpherical=coordinatesCartesian
     call Assert('cartesian to spherical'  ,coordinatesSpherical  %position,spherical  (:,i),absTol=1.0d-6)
     ! Cartesian to cylindrical.
     coordinatesCartesian  =cartesian  (:,i)
     coordinatesCylindrical=coordinatesCartesian
     call Assert('cartesian to cylindrical',coordinatesCylindrical%position,cylindrical(:,i),absTol=1.0d-6)
     ! Spherical to cylindrical.
     coordinatesSpherical  =spherical  (:,i)
     coordinatesCylindrical=coordinatesSpherical
     call Assert('spherical to cylindrical',coordinatesCylindrical%position,cylindrical(:,i),absTol=1.0d-6)
     ! Cylindrical to spherical.
     coordinatesCylindrical=cylindrical(:,i)
     coordinatesSpherical  =coordinatesCylindrical
     call Assert('cylindrical to spherical',coordinatesSpherical  %position,spherical  (:,i),absTol=1.0d-6)
     call Unit_Tests_End_Group()
  end do

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Test_Coordinates
