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
Contains a program to coordinate system functions.
!!}

program Test_Coordinate_Systems
  !!{
  Tests of coordinate system functions.
  !!}
  use :: Coordinate_Systems      , only : Coordinates_Cartesian_To_Cylindrical, Coordinates_Cartesian_To_Spherical, Coordinates_Cylindrical_To_Spherical, Coordinates_Spherical_To_Cylindrical
  use :: Display                 , only : displayVerbositySet                 , verbosityLevelStandard
  use :: Numerical_Constants_Math, only : Pi
  use :: Unit_Tests              , only : Assert                              , Unit_Tests_Begin_Group            , Unit_Tests_End_Group                , Unit_Tests_Finish
  implicit none
  double precision        , dimension(3,8) :: cartesian  =reshape([0.0d0,0.0d0,0.0d0,1.0d0,1.0d0,1.0d0,1.0d0,0.0d0,0.0d0,0.0d0,1.0d0,0.0d0,0.0d0,0.0d0,1.0d0,-1.0d0,0.0d0,0.0d0,0.0d0,-1.0d0,0.0d0,0.0d0,0.0d0,-1.0d0],shape(cartesian))
  double precision        , dimension(3,8) :: spherical  =reshape([0.0d0,0.0d0,0.0d0,sqrt(3.0d0),0.955316618d0,Pi/4.0d0,1.0d0,Pi/2.0d0,0.0d0,1.0d0,Pi/2.0d0,Pi/2.0d0,1.0d0,0.0d0,0.0d0,1.0d0,Pi/2.0d0,Pi,1.0d0,Pi/2.0d0,-0.5d0*Pi,1.0d0,Pi,0.0d0],shape(spherical))
  double precision        , dimension(3,8) :: cylindrical=reshape([0.0d0,0.0d0,0.0d0,sqrt(2.0d0),Pi/4.0d0,1.0d0,1.0d0,0.0d0,0.0d0,1.0d0,Pi/2.0d0,0.0d0,0.0d0,0.0d0,1.0d0,1.0d0,Pi,0.0d0,1.0d0,-0.5d0*Pi,0.0d0,0.0d0,0.0d0,-1.0d0],shape(cylindrical))
  integer                                  :: i
  character       (len=43)                 :: groupName

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Coordinate systems")

  do i=1,size(cartesian,dim=2)
     write (groupName,'(a,3(f4.1,a))') 'convert point [',cartesian(1,i),',',cartesian(2,i),',',cartesian(3,i),'] (cartesian):'
     call Unit_Tests_Begin_Group(groupName)
     call Assert('cartesian to spherical'  ,Coordinates_Cartesian_To_Spherical  (cartesian  (:,i)),spherical  (:,i),absTol=1.0d-6)
     call Assert('cartesian to cylindrical',Coordinates_Cartesian_To_Cylindrical(cartesian  (:,i)),cylindrical(:,i),absTol=1.0d-6)
     call Assert('spherical to cylindrical',Coordinates_Spherical_To_Cylindrical(spherical  (:,i)),cylindrical(:,i),absTol=1.0d-6)
     call Assert('cylindrical to spherical',Coordinates_Cylindrical_To_Spherical(cylindrical(:,i)),spherical  (:,i),absTol=1.0d-6)
     call Unit_Tests_End_Group()
  end do

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Test_Coordinate_Systems
