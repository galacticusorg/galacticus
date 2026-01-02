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
Contains a program to test differentiation functions.
!!}

program Test_Differentiation
  !!{
  Tests that numerical differentiation functions work.
  !!}
  use :: Display                       , only : displayVerbositySet, verbosityLevelStandard
  use :: Numerical_Differentiation     , only : differentiator
  use :: Test_Differentiation_Functions, only : function1
  use :: Unit_Tests                    , only : Assert             , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  integer                         , parameter            :: countStep      =10
  double precision                , dimension(countStep) :: x              =[0.0d0,0.1d0,0.3d0,0.5d0,0.7d0,0.9d0,1.1d0,1.3d0,1.5d0,1.7d0], y
  double precision                , parameter            :: stepSize       =0.1d0
  type            (differentiator)                       :: differentiator_
  integer                                                :: i

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Numerical differentiation")
  differentiator_=differentiator(function1)
  do i=1,countStep
     y(i)=differentiator_%derivative(x(i),stepSize)
  end do
  call Assert("d/dx(sin(x))=cos(x)",y,cos(x),relTol=1.0d-6)
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()

end program Test_Differentiation
