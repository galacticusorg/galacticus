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
Contains a program to test integration routines.
!!}

program Test_Integration
  !!{
  Tests that numerical integration routines work.
  !!}
  use :: Display                   , only : displayVerbositySet, verbosityLevelStandard
  use :: Numerical_Constants_Math  , only : Pi
  use :: Numerical_Integration     , only : integrator
  use :: Test_Integration_Functions, only : Integrand1         , Integrand2            , Integrand3          , Integrand4
  use :: Unit_Tests                , only : Assert             , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  double precision                          :: integral
  type            (integrator), allocatable :: integrator_

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)

  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Numerical integration")

  ! Test simple integrations.
  allocate(integrator_)
  integrator_=integrator(Integrand1,toleranceRelative=1.0d-6)
  integral=integrator_%integrate(0.0d0,1.0d0)
  deallocate(integrator_)
  call Assert("integrate f(x)=x          from x=0……1"          ,integral,0.5d0             ,relTol=1.0d-6)

  allocate(integrator_)
  integrator_=integrator(Integrand2,toleranceAbsolute=1.0d-6)
  integral=integrator_%integrate(0.0d0,2.0d0*Pi)
  deallocate(integrator_)
  call Assert("integrate f(x)=sin(x)     from x=0…2π"          ,integral,0.0d0             ,absTol=1.0d-6)

  allocate(integrator_)
  integrator_=integrator(Integrand3,toleranceRelative=1.0d-6)
  integral=integrator_%integrate(0.0d0,10.0d0)
  deallocate(integrator_)
  call Assert("integrate f(x)=1/√x       from x=0…10"          ,integral,2.0d0*sqrt(10.0d0),relTol=1.0d-6)

  ! Test 2D integrations.
  allocate(integrator_)
  integrator_=integrator(Integrand4,toleranceRelative=1.0d-6)
  integral=integrator_%integrate(0.0d0,2.0d0*Pi)
  deallocate(integrator_)
  call Assert("integrate f(x,y)=y·cos(x) from x=0…2π and y=0…x",integral,2.0d0*Pi          ,relTol=1.0d-6)

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()

end program Test_Integration
