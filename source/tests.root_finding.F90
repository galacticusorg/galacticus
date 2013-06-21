!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a program to test root finding routines.

program Test_Root_Finding
  !% Tests that routine finding routines work.
  use Unit_Tests
  use ISO_Varying_String
  use FGSL
  use Root_Finder
  use Test_Root_Finding_Functions
  use, intrinsic :: ISO_C_Binding
  implicit none
  type            (rootFinder)               :: finder         
  double precision                           :: xGuess, xRoot  
  double precision            , dimension(2) :: xRange         
  
  ! Begin unit tests.                                                          
  call Unit_Tests_Begin_Group("Root finding")

  ! Test root finding.
  xRange=[-1.0d0,+1.0d0]
  call finder%tolerance(1.0d-6,1.0d-6)
  call finder%rootFunction(Root_Function_1)
  xRoot=finder%find(rootRange=xRange)
  call Assert('root of f(x)=x',xRoot,0.0d0,absTol=1.0d-6,relTol=1.0d-6)

  xRange=[-1.0d0,+1.0d0]
  call finder%tolerance(1.0d-6,1.0d-6)
  call finder%rootFunction(Root_Function_2)
  xRoot=finder%find(rootRange=xRange)
  call Assert('root of f(x)=x² - 5x + 1 in range -1 < x <  1',xRoot,0.5d0*(5.0d0-sqrt(21.0d0)),absTol=1.0d-6,relTol=1.0d-6)

  xRange=[2.0d0,10.0d0]
  call finder%tolerance(1.0d-6,1.0d-6)
  call finder%rootFunction(Root_Function_2)
  xRoot=finder%find(rootRange=xRange)
  call Assert('root of f(x)=x² - 5x + 1 in range  2 < x < 10',xRoot,0.5d0*(5.0d0+sqrt(21.0d0)),absTol=1.0d-6,relTol=1.0d-6)
  
  xRange=[-1.0d0,+1.0d0]
  call finder%tolerance(1.0d-6,1.0d-6)
  call finder%rootFunction(Root_Function_3)
  xRoot=finder%find(rootRange=xRange)
  call Assert('root of f(x)=x × exp(-x) + 1',xRoot,-0.567143d0,absTol=1.0d-6,relTol=1.0d-6)

  ! Test with root bracketing.
  xGuess=0.0d0
  call finder%tolerance(1.0d-6,1.0d-6)
  call finder%rangeExpand(rangeExpandUpward=0.1d0,rangeExpandDownward=-0.1d0,rangeExpandType=rangeExpandAdditive)
  call finder%rootFunction(Root_Function_3)
  xRoot=finder%find(rootGuess=xGuess)
  call Assert('root of f(x)=x × exp(-x) + 1; with bracketing',xRoot,-0.567143d0,absTol=1.0d-6,relTol=1.0d-6)

  ! Test with root bracketing and limit.
  xGuess=0.0d0
  call finder%tolerance(1.0d-6,1.0d-6)
  call finder%rangeExpand(rangeExpandUpward=0.1d0,rangeExpandDownward=-0.1d0,rangeExpandType=rangeExpandAdditive,rangeUpwardLimit=1.0d0,rangeDownwardLimit=-5.0d0)
  call finder%rootFunction(Root_Function_3)
  xRoot=finder%find(rootGuess=xGuess)
  call Assert('root of f(x)=x × exp(-x) + 1; with bracketing + limit',xRoot,-0.567143d0,absTol=1.0d-6,relTol=1.0d-6)
  
  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Test_Root_Finding
