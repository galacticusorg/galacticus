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
Contains a program which tests functionality of the event hook infrastructure.
!!}

program Tests_Events_Hooks
  !!{
  Tests the event hooks infrastructure.
  !!}
  use :: Display                   , only : displayVerbositySet       , verbosityLevelStandard
  use :: Events_Hooks              , only : eventsHooksInitialize
  use :: Tests_Event_Hook_Functions, only : hookedFunctionCalled      , hookedFunctionOrder1After2 , hookedFunctionOrder1After4, hookedFunctionOrder2After3, &
          &                                 hookedFunctionOrder2After4, hookedFunctionOrder3Before1, testEventHooksInitialize
  use :: Unit_Tests                , only : Assert                    , Unit_Tests_Begin_Group     , Unit_Tests_End_Group      , Unit_Tests_Finish
  implicit none
  integer :: testValue
  
  call displayVerbositySet(verbosityLevelStandard)
  call eventsHooksInitialize         (                 )
  call testEventHooksInitialize      (                 )
  call Unit_Tests_Begin_Group        ("Event hooks"    )
  testValue=99
  !![
  <eventHook name="testEvent">
   <callWith>testValue</callWith>
  </eventHook>
  !!]
  call Assert("Hooked functions were all called" ,hookedFunctionCalled       ,spread(.true.,1,size(hookedFunctionCalled)))
  call Assert("Hooked function 2 called after 3" ,hookedFunctionOrder2After3 ,       .true.                              )
  call Assert("Hooked function 2 called after 4" ,hookedFunctionOrder2After4 ,       .true.                              )
  call Assert("Hooked function 3 called before 1",hookedFunctionOrder3Before1,       .true.                              )
  call Assert("Hooked function 1 called after 2" ,hookedFunctionOrder1After2 ,       .true.                              )
  call Assert("Hooked function 1 called after 4" ,hookedFunctionOrder1After4 ,       .true.                              )
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish   ()
end program Tests_Events_Hooks
