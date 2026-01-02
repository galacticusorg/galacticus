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
Contains a module of functions for testing the functionality of the event hook infrastructure.
!!}

module Tests_Event_Hook_Functions
  !!{
  Module providing functions used in testing the functionality of the event hook infrastructure.
  !!}
  private
  public :: testEventHooksInitialize

  ! Records of which functions were actually called.
  logical, public, dimension(4) :: hookedFunctionCalled       =.false.
  logical, public               :: hookedFunctionOrder2After3 =.false.
  logical, public               :: hookedFunctionOrder2After4 =.false.
  logical, public               :: hookedFunctionOrder3Before1=.false.
  logical, public               :: hookedFunctionOrder1After2 =.false.
  logical, public               :: hookedFunctionOrder1After4 =.false.
  
contains

  subroutine testEventHooksInitialize()
    !!{
    Initialize all hooks into the test event.
    !!}
    use :: Events_Hooks, only : testEventEvent, openMPThreadBindingAtLevel, dependencyExact, dependencyRegEx, dependencyDirectionAfter, dependencyDirectionBefore
    implicit none
    integer                                :: dummySelf
    type   (dependencyRegEx), dimension(1) :: dependencies

    dependencies(1)=dependencyRegEx(dependencyDirectionAfter,'^function[24]$')
    call testEventEvent%attach(dummySelf,hookedFunction1,openMPThreadBindingAtLevel,label="function1",dependencies= dependencies                                                                                                 )
    call testEventEvent%attach(dummySelf,hookedFunction2,openMPThreadBindingAtLevel,label="function2",dependencies=[dependencyExact(dependencyDirectionAfter ,'function3'),dependencyExact(dependencyDirectionAfter,'function4')])
    call testEventEvent%attach(dummySelf,hookedFunction3,openMPThreadBindingAtLevel,label="function3",dependencies=[dependencyExact(dependencyDirectionBefore,'function1')                                                      ])    
    call testEventEvent%attach(dummySelf,hookedFunction4,openMPThreadBindingAtLevel,label="function4"                                                                                                                            )
    return
  end subroutine testEventHooksInitialize

  subroutine hookedFunction1(self,testValue)
    !!{
    A function used in testing the event hook infrastructure.
    !!}
    implicit none
    class  (*), intent(inout) :: self
    integer   , intent(in   ) :: testValue
    !$GLC attributes unused :: self, testValue

    hookedFunctionCalled(1)=.true.
    hookedFunctionOrder1After2=hookedFunctionCalled(2)
    hookedFunctionOrder1After4=hookedFunctionCalled(4)
    return
  end subroutine hookedFunction1
  
  subroutine hookedFunction2(self,testValue)
    !!{
    A function used in testing the event hook infrastructure.
    !!}
    implicit none
    class  (*), intent(inout) :: self
    integer   , intent(in   ) :: testValue
    !$GLC attributes unused :: self, testValue

    hookedFunctionCalled(2)=.true.
    hookedFunctionOrder2After3=hookedFunctionCalled(3)
    hookedFunctionOrder2After4=hookedFunctionCalled(4)
    return
  end subroutine hookedFunction2

  subroutine hookedFunction3(self,testValue)
    !!{
    A function used in testing the event hook infrastructure.
    !!}
    implicit none
    class  (*), intent(inout) :: self
    integer   , intent(in   ) :: testValue
    !$GLC attributes unused :: self, testValue

    hookedFunctionCalled(3)=.true.
    hookedFunctionOrder3Before1=.not.hookedFunctionCalled(1)
    return
  end subroutine hookedFunction3

  subroutine hookedFunction4(self,testValue)
    !!{
    A function used in testing the event hook infrastructure.
    !!}
    implicit none
    class  (*), intent(inout) :: self
    integer   , intent(in   ) :: testValue
    !$GLC attributes unused :: self, testValue

    hookedFunctionCalled(4)=.true.
    return
  end subroutine hookedFunction4

end module Tests_Event_Hook_Functions
