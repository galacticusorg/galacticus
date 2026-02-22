!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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
Contains a module that provides a wrapper class for testing the {\normalfont \ttfamily resourceManager} class.
!!}

module Test_Resource_Manager_Wrapper
  !!{
  Provides a wrapper class for testing the {\normalfont \ttfamily resourceManager} class.
  !!}
  use :: Resource_Manager, only : resourceManager
  private
  public :: resourceHolder, countDestructions
  
  type :: sharedResource
     !!{
     A dummy type to act as a shared resource.
     !!}
     private
     integer :: dummy
  end type sharedResource

  type :: resourceHolder
     !!{
     A dummy type to act as a class that wants to retain a shared resource.
     !!}
     type(sharedResource ), pointer :: sharedObject => null()
     type(resourceManager)          :: manager
   contains
     final :: resourceHolderDestruct
  end type resourceHolder

  interface resourceHolder
     !!{
     Constructors for the dummy resource holder class.
     !!}
     module procedure resourceHolderConstructFirst
     module procedure resourceHolderConstructSecond
  end interface resourceHolder

  ! Counter for the number of calls to our destructor.
  integer :: countDestructions=0
  
contains

  function resourceHolderConstructFirst() result(self)
    !!{
    First-level constructor for the dummy resource holder class.
    !!}
    use :: Display, only : displayMessage
    implicit none
    type(resourceHolder) :: self

    call displayMessage("resourceHolder: construct (1st phase; start)")
    self=resourceHolder(3)
    call displayMessage("resourceHolder: construct (1st phase; end)")
    return
  end function resourceHolderConstructFirst

  function resourceHolderConstructSecond(i) result(self)
    !!{
    Second-level constructor for the dummy resource holder class.
    !!}
    use :: Display, only : displayMessage
    implicit none
    type   (resourceHolder)                :: self
    integer                , intent(in   ) :: i
    class  (*             ), pointer       :: actual_
    !$GLC attributes unused :: i

    call displayMessage("resourceHolder: construct (2nd phase; start)")
    allocate(self%sharedObject)
    self%sharedObject=sharedResource(1)
    !![
    <workaround type="gfortran" PR="105807" url="https:&#x2F;&#x2F;gcc.gnu.org&#x2F;bugzilla&#x2F;show_bug.cgi=105807">
    <description>ICE when passing a derived type component to a class(*) function argument.</description>
    </workaround>
    !!]
    actual_      => self%sharedObject
    self%manager =  resourceManager(actual_)
    call displayMessage("resourceHolder: construct (2nd phase; end)")
    return
  end function resourceHolderConstructSecond

  subroutine resourceHolderDestruct(self)
    !!{
    Destructor for the dummy resource holder class.
    !!}
    use :: Display, only : displayMessage
    implicit none
    type(resourceHolder), intent(inout) :: self

    call displayMessage("resourceHolder: destruct")
    countDestructions=countDestructions+1
    return
  end subroutine resourceHolderDestruct

end module Test_Resource_Manager_Wrapper
