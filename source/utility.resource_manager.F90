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
Contains a module which implements a class that manages shared resources.
!!}

module Resource_Manager
  !!{
  Implements a class that manages shared resources (typically pointers to objects shared by multiple other objects) via a
  reference counting approach and destructs them when no more references exist. Similar in approach to a
  \href{https://en.cppreference.com/w/cpp/memory/shared_ptr}{\normalfont \ttfamily shared\_ptr} in C++.
  !!}
  !$ use :: OMP_Lib, only : omp_lock_kind
  private
  public :: resourceManager, resourceManagerForceReportOn, resourceManagerForceReportOff

  type :: resourceManager
     !!{
     A class that manages shared resources (typically pointers to objects shared by multiple other objects) via a reference
     counting approach and destructs them when no more references exist. Similar in approach to a
     \href{https://en.cppreference.com/w/cpp/memory/shared_ptr}{\normalfont \ttfamily shared\_ptr} in C++.  
     !!}
     class     (*            ), pointer :: resource  => null()
     integer                  , pointer :: counter   => null()
     !$ integer(omp_lock_kind), pointer :: lock      => null()
     logical                            :: reportOn_ =  .false.
   contains
     !![
     <methods>
        <method method="assignment(=)" description="Assign the reference manager, incrementing the reference count of the managed resource."/>
        <method method="release"       description="Release the managed object."                                                            />
        <method method="count"         description="Return the current reference count to the managed object."                              />
        <method method="reportOn"      description="Report on changes to the reference count to the managed object."                        />
     </methods>
     !!]
     final     :: resourceManagerDestructor
     procedure :: release                   => resourceManagerRelease
     procedure :: count                     => resourceManagerCount
     procedure :: reportOn                  => resourceManagerReportOn
     procedure :: resourceManagerAssign
     generic   :: assignment(=)             => resourceManagerAssign
  end type resourceManager

  interface resourceManager
     !!{
     Constructors for the {\normalfont \ttfamily resourceManager} class.
     !!}
     module procedure resourceManagerConstructor
  end interface resourceManager

  ! Option to force reporting of resource state.
  logical :: forceReportOn=.false.
  !$omp threadprivate(forceReportOn)
  
contains

  subroutine resourceManagerForceReportOn()
    !!{
    Force all {\normalfont \ttfamily resourceManager} objects to be created with reporting enabled.
    !!}
    implicit none

    forceReportOn=.true.
    return
  end subroutine resourceManagerForceReportOn

  subroutine resourceManagerForceReportOff()
    !!{
    Cease forcing all {\normalfont \ttfamily resourceManager} objects to be created with reporting enabled.
    !!}
    implicit none

    forceReportOn=.false.
    return
  end subroutine resourceManagerForceReportOff

  function resourceManagerConstructor(resource,reportOn,threadSafe) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily resourceManager} class. This should be called with a pointer to the resource to
    manage after it is first created.
    !!}
    use :: Display           , only : displayMessage
    use :: String_Handling   , only : operator(//)
    use :: ISO_Varying_String, only : operator(//)  , var_str
    implicit none
    type   (resourceManager)                          :: self
    class  (*              ), intent(in   ), pointer  :: resource
    logical                 , intent(in   ), optional :: reportOn, threadSafe
    !![
    <optionalArgument name="reportOn"   defaultsTo=".false." />
    <optionalArgument name="threadSafe" defaultsTo=".false." />
    !!]
    
    ! Retain a pointer to the shared resource.
    self%resource => resource
    ! Create a shared counter for the resource and initialize the reference count to 1.
    allocate(self%counter)
    self%counter=1
    ! If thread safety is requested, create an OpenMP lock.
    if (threadSafe_) then
       !$ allocate(self%lock)
       !$ call OMP_Init_Lock(self%lock)
    end if
    ! Set reporting state.
    self%reportOn_=reportOn_ .or. forceReportOn
    if (self%reportOn_) then
       call displayMessage(var_str('report on managed resource [loc:')//loc(self%resource)//' ] references - count = '//self%counter)
       call backtrace()
    end if
    return
  end function resourceManagerConstructor

  subroutine resourceManagerDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily resourceManager} class.
    !!}
    implicit none
    type(resourceManager), intent(inout) :: self

    call self%release()
    return
  end subroutine resourceManagerDestructor

  subroutine resourceManagerAssign(to,from)
    !!{
    Assign a {\normalfont \ttfamily resourceManager} object.
    !!}
    !$ use :: OMP_Lib           , only : OMP_Set_Lock  , OMP_Unset_Lock
    use    :: Display           , only : displayMessage
    use    :: String_Handling   , only : operator(//)
    use    :: ISO_Varying_String, only : operator(//)  , var_str
    implicit none
    class(resourceManager), intent(  out) :: to
    class(resourceManager), intent(in   ) :: from

    if (associated(from%counter)) then
       !$ if (associated(from%lock)) call OMP_Set_Lock(from%lock)
       ! Copy pointers to the shared resource and shared counter.
       to   %resource  => from%resource
       to   %counter   => from%counter
       !$ to%lock      => from%lock
       ! Increment the reference count to our shared object.
       to%counter   =  to  %counter  +1
       to%reportOn_ =  from%reportOn_
       if (to%reportOn_) then
          call displayMessage(var_str('increment managed resource [loc:')//loc(to%resource)//' ] references - count = '//to%counter)
          call backtrace()
       end if
       !$ if (associated(from%lock)) call OMP_Unset_Lock(from%lock)
    else
       ! No resource to manage - set null pointers.
       to   %resource => null()
       to   %counter  => null()
       !$ to%lock     => null()
    end if
    return
  end subroutine resourceManagerAssign
    
  subroutine resourceManagerRelease(self)
    !!{
    Release the managed resource.
    !!}
    !$ use :: OMP_Lib           , only : OMP_Set_Lock  , OMP_Unset_Lock
    use    :: Display           , only : displayMessage
    use    :: String_Handling   , only : operator(//)
    use    :: ISO_Varying_String, only : operator(//)  , var_str
    implicit none
    class(resourceManager), intent(inout) :: self
    
    ! If no counter has been created, then we have no resource that we are managing. We can simply return.
    if (.not.associated(self%counter)) return
    ! Lock if required.
    !$ if (associated(self%lock)) call OMP_Set_Lock(self%lock)
    ! Decrement the reference count to our shared resource.
    self%counter=self%counter-1
    if (self%reportOn_) then
       call displayMessage(var_str('decrement managed resource [loc:')//loc(self%resource)//' ] references - count = '//self%counter)
       call backtrace()
    end if
    ! If no more references to the shared resource exist we can destroy it (and destroy the shared counter also).
    if (self%counter == 0) then
       deallocate(self%resource)
       deallocate(self%counter )
       !$ if (associated(self%lock)) then
       !$    call OMP_Unset_Lock(self%lock)
       !$    deallocate(self%lock)
       !$    nullify   (self%lock)
       !$ end if
    end if
    ! Unlock if needed.
    !$ if (associated(self%lock)) call OMP_Unset_Lock(self%lock)
    ! Nullify pointers to avoid any dangling pointer issues.
    nullify   (self%counter )
    nullify   (self%resource)
    !$ nullify(self%lock    )
    return
  end subroutine resourceManagerRelease

  integer function resourceManagerCount(self)
    !!{
    Return the current reference count to the managed object.
    !!}
    !$ use :: OMP_Lib, only : OMP_Set_Lock, OMP_Unset_Lock
    implicit none
    class(resourceManager), intent(in   ) :: self

    if (associated(self%counter)) then
       !$ if (associated(self%lock)) call OMP_Set_Lock  (self%lock)
       resourceManagerCount=self%counter
       !$ if (associated(self%lock)) call OMP_Unset_Lock(self%lock)
    else
       resourceManagerCount=0
    end if
    return
  end function resourceManagerCount

  subroutine resourceManagerReportOn(self)
    !!{
    Report on the managed resource.
    !!}
    use :: Display           , only : displayMessage
    use :: String_Handling   , only : operator(//)
    use :: ISO_Varying_String, only : operator(//)  , var_str
    implicit none
    class(resourceManager), intent(inout) :: self

    self%reportOn_=.true.
    call displayMessage(var_str('report on managed resource [loc:')//loc(self%resource)//' ] references - count = '//self%counter)
    call backtrace()
    return
  end subroutine resourceManagerReportOn

end module Resource_Manager
