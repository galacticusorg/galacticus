!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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
Contains module which handles hooking of objects into events.
!!}

module Events_Hooks
  !!{
  Handles hooking of object function class into events.
  !!}
  use :: Locks              , only : ompReadWriteLock
  use :: Regular_Expressions, only : regEx
  private
  public :: hook                 , hookUnspecified, dependencyExact, dependencyRegEx

  type :: dependency
     !!{
     Base class for event hook dependencies.
     !!}
     private
     integer            :: direction
     character(len=128) :: label
  end type dependency

  type, extends(dependency) :: dependencyExact
     !!{
     Type for exactly matching event hook dependencies.
     !!}
  end type dependencyExact

  type, extends(dependency) :: dependencyRegEx
     !!{
     Type for matching event hook dependencies via a regular expression.
     !!}
     private
     type(regEx) :: regEx_
  end type dependencyRegEx

  interface dependencyExact
     !!{
     Constructor for exactly matching event hook dependencies.
     !!}
     module procedure dependencyExactConstructor
  end interface dependencyExact

  interface dependencyRegEx
     !!{
     Constructor for matching event hook dependencies via a regular expression.
     !!}
     module procedure dependencyRegExConstructor
  end interface dependencyRegEx

  !![
  <enumeration>
   <name>dependencyDirection</name>
   <description>Used to specify direction of event hook dependencies.</description>
   <visibility>public</visibility>
   <entry label="before"/>
   <entry label="after" />
  </enumeration>
  !!]
  
  type :: hook
     !!{
     Base class for individual hooked function calls. Stores the object to be passed as the first argument to the function.
     !!}
     class    (*             ), pointer                   :: object_             => null()
     class    (hook          ), pointer                   :: next                => null()
     integer                                              :: openMPThreadBinding          , openMPLevel
     integer                  , dimension(:), allocatable :: openMPThread
     character(len=128       )                            :: label
     class    (dependency    ), dimension(:), allocatable :: dependencies
  end type hook

  type, extends(hook) :: hookUnspecified
     !!{
     Class for hooked function calls with unspecified interfaces.
     !!}
     procedure(), pointer, nopass :: function_ => null()
  end type hookUnspecified

  type :: hookList
     !!{
     List of pointers to hooks.
     !!}
     class(hook), pointer :: hook_
  end type hookList
  
  type :: eventHook
     !!{
     Class used to define a set of hooked function calls for a given event.
     !!}
     private
     integer                               :: count_       =  0
     !$ type   (ompReadWriteLock)          :: lock_
     !$ logical                            :: initialized_ =  .false.
     class     (hook            ), pointer :: first_       => null()
#ifdef OMPPROFILE
     !$ double precision                   :: waitTimeRead = 0.0d0   , waitTimeWrite=0.0d0
#endif
   contains
     !![
     <methods>
       <method description="Return a count of the number of hooks into this event." method="count" />
       <method description="Return a pointer to the first hook into this event." method="first" />
       <method description="Initialize the event." method="initialize" />
       <method description="Lock the event." method="lock" />
       <method description="Unlock the event." method="unlock" />
       <method description="Reorder hooked functions to resolved any dependencies." method="resolveDependencies" />
     </methods>
     !!]
     procedure :: count               => eventHookCount
     procedure :: first               => eventHookFirst
     procedure :: initialize          => eventHookInitialize
     procedure :: lock                => eventHookLock
     procedure :: unlock              => eventHookUnlock
     procedure :: resolveDependencies => eventHookResolveDependencies
  end type eventHook

  type, extends(eventHook) :: eventHookUnspecified
     !!{
     Class used to define a set of hooked function calls for a given event.
     !!}
     private
   contains
     !![
     <methods>
       <method description="Attach a hook to the event." method="attach" />
       <method description="Return true if the object is attached to this event." method="isAttached" />
       <method description="Detach a hook from the event." method="detach" />
     </methods>
     !!]
     procedure :: attach     => eventHookUnspecifiedAttach
     procedure :: isAttached => eventHookUnspecifiedIsAttached
     procedure :: detach     => eventHookUnspecifiedDetach
  end type eventHookUnspecified

  ! The following is an enumeration of ways in which hooked functions can be bound to OpenMP threads. The behavior is as follows:
  !
  !      none: The hooked function is not bound to any OpenMP thread - it will always be called whenever the event is triggered,
  !            irrespective of which OpenMP thread triggered the event.
  !   atLevel: The hooked function is only called if: a) the current OpenMP level matches the level at which the function was
  !            hooked, and; b) the OpenMP thread number of the event triggering thread matches that of the thread which hooked the
  !            function at all OpenMP levels.
  ! allLevels: The hooked function is only called if: a) the current OpenMP level is equal to or greater than the level at which
  !            the function was hooked, and; b) the OpenMP thread number of the event triggering thread matches that of the thread
  !            which hooked the function at all OpenMP levels, where the thread number of the hooked function at OpenMP levels
  !            above that at which was hooked is taken to be equal to the thread number of the highest OpenMP level when the
  !            function was hooked.
  !
  ! Ideally there would be no need for the "allLevels" option, but currently some default functionClass objects are used - for
  ! thread-0 the default at levels greater than 0 is the exact same object as that at level-0 - in many cases
  ! (e.g. calculationResets) we do want to call this function in the case of the event being triggered even though it exists at a
  ! lower level. When default functionClass objects are removed this option should be deprecated.
  !![
  <enumeration>
   <name>openMPThreadBinding</name>
   <description>Used to specify how hooked functions are bound to OpenMP threads.</description>
   <visibility>public</visibility>
   <entry label="none"     />
   <entry label="atLevel"  />
   <entry label="allLevels"/>
  </enumeration>
  !!]

  !![
  <eventHookManager/>
  !!]

contains

  subroutine eventHookInitialize(self)
    !!{
    Initialize the OpenMP lock in an event object.
    !!}
    class(eventHook), intent(inout) :: self

    !$ if (.not.self%initialized_) then
    !$   self%lock_=ompReadWriteLock()
    !$   self%initialized_=.true.
    !$ end if
    return
  end subroutine eventHookInitialize

  subroutine eventHookLock(self,writeLock)
    !!{
    Lock the event to avoid race conditions between OpenMP threads.
    !!}
#ifdef OMPPROFILE
    !$ use :: OMP_Lib, only : OMP_Get_WTime
#endif
    implicit none
#ifdef OMPPROFILE
    double precision                            :: ompProfileTimeWaitStart, ompProfileTimeWaitEnd
#endif
    class  (eventHook), intent(inout)           :: self
    logical           , intent(in   ), optional :: writeLock
    !![
    <optionalArgument name="writeLock" defaultsTo=".true."/>
    !!]

    if (writeLock_) then
#ifdef OMPPROFILE
       !$ ompProfileTimeWaitStart=OMP_Get_WTime()
#endif
       !$ call self%lock_%setWrite(haveReadLock=.false.)
#ifdef OMPPROFILE
       !$ ompProfileTimeWaitEnd=OMP_Get_WTime()
       !$ ompProfileTimeWaitEnd=ompProfileTimeWaitEnd-ompProfileTimeWaitStart
       !$omp atomic
       !$ self%waitTimeWrite=self%waitTimeWrite+ompProfileTimeWaitEnd
#endif
    else
#ifdef OMPPROFILE
       !$ ompProfileTimeWaitStart=OMP_Get_WTime()
#endif
       !$ call self%lock_%setRead (                    )
#ifdef OMPPROFILE
       !$ ompProfileTimeWaitEnd=OMP_Get_WTime()
       !$ ompProfileTimeWaitEnd=ompProfileTimeWaitEnd-ompProfileTimeWaitStart
       !$omp atomic
       !$ self%waitTimeRead =self%waitTimeRead +ompProfileTimeWaitEnd
#endif
    end if
    return
  end subroutine eventHookLock

  subroutine eventHookUnlock(self,writeLock)
    !!{
    Unlock the event to avoid race conditions between OpenMP threads.
    !!}
    implicit none
    class  (eventHook), intent(inout)           :: self
    logical           , intent(in   ), optional :: writeLock
    !![
    <optionalArgument name="writeLock" defaultsTo=".true."/>
    !!]

    if (writeLock_) then
       !$ call self%lock_%unsetWrite(haveReadLock=.false.)
    else
       !$ call self%lock_%unsetRead (                    )
    end if
    return
  end subroutine eventHookUnlock
  
  subroutine eventHookUnspecifiedAttach(self,object_,function_,openMPThreadBinding,label,dependencies)
    !!{
    Attach an object to an event hook.
    !!}
    use    :: Galacticus_Error, only : Galacticus_Error_Report
    !$ use :: OMP_Lib         , only : OMP_Get_Ancestor_Thread_Num, OMP_Get_Level
    implicit none
    class     (eventHookUnspecified), intent(inout)                         :: self
    class     (*                   ), intent(in   ), target                 :: object_
    integer                         , intent(in   ), optional               :: openMPThreadBinding
    character (len=*               ), intent(in   ), optional               :: label
    class     (dependency          ), intent(in   ), optional, dimension(:) :: dependencies
    procedure (                    )                                        :: function_
    class     (hook                )                         , pointer      :: hook_
    !$ integer                                                              :: i
    !![
    <optionalArgument name="openMPThreadBinding" defaultsTo="openMPThreadBindingNone" />
    !!]

    ! Lock the object.
    !$ if (.not.self%initialized_) call Galacticus_Error_Report('event has not been initialized'//{introspection:location})
    call self%lock()
    ! Allocate the next entry in our list of hooks.
    if (associated(self%first_)) then
       hook_ => self%first_
       do while (associated(hook_%next))
          hook_ => hook_%next
       end do
       allocate(hookUnspecified :: hook_%next )
       hook_ => hook_%next
    else
       allocate(hookUnspecified :: self%first_)
       hook_ => self%first_
    end if
    ! Create the new hook.
    select type (hook_)
    type is (hookUnspecified)
       hook_%object_             => object_
       hook_%function_           => function_
       hook_%openMPThreadBinding =  openMPThreadBinding_
       if (present(label)) then
          hook_%label=label
       else
          hook_%label=""
       end if
       !$ if (hook_%openMPThreadBinding == openMPThreadBindingAtLevel .or. hook_%openMPThreadBinding == openMPThreadBindingAllLevels) then
       !$    hook_%openMPLevel=OMP_Get_Level()
       !$    allocate(hook_%openMPThread(0:hook_%openMPLevel))
       !$    do i=0,hook_%openMPLevel
       !$       hook_%openMPThread(i)=OMP_Get_Ancestor_Thread_Num(i)
       !$    end do
       !$ end if
    end select
    ! Increment the count of hooks into this event and resolve dependencies.
    self%count_=self%count_+1
    call self%resolveDependencies(hook_,dependencies)
    call self%unlock             (                  )
    return
  end subroutine eventHookUnspecifiedAttach

  subroutine eventHookResolveDependencies(self,hookNew,dependencies)
    !!{
    Resolve dependencies between hooked function calls.
    !!}
    use :: Galacticus_Error   , only : Galacticus_Error_Report, errorStatusSuccess
    use :: Sorting_Topological, only : Sort_Topological
    implicit none
    class    (eventHook ), intent(inout)                           :: self
    class    (hook      ), intent(inout)                           :: hookNew
    class    (dependency), intent(in   ), dimension(:  ), optional :: dependencies
    integer              , allocatable  , dimension(:  )           :: order
    integer              , allocatable  , dimension(:,:)           :: dependentIndices, dependentIndicesTmp
    type     (hookList  ), allocatable  , dimension(:  )           :: hooksUnordered  , hooksOrdered
    class    (hook      )               , pointer                  :: hook_           , hook__
    integer                                                        :: i               , j                  , &
         &                                                            k               , dependencyCount    , &
         &                                                            countOrdered    , status             , &
         &                                                            l
    logical                                                        :: matches
    
    ! Add dependencies to the new hooked function.
    if (present(dependencies)) then
       allocate(hookNew%dependencies(size(dependencies)),mold=dependencies)
       hookNew%dependencies=dependencies
       do i=1,size(dependencies)
          select type (dependency_ => hookNew%dependencies(i))
          type is (dependencyRegEx)
             dependency_%regEx_=regEx(dependency_%label)
          end select
       end do
    end if
    ! Build the dependency array.
    allocate(dependentIndices(1,2))
    hook_           => self%first_
    i               =  0
    dependencyCount =  0
    do while (associated(hook_))
       i=i+1
       if (allocated(hook_%dependencies)) then
          do k=1,size(hook_%dependencies)
             j    = 0
             hook__ => self%first_
             do while (associated(hook__))
                j      =j      +1
                matches=.false.
                select type (dependency_ => hook_%dependencies(k))
                type is (dependencyExact)
                   ! Exact match dependency.
                   matches=dependency_%label  ==      hook__%label
                type is (dependencyRegEx )
                   ! Regular expression dependency.
                   matches=dependency_%regEx_%matches(hook__%label)
                class default
                   call Galacticus_Error_Report('unknown dependency'//{introspection:location})
                end select
                if (matches) then                   
                   dependencyCount=dependencyCount+1
                   if (dependencyCount > size(dependentIndices,dim=1)) then
                      call move_alloc(dependentIndices,dependentIndicesTmp)
                      allocate(dependentIndices(2*size(dependentIndicesTmp,dim=1),2))
                      do l=1,size(dependentIndicesTmp,dim=1)
                         dependentIndices(l,:)=dependentIndicesTmp(l,:)
                      end do
                      deallocate(dependentIndicesTmp)
                   end if
                   select case (hook_%dependencies(k)%direction)
                   case (dependencyDirectionBefore)
                      dependentIndices(dependencyCount,:)=[j,i]
                   case (dependencyDirectionAfter)
                      dependentIndices(dependencyCount,:)=[i,j]
                   case default
                      call Galacticus_Error_Report('unknown dependency direction'//{introspection:location})
                   end select
                end if
                hook__ => hook__%next
             end do
           end do
       end if
       hook_ => hook_%next
    end do
    ! Generate an ordering which satisfies all dependencies.
    allocate(order(self%count_))
    call Sort_Topological(self%count_,dependencyCount,dependentIndices(1:dependencyCount,:),order,countOrdered,status)
    if (status /= errorStatusSuccess) call Galacticus_Error_Report('unable to resolve hooked function dependencies'//{introspection:location})
    ! Build an array of pointers to our hooks with this ordering.
    allocate(hooksUnordered(self%count_))
    allocate(hooksOrdered  (self%count_))
    hook_ => self%first_
    i     =  0
    do while (associated(hook_))
       i=i+1
       hooksUnordered(i)%hook_ => hook_
       hook_ => hook_%next
    end do
    do i=1,self%count_
       hooksOrdered(i)%hook_ => hooksUnordered(order(i))%hook_
    end do
    self%first_ => hooksOrdered(1)%hook_
    do i=1,self%count_-1
       hooksOrdered(i)%hook_%next => hooksOrdered(i+1)%hook_
    end do
    hooksOrdered(self%count_)%hook_%next => null()
    ! Clean up.
    deallocate(dependentIndices)
    deallocate(order           )
    deallocate(hooksUnordered  )
    deallocate(hooksOrdered    )
    nullify   (hook_           )
    nullify   (hook__          )
    return
  end subroutine eventHookResolveDependencies
  
  logical function eventHookUnspecifiedIsAttached(self,object_,function_)
    !!{
    Return true if an object is attached to an event hook.
    !!}
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class    (eventHookUnspecified), intent(inout)          :: self
    class    (*                   ), intent(in   ), target  :: object_
    procedure(                    )                         :: function_
    class    (hook                )               , pointer :: hook_

    ! Lock the object.
    !$ if (.not.self%initialized_) call Galacticus_Error_Report('event has not been initialized'//{introspection:location})
    call self%lock(writeLock=.false.)
    if (associated(self%first_)) then
       hook_ => self%first_
       do while (associated(hook_))
          select type (hook_)
          type is (hookUnspecified)
             if (associated(hook_%object_,object_).and.associated(hook_%function_,function_)) then
                eventHookUnspecifiedIsAttached=.true.
                call self%unlock(writeLock=.false.)
                return
             end if
          end select
          hook_ => hook_%next
       end do
    end if
    eventHookUnspecifiedIsAttached=.false.
    call self%unlock(writeLock=.false.)
    return
  end function eventHookUnspecifiedIsAttached

  subroutine eventHookUnspecifiedDetach(self,object_,function_)
    !!{
    Attach an object to an event hook.
    !!}
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class    (eventHookUnspecified), intent(inout)          :: self
    class    (*                   ), intent(in   ), target  :: object_
    procedure(                    )                         :: function_
    class    (hook                )               , pointer :: hook_    , hookPrevious_

    ! Lock the object.
    !$ if (.not.self%initialized_) call Galacticus_Error_Report('event has not been initialized'//{introspection:location})
    call self%lock()
    if (associated(self%first_)) then
       hookPrevious_ => null()
       hook_         => self%first_
       do while (associated(hook_))
          select type (hook_)
          type is (hookUnspecified)
             if (associated(hook_%object_,object_).and.associated(hook_%function_,function_)) then
                self%count_=self%count_-1
                if (associated(hookPrevious_)) then
                   hookPrevious_%next   => hook_%next
                else
                   self         %first_ => hook_%next
                end if
                deallocate(hook_)
                call self%unlock()
                return
             end if
          end select
          hookPrevious_ => hook_
          hook_         => hook_%next
       end do
    end if
    call Galacticus_Error_Report('object/function not attached to this event'//{introspection:location})
    call self%unlock()
    return
  end subroutine eventHookUnspecifiedDetach

  integer function eventHookCount(self)
    !!{
    Return a count of the number of hooks into this event.
    !!}
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class(eventHook), intent(inout):: self

    !$ if (.not.self%initialized_) call Galacticus_Error_Report('event has not been initialized'//{introspection:location})
    call self%lock(writeLock=.false.)
    eventHookCount=self%count_
    call self%unlock(writeLock=.false.)
    return
  end function eventHookCount

  function eventHookFirst(self)
    !!{
    Return a pointer to the first hook into this event.
    !!}
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class(hook     ), pointer       :: eventHookFirst
    class(eventHook), intent(inout) :: self

    !$ if (.not.self%initialized_) call Galacticus_Error_Report('event has not been initialized'//{introspection:location})
    eventHookFirst => self%first_
    return
  end function eventHookFirst

  function dependencyExactConstructor(direction,label) result(self)
    !!{
    Constructor for an exact dependency.
    !!}
    implicit none
    type     (dependencyExact) :: self
    integer                    :: direction
    character(len=*          ) :: label
    !![
    <constructorAssign variables="direction, label"/>
    !!]

    return
  end function dependencyExactConstructor
  
  function dependencyRegExConstructor(direction,label) result(self)
    !!{
    Constructor for a regular expression dependency.
    !!}
    implicit none
    type     (dependencyRegEx) :: self
    integer                    :: direction
    character(len=*          ) :: label
    !![
    <constructorAssign variables="direction, label"/>
    !!]

    self%regEx_=regEx(label)
    return
  end function dependencyRegExConstructor

  !![
  <hdfPreCloseTask>
   <unitName>eventsHooksWaitTimes</unitName>
  </hdfPreCloseTask>
  !!]

end module Events_Hooks
