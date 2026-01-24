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
Contains module which handles hooking of objects into events.
!!}

module Events_Hooks
  !!{
  Handles hooking of object function class into events.
  !!}
  use :: Regular_Expressions, only : regEx
  use :: Locks              , only : ompLock, ompReadWriteLock
  private
  public :: hook                 , hookUnspecified        , dependencyExact              , dependencyRegEx, &
       &    eventsHooksInitialize, eventsHooksFutureThread, eventsHooksAtLevelToAllLevels

  !![
  <enumeration>
   <name>dependencyDirection</name>
   <description>Used to specify direction of event hook dependencies.</description>
   <visibility>public</visibility>
   <entry label="before"/>
   <entry label="after" />
  </enumeration>
  !!]
  
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

  type :: dependency
     !!{
     Base class for event hook dependencies.
     !!}
     private
     type     (enumerationDependencyDirectionType) :: direction
     character(len=128                           ) :: label
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
   contains
     procedure ::                  dependencyRegExAssign
     generic   :: assignment(=) => dependencyRegExAssign
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

  type :: hook
     !!{
     Base class for individual hooked function calls. Stores the object to be passed as the first argument to the function.
     !!}
     class    (*                                 ), pointer                   :: object_             => null()
     type     (enumerationOpenMPThreadBindingType)                            :: openMPThreadBinding
     integer                                                                  :: openMPLevel                  , eventID
     integer                                      , dimension(:), allocatable :: openMPThread
     character(len=128                           )                            :: label
     class    (dependency                        ), dimension(:), allocatable :: dependencies
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
     class(hook), pointer :: hook_ => null()
  end type hookList
  
  type :: eventHook
     !!{
     Class used to define a set of hooked function calls for a given event.
     !!}
     private
     logical                                                                  :: isGlobal    =.false.
     !$ type            (ompReadWriteLock)                                    :: lock_
#ifdef OMPPROFILE
     !$ double precision                                                      :: waitTimeRead=0.0d0  , waitTimeWrite=0.0d0
#endif
     character          (len=128         )                                    :: label
     integer                                                                  :: count_      =0
     type               (hookList        ), allocatable, dimension(:), public :: hooks_
   contains
     !![
     <methods>
       <method description="Return a count of the number of hooks into this event."  method="count"              />
       <method description="Reorder hooked functions to resolved any dependencies."  method="resolveDependencies"/>
       <method description="Filter events to match the current OpenMP thread/level." method="filter"             />
       <method description="Lock the event."                                         method="lock"               />
       <method description="Unlock the event."                                       method="unlock"             />
     </methods>
     !!]
     procedure :: count               => eventHookCount
     procedure :: resolveDependencies => eventHookResolveDependencies
     procedure :: filter              => eventHookFilter
     procedure :: lock                => eventHookLock
     procedure :: unlock              => eventHookUnlock
  end type eventHook

  type :: eventHookList
     !!{
     List of event hooks.
     !!}
     class(eventHook    ), allocatable :: eventHook_
     class(eventHookList), pointer     :: next       => null()
  end type eventHookList
  
  type, extends(eventHook) :: eventHookUnspecified
     !!{
     Class used to define a set of hooked function calls for a given event.
     !!}
     private
   contains
     !![
     <methods>
       <method description="Attach a hook to the event."                          method="attach"    />
       <method description="Return true if the object is attached to this event." method="isAttached"/>
       <method description="Detach a hook from the event."                        method="detach"    />
     </methods>
     !!]
     procedure :: attach     => eventHookUnspecifiedAttach
     procedure :: isAttached => eventHookUnspecifiedIsAttached
     procedure :: detach     => eventHookUnspecifiedDetach
  end type eventHookUnspecified

  ! Lock used to guard shared memory used for copyin/out operations.
  type   (ompLock) :: copyLock

  ! Globally-unique ID for events.
  integer          :: eventID           = 0

  ! Future-thread number, used for setting up events that attach to yet-to-be-created threads.
  integer          :: futureThread_      =-1
  !$omp threadprivate(futureThread_)

  ! State controlling whether "atLevel" attachment should be promoted to "allLevels" attachment. This is needed for objects
  ! created in the master thread. Once all objects are deepCopied from the master thread this should no longer be needed.
  logical          :: atLevelToAllLevels_=.false.
  !$omp threadprivate(atLevelToAllLevels_)
  
  !![
  <eventHookManager/>
  !!]

contains

  subroutine dependencyRegExAssign(self,from)
    !!{
    Perform assignment of reg-ex dependencies.
    !!}
    implicit none
    class(dependencyRegEx), intent(inout) :: self
    class(dependencyRegEx), intent(in   ) :: from

    self%direction=from%direction
    self%label    =from%label
    self%regEx_   =from%regEx_
    return
  end subroutine dependencyRegExAssign

  subroutine eventsHooksFutureThread(futureThread)
    !!{
    Set the future thread to which events will attach. If no argument is given, future threads are disabled.
    !!}
    implicit none
    integer, intent(in   ), optional :: futureThread

    if (present(futureThread)) then
       futureThread_=futureThread
    else
       futureThread_=-1
    end if
    return
  end subroutine eventsHooksFutureThread

  subroutine eventsHooksAtLevelToAllLevels(atLevelToAllLevels)
    !!{
    Set the promotion state for events attaching ``at-level''.
    !!}
    implicit none
    logical, intent(in   ) :: atLevelToAllLevels

    atLevelToAllLevels_=atLevelToAllLevels
    return
  end subroutine eventsHooksAtLevelToAllLevels

  subroutine eventHookUnspecifiedAttach(self,object_,function_,openMPThreadBinding,label,dependencies)
    !!{
    Attach an object to an event hook.
    !!}
    use    :: Display           , only : displayMessage             , verbosityLevelInfo
    use    :: Error             , only : Error_Report
    use    :: ISO_Varying_String, only : varying_string             , var_str           , assignment(=), operator(//)
    use    :: String_Handling   , only : operator(//)
    !$ use :: OMP_Lib           , only : OMP_Get_Ancestor_Thread_Num, OMP_Get_Level
    implicit none
    class     (eventHookUnspecified              ), intent(inout)                            :: self
    class     (*                                 ), intent(in   ), target                    :: object_
    type      (enumerationOpenMPThreadBindingType), intent(in   ), optional                  :: openMPThreadBinding
    character (len=*                             ), intent(in   ), optional                  :: label
    class     (dependency                        ), intent(in   ), optional   , dimension(:) :: dependencies
    procedure (                                  )                                           :: function_
    type      (hookList                          )               , allocatable, dimension(:) :: hooksTmp
    type      (hookUnspecified                   )                            , pointer      :: hook_
    type      (varying_string                    )                                           :: threadLabel
    integer                                                                                  :: ompLevelEffective
    !$ integer                                                                               :: i
    !![
    <optionalArgument name="openMPThreadBinding" defaultsTo="openMPThreadBindingNone" />
    !!]

    ! Validate the thread binding model.
    if (self%isGlobal) then
       if (openMPThreadBinding_ /= openMPThreadBindingNone) call Error_Report("global event hooks permit only 'openMPThreadBindingNone'"         //{introspection:location})
    else
       if (openMPThreadBinding_ == openMPThreadBindingNone) call Error_Report("threadprivate event hooks do not permit 'openMPThreadBindingNone'"//{introspection:location})
    end if
    ! Check if atLevel attachment should be promoted.
    if (atLevelToAllLevels_ .and. openMPThreadBinding_ == openMPThreadBindingAtLevel) &
         openMPThreadBinding_=openMPThreadBindingAllLevels
    ! Lock the object.
    !$ if (self%isGlobal) call self%lock()
    ! Resize the array of hooks.
    if (allocated(self%hooks_)) then
       call move_alloc(self%hooks_,hooksTmp)
       allocate(self%hooks_(self%count_+1))
       self%hooks_(1:self%count_)=hooksTmp
       deallocate(hooksTmp)
    else
       allocate(self%hooks_(1))
    end if
    ! Create the new hook.
    allocate(hook_)
    hook_%object_             => object_
    hook_%function_           => function_
    hook_%openMPThreadBinding =  openMPThreadBinding_
    if (present(label)) then
       hook_%label=label
    else
       hook_%label=""
    end if
    !$omp atomic
    eventID            =eventID+1
    hook_      %eventID=eventID
    threadLabel        =""
    !$ threadLabel=" from thread "
    !$ ompLevelEffective=OMP_Get_Level()
    !$ if (futureThread_ /= -1) ompLevelEffective=ompLevelEffective+1
    !$ hook_%openMPLevel=ompLevelEffective
    !$ allocate(hook_%openMPThread(0:hook_%openMPLevel))
    !$ do i=0,hook_%openMPLevel
    !$    if (i == hook_%openMPLevel .and. futureThread_ /= -1) then
    !$      hook_%openMPThread(i)=futureThread_
    !$    else
    !$      hook_%openMPThread(i)=OMP_Get_Ancestor_Thread_Num(i)
    !$    end if
    !$    if (i > 0) threadLabel=threadLabel//" -> "
    !$    threadLabel=threadLabel//hook_%openMPThread(i)
    !$ end do
    ! Insert the hook into the list.
    self%hooks_(self%count_+1)%hook_ => hook_
    ! Increment the count of hooks into this event and resolve dependencies.
    self%count_=self%count_+1
    call self%resolveDependencies(hook_,dependencies)
    ! Report
    call displayMessage(var_str("attaching '")//trim(hook_%label)//"' ["//hook_%eventID//"] to event "//trim(self%label)//threadLabel//" [count="//self%count_//"]",verbosityLevelInfo)
    !$ if (self%isGlobal) call self%unlock()
    return
  end subroutine eventHookUnspecifiedAttach

  subroutine eventHookResolveDependencies(self,hookNew,dependencies)
    !!{
    Resolve dependencies between hooked function calls.
    !!}
    use :: Error              , only : Error_Report    , errorStatusSuccess
    use :: Sorting_Topological, only : Sort_Topological
    implicit none
    class    (eventHook ), intent(inout)                           :: self
    class    (hook      ), intent(inout)                           :: hookNew
    class    (dependency), intent(in   ), dimension(:  ), optional :: dependencies
    integer              , allocatable  , dimension(:  )           :: order
    integer              , allocatable  , dimension(:,:)           :: dependentIndices, dependentIndicesTmp
    type     (hookList  ), allocatable  , dimension(:  )           :: hooksOrdered
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
    dependencyCount =  0
    do i=1,self%count_
       if (allocated(self%hooks_(i)%hook_%dependencies)) then
          do k=1,size(self%hooks_(i)%hook_%dependencies)
             j    = 0
             do j=1,self%count_
                matches=.false.
                select type (dependency_ => self%hooks_(i)%hook_%dependencies(k))
                type is (dependencyExact)
                   ! Exact match dependency.
                   matches=dependency_%label  ==      self%hooks_(j)%hook_%label
                type is (dependencyRegEx )
                   ! Regular expression dependency.
                   matches=dependency_%regEx_%matches(self%hooks_(j)%hook_%label)
                class default
                   call Error_Report('unknown dependency'//{introspection:location})
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
                   select case (self%hooks_(i)%hook_%dependencies(k)%direction%ID)
                   case (dependencyDirectionBefore%ID)
                      dependentIndices(dependencyCount,:)=[j,i]
                   case (dependencyDirectionAfter %ID)
                      dependentIndices(dependencyCount,:)=[i,j]
                   case default
                      call Error_Report('unknown dependency direction'//{introspection:location})
                   end select
                end if
             end do
           end do
       end if
    end do
    ! If there are dependencies present then generate an ordering which satisfies all dependencies.
    if (dependencyCount > 0) then
       allocate(order(self%count_))
       call Sort_Topological(self%count_,dependencyCount,dependentIndices(1:dependencyCount,:),order,countOrdered,status)
       if (status /= errorStatusSuccess) call Error_Report('unable to resolve hooked function dependencies'//{introspection:location})
       ! Build an array of pointers to our hooks with this ordering.
       allocate(hooksOrdered(self%count_))
       do i=1,self%count_
          hooksOrdered(i)%hook_ => self%hooks_(order(i))%hook_
       end do
       deallocate(self%hooks_)
       call move_alloc(hooksOrdered,self%hooks_)
    end if
    return
  end subroutine eventHookResolveDependencies
  
  logical function eventHookUnspecifiedIsAttached(self,object_,function_)
    !!{
    Return true if an object is attached to an event hook.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class    (eventHookUnspecified), intent(inout)          :: self
    class    (*                   ), intent(in   ), target  :: object_
    procedure(                    )                         :: function_
    integer                                                 :: i
    
    !$ if (self%isGlobal) call self%lock(writeLock=.false.)
    if (allocated(self%hooks_)) then
       do i=1,self%count_
          select type (hook_ => self%hooks_(i)%hook_)
          type is (hookUnspecified)
             if (associated(hook_%object_,object_).and.associated(hook_%function_,function_)) then
                eventHookUnspecifiedIsAttached=.true.
                !$ if (self%isGlobal) call self%unlock(writeLock=.false.)
                return
             end if
          end select
       end do
    end if
    eventHookUnspecifiedIsAttached=.false.
    !$ if (self%isGlobal) call self%unlock(writeLock=.false.)
    return
  end function eventHookUnspecifiedIsAttached

  subroutine eventHookUnspecifiedDetach(self,object_,function_)
    !!{
    Attach an object to an event hook.
    !!}
    use    :: Display           , only : displayMessage             , verbosityLevelInfo
    use    :: Error             , only : Error_Report
    use    :: ISO_Varying_String, only : varying_string             , var_str           , assignment(=), operator(//)
    use    :: String_Handling   , only : operator(//)
    !$ use :: OMP_Lib           , only : OMP_Get_Ancestor_Thread_Num, OMP_Get_Level
    implicit none
    class    (eventHookUnspecified), intent(inout)               :: self
    class    (*                   ), intent(in   ), target       :: object_
    procedure(                    )                              :: function_
    type     (hookList            ), allocatable  , dimension(:) :: hooksTmp
    type     (varying_string      )                              :: threadLabel
    integer                                                      :: i          , j
    
    !$ if (self%isGlobal) call self%lock()
    if (allocated(self%hooks_)) then
       do i=1,self%count_
          select type (hook_ => self%hooks_(i)%hook_)
          type is (hookUnspecified)
             if (associated(hook_%object_,object_).and.associated(hook_%function_,function_)) then
                ! Report
                threadLabel   =""
                !$ threadLabel=" from thread "
                !$ do j=0,OMP_Get_Level()
                !$    if (j > 0) threadLabel=threadLabel//" -> "
                !$    threadLabel=threadLabel//OMP_Get_Ancestor_Thread_Num(j)
                !$ end do
                call displayMessage(var_str("detaching '")//trim(self%hooks_(i)%hook_%label)//"' ["//self%hooks_(i)%hook_%eventID//"] from event"//trim(self%label)//threadLabel//" [count="//self%count_//"]",verbosityLevelInfo)
                deallocate(self%hooks_(i)%hook_)
                if (self%count_ > 1) then
                   call move_alloc(self%hooks_,hooksTmp)
                   allocate(self%hooks_(self%count_-1))
                   if (i >           1) self%hooks_(1:          i-1)=hooksTmp(1  :          i-1)
                   if (i < self%count_) self%hooks_(i:self%count_-1)=hooksTmp(i+1:self%count_  )
                   deallocate(hooksTmp)
                else
                   deallocate(self%hooks_)
                end if
                self%count_=self%count_-1
               !$ if (self%isGlobal) call self%unlock()
                return
             end if
          end select
       end do
    end if
    call Error_Report('object/function not attached to this event'//{introspection:location})
    !$ if (self%isGlobal) call self%unlock()
    return
  end subroutine eventHookUnspecifiedDetach

  integer function eventHookCount(self)
    !!{
    Return a count of the number of hooks into this event.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(eventHook), intent(inout):: self

    !$ if (self%isGlobal) call self%  lock(writeLock=.false.)
    eventHookCount=self%count_
    !$ if (self%isGlobal) call self%unlock(writeLock=.false.)
    return
  end function eventHookCount

  subroutine eventHookFilter(self)
    !!{
    Filter hooked functions for the current OpenMP thread/level.
    !!}
    !$ use :: Error  , only : Error_Report
    !$ use :: OMP_Lib, only : OMP_Get_Ancestor_Thread_Num, OMP_Get_Level
    implicit none
    class  (eventHook), intent(inout)               :: self
    logical           , allocatable  , dimension(:) :: functionActive_
    integer           , allocatable  , dimension(:) :: ompAncestorThreadNum_
    type   (hookList ), allocatable  , dimension(:) :: hooksTmp
    integer                                         :: i                    , j               , &
         &                                             ompLevel_            , ompLevelCurrent_
    
    ! Nothing to do if no hooks are attached.
    if (.not.allocated(self%hooks_)) return
    !$ ! Get the current OMP level and ancestor thread numbers.
    !$ ompLevelCurrent_=OMP_Get_Level()
    !$ allocate(ompAncestorThreadNum_(0:ompLevelCurrent_))
    !$ do ompLevel_=0,ompLevelCurrent_
    !$    ompAncestorThreadNum_(ompLevel_)=OMP_Get_Ancestor_Thread_Num(ompLevel_)
    !$ end do
    !$ allocate(functionActive_(self%count_))
    !$ ! Examine each hooked function to see if it is active in this OpenMP thread.
    !$ do i=1,self%count_
    !$    select case (self%hooks_(i)%hook_%openMPThreadBinding%ID)
    !$    case (openMPThreadBindingNone%ID)
    !$       ! Not bound to any OpenMP thread, so always call.
    !$       functionActive_(i)=.true.
    !$    case (openMPThreadBindingAtLevel%ID)
    !$       ! Binds at the OpenMP level - check levels match, and that this hooked object matches the OpenMP thread number across all levels.
    !$       if (self%hooks_(i)%hook_%openMPLevel == ompLevelCurrent_) then
    !$          functionActive_(i)=.true.
    !$          do ompLevel_=self%hooks_(i)%hook_%openMPLevel,0,-1
    !$             if (self%hooks_(i)%hook_%openMPThread(ompLevel_) /= ompAncestorThreadNum_(ompLevel_)) then
    !$                functionActive_(i)=.false.
    !$                exit
    !$             end if
    !$          end do
    !$       else
    !$          functionActive_(i)=.false.
    !$       end if
    !$    case (openMPThreadBindingAllLevels%ID)
    !$       ! Binds at all levels at or above the level of the hooked object - check this condition is met, and that the hooked object matches the OpenMP thread number across all levels.
    !$       if (self%hooks_(i)%hook_%openMPLevel <= ompLevelCurrent_) then
    !$          functionActive_(i)=.true.
    !$          do ompLevel_=ompLevelCurrent_,0,-1
    !$             if (self%hooks_(i)%hook_%openMPThread(min(ompLevel_,self%hooks_(i)%hook_%openMPLevel)) /= ompAncestorThreadNum_(ompLevel_)) then
    !$                functionActive_(i)=.false.
    !$                exit
    !$             end if
    !$          end do
    !$       else
    !$          functionActive_(i)=.false.
    !$       end if
    !$    case default
    !$       functionActive_(i)=.false.
    !$       call Error_Report('unknown OpenMP binding'//{introspection:location})
    !$    end select
    !$ end do
    !$ ! Create a filtered list of hooks.
    !$ if (count(functionActive_) > 0) then
    !$    allocate(hooksTmp(count(functionActive_)))
    !$    j=0
    !$    do i=1,self%count_
    !$       if (functionActive_(i)) then
    !$          j=j+1
    !$          hooksTmp(j)=self%hooks_(i)
    !$       end if
    !$    end do
    !$    deallocate(self%hooks_)
    !$    call move_alloc(hooksTmp,self%hooks_)
    !$ else
    !$    deallocate(self%hooks_)
    !$ end if
    !$ self%count_=count(functionActive_)
    return
  end subroutine eventHookFilter
  
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
  
  function dependencyExactConstructor(direction,label) result(self)
    !!{
    Constructor for an exact dependency.
    !!}
    implicit none
    type     (dependencyExact                   ) :: self
    type     (enumerationDependencyDirectionType) :: direction
    character(len=*                             ) :: label
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
    type     (dependencyRegEx                   ) :: self
    type     (enumerationDependencyDirectionType) :: direction
    character(len=*                             ) :: label
    !![
    <constructorAssign variables="direction, label"/>
    !!]

    self%regEx_=regEx(label)
    return
  end function dependencyRegExConstructor

  !![
  <outputFileClose function="eventsHooksWaitTimes"/>
  !!]

end module Events_Hooks
