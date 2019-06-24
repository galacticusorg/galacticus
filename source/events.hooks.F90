!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

!% Contains module which handles hooking of objects into events.

module Events_Hooks
  !% Handles hooking of object function class into events.
  !$ use OMP_Lib
  private
  public :: hook, hookUnspecified

  type :: hook
     !% Base class for individual hooked function calls. Stores the object to be passed as the first argument to the function.
     class  (*   ), pointer                   :: object_             => null()
     class  (hook), pointer                   :: next                => null()
     integer                                  :: openMPThreadBinding          , openMPLevel
     integer      , dimension(:), allocatable :: openMPThread
  end type hook

  type, extends(hook) :: hookUnspecified
     !% Class for hooked function calls with unspecified interfaces.
     procedure(), pointer, nopass :: function_ => null()
  end type hookUnspecified
  
  type :: eventHook
     !% Class used to define a set of hooked function calls for a given event.
     private
     integer                            :: count_       =  0
     !$ integer(omp_lock_kind)          :: lock_
     !$ logical                         :: initialized_ =  .false.
     class     (hook         ), pointer :: first_       => null()
   contains
     !@ <objectMethods>
     !@   <object>eventHook</object>
     !@   <objectMethod>
     !@     <method>count</method>
     !@     <type>\intzero</type>
     !@     <arguments></arguments>
     !@     <description>Return a count of the number of hooks into this event.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>first</method>
     !@     <type></type>
     !@     <arguments>\textcolor{red}{\textless *type(hook)\textgreater}</arguments>
     !@     <description>Return a pointer to the first hook into this event.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>initialize</method>
     !@     <type></type>
     !@     <arguments></arguments>
     !@     <description>Initialize the event.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>lock</method>
     !@     <type>void</type>
     !@     <arguments></arguments>
     !@     <description>Lock the event.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>unlock</method>
     !@     <type>void</type>
     !@     <arguments></arguments>
     !@     <description>Unlock the event.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     final     ::               eventHookDestructor
     procedure :: count      => eventHookCount
     procedure :: first      => eventHookFirst
     procedure :: initialize => eventHookInitialize
     procedure :: lock       => eventHookLock
     procedure :: unlock     => eventHookUnlock
  end type eventHook
  
  type, extends(eventHook) :: eventHookUnspecified
     !% Class used to define a set of hooked function calls for a given event.
     private
   contains
     !@ <objectMethods>
     !@   <object>eventHookUnspecified</object>
     !@   <objectMethod>
     !@     <method>attach</method>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless class(*)\textgreater} *object\_\argin, \textcolor{red}{\textless procedure()\textgreater} *function\_\argin</arguments>
     !@     <description>Attach a hook to the event.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>detach</method>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless class(*)\textgreater} *object\_\argin, \textcolor{red}{\textless procedure()\textgreater} *function\_\argin</arguments>
     !@     <description>Detach a hook from the event.</description>
     !@   </objectMethod>
     !@ </objectMethods>
      procedure :: attach => eventHookUnspecifiedAttach
      procedure :: detach => eventHookUnspecifiedDetach
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
  !# <enumeration>
  !#  <name>openMPThreadBinding</name>
  !#  <description>Used to specify how hooked functions are bound to OpenMP threads.</description>
  !#  <visibility>public</visibility>
  !#  <entry label="none"     />
  !#  <entry label="atLevel"  />
  !#  <entry label="allLevels"/>
  !# </enumeration>

  !# <eventHookManager/>
  
contains

  subroutine eventHookInitialize(self)
    !% Initialize the OpenMP lock in an event object.
    class(eventHook), intent(inout) :: self

    !$ if (.not.self%initialized_) then
    !$   call OMP_Init_Lock(self%lock_)
    !$   self%initialized_=.true.
    !$ end if
    return
  end subroutine eventHookInitialize
  
  subroutine eventHookDestructor(self)
    !% Destructor for event hook class.
    type(eventHook), intent(inout) :: self

    !$ if (self%initialized_) call OMP_Destroy_Lock(self%lock_)
    return
  end subroutine eventHookDestructor

  subroutine eventHookLock(self)
    !% Lock the event to avoid race conditions between OpenMP threads.
    !$ use OMP_Lib
    implicit none
    class(eventHook), intent(inout) :: self

    !$ call OMP_Set_Lock(self%lock_)
    return
  end subroutine eventHookLock

  subroutine eventHookUnlock(self)
    !% Unlock the event to avoid race conditions between OpenMP threads.
    !$ use OMP_Lib
    implicit none
    class(eventHook), intent(inout) :: self

    !$ call OMP_Unset_Lock(self%lock_)
    return
  end subroutine eventHookUnlock
  
  subroutine eventHookUnspecifiedAttach(self,object_,function_,openMPThreadBinding)
    !% Attach an object to an event hook.
    !$ use OMP_Lib
    use Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class    (eventHookUnspecified), intent(inout)           :: self
    class    (*                   ), intent(in   ), target   :: object_
    integer                        , intent(in   ), optional :: openMPThreadBinding
    procedure(                    )                          :: function_
    class    (hook                )               , pointer  :: hook_
    integer                                                  :: i
    !# <optionalArgument name="openMPThreadBinding" defaultsTo="openMPThreadBindingNone" />

    ! Lock the object.
    !$ if (.not.self%initialized_) call Galacticus_Error_Report('event has not been initialized'//{introspection:location})
    !$ call OMP_Set_Lock(self%lock_)
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
       if (hook_%openMPThreadBinding == openMPThreadBindingAtLevel .or. hook_%openMPThreadBinding == openMPThreadBindingAllLevels) then
          hook_%openMPLevel=OMP_Get_Level()
          allocate(hook_%openMPThread(0:hook_%openMPLevel))
          do i=0,hook_%openMPLevel
             hook_%openMPThread(i)=OMP_Get_Ancestor_Thread_Num(i)
          end do
       end if
    end select
    ! Increment the count of hooks into this event.
    self%count_=self%count_+1
    !$ call OMP_Unset_Lock(self%lock_)
    return
  end subroutine eventHookUnspecifiedAttach

  subroutine eventHookUnspecifiedDetach(self,object_,function_)
    !% Attach an object to an event hook.
    use Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class    (eventHookUnspecified), intent(inout)          :: self
    class    (*                   ), intent(in   ), target  :: object_
    procedure(                    )                         :: function_
    class    (hook                )               , pointer :: hook_    , hookPrevious_
    
    ! Lock the object.
    !$ if (.not.self%initialized_) call Galacticus_Error_Report('event has not been initialized'//{introspection:location})
    !$ call OMP_Set_Lock(self%lock_)
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
                !$ call OMP_Unset_Lock(self%lock_)
                return
             end if
          end select
          hookPrevious_ => hook_
          hook_         => hook_%next
       end do
    end if
    call Galacticus_Error_Report('object/function not attached to this event'//{introspection:location})
    !$ call OMP_Unset_Lock(self%lock_)
    return
  end subroutine eventHookUnspecifiedDetach

  integer function eventHookCount(self)
    !% Return a count of the number of hooks into this event.
    use Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class(eventHook), intent(inout):: self

    !$ if (.not.self%initialized_) call Galacticus_Error_Report('event has not been initialized'//{introspection:location})
    !$ call OMP_Set_Lock(self%lock_)
    eventHookCount=self%count_
    !$ call OMP_Unset_Lock(self%lock_)
    return
  end function eventHookCount

  function eventHookFirst(self)
    !% Return a pointer to the first hook into this event.
    use Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class(hook     ), pointer      :: eventHookFirst
    class(eventHook), intent(inout):: self

    !$ if (.not.self%initialized_) call Galacticus_Error_Report('event has not been initialized'//{introspection:location})
    eventHookFirst => self%first_
    return
  end function eventHookFirst
  
end module Events_Hooks
