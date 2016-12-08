!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
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

  !% Contains a module which implements advanced locks.

module Locks
  !% Provides advanced locks.
  use OMP_Lib
  implicit none
  private
  public :: ompReadWriteLock

  type :: ompReadWriteLock
     !% OpenMP lock type which supports read/write locking.
     integer(omp_lock_kind), allocatable, dimension(:) :: locks
   contains
     !@ <objectMethods>
     !@   <object>ompReadWriteLock</object>
     !@   <objectMethod>
     !@     <method>setRead</method>
     !@     <type>\void</type>
     !@     <arguments></arguments>
     !@     <description>Obtain a read (non-blocking) lock on the object.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>unsetRead</method>
     !@     <type>\void</type>
     !@     <arguments></arguments>
     !@     <description>Release a read (non-blocking) lock on the object.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>setWrite</method>
     !@     <type>\void</type>
     !@     <arguments>\logicalzero\ [haveReadLock]\argin</arguments>
     !@     <description>Obtain a write (blocking) lock on the object. The lock will block until all other read/write locks on the object are released and while held will prevent any read locks from being obtained. If the thread requesting the write lock already has a read lock it should set {\normalfont \ttfamily haveReadLock=.true.} when calling this function.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>unsetWrite</method>
     !@     <type>\void</type>
     !@     <arguments>\logicalzero\ [haveReadLock]\argin</arguments>
     !@     <description>Release a write (blocking) lock on the object. If the thread releasing the write lock already had a read lock it should set {\normalfont \ttfamily haveReadLock=.true.} when calling this function to ensure that read locked is retained.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     final     ::               ompReadWriteLockDestructor
     procedure :: setRead    => ompReadWriteLockSetRead
     procedure :: unsetRead  => ompReadWriteLockUnsetRead
     procedure :: setWrite   => ompReadWriteLockSetWrite
     procedure :: unsetWrite => ompReadWriteLockUnsetWrite
  end type ompReadWriteLock

  interface ompReadWriteLock
     !% Interface to constructors for OpenMP read/write locks.
     module procedure :: ompReadWriteLockConstructor
  end interface ompReadWriteLock

contains
  
  function ompReadWriteLockConstructor() result (self)
    !% Constructor for OpenMP read/write lock objects.
    implicit none
    type   (ompReadWriteLock) :: self
    integer                   :: i

    ! Allocate a lock for each thread in the current scope.
    allocate(self%locks(0:omp_get_num_threads()-1))
    ! Initialize each lock.
    do i=0,ubound(self%locks,dim=1)
       call OMP_Init_Lock(self%locks(i))
    end do
    return
  end function ompReadWriteLockConstructor

  subroutine ompReadWriteLockDestructor(self)
    !% Destructor for OpenMP read/write lock objects.
    implicit none
    type   (ompReadWriteLock), intent(inout) :: self
    integer                                  :: i

    ! Destroy each lock.
    do i=0,ubound(self%locks,dim=1)
       call OMP_Destroy_Lock(self%locks(i))
    end do
    ! Deallocate the locks.
    deallocate(self%locks)
    return
  end subroutine ompReadWriteLockDestructor

  subroutine ompReadWriteLockSetRead(self)
    !% Get a read lock on an OpenMP read/write lock objects.
    implicit none
    class(ompReadWriteLock), intent(inout) :: self

    call OMP_Set_Lock(self%locks(omp_get_thread_num()))
    return
   end subroutine ompReadWriteLockSetRead

  subroutine ompReadWriteLockUnsetRead(self)
    !% Release a read lock on an OpenMP read/write lock objects.
    implicit none
    class(ompReadWriteLock), intent(inout) :: self

    call OMP_Unset_Lock(self%locks(omp_get_thread_num()))
    return
  end subroutine ompReadWriteLockUnsetRead
  
  subroutine ompReadWriteLockSetWrite(self,haveReadLock)
    !% Get a write lock on an OpenMP read/write lock objects.
    implicit none
    class  (ompReadWriteLock), intent(inout)               :: self
    logical                  , intent(in   ), optional     :: haveReadLock
    integer                                                :: i
    !# <optionalArgument name="haveReadLock" defaultsTo=".true." />

    ! If we have a read lock, release it to avoid deadlocks.
    if (haveReadLock_) call self%unsetRead()
    ! We must obtain all locks in sequence to avoid the possibility of deadlocking against another thread attempting to obtain a
    ! write lock.
    do i=0,ubound(self%locks,dim=1)
       call OMP_Set_Lock(self%locks(i))
    end do
    return
  end subroutine ompReadWriteLockSetWrite

  subroutine ompReadWriteLockUnsetWrite(self,haveReadLock)
    !% Release a write lock on an OpenMP read/write lock objects.
    implicit none
    class  (ompReadWriteLock), intent(inout)           :: self
    logical                  , intent(in   ), optional :: haveReadLock
    integer                                            :: i
    !# <optionalArgument name="haveReadLock" defaultsTo=".true." />

    do i=0,ubound(self%locks,dim=1)
       call OMP_Unset_Lock(self%locks(i))
    end do
    ! Reobtain a read lock if we had one previously.
    if (haveReadLock_) call self%setRead()
    return
  end subroutine ompReadWriteLockUnsetWrite

end module Locks
