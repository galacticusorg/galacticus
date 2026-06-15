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

!!{RST
Contains a module which implements a reusable pool of reference-counted :galacticus-class:`functionClass` objects.
!!}

module Object_Pools
  !!{RST
  Implements a reusable pool of reference-counted :galacticus-class:`functionClass` objects, so that objects which would otherwise be created and destroyed on every call (e.g. per-:term:`node` mass distributions) can instead be re-used. An object held in a pool is considered available for re-use when its reference count has fallen to :math:`1`, i.e. the only remaining reference to it is that held by the pool itself.
  !!}
  use :: Function_Classes, only : functionClass
  implicit none
  private
  public :: objectPool

  type :: poolSlot
     !!{RST
     A single slot in an :galacticus-class:`objectPool`, holding a pointer to a pooled object.
     !!}
     class(functionClass), pointer :: object_ => null()
  end type poolSlot

  type :: objectPool
     !!{RST
     A pool of reference-counted :galacticus-class:`functionClass` objects.
     !!}
     type(poolSlot), allocatable, dimension(:) :: slots
   contains
     !![
     <methods docformat="rst">
       <method method="acquire" description="Return the index of an available pool slot, growing the pool if necessary."/>
       <method method="destroy" description="Release all pooled objects."                                               />
     </methods>
     !!]
     procedure :: acquire => objectPoolAcquire
     procedure :: destroy => objectPoolDestroy
  end type objectPool

contains

  subroutine objectPoolAcquire(self,index,reused)
    !!{RST
    Return the  index of a pool slot for the caller to use. If a slot holding an available object (reference count equal to :math:`1`, i.e. held only by the pool) is found then  reused is set to  .true. and the caller should re-initialize the object in  self%slots(index)%object_. Otherwise the pool is grown by one slot,  reused is set to  .false., and the caller should allocate and construct a new object into that slot.
    !!}
    implicit none
    class  (objectPool), intent(inout)               :: self
    integer            , intent(  out)               :: index
    logical            , intent(  out)               :: reused
    type   (poolSlot  ), allocatable  , dimension(:) :: slotsTmp
    integer                                          :: i

    if (allocated(self%slots)) then
       ! Search for an available slot, i.e. one whose object is referenced only by the pool.
       do i=1,size(self%slots)
          if (self%slots(i)%object_%referenceCount == 1) then
             index =i
             reused=.true.
             return
          end if
       end do
       ! No slot was available - grow the pool by one slot.
       call move_alloc(self%slots,slotsTmp)
       allocate(self%slots(size(slotsTmp)+1))
       do i=1,size(slotsTmp)
          self%slots(i)%object_ => slotsTmp(i)%object_
       end do
       index=size(self%slots)
    else
       ! The pool is empty - create the first slot.
       allocate(self%slots(1))
       index=1
    end if
    reused=.false.
    return
  end subroutine objectPoolAcquire

  subroutine objectPoolDestroy(self)
    !!{RST
    Release all objects held by the pool. This should be called from the destructor of the object that owns the pool.
    !!}
    implicit none
    class(objectPool   ), intent(inout) :: self
    class(functionClass), pointer       :: object_
    integer                             :: i

    if (.not.allocated(self%slots)) return
    do i=1,size(self%slots)
       object_ => self%slots(i)%object_
       !![
       <objectDestructor name="object_"/>
       !!]
       self%slots(i)%object_ => null()
    end do
    deallocate(self%slots)
    return
  end subroutine objectPoolDestroy

end module Object_Pools
