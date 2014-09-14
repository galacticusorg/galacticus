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

!% Contains a module which wraps the \href{http://www.cs.umd.edu/~mount/ANN/}{ANN} (Approximate Nearest Neighbor) library.

! Specify an explicit dependence on the ANN.o object file.
!: ./work/build/ANN.o

module Nearest_Neighbors
  !% Wraps the \href{http://www.cs.umd.edu/~mount/ANN/}{ANN} (Approximate Nearest Neighbor) library.
  use, intrinsic :: ISO_C_Binding
  private
  public :: nearestNeighbors, nearestNeighborsClose

  type :: nearestNeighbors
     !% Wrapper object for nearest neighbor searching.
     type(c_ptr) :: ANNkd_tree
   contains
     final     ::                      nearestNeighborsDestructor
     procedure :: search            => nearestNeighborsSearch
     procedure :: searchFixedRadius => nearestNeighborsSearchFixedRadius
  end type nearestNeighbors

  interface nearestNeighbors
     module procedure :: nearestNeighborsConstructor
  end interface nearestNeighbors

  interface
     function nearestNeighborsConstructorC(n,d,pa) bind(c,name='nearestNeighborsConstructorC')
       !% Template for a C function that constructs an {\tt ANNkd\_tree}.
       import
       type   (c_ptr   )                              :: nearestNeighborsConstructorC
       integer(c_int   ), value       , intent(in   ) :: n, d
       real   (c_double), dimension(*), intent(in   ) :: pa
     end function nearestNeighborsConstructorC
  end interface

  interface
    subroutine nearestNeighborsSearchC(ANN,point,neighborCount,tolerance,neighborIndex,neighborDistance) bind(c,name='nearestNeighborsSearchC')
       !% Template for a C function that searches for nearest neighbors.
      import
      type            (c_ptr), value       , intent(in   ) :: ANN
      double precision       , dimension(*), intent(in   ) :: point
      integer                , value       , intent(in   ) :: neighborCount
      double precision       , value       , intent(in   ) :: tolerance
      integer                , dimension(*), intent(  out) :: neighborIndex
      double precision       , dimension(*), intent(  out) :: neighborDistance
    end subroutine nearestNeighborsSearchC
  end interface

  interface
     integer function nearestNeighborsSearchFixedRadiusC(ANN,point,radiusSquared,neighborCount,neighborIndex,neighborDistance,tolerance) bind(c,name='nearestNeighborsSearchFixedRadiusC')
       !% Template for a C function that searches for nearest neighbors.
      import
      type            (c_ptr), value       , intent(in   ) :: ANN
      double precision       , dimension(*), intent(in   ) :: point
      double precision       , value       , intent(in   ) :: radiusSquared
      integer                , value       , intent(in   ) :: neighborCount
      integer                , dimension(*), intent(  out) :: neighborIndex
      double precision       , dimension(*), intent(  out) :: neighborDistance
      double precision       , value       , intent(in   ) :: tolerance
    end function nearestNeighborsSearchFixedRadiusC
  end interface

  interface
     subroutine nearestNeighborsDestructorC(ANN) bind(c,name='nearestNeighborsDestructorC')
       !% Template for a C function that destroys an {\tt ANNkd\_tree}.
       import
       type(c_ptr), intent(in   ), value :: ANN
     end subroutine nearestNeighborsDestructorC
  end interface
  
  interface
     subroutine nearestNeighborsCloseC() bind(c,name='nearestNeighborsCloseC')
       !% Template for a C function that closes the ANN library.
       import
     end subroutine nearestNeighborsCloseC
  end interface
  
contains
  
  function nearestNeighborsConstructor(points)
    !% Constructs a nearest neighor search object.
    implicit none
    type            (nearestNeighbors)                                :: nearestNeighborsConstructor
    double precision                  , intent(in   ), dimension(:,:) :: points
    
    nearestNeighborsConstructor%ANNkd_tree=nearestNeighborsConstructorC(size(points,dim=1),size(points,dim=2),points)
    return
  end function nearestNeighborsConstructor

  subroutine nearestNeighborsDestructor(self)
    !% Destroys a nearest neighbor search object.
    implicit none
    type(nearestNeighbors), intent(inout) :: self

    call nearestNeighborsDestructorC(self%ANNkd_tree)
    return
  end subroutine nearestNeighborsDestructor

  subroutine nearestNeighborsClose()
    !% Closes the ANN (Approximate Nearest Neighbor) library.
    implicit none

    call nearestNeighborsCloseC()
    return
  end subroutine nearestNeighborsClose

  subroutine nearestNeighborsSearch(self,point,neighborCount,tolerance,neighborIndex,neighborDistance)
    !% Return indices and distances to the (approximate) nearest neighbors.
    implicit none
    class           (nearestNeighbors)              , intent(inout) :: self
    double precision                  , dimension(:), intent(in   ) :: point
    integer                                         , intent(in   ) :: neighborCount
    double precision                                , intent(in   ) :: tolerance
    integer                           , dimension(:), intent(  out) :: neighborIndex
    double precision                  , dimension(:), intent(  out) :: neighborDistance
    
    ! Perform the search.
    call nearestNeighborsSearchC(self%ANNkd_tree,point,neighborCount,tolerance,neighborIndex,neighborDistance)
    ! Adjust indices to Fortran standard.
    neighborIndex=neighborIndex+1
    return
  end subroutine nearestNeighborsSearch

  subroutine nearestNeighborsSearchFixedRadius(self,point,radius,tolerance,neighborCount,neighborIndex,neighborDistance)
    !% Return indices and distances to all neighbors within a given {\tt radius}.
    use Memory_Management
    implicit none
    class           (nearestNeighbors)                           , intent(inout) :: self
    double precision                  , dimension(:)             , intent(in   ) :: point
    double precision                                             , intent(in   ) :: radius
    double precision                                             , intent(in   ) :: tolerance
    integer                           , dimension(:), allocatable, intent(inout) :: neighborIndex
    double precision                  , dimension(:), allocatable, intent(inout) :: neighborDistance
    integer                                                      , intent(  out) :: neighborCount
    double precision                                                             :: radiusSquared
    integer                                                                      :: arraySize

    ! Get the squared radius.    
    radiusSquared=radius**2
    ! Call once with current array sizes.
    if (allocated(neighborIndex)) then
       arraySize=size(neighborIndex)
    else
       arraySize=0
    end if
    neighborCount=nearestNeighborsSearchFixedRadiusC(self%ANNkd_tree,point,radiusSquared,arraySize,neighborIndex,neighborDistance,tolerance)    
    ! Resize arrays if necessary.
    if (neighborCount > arraySize) then
       if (allocated(neighborIndex   )) call Dealloc_Array(neighborIndex   )
       if (allocated(neighborDistance)) call Dealloc_Array(neighborDistance)
       call Alloc_Array(neighborIndex   ,[neighborCount])
       call Alloc_Array(neighborDistance,[neighborCount])
       ! Call again to get all neighbors.
       neighborCount=nearestNeighborsSearchFixedRadiusC(self%ANNkd_tree,point,radiusSquared,neighborCount,neighborIndex,neighborDistance,tolerance)
    end if
    ! Adjust indices to Fortran standard.
    if (neighborCount > 0) neighborIndex=neighborIndex+1
    return
  end subroutine nearestNeighborsSearchFixedRadius

end module Nearest_Neighbors
