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
Contains a module which wraps the \href{http://www.cs.umd.edu/~mount/ANN/}{ANN} (Approximate Nearest Neighbor) library.
!!}

! Specify an explicit dependence on the ANN.o object file.
!: $(BUILDPATH)/ANN.o

module Nearest_Neighbors
  !!{
  Wraps the \href{http://www.cs.umd.edu/~mount/ANN/}{ANN} (Approximate Nearest Neighbor) library.
  !!}
  use, intrinsic :: ISO_C_Binding, only : C_Null_Ptr, c_ptr, c_int, c_double
  private
  public :: nearestNeighbors, nearestNeighborsClose

  type :: nearestNeighbors
     !!{
     Wrapper object for nearest neighbor searching.
     !!}
     type(c_ptr) :: ANNkd_tree=C_Null_Ptr
   contains
     !![
     <methods>
       <method description="Compute indices and distances to the approximate nearest neighbors." method="search" />
       <method description="Compute indices and distances to the approximate nearest neighbors." method="searchFixedRadius" />
     </methods>
     !!]
     final     ::                      nearestNeighborsDestructor
     procedure :: search            => nearestNeighborsSearch
     procedure :: searchFixedRadius => nearestNeighborsSearchFixedRadius
  end type nearestNeighbors

  interface nearestNeighbors
     module procedure :: nearestNeighborsConstructor
  end interface nearestNeighbors

#ifdef ANNAVAIL
  interface
     function nearestNeighborsConstructorC(n,d,pa) bind(c,name='nearestNeighborsConstructorC')
       !!{
       Template for a C function that constructs an {\normalfont \ttfamily ANNkd\_tree}.
       !!}
       import
       type   (c_ptr   )                              :: nearestNeighborsConstructorC
       integer(c_int   ), value       , intent(in   ) :: n, d
       real   (c_double), dimension(*), intent(in   ) :: pa
     end function nearestNeighborsConstructorC
  end interface

  interface
    subroutine nearestNeighborsSearchC(ANN,point,neighborCount,tolerance,neighborIndex,neighborDistance) bind(c,name='nearestNeighborsSearchC')
       !!{
       Template for a C function that searches for nearest neighbors.
       !!}
      import
      type   (c_ptr   ), value       , intent(in   ) :: ANN
      real   (c_double), dimension(*), intent(in   ) :: point
      integer(c_int   ), value       , intent(in   ) :: neighborCount
      real   (c_double), value       , intent(in   ) :: tolerance
      integer(c_int   ), dimension(*), intent(  out) :: neighborIndex
      real   (c_double), dimension(*), intent(  out) :: neighborDistance
    end subroutine nearestNeighborsSearchC
  end interface

  interface
     function nearestNeighborsSearchFixedRadiusC(ANN,point,radiusSquared,neighborCount,neighborIndex,neighborDistance,tolerance) bind(c,name='nearestNeighborsSearchFixedRadiusC')
       !!{
       Template for a C function that searches for nearest neighbors.
       !!}
       import
       integer(c_int   )                              :: nearestNeighborsSearchFixedRadiusC
       type   (c_ptr   ), value       , intent(in   ) :: ANN
       real   (c_double), dimension(*), intent(in   ) :: point
       real   (c_double), value       , intent(in   ) :: radiusSquared
       integer(c_int   ), value       , intent(in   ) :: neighborCount
       integer(c_int   ), dimension(*), intent(  out) :: neighborIndex
       real   (c_double), dimension(*), intent(  out) :: neighborDistance
       real   (c_double), value       , intent(in   ) :: tolerance
     end function nearestNeighborsSearchFixedRadiusC
  end interface

  interface
     subroutine nearestNeighborsDestructorC(ANN) bind(c,name='nearestNeighborsDestructorC')
       !!{
       Template for a C function that destroys an {\normalfont \ttfamily ANNkd\_tree}.
       !!}
       import
       type(c_ptr), intent(in   ), value :: ANN
     end subroutine nearestNeighborsDestructorC
  end interface

  interface
     subroutine nearestNeighborsCloseC() bind(c,name='nearestNeighborsCloseC')
       !!{
       Template for a C function that closes the ANN library.
       !!}
       import
     end subroutine nearestNeighborsCloseC
  end interface
#endif

contains

  function nearestNeighborsConstructor(points) result(self)
    !!{
    Constructs a nearest neighbor search object.
    !!}
#ifndef ANNAVAIL
    use :: Error, only : Error_Report
#endif
    implicit none
    type            (nearestNeighbors)                                :: self
    double precision                  , intent(in   ), dimension(:,:) :: points

#ifdef ANNAVAIL
    self%ANNkd_tree=nearestNeighborsConstructorC(size(points,dim=1),size(points,dim=2),points)
#else
    !$GLC attributes unused :: points
    self%ANNkd_tree=C_Null_Ptr
    call Error_Report('ANN library is required but was not found'//{introspection:location})
#endif
    return
  end function nearestNeighborsConstructor

  subroutine nearestNeighborsDestructor(self)
    !!{
    Destroys a nearest neighbor search object.
    !!}
#ifdef ANNAVAIL
    use, intrinsic :: ISO_C_Binding, only : c_associated
#else
    use            :: Error        , only : Error_Report
#endif
    implicit none
    type(nearestNeighbors), intent(inout) :: self

#ifdef ANNAVAIL
    if (c_associated(self%ANNkd_tree)) call nearestNeighborsDestructorC(self%ANNkd_tree)
#else
    !$GLC attributes unused :: self
    call Error_Report('ANN library is required but was not found'//{introspection:location})
#endif
    return
  end subroutine nearestNeighborsDestructor

  subroutine nearestNeighborsClose()
    !!{
    Closes the ANN (Approximate Nearest Neighbor) library.
    !!}
#ifndef ANNAVAIL
    use :: Error, only : Error_Report
#endif
    implicit none

#ifdef ANNAVAIL
    call nearestNeighborsCloseC()
#else
    call Error_Report('ANN library is required but was not found'//{introspection:location})
#endif
    return
  end subroutine nearestNeighborsClose

  subroutine nearestNeighborsSearch(self,point,neighborCount,tolerance,neighborIndex,neighborDistance)
    !!{
    Return indices and distances to the (approximate) nearest neighbors.
    !!}
#ifndef ANNAVAIL
    use :: Error, only : Error_Report
#endif
    implicit none
    class           (nearestNeighbors)              , intent(inout) :: self
    double precision                  , dimension(:), intent(in   ) :: point
    integer                                         , intent(in   ) :: neighborCount
    double precision                                , intent(in   ) :: tolerance
    integer                           , dimension(:), intent(  out) :: neighborIndex
    double precision                  , dimension(:), intent(  out) :: neighborDistance

#ifdef ANNAVAIL
    ! Perform the search.
    call nearestNeighborsSearchC(self%ANNkd_tree,point,neighborCount,tolerance,neighborIndex,neighborDistance)
    ! Adjust indices to Fortran standard.
    neighborIndex=neighborIndex+1
#else
    !$GLC attributes unused :: self, point, neighborCount, tolerance
    neighborIndex   =0
    neighborDistance=0.0d0
    call Error_Report('ANN library is required but was not found'//{introspection:location})
#endif
    return
  end subroutine nearestNeighborsSearch

  subroutine nearestNeighborsSearchFixedRadius(self,point,radius,tolerance,neighborCount,neighborIndex,neighborDistance)
    !!{
    Return indices and distances to all neighbors within a given {\normalfont \ttfamily radius}.
    !!}
#ifndef ANNAVAIL
    use :: Error            , only : Error_Report
#endif
    implicit none
    class           (nearestNeighbors)                           , intent(inout)           :: self
    double precision                  , dimension(:)             , intent(in   )           :: point
    double precision                                             , intent(in   )           :: radius
    double precision                                             , intent(in   )           :: tolerance
    integer                           , dimension(:), allocatable, intent(inout), optional :: neighborIndex
    double precision                  , dimension(:), allocatable, intent(inout), optional :: neighborDistance
    integer                                                      , intent(  out)           :: neighborCount
#ifdef ANNAVAIL
    double precision                                                                       :: radiusSquared
    integer                           , dimension(0)                                       :: neighborIndex0
    double precision                  , dimension(0)                                       :: neighborDistance0
    integer                                                                                :: arraySize

    ! Get the squared radius.
    radiusSquared=radius**2
    ! Call once with current array sizes.
    arraySize=0
    if (present(neighborIndex)) then
       if (allocated(neighborIndex)) arraySize=size(neighborIndex)
       neighborCount=nearestNeighborsSearchFixedRadiusC(self%ANNkd_tree,point,radiusSquared,arraySize,neighborIndex,neighborDistance,tolerance)
       ! Resize arrays if necessary.
       if (neighborCount > arraySize) then
          if (allocated(neighborIndex   )) deallocate(neighborIndex   )
          if (allocated(neighborDistance)) deallocate(neighborDistance)
          allocate(neighborIndex   (neighborCount))
          allocate(neighborDistance(neighborCount))
          ! Call again to get all neighbors.
          neighborCount=nearestNeighborsSearchFixedRadiusC(self%ANNkd_tree,point,radiusSquared,neighborCount,neighborIndex,neighborDistance,tolerance)
       end if
    else
       neighborCount=nearestNeighborsSearchFixedRadiusC(self%ANNkd_tree,point,radiusSquared,arraySize,neighborIndex0,neighborDistance0,tolerance)
    end if
    ! Adjust indices to Fortran standard.
    if (neighborCount > 0 .and. present(neighborIndex)) neighborIndex=neighborIndex+1
#else
    !$GLC attributes unused :: neighborDistance, neighborIndex, point, radius, self, tolerance
    neighborCount=0
    call Error_Report('ANN library is required but was not found'//{introspection:location})
#endif
    return
  end subroutine nearestNeighborsSearchFixedRadius

end module Nearest_Neighbors
