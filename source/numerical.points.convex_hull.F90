!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
Contains a module which provieds tools for working with convex hulls of sets of points.
!!}

! Specify an explicit dependence on the qhull.o object file.
!: $(BUILDPATH)/qhull.o

module Points_Convex_Hull
  !!{
  Provide tools for working with convex hulls of sets of points.
  !!}
  use, intrinsic :: ISO_C_Binding, only : c_ptr , c_null_ptr, c_double, c_int, &
       &                                  c_long, c_bool
  implicit none
  private
  public :: convexHull

  type :: convexHull
     !!{
     Type used to wrap \href{http://www.qhull.org}{qhull} convex hull objects.
     !!}
     type(c_ptr) :: qhull=c_null_ptr
   contains
     !![
     <methods>
       <method method="volume"        description="Return the volume of the convex hull."                    />
       <method method="pointIsInHull" description="Return true if the given point is inside the convex hull."/>
     </methods>
     !!]
     final     ::                  convexHullDestructor
     procedure :: volume        => convexHullVolume
     procedure :: pointIsInHull => convexHullPointIsInHull
  end type convexHull

  interface convexHull
     !!{
     Interface to constructors for convex hull objects.
     !!}
     module procedure convexHullConstructor
  end interface convexHull
  
#ifdef QHULLAVAIL
  interface
     function convexHullConstructorC(n,points,status) bind(c,name='convexHullConstructorC')
       !!{
       Interface to the C++ convex hull constructor.
       !!}
       import c_ptr, c_double, c_long, c_int
       type   (c_ptr   )                 :: convexHullConstructorC
       integer(c_long  ), value          :: n
       real   (c_double), dimension(:,:) :: points
       integer(c_int   )                 :: status
     end function convexHullConstructorC
     subroutine convexHullDestructorC(qhull) bind(c,name='convexHullDestructorC')
       !!{
       Interface to the C++ convex hull destructor.
       !!}
       import c_ptr
       type(c_ptr), value :: qhull
     end subroutine convexHullDestructorC
     function convexHullVolumeC(qhull) bind(c,name='convexHullVolumeC')
       !!{
       Interface to the C++ convex hull volume function.
       !!}
       import c_ptr, c_double
       real(c_double)        :: convexHullVolumeC
       type(c_ptr   ), value :: qhull
     end function convexHullVolumeC
     function convexHullPointIsInHullC(qhull,point) bind(c,name='convexHullPointIsInHullC')
       !!{
       Interface to the C++ convex hull ``point is in hull'' function.
       !!}
       import c_ptr, c_double, c_bool
       logical(c_bool  )               :: convexHullPointIsInHullC
       real   (c_double), dimension(3) :: point
       type   (c_ptr   ), value        :: qhull
     end function convexHullPointIsInHullC
  end interface
#endif
  
contains

  function convexHullConstructor(points) result(self)
    !!{
    Constructor for {\normalfont \ttfamily convexHull} objects.
    !!}
#ifdef QHULLAVAIL
    use :: Error, only : Error_Report
#endif
    implicit none
    type            (convexHull)                                :: self
    double precision            , dimension(:,:), intent(in   ) :: points    
#ifdef QHULLAVAIL
    integer                                                     :: status

    if (size(points,dim=1) /= 3) call Error_Report('3D points are required'//{introspection:location})
    self%qhull=convexHullConstructorC(size(points,dim=2,kind=c_long),points,status)
    if (status /= 0) call Error_Report('convex hull construction failed'//{introspection:location})
#else
    !$GLC attributes unused :: points
    self%qhull=c_null_ptr
#endif
    return
  end function convexHullConstructor

  subroutine convexHullDestructor(self)
    !!{
    Destructor for {\normalfont \ttfamily convexHull} objects.
    !!}
#ifdef QHULLAVAIL
    use, intrinsic :: ISO_C_Binding, only : c_associated
#endif
    implicit none
    type(convexHull), intent(inout) :: self
    
#ifdef QHULLAVAIL
    if (c_associated(self%qhull)) call convexHullDestructorC(self%qhull)
#else
    !$GLC attributes unused :: self
#endif
    return
  end subroutine convexHullDestructor

  double precision function convexHullVolume(self)
    !!{
    Return the volume of a convex hull.
    !!}
#ifndef QHULLAVAIL
    use :: Error, only : Error_Report
#endif
    implicit none
    class(convexHull), intent(inout) :: self

#ifdef QHULLAVAIL
    convexHullVolume=convexHullVolumeC(self%qhull)
#else
    !$GLC attributes unused :: self
    convexHullVolume=0.0d0
    call Error_Report('qhull library is unavailable'//{introspection:location})
#endif
    return
  end function convexHullVolume
  
  logical function convexHullPointIsInHull(self,point)
    !!{
    Return the volume of a convex hull.
    !!}
#ifndef QHULLAVAIL
    use :: Error, only : Error_Report
#endif
    implicit none
    class(convexHull), intent(inout) :: self
    real (c_double  ), dimension(3)  :: point

#ifdef QHULLAVAIL
    convexHullPointIsInHull=convexHullPointIsInHullC(self%qhull,point)
#else
    !$GLC attributes unused :: self, point
    convexHullPointIsInHull=.false.
    call Error_Report('qhull library is unavailable'//{introspection:location})
#endif
    return
  end function convexHullPointIsInHull
  
end module Points_Convex_Hull
