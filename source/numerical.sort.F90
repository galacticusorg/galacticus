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
Contains a module which implements sorting sequences.
!!}

! Add dependency on GSL library.
!; gsl

module Sorting
  !!{
  Implements sorting.
  !!}
  use, intrinsic :: ISO_C_Binding     , only : c_ptr, c_funptr, c_size_t, c_int
  use            :: Kind_Numbers      , only : kind_int8
  use            :: ISO_Varying_String, only : varying_string
  implicit none
  private
  public :: sortByIndex, sort, sortIndex, sortSmallest, sortLargest, sortSmallestIndex, sortLargestIndex

  !![
  <generic identifier="Type">
   <instance label="integer"  intrinsic="integer"              exemplar="0"          />
   <instance label="integer8" intrinsic="integer(kind_int8)"   exemplar="0_kind_int8"/>
   <instance label="double"   intrinsic="double precision"     exemplar="0.0d0"      />
   <instance label="varstr"   intrinsic="type(varying_string)" exemplar="var_str('')"/>
  </generic>
  !!]

  interface sort
     !!{
     Generic interface to in-place sort functions.
     !!}
     module procedure sort{Type¦label}
     module procedure sortBoth{Type¦label}
  end interface sort
    
  interface sortIndex
     !!{
     Generic interface to sort functions that return the sort index.
     !!}
     module procedure sortIndex{Type¦label}
  end interface sortIndex
    
  interface sortByIndex
     !!{
     Generic interface to in-place sort functions using a supplied index.
     !!}
     module procedure sortByIndex{Type¦label}
  end interface sortByIndex

  interface
     subroutine gsl_heapsort(array,count_,size_,compare) bind(c,name='gsl_heapsort')
       !!{
       Template for the GSL heapsort function.
       !!}
       import
       type   (c_ptr   ), value :: array
       integer(c_size_t), value :: count_ , size_
       type   (c_funptr), value :: compare
     end subroutine gsl_heapsort

     function gsl_heapsort_index(order,array,count_,size_,compare) bind(c,name='gsl_heapsort_index')
       !!{
       Template for the GSL heapsort index function.
       !!}
       import
       integer(c_int   )        :: gsl_heapsort_index
       type   (c_ptr   ), value :: order
       type   (c_ptr   ), value :: array
       integer(c_size_t), value :: count_            , size_
       type   (c_funptr), value :: compare
     end function gsl_heapsort_index

     function gsl_sort_smallest(dest,k,src,stride,n) bind(c,name='gsl_sort_smallest')
       !!{
       Template for the GSL sort smallest function.
       !!}
       import
       integer(c_int   )        :: gsl_sort_smallest
       type   (c_ptr   ), value :: dest
       integer(c_size_t), value :: k
       type   (c_ptr   ), value :: src
       integer(c_size_t), value :: stride           , n
     end function gsl_sort_smallest

     function gsl_sort_largest(dest,k,src,stride,n) bind(c,name='gsl_sort_largest')
       !!{
       Template for the GSL sort largest function.
       !!}
       import
       integer(c_int   )        :: gsl_sort_largest
       type   (c_ptr   ), value :: dest
       integer(c_size_t), value :: k
       type   (c_ptr   ), value :: src
       integer(c_size_t), value :: stride          , n
     end function gsl_sort_largest

     function gsl_sort_smallest_index(p,k,src,stride,n) bind(c,name='gsl_sort_smallest_index')
       !!{
       Template for the GSL sort smallest index function.
       !!}
       import
       integer(c_int   )        :: gsl_sort_smallest_index
       type   (c_ptr   ), value :: p
       integer(c_size_t), value :: k
       type   (c_ptr   ), value :: src
       integer(c_size_t), value :: stride                 , n
     end function gsl_sort_smallest_index

     function gsl_sort_largest_index(p,k,src,stride,n) bind(c,name='gsl_sort_largest_index')
       !!{
       Template for the GSL sort largest index function.
       !!}
       import
       integer(c_int   )        :: gsl_sort_largest_index
       type   (c_ptr   ), value :: p
       integer(c_size_t), value :: k
       type   (c_ptr   ), value :: src
       integer(c_size_t), value :: stride                , n
     end function gsl_sort_largest_index
  end interface
  
contains

  subroutine sort{Type¦label}(array)
    !!{
    Given an unsorted {\normalfont \ttfamily array}, sorts it in place.
    !!}
    use                         , intrinsic :: ISO_C_Binding     , only : c_loc  , c_funloc
    {Type¦match¦^(varstr)$¦  use            :: ISO_Varying_String, only : var_str¦}
    implicit none
    {Type¦intrinsic}, dimension(:), intent(inout), target :: array

    call GSL_HeapSort(c_loc(array),size(array,kind=c_size_t),sizeof({Type¦exemplar}),c_funloc(compare{Type¦label}))    
    return
  end subroutine sort{Type¦label}
  
  subroutine sortBoth{Type¦label}(array,array2)
    !!{
    Given an unsorted double precision {\normalfont \ttfamily array}, sorts it in place while also rearranging {\normalfont \ttfamily array2} in the same way.
    !!}
    implicit none
    {Type¦intrinsic}                , dimension(:                        ), intent(inout) :: array, array2
    integer         (kind=c_size_t ), dimension(size(array,kind=c_size_t))                :: order
    {Type¦intrinsic}                , dimension(size(array,kind=c_size_t))                :: tmp
    integer         (kind=c_size_t )                                                      :: i

    order=sortIndex(array)
    forall(i=1:size(array,kind=c_size_t))
       tmp(i)=array(order(i))
    end forall
    array=tmp
    forall(i=1:size(array,kind=c_size_t))
       tmp(i)=array2(order(i))
    end forall
    array2=tmp
    return
  end subroutine sortBoth{Type¦label}

  function sortIndex{Type¦label}(array) result(order)
    !!{
    Given an unsorted {\normalfont \ttfamily array}, return the sort index.
    !!}
    use                         , intrinsic :: ISO_C_Binding     , only : c_loc  , c_funloc
    {Type¦match¦^(varstr)$¦  use            :: ISO_Varying_String, only : var_str¦}
    implicit none
    {Type¦intrinsic} , dimension(:                        ), intent(in   ), target :: array
    integer(c_size_t), dimension(size(array,kind=c_size_t)), target                :: order
    integer(c_int   )                                                              :: status

    status=GSL_HeapSort_Index(c_loc(order),c_loc(array),size(array,kind=c_size_t),sizeof({Type¦exemplar}),c_funloc(compare{Type¦label}))    
    order=order+1_c_size_t
    return
  end function sortIndex{Type¦label}

  subroutine sortByIndex{Type¦label}(array,index)
    !!{
    Given an {\normalfont \ttfamily array}, sort it in place using the supplied index.
    !!}
    implicit none
    {Type¦intrinsic}               , dimension(:          ), intent(inout) :: array
    integer         (kind=c_size_t), dimension(:          ), intent(in   ) :: index
    {Type¦intrinsic}               , dimension(size(array))                :: arrayTmp
    integer         (kind=c_size_t)                                        :: i

    forall(i=1:size(array))
       arrayTmp(i)=array(index(i))
    end forall
    array=arrayTmp
    return
  end subroutine sortByIndex{Type¦label}

  function compare{Type¦label}(x,y) bind(c)
    !!{
    Comparison function for sorting.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_double, c_f_pointer, c_int, c_ptr
    {Type¦match¦^(varstr)$¦use ISO_Varying_String, only : operator(<), operator(>)¦}    
    integer(kind=c_int   )          :: compare{Type¦label}
    type   (c_ptr        ), value   :: x                  , y
    {Type¦intrinsic}      , pointer :: rx                 , ry

    call c_f_pointer(x,rx)
    call c_f_pointer(y,ry)
    if      (rx > ry) then
       compare{Type¦label}=+1
    else if (rx < ry) then
       compare{Type¦label}=-1
    else
       compare{Type¦label}= 0
    end if
    return
  end function compare{Type¦label}

  function sortSmallest(array,k,mask) result(arraySort)
    !!{
    Given an {\normalfont \ttfamily array}, find the k smallest elements.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_loc
    implicit none
    double precision          , dimension(:), target  , intent(in   ) :: array
    integer         (c_size_t),                         intent(in   ) :: k
    logical                   , dimension(:), optional, intent(in   ) :: mask
    double precision          , dimension(k), target                  :: arraySort
    double precision          , dimension(:), target  , allocatable   :: arrayMasked
    integer         (c_int   )                                        :: status

    if (present(mask)) then
       allocate(arrayMasked(size(array)))
       where(mask)
          arrayMasked=array
       elsewhere
          arrayMasked=huge(0.0d0)
       end where
       status=GSL_Sort_Smallest(c_loc(arraySort),k,c_loc(arrayMasked),1_c_size_t,size(array,kind=c_size_t))
       deallocate(arrayMasked)
    else
       status=GSL_Sort_Smallest(c_loc(arraySort),k,c_loc(array      ),1_c_size_t,size(array,kind=c_size_t))
    end if
    return
  end function sortSmallest

  function sortLargest(array,k,mask) result(arraySort)
    !!{
    Given an {\normalfont \ttfamily array}, find the k largest elements.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_loc
    implicit none
    double precision          , dimension(:), target  , intent(in   ) :: array
    integer         (c_size_t),                         intent(in   ) :: k
    logical                   , dimension(:), optional, intent(in   ) :: mask
    double precision          , dimension(k), target                  :: arraySort
    double precision          , dimension(:), target  , allocatable   :: arrayMasked
    integer         (c_int   )                                        :: status

    if (present(mask)) then
       allocate(arrayMasked(size(array)))
       where(mask)
          arrayMasked=array
       elsewhere
          arrayMasked=-huge(0.0d0)
       end where
       status=GSL_Sort_Largest(c_loc(arraySort),k,c_loc(arrayMasked),1_c_size_t,size(array,kind=c_size_t))
       deallocate(arrayMasked)
    else
       status=GSL_Sort_Largest(c_loc(arraySort),k,c_loc(array      ),1_c_size_t,size(array,kind=c_size_t))
    end if
    return
  end function sortLargest

  function sortSmallestIndex(array,k,mask) result(index)
    !!{
    Given an {\normalfont \ttfamily array}, find the indices of k smallest elements.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_loc
    implicit none
    double precision          , dimension(:), target  , intent(in   ) :: array
    integer         (c_size_t),                         intent(in   ) :: k
    logical                   , dimension(:), optional, intent(in   ) :: mask
    integer         (c_size_t), dimension(k), target                  :: index
    double precision          , dimension(:), target  , allocatable   :: arrayMasked
    integer         (c_int   )                                        :: status

    if (present(mask)) then
       allocate(arrayMasked(size(array)))
       where(mask)
          arrayMasked=array
       elsewhere
          arrayMasked=huge(0.0d0)
       end where
       status=GSL_Sort_Smallest_Index(c_loc(index),k,c_loc(arrayMasked),1_c_size_t,size(array,kind=c_size_t))
       deallocate(arrayMasked)
    else
       status=GSL_Sort_Smallest_Index(c_loc(index),k,c_loc(array      ),1_c_size_t,size(array,kind=c_size_t))
    end if
    index=index+1_c_size_t
    return
  end function sortSmallestIndex

  function sortLargestIndex(array,k,mask) result(index)
    !!{
     Given an {\normalfont \ttfamily array}, find the indices of k largest elements.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_loc
    implicit none
    double precision          , dimension(:), target  , intent(in   ) :: array
    integer         (c_size_t),                         intent(in   ) :: k
    logical                   , dimension(:), optional, intent(in   ) :: mask
    integer         (c_size_t), dimension(k), target                  :: index
    double precision          , dimension(:), target  , allocatable   :: arrayMasked
    integer         (c_int   )                                        :: status

    if (present(mask)) then
       allocate(arrayMasked(size(array)))
       where(mask)
          arrayMasked=array
       elsewhere
          arrayMasked=-huge(0.0d0)
       end where
       status=GSL_Sort_Largest_Index(c_loc(index),k,c_loc(arrayMasked),1_c_size_t,size(array,kind=c_size_t))
       deallocate(arrayMasked)
    else
       status=GSL_Sort_Largest_Index(c_loc(index),k,c_loc(array      ),1_c_size_t,size(array,kind=c_size_t))
    end if
    index=index+1_c_size_t
    return
  end function sortLargestIndex

end module Sorting
