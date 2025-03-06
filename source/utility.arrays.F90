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
Contains a module which implements useful operations on arrays.
!!}

module Array_Utilities
  !!{
  Contains routines which implement useful operations on arrays.
  !!}
  implicit none
  private
  public :: Array_Reverse, Array_Cumulate, Array_Is_Monotonic      , Array_Is_Uniform, &
       &    Array_Which  , Array_Index   , operator(.intersection.), slice5Dto3D     , &
       &    slice5Dto2D

  interface operator(.intersection.)
     module procedure Array_Intersection_Varying_String
  end interface operator(.intersection.)

  interface Array_Reverse
     !!{
     Interface to generic routines which reverse the direction of an array.
     !!}
     module procedure Array_Reverse_Real
     module procedure Array_Reverse_Double
     module procedure Array_Reverse_SizeT
  end interface Array_Reverse

  interface Array_Cumulate
     !!{
     Interface to generic routines which cumulate values in an array.
     !!}
     module procedure Array_Cumulate_Double
  end interface Array_Cumulate

  interface Array_Is_Monotonic
     !!{
     Interface to generic routines which check if an array is monotonic.
     !!}
     module procedure Array_Is_Monotonic_Integer8
     module procedure Array_Is_Monotonic_Double
  end interface Array_Is_Monotonic

  interface Array_Index
     !!{
     Interface to generic routines which return a subset of an array given indices into the array.
     !!}
     module procedure Array_Index_Integer8
     module procedure Array_Index_Integer
     module procedure Array_Index_Double
     module procedure Array_Index_Double_2D
  end interface Array_Index

  ! Types of direction for monotonic arrays.
  integer, parameter, public :: directionDecreasing=-1
  integer, parameter, public :: directionIncreasing=+1

contains

  function Array_Reverse_SizeT(array) result (reversedArray)
    !!{
    Reverses the direction of a real array.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    implicit none
    integer(c_size_t), intent(in   )          :: array        (:)
    integer(c_size_t), dimension(size(array)) :: reversedArray
    integer                                   :: i

    forall (i=1:size(array))
       reversedArray(i)=array(size(array)+1-i)
    end forall
    return
  end function Array_Reverse_SizeT

  function Array_Reverse_Real(array) result (reversedArray)
    !!{
    Reverses the direction of a real array.
    !!}
    implicit none
    real   , intent(in   )          :: array        (:)
    real   , dimension(size(array)) :: reversedArray
    integer                         :: i

    forall (i=1:size(array))
       reversedArray(i)=array(size(array)+1-i)
    end forall
    return
  end function Array_Reverse_Real

  function Array_Reverse_Double(array) result (reversedArray)
    !!{
    Reverses the direction of a double precision array.
    !!}
    implicit none
    double precision, intent(in   )          :: array        (:)
    double precision, dimension(size(array)) :: reversedArray
    integer                                  :: i

    forall (i=1:size(array))
       reversedArray(i)=array(size(array)+1-i)
    end forall
    return
  end function Array_Reverse_Double

  function Array_Cumulate_Double(array) result (cumulatedArray)
    !!{
    Cumulates values in a double precision array.
    !!}
    implicit none
    double precision, intent(in   )          :: array         (:)
    double precision, dimension(size(array)) :: cumulatedArray
    integer                                  :: i

    cumulatedArray(1)=array(1)
    if (size(array) > 1) then
       do i=2,size(array)
          cumulatedArray(i)=cumulatedArray(i-1)+array(i)
       end do
    end if
    return
  end function Array_Cumulate_Double

  logical function Array_Is_Monotonic_Double(array,direction,allowEqual)
    !!{
    Checks if a double precision array is monotonic.
    !!}
    implicit none
    double precision, intent(in   )           :: array           (:)
    integer         , intent(in   ), optional :: direction
    logical         , intent(in   ), optional :: allowEqual
    integer                                   :: i
    logical                                   :: allowEqualActual   , arrayIsFlat, isIncreasing

    ! Single element arrays count as monotonic.
    if (size(array) <= 1) then
       Array_Is_Monotonic_Double=.true.
       return
    end if

    ! Determine if equal points are allowed.
    if (present(allowEqual)) then
       allowEqualActual=allowEqual
    else
       allowEqualActual=.false.
    end if

    ! Determine if the array is increasing or decreasing at the start.
    if (allowEqualActual) then
       ! Equal entries are allowed, so scan the array to find the first element which is not equal to the first element.
       i=2                ! Begin with the second element.
       arrayIsFlat=.true. ! Assume that the array is flat (all elements the same) until we can prove otherwise.
       do while (i <= size(array))
          ! Check for a difference in the element compared to the first element.
          if (array(i) /= array(1)) then
             ! Element does differ: 1) determine in what sense it differs; 2) flag that the array is not flat; 3) exit the loop.
             isIncreasing=array(i) > array(1)
             arrayIsFlat=.false.
             exit
          end if
          i=i+1
       end do
       ! If no element differing from the first was found, then the array is flat and so is deemed to be monotonic in this case
       ! (since equal entries are allowed).
       if (arrayIsFlat) then
          Array_Is_Monotonic_Double=.true.
          return
       end if
    else
       if (array(2) == array(1)) then
          Array_Is_Monotonic_Double=.false.
          return
       end if
       isIncreasing=array(2) > array(1)
    end if

    ! Check direction is correct if this was specified.
    if (present(direction)) then
       if (   (direction == directionIncreasing .and. .not.isIncreasing) .or.        &
            & (direction == directionDecreasing .and.      isIncreasing)      ) then
          Array_Is_Monotonic_Double=.false.
          return
       end if
    end if

    ! Check elements are monotonic. We will exit immediately on finding any non-monotonicity, so set the result to false.
    Array_Is_Monotonic_Double=.false.
    do i=2,size(array)
       select case (isIncreasing)
       case (.true.)
          select case (allowEqualActual)
          case (.false.)
             if (array(i) <= array(i-1)) return
          case (.true.)
             if (array(i) <  array(i-1)) return
          end select
       case (.false.)
          select case (allowEqualActual)
          case (.false.)
             if (array(i) >= array(i-1)) return
          case (.true.)
             if (array(i) >  array(i-1)) return
          end select
     end select
    end do
    ! No elements failed the test, so the array must be monotonic.
    Array_Is_Monotonic_Double=.true.
    return
  end function Array_Is_Monotonic_Double

  subroutine Array_Which(mask,indices)
    !!{
    Return an array of indices for which {\normalfont \ttfamily mask} is true.
    !!}
    use :: Error, only : Error_Report
    implicit none
    logical, intent(in   ) :: mask   (:)
    integer, intent(  out) :: indices(:)
    integer                :: index     , matchCount

    matchCount=0
    indices   =0
    do index=1,size(mask)
       if (mask(index)) then
          matchCount=matchCount+1
          if (matchCount > size(indices)) call Error_Report('indices array is too small'//{introspection:location})
          indices(matchCount)=index
       end if
    end do
    return
  end subroutine Array_Which

  function Array_Index_Double(array,indices) result (arraySubset)
    !!{
    Return a subset of a double precision array given a set of indices into the array.
    !!}
    implicit none
    double precision, dimension(:)            , intent(in   ) :: array
    integer         , dimension(:)            , intent(in   ) :: indices
    double precision, dimension(size(indices))                :: arraySubset
    integer                                                   :: i

    forall(i=1:size(indices))
       arraySubset(i)=array(indices(i))
    end forall
    return
  end function Array_Index_Double

  function Array_Index_Integer(array,indices) result (arraySubset)
    !!{
    Return a subset of an integer array given a set of indices into the array.
    !!}
    implicit none
    integer, dimension(:)            , intent(in   ) :: array
    integer, dimension(:)            , intent(in   ) :: indices
    integer, dimension(size(indices))                :: arraySubset
    integer                                          :: i

    forall(i=1:size(indices))
       arraySubset(i)=array(indices(i))
    end forall
    return
  end function Array_Index_Integer

  function Array_Index_Integer8(array,indices) result (arraySubset)
    !!{
    Return a subset of an integer array given a set of indices into the array.
    !!}
    use :: Kind_Numbers, only : kind_int8
    implicit none
    integer(kind=kind_int8), dimension(:)            , intent(in   ) :: array
    integer                , dimension(:)            , intent(in   ) :: indices
    integer(kind=kind_int8), dimension(size(indices))                :: arraySubset
    integer                                                          :: i

    forall(i=1:size(indices))
       arraySubset(i)=array(indices(i))
    end forall
    return
  end function Array_Index_Integer8

  function Array_Index_Double_2D(array,indices,indexOn) result (arraySubset)
    !!{
    Return a subset of a 2D double precision array given a set of indices into the array.
    !!}
    use :: Error, only : Error_Report
    implicit none
    double precision             , dimension(:,:), intent(in   )           :: array
    integer                      , dimension(:  ), intent(in   )           :: indices
    integer                                      , intent(in   ), optional :: indexOn
    double precision, allocatable, dimension(:,:)                          :: arraySubset
    integer                                                                :: i          , indexOnActual

    indexOnActual=2
    if (present(indexOn)) then
       if (indexOn < 1 .or. indexOn > 2) call Error_Report('1≤indexOn≤2'//{introspection:location})
       indexOnActual=indexOn
    end if
    select case (indexOnActual)
    case (1)
       allocate(arraySubset(size(indices),size(array,dim=2)))
       forall(i=1:size(indices))
          arraySubset(i,:)=array(indices(i),:)
       end forall
    case (2)
       allocate(arraySubset(size(array,dim=1),size(indices)))
       forall(i=1:size(indices))
          arraySubset(:,i)=array(:,indices(i))
       end forall
    end select
    return
  end function Array_Index_Double_2D

  logical function Array_Is_Monotonic_Integer8(array,direction,allowEqual)
    !!{
    Checks if an integer array is monotonic.
    !!}
    use :: Kind_Numbers, only : kind_int8
    implicit none
    integer(kind=kind_int8), intent(in   )           :: array           (:)
    integer                , intent(in   ), optional :: direction
    logical                , intent(in   ), optional :: allowEqual
    integer                                          :: i
    logical                                          :: allowEqualActual   , arrayIsFlat, isIncreasing

    ! Single element arrays count as monotonic.
    if (size(array) <= 1) then
       Array_Is_Monotonic_Integer8=.true.
       return
    end if

    ! Determine if equal points are allowed.
    if (present(allowEqual)) then
       allowEqualActual=allowEqual
    else
       allowEqualActual=.false.
    end if

    ! Determine if the array is increasing or decreasing at the start.
    if (allowEqualActual) then
       ! Equal entries are allowed, so scan the array to find the first element which is not equal to the first element.
       i=2                ! Begin with the second element.
       arrayIsFlat=.true. ! Assume that the array is flat (all elements the same) until we can prove otherwise.
       do while (i <= size(array))
          ! Check for a difference in the element compared to the first element.
          if (array(i) /= array(1)) then
             ! Element does differ: 1) determine in what sense it differs; 2) flag that the array is not flat; 3) exit the loop.
             isIncreasing=array(i) > array(1)
             arrayIsFlat=.false.
             exit
          end if
          i=i+1
       end do
       ! If no element differing from the first was found, then the array is flat and so is deemed to be monotonic in this case
       ! (since equal entries are allowed).
       if (arrayIsFlat) then
          Array_Is_Monotonic_Integer8=.true.
          return
       end if
    else
       if (array(2) == array(1)) then
          Array_Is_Monotonic_Integer8=.false.
          return
       end if
       isIncreasing=array(2) > array(1)
    end if

    ! Check direction is correct if this was specified.
    if (present(direction)) then
       if (   (direction == directionIncreasing .and. .not.isIncreasing) .or.        &
            & (direction == directionDecreasing .and.      isIncreasing)      ) then
          Array_Is_Monotonic_Integer8=.false.
          return
       end if
    end if

    ! Check elements are monotonic. We will exit immediately on finding any non-monotonicity, so set the result to false.
    Array_Is_Monotonic_Integer8=.false.
    do i=2,size(array)
       select case (isIncreasing)
       case (.true.)
          select case (allowEqualActual)
          case (.false.)
             if (array(i) <= array(i-1)) return
          case (.true.)
             if (array(i) <  array(i-1)) return
          end select
       case (.false.)
          select case (allowEqualActual)
          case (.false.)
             if (array(i) >= array(i-1)) return
          case (.true.)
             if (array(i) >  array(i-1)) return
          end select
     end select
    end do
    ! No elements failed the test, so the array must be monotonic.
    Array_Is_Monotonic_Integer8=.true.
    return
  end function Array_Is_Monotonic_Integer8

  function Array_Intersection_Varying_String(a,b)
    use :: ISO_Varying_String, only : operator(==), varying_string
    implicit none
    type   (varying_string), allocatable  , dimension(:) :: Array_Intersection_Varying_String
    type   (varying_string), intent(in   ), dimension(:) :: a, b
    integer                                              :: i, c

    c=0
    do i=1,size(a)
       if (any(b == a(i)) .and. count(a(1:i) == a(i)) == 1) c=c+1
    end do
    allocate(Array_Intersection_Varying_String(c))
    if (c > 0) then
       c=0
       do i=1,size(a)
          if (any(b == a(i)) .and. count(a(1:i) == a(i)) == 1) then
             c=c+1
             Array_Intersection_Varying_String(c)=a(i)
          end if
       end do
    end if
    return
  end function Array_Intersection_Varying_String

  logical function Array_Is_Uniform(array,tolerance,logarithmic)
    !!{
    Return true if an array is uniformly distributed (optionally in the logarithm of its
    values) to the given tolerance.
    !!}
    use :: Numerical_Comparison, only : Values_Agree
    implicit none
    double precision, intent(in   ), dimension(:) :: array
    double precision, intent(in   )               :: tolerance
    logical         , intent(in   ), optional     :: logarithmic
    double precision                              :: increment  , incrementExpected
    integer                                       :: i
    !![
    <optionalArgument name="logarithmic" defaultsTo=".false." />
    !!]

    Array_Is_Uniform=.true.
    if (size(array) <= 1) return
    if (logarithmic_) then
       incrementExpected=log(array(size(array))/array(1))/dble(size(array)-1)
    else
       incrementExpected=   (array(size(array))-array(1))/dble(size(array)-1)
    end if
    do i=2,size(array)
       if (logarithmic_) then
          increment=log(array(i)/array(i-1))
       else
          increment=    array(i)-array(i-1)
       end if
       Array_Is_Uniform=Values_Agree(increment,incrementExpected,absTol=tolerance)
       if (.not.Array_Is_Uniform) return
    end do
    return
  end function Array_Is_Uniform

  function slice5Dto3D(array,dimension_,index_) result(slice)
    !!{
    Return a 3D slice of a 5D array.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    use            :: Error        , only : Error_Report
    implicit none
    double precision          , allocatable  , dimension(:,:,:    ) :: slice
    double precision          , intent(in   ), dimension(:,:,:,:,:) :: array
    integer         (c_size_t), intent(in   ), dimension(2        ) :: dimension_  , index_
    integer         (c_size_t)               , dimension(5        ) :: permutation , shape_
    double precision          , allocatable  , dimension(:,:,:,:,:) :: arrayOrdered
    integer                                                         :: i           , j

    ! Validate inputs.
    if (any(dimension_ < 1 .or. dimension_ > 5)) call Error_Report('dimension is out of range'//{introspection:location})
    do i=1,1
       if (any(dimension_(i) == dimension_(i+1:2))) call Error_Report('duplicated dimension'//{introspection:location})
    end do
    ! Permute the array so the dimensions to remove are the final two.
    j=0
    do i=1,5
       if (.not.any(dimension_ == i)) then
          j             =j+1
          permutation(j)=               i
          shape_     (j)=size(array,dim=i)
       end if
    end do
    permutation(4)=               dimension_(1)
    permutation(5)=               dimension_(2)
    shape_     (4)=size(array,dim=dimension_(1))
    shape_     (5)=size(array,dim=dimension_(2))
    arrayOrdered=reshape(array,shape_,order=permutation)
    ! Slice the final two dimensions of the re-ordered array to get out request slice.
    !![
    <allocate variable="slice" shape="shape_"/>
    !!]
    slice=arrayOrdered(:,:,:,index_(1),index_(2))
    return
  end function slice5Dto3D

  function slice5Dto2D(array,dimension_,index_) result(slice)
    !!{
    Return a 2D slice of a 5D array.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    use            :: Error        , only : Error_Report
    implicit none
    double precision          , allocatable  , dimension(:,:      ) :: slice
    double precision          , intent(in   ), dimension(:,:,:,:,:) :: array
    integer         (c_size_t), intent(in   ), dimension(3        ) :: dimension_  , index_
    integer         (c_size_t)               , dimension(5        ) :: permutation , shape_
    double precision          , allocatable  , dimension(:,:,:,:,:) :: arrayOrdered
    integer                                                         :: i           , j

    ! Validate inputs.
    if (any(dimension_ < 1 .or. dimension_ > 5)) call Error_Report('dimension is out of range'//{introspection:location})
    do i=1,2
       if (any(dimension_(i) == dimension_(i+1:3))) call Error_Report('duplicated dimension'//{introspection:location})
    end do
    ! Permute the array so the dimensions to remove are the final three.
    j=0
    do i=1,5
       if (.not.any(dimension_ == i)) then
          j             =j+1
          permutation(j)=               i
          shape_     (j)=size(array,dim=i)
       end if
    end do
    permutation(3)=               dimension_(1)
    permutation(4)=               dimension_(2)
    permutation(5)=               dimension_(3)
    shape_     (3)=size(array,dim=dimension_(1))
    shape_     (4)=size(array,dim=dimension_(2))
    shape_     (5)=size(array,dim=dimension_(3))
    arrayOrdered=reshape(array,shape_,order=permutation)
    ! Slice the final three dimensions of the re-ordered array to get out request slice.
    !![
    <allocate variable="slice" shape="shape_"/>
    !!]
    slice=arrayOrdered(:,:,index_(1),index_(2),index_(3))
    return
  end function slice5Dto2D
  
end module Array_Utilities
