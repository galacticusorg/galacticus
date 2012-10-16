!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


!% Contains a module which implements useful operations on arrays.

module Array_Utilities
  !% Contains routines which implement useful operations on arrays.
  implicit none
  private
  public :: Array_Reverse, Array_Cumulate, Array_Is_Monotonic, Array_Which, Array_Index

  interface Array_Reverse
     !% Interface to generic routines which reverse the direction of an array.
     module procedure Array_Reverse_Real
     module procedure Array_Reverse_Double
  end interface

  interface Array_Cumulate
     !% Interface to generic routines which cumulate values in an array.
     module procedure Array_Cumulate_Double
  end interface

  interface Array_Is_Monotonic
     !% Interface to generic routines which check if an array is monotonic.
     module procedure Array_Is_Monotonic_Integer8
     module procedure Array_Is_Monotonic_Double
  end interface

  interface Array_Index
     !% Interface to generic routines which return a subset of an array given indices into the array.
     module procedure Array_Index_Integer8
     module procedure Array_Index_Integer
     module procedure Array_Index_Double
     module procedure Array_Index_Double_2D
  end interface Array_Index

  ! Types of direction for monotonic arrays.
  integer, parameter, public :: directionDecreasing=-1
  integer, parameter, public :: directionIncreasing= 1

contains

  function Array_Reverse_Real(array) result (reversedArray)
    !% Reverses the direction of a real array.
    implicit none
    real, intent(in)             :: array(:)
    real, dimension(size(array)) :: reversedArray
    integer                      :: i
    
    forall (i=1:size(array))
       reversedArray(i)=array(size(array)+1-i)
    end forall
    return
  end function Array_Reverse_Real

  function Array_Reverse_Double(array) result (reversedArray)
    !% Reverses the direction of a double precision array.
    implicit none
    double precision, intent(in)             :: array(:)
    double precision, dimension(size(array)) :: reversedArray
    integer                                  :: i
    
    forall (i=1:size(array))
       reversedArray(i)=array(size(array)+1-i)
    end forall
    return
  end function Array_Reverse_Double

  function Array_Cumulate_Double(array) result (cumulatedArray)
    !% Cumulates values in a double precision array.
    implicit none
    double precision, intent(in)             :: array(:)
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
    !% Checks if a double precision array is monotonic.
    implicit none
    double precision, intent(in)           :: array(:)
    integer,          intent(in), optional :: direction
    logical,          intent(in), optional :: allowEqual
    integer                                :: i
    logical                                :: isIncreasing,allowEqualActual,arrayIsFlat

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
    !% Return an array of indices for which {\tt mask} is true.
    use Galacticus_Error
    implicit none
    logical, intent(in)  :: mask(:)
    integer, intent(out) :: indices(:)
    integer              :: index,matchCount
    
    matchCount=0
    indices   =0
    do index=1,size(mask)
       if (mask(index)) then
          matchCount=matchCount+1
          if (matchCount > size(indices)) call Galacticus_Error_Report('Array_Which','indices array is too small')
          indices(matchCount)=index
       end if
    end do
    return
  end subroutine Array_Which
  
  function Array_Index_Double(array,indices) result (arraySubset)
    !% Return a subset of a double precision array given a set of indices into the array.
    implicit none
    double precision, dimension(:),            intent(in) :: array
    integer,          dimension(:),            intent(in) :: indices
    double precision, dimension(size(indices))            :: arraySubset
    integer                                               :: i

    forall(i=1:size(indices))
       arraySubset(i)=array(indices(i))
    end forall
    return
  end function Array_Index_Double

  function Array_Index_Integer(array,indices) result (arraySubset)
    !% Return a subset of an integer array given a set of indices into the array.
    implicit none
    integer, dimension(:),            intent(in) :: array
    integer, dimension(:),            intent(in) :: indices
    integer, dimension(size(indices))            :: arraySubset
    integer                                      :: i

    forall(i=1:size(indices))
       arraySubset(i)=array(indices(i))
    end forall
    return
  end function Array_Index_Integer

  function Array_Index_Integer8(array,indices) result (arraySubset)
    !% Return a subset of an integer array given a set of indices into the array.
    use Kind_Numbers
    implicit none
    integer(kind=kind_int8), dimension(:),            intent(in) :: array
    integer,                 dimension(:),            intent(in) :: indices
    integer(kind=kind_int8), dimension(size(indices))            :: arraySubset
    integer                                                      :: i

    forall(i=1:size(indices))
       arraySubset(i)=array(indices(i))
    end forall
    return
  end function Array_Index_Integer8

  function Array_Index_Double_2D(array,indices,indexOn) result (arraySubset)
    !% Return a subset of a 2D double precision array given a set of indices into the array.
    use Galacticus_Error
    implicit none
    double precision, dimension(:,:), intent(in   ) :: array
    integer,          dimension(:  ), intent(in   ) :: indices
    integer,          optional      , intent(in   ) :: indexOn
    double precision, dimension(:,:), allocatable   :: arraySubset
    integer                                         :: i,indexOnActual

    indexOnActual=2
    if (present(indexOn)) then
       if (indexOn < 1 .or. indexOn > 2) call Galacticus_Error_Report('Array_Index_Double_2D','1≤indexOn≤2')
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
    !% Checks if an integer array is monotonic.
    use Kind_Numbers
    implicit none
    integer(kind=kind_int8), intent(in)           :: array(:)
    integer,                 intent(in), optional :: direction
    logical,                 intent(in), optional :: allowEqual
    integer                                       :: i
    logical                                       :: isIncreasing,allowEqualActual,arrayIsFlat

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

end module Array_Utilities
