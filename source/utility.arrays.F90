!! Copyright 2009, 2010, Andrew Benson <abenson@caltech.edu>
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
  private
  public :: Array_Reverse, Array_Cumulate, Array_Is_Monotonic

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
     module procedure Array_Is_Monotonic_Double
  end interface

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
    logical                                :: isIncreasing,allowEqualActual

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

    ! Determine if the array is increasing of decreasing in the first two elements.
    isIncreasing=array(size(array)) > array(1)

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

end module Array_Utilities
