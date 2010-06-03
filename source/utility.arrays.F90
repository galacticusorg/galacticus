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






!% Contains a module which implements useful operations on arrays.

module Array_Utilities
  !% Contains routines which implement useful operations on arrays.
  private
  public :: Array_Reverse

  interface Array_Reverse
     !% Interface to generic routines which reverse the direction of an array.
     module procedure Array_Reverse_Real
     module procedure Array_Reverse_Double
  end interface

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

end module Array_Utilities
