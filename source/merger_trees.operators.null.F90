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

!% Contains a module which implements a null operator on merger trees.

  !# <mergerTreeOperator name="mergerTreeOperatorNull">
  !#  <description>Provides a null operator on merger trees.</description>
  !# </mergerTreeOperator>
  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorNull
     !% A null merger tree operator class.
     private
   contains
     final     ::            nullDestructor
     procedure :: operate => nullOperate
  end type mergerTreeOperatorNull

  interface mergerTreeOperatorNull
     !% Constructors for the null merger tree operator class.
     module procedure nullConstructorParameters
     module procedure nullConstructorInternal
  end interface mergerTreeOperatorNull

contains

  function nullConstructorParameters(parameters)
    !% Constructor for the null merger tree operator class which takes a parameter set as input.
    use Input_Parameters2
    implicit none
    type(mergerTreeOperatorNull)                :: nullConstructorParameters
    type(inputParameters       ), intent(inout) :: parameters
    
    return
  end function nullConstructorParameters

  function nullConstructorInternal()
    !% Internal constructor for the null merger tree operator class.
    implicit none
    type(mergerTreeOperatorNull) :: nullConstructorInternal

    return
  end function nullConstructorInternal

  elemental subroutine nullDestructor(self)
    !% Destructor for the merger tree operator function class.
    implicit none
    type(mergerTreeOperatorNull), intent(inout) :: self

    ! Nothing to do.
    return
  end subroutine nullDestructor

  subroutine nullOperate(self,tree)
    !% Perform a null operation on a merger tree.
    implicit none
    class(mergerTreeOperatorNull), intent(inout)         :: self
    type (mergerTree            ), intent(inout), target :: tree

    ! Nothing to do.
    return
  end subroutine nullOperate
