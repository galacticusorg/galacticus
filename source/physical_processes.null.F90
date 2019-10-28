!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

!% Contains a module which implements a null physical process class.

  !# <physicalProcess name="physicalProcessNull">
  !#  <description>A null physical process class.</description>
  !# </physicalProcess>
  type, extends(physicalProcessClass) :: physicalProcessNull
     !% A null physical process class.
     private
   contains
     procedure :: type => nullNodePromote
  end type physicalProcessNull

  interface physicalProcessNull
     !% Constructors for the ``null'' physical process class.
     module procedure nullConstructorParameters
  end interface physicalProcessNull

contains

  function nullConstructorParameters(parameters) result(self)
    !% Constructor for the ``null'' physical process class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(physicalProcessNull)                :: self
    type(inputParameters    ), intent(inout) :: parameters
    !GCC$ attributes unused :: parameters

    self=physicalProcessNull()
    return
  end function nullConstructorParameters

  subroutine nullNodePromotion(self,node)
    !% Act on node promotion.
    implicit none
    class(physicalProcessNull), intent(inout) :: self
    type (treeNode           ), intent(inout) :: node
    !GCC$ attributes unused :: self, node

    ! This is a null process - do nothing.
    return
  end subroutine nullNodePromotion
