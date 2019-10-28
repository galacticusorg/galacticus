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

!% Contains a module which implements a physical process class that shifts node indices at node promotion.

  !# <physicalProcess name="physicalProcessIndexShift">
  !#  <description>A  physical process class that shifts node indices at node promotion.</description>
  !# </physicalProcess>
  type, extends(physicalProcessClass) :: physicalProcessIndexShift
     !% A  physical process class that shifts node indices at node promotion.
     private
   contains
     procedure :: nodePromote => indexShiftNodePromote
  end type physicalProcessIndexShift

  interface physicalProcessIndexShift
     !% Constructors for the {\normalfont \ttfamily indexShift} physical process class.
     module procedure indexShiftConstructorParameters
  end interface physicalProcessIndexShift

contains

  function indexShiftConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily indexShift} physical process class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(physicalProcessIndexShift)                :: self
    type(inputParameters          ), intent(inout) :: parameters
    !GCC$ attributes unused :: parameters
    
    self=physicalProcessIndexShift()
    return
  end function indexShiftConstructorParameters

  subroutine indexShiftNodePromote(self,node)
    !% Act on node promotion.
    implicit none
    class(physicalProcessIndexShift), intent(inout) :: self
    type (treeNode                 ), intent(inout) :: node
    !GCC$ attributes unused :: self

    ! Shift the index from our node to the parent node.
    call node%parent%indexSet(node%index())
    return
  end subroutine indexShiftNodePromote
