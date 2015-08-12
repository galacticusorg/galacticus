!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements a galactic filter which always passes.

  !# <galacticFilter name="galacticFilterAlways">
  !#  <description>A galactic filter which always passes. (Used mostly for testing purposes.)</description>
  !# </galacticFilter>
  type, extends(galacticFilterClass) :: galacticFilterAlways
     !% A galactic filter class which always passes.
     private
   contains
     final     ::           alwaysDestructor
     procedure :: passes => alwaysPasses
  end type galacticFilterAlways

  interface galacticFilterAlways
     !% Constructors for the ``always'' galactic filter class.
     module procedure alwaysConstructorParameters
     module procedure alwaysConstructorInternal
  end interface galacticFilterAlways

contains

  function alwaysConstructorParameters(parameters)
    !% Constructor for the ``always'' galactic filter class which takes a parameter set as input.
    use Input_Parameters2
    implicit none
    type(galacticFilterAlways)                :: alwaysConstructorParameters
    type(inputParameters     ), intent(in   ) :: parameters

    return
  end function alwaysConstructorParameters

  function alwaysConstructorInternal()
    !% Internal constructor for the ``always'' galactic filter class.
    implicit none
    type(galacticFilterAlways) :: alwaysConstructorInternal

    return
  end function alwaysConstructorInternal

  elemental subroutine alwaysDestructor(self)
    !% Destructor for the ``always'' galactic filter class.
    implicit none
    type(galacticFilterAlways), intent(inout) :: self

    ! Nothing to do.
    return
  end subroutine alwaysDestructor

  logical function alwaysPasses(self,node)
    !% Implement an always-pass galactic filter.
    implicit none
    class(galacticFilterAlways), intent(inout) :: self
    type (treeNode            ), intent(inout) :: node
    
    alwaysPasses=.true.
    return
  end function alwaysPasses
