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

  !!{RST
  Implementation of an "null" solver for galactic structure.
  !!}

  !![
  <galacticStructureSolver name="galacticStructureSolverNull" docformat="rst">
   <description>
   A no-op galactic structure solver that performs no radius solving, useful as a placeholder when galactic structure calculations are not needed or as a baseline for testing.
   </description>
  </galacticStructureSolver>
  !!]
  type, extends(galacticStructureSolverHooked) :: galacticStructureSolverNull
     !!{RST
     Implementation of an "null" solver for galactic structure.
     !!}
     private
   contains
     final     ::          nullDestructor
     procedure :: solve => nullSolve
  end type galacticStructureSolverNull

  interface galacticStructureSolverNull
     !!{RST
     Constructors for the :galacticus-class:`galacticStructureSolverNull` galactic structure solver class.
     !!}
     module procedure nullConstructorParameters
  end interface galacticStructureSolverNull

contains

  function nullConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`galacticStructureSolverNull` galactic structure solver class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(galacticStructureSolverNull)                :: self
    type(inputParameters            ), intent(inout) :: parameters

    self=galacticStructureSolverNull()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function nullConstructorParameters

  subroutine nullDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`galacticStructureSolverNull` galactic structure solver class.
    !!}
    implicit none
    type(galacticStructureSolverNull), intent(inout) :: self

    call self%detachHooks()
    return
  end subroutine nullDestructor

  subroutine nullSolve(self,node,plausibilityOnly)
    !!{RST
    Solve for the structure of galactic components.
    !!}
    use :: Galactic_Structure_Radius_Solver_Utilities, only : radiusSolverPlausibilities
    implicit none
    class  (galacticStructureSolverNull), intent(inout)           :: self
    type   (treeNode                   ), intent(inout), target   :: node
    logical                             , intent(in   ), optional :: plausibilityOnly
    !$GLC attributes unused :: self, plausibilityOnly

    node%isPhysicallyPlausible=.true.
    node%isSolvable           =.true.
    call radiusSolverPlausibilities(node)
    return
  end subroutine nullSolve
