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
  Implementation of an ``null'' solver for galactic structure.
  !!}

  !![
  <galacticStructureSolver name="galacticStructureSolverNull">
   <description>An ``null'' solver for galactic structure.</description>
  </galacticStructureSolver>
  !!]
  type, extends(galacticStructureSolverClass) :: galacticStructureSolverNull
     !!{
     Implementation of an ``null'' solver for galactic structure.
     !!}
     private
   contains
     final     ::             nullDestructor
     procedure :: solve    => nullSolve
     procedure :: autoHook => nullAutoHook
  end type galacticStructureSolverNull

  interface galacticStructureSolverNull
     !!{
     Constructors for the \refClass{galacticStructureSolverNull} galactic structure solver class.
     !!}
     module procedure nullConstructorParameters
  end interface galacticStructureSolverNull

contains

  function nullConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{galacticStructureSolverNull} galactic structure solver class which takes a
    parameter set as input.
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

  subroutine nullAutoHook(self)
    !!{
    Attach to various event hooks.
    !!}
    use :: Events_Hooks, only : nodePromotionEvent  , openMPThreadBindingAtLevel, postEvolveEvent, preDerivativeEvent, &
          &                     satelliteMergerEvent
    implicit none
    class(galacticStructureSolverNull), intent(inout) :: self

    call   preDerivativeEvent%attach(self,nullSolvePreDeriativeHook,openMPThreadBindingAtLevel,label='galacticStructureSolverNull')
    call      postEvolveEvent%attach(self,nullSolveHook            ,openMPThreadBindingAtLevel,label='galacticStructureSolverNull')
    call satelliteMergerEvent%attach(self,nullSolveHook            ,openMPThreadBindingAtLevel,label='galacticStructureSolverNull')
    call   nodePromotionEvent%attach(self,nullSolveHook            ,openMPThreadBindingAtLevel,label='galacticStructureSolverNull')
    return
  end subroutine nullAutoHook

  subroutine nullDestructor(self)
    !!{
    Destructor for the \refClass{galacticStructureSolverNull} galactic structure solver class.
    !!}
    use :: Events_Hooks, only : nodePromotionEvent, postEvolveEvent, preDerivativeEvent, satelliteMergerEvent
    implicit none
    type(galacticStructureSolverNull), intent(inout) :: self

    if (  preDerivativeEvent%isAttached(self,nullSolvePreDeriativeHook)) call   preDerivativeEvent%detach(self,nullSolvePreDeriativeHook)
    if (     postEvolveEvent%isAttached(self,nullSolveHook            )) call      postEvolveEvent%detach(self,nullSolveHook            )
    if (satelliteMergerEvent%isAttached(self,nullSolveHook            )) call satelliteMergerEvent%detach(self,nullSolveHook            )
    if (  nodePromotionEvent%isAttached(self,nullSolveHook            )) call   nodePromotionEvent%detach(self,nullSolveHook            )
    return
  end subroutine nullDestructor

  subroutine nullSolveHook(self,node)
    !!{
    Hookable wrapper around the solver.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(*       ), intent(inout)         :: self
    type (treeNode), intent(inout), target :: node

    select type (self)
    type is (galacticStructureSolverNull)
       call self%solve(node)
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine nullSolveHook

  subroutine nullSolvePreDeriativeHook(self,node,propertyType)
    !!{
    Hookable wrapper around the solver.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class  (*       ), intent(inout)         :: self
    type   (treeNode), intent(inout), target :: node
    integer          , intent(in   )         :: propertyType
    !$GLC attributes unused :: propertyType

    select type (self)
    type is (galacticStructureSolverNull)
       call self%solve(node)
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine nullSolvePreDeriativeHook

  subroutine nullSolve(self,node,plausibilityOnly)
    !!{
    Solve for the structure of galactic components.
    !!}
    include 'galactic_structure.radius_solver.plausible.modules.inc'
    implicit none
    class  (galacticStructureSolverNull), intent(inout)           :: self
    type   (treeNode                   ), intent(inout), target   :: node
    logical                             , intent(in   ), optional :: plausibilityOnly
    !$GLC attributes unused :: self, plausibilityOnly

    node%isPhysicallyPlausible=.true.
    node%isSolvable           =.true.
    include 'galactic_structure.radius_solver.plausible.inc'
    return
  end subroutine nullSolve
