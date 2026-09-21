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

  !+    Contributions to this file made by: Claude.

  !!{RST
  Implements an abstract galactic structure solver class which solves for structure in response to node events.
  !!}

  !![
  <galacticStructureSolver name="galacticStructureSolverHooked" abstract="yes" docformat="rst">
   <description>
   An abstract galactic structure solver class which solves for the structure of a node before its derivatives are computed, after it has evolved, after a satellite merges into it, and after it is promoted. This is an abstract class---the solver itself must be provided by a concrete class, which may also override ``solvePreDerivative`` to change how the solver responds to the pre-derivative event.
   </description>
  </galacticStructureSolver>
  !!]
  type, abstract, extends(galacticStructureSolverClass) :: galacticStructureSolverHooked
     !!{RST
     An abstract galactic structure solver class which solves for structure in response to node events.
     !!}
     private
   contains
     !![
     <methods docformat="rst">
       <method description="Solve for the structure of the given node in response to the pre-derivative event, for properties of the given type." method="solvePreDerivative" />
       <method description="Detach from the node events. Must be called by the destructor of each concrete class."                                   method="detachHooks"        />
     </methods>
     !!]
     procedure :: autoHook           => hookedAutoHook
     procedure :: detachHooks        => hookedDetachHooks
     procedure :: solvePreDerivative => hookedSolvePreDerivative
  end type galacticStructureSolverHooked

contains

  subroutine hookedAutoHook(self)
    !!{RST
    Attach to the node events.
    !!}
    use :: Events_Hooks      , only : dependencyDirectionAfter, dependencyRegEx     , nodePromotionEvent, openMPThreadBindingAtLevel, &
          &                           postEvolveEvent         , preDerivativeEvent  , satelliteMergerEvent
    use :: ISO_Varying_String, only : char
    implicit none
    class(galacticStructureSolverHooked), intent(inout) :: self
    type (dependencyRegEx              ), dimension(1)  :: dependencies

    dependencies(1)=dependencyRegEx(dependencyDirectionAfter,'^nodeComponent')
    call   preDerivativeEvent%attach(self,hookedSolvePreDerivativeHook,openMPThreadBindingAtLevel,label=char(self%objectType())                          )
    call      postEvolveEvent%attach(self,hookedSolveHook             ,openMPThreadBindingAtLevel,label=char(self%objectType()),dependencies=dependencies)
    call satelliteMergerEvent%attach(self,hookedSolveHook             ,openMPThreadBindingAtLevel,label=char(self%objectType()),dependencies=dependencies)
    call   nodePromotionEvent%attach(self,hookedSolveHook             ,openMPThreadBindingAtLevel,label=char(self%objectType()),dependencies=dependencies)
    return
  end subroutine hookedAutoHook

  subroutine hookedDetachHooks(self)
    !!{RST
    Detach from the node events.
    !!}
    use :: Events_Hooks, only : nodePromotionEvent, postEvolveEvent, preDerivativeEvent, satelliteMergerEvent
    implicit none
    class(galacticStructureSolverHooked), intent(inout) :: self

    if (  preDerivativeEvent%isAttached(self,hookedSolvePreDerivativeHook)) call   preDerivativeEvent%detach(self,hookedSolvePreDerivativeHook)
    if (     postEvolveEvent%isAttached(self,hookedSolveHook             )) call      postEvolveEvent%detach(self,hookedSolveHook             )
    if (satelliteMergerEvent%isAttached(self,hookedSolveHook             )) call satelliteMergerEvent%detach(self,hookedSolveHook             )
    if (  nodePromotionEvent%isAttached(self,hookedSolveHook             )) call   nodePromotionEvent%detach(self,hookedSolveHook             )
    return
  end subroutine hookedDetachHooks

  subroutine hookedSolvePreDerivative(self,node,propertyType)
    !!{RST
    Solve for the structure of the given node in response to the pre-derivative event.
    !!}
    implicit none
    class  (galacticStructureSolverHooked), intent(inout)         :: self
    type   (treeNode                     ), intent(inout), target :: node
    integer                               , intent(in   )         :: propertyType
    !$GLC attributes unused :: propertyType

    call self%solve(node)
    return
  end subroutine hookedSolvePreDerivative

  subroutine hookedSolveHook(self,node)
    !!{RST
    Hookable wrapper around the solver.
    !!}
    use :: Error             , only : Error_Report
    use :: ISO_Varying_String, only : char
    use :: Function_Classes  , only : functionClass
    implicit none
    class(*       ), intent(inout)         :: self
    type (treeNode), intent(inout), target :: node

    select type (self)
    class is (galacticStructureSolverHooked)
       call self%solve(node)
    class is (functionClass)
       call Error_Report('object is not of [galacticStructureSolverHooked] class, but of ['//char(self%objectType())//'] class'//{introspection:location})
    class default
       call Error_Report('object is not of [galacticStructureSolverHooked] class'//{introspection:location})
    end select
    return
  end subroutine hookedSolveHook

  subroutine hookedSolvePreDerivativeHook(self,node,propertyType)
    !!{RST
    Hookable wrapper around the solver for pre-derivative events.
    !!}
    use :: Error             , only : Error_Report
    use :: ISO_Varying_String, only : char
    use :: Function_Classes  , only : functionClass
    implicit none
    class  (*       ), intent(inout)         :: self
    type   (treeNode), intent(inout), target :: node
    integer          , intent(in   )         :: propertyType

    select type (self)
    class is (galacticStructureSolverHooked)
       call self%solvePreDerivative(node,propertyType)
    class is (functionClass)
       call Error_Report('object is not of [galacticStructureSolverHooked] class, but of ['//char(self%objectType())//'] class'//{introspection:location})
    class default
       call Error_Report('object is not of [galacticStructureSolverHooked] class'//{introspection:location})
    end select
    return
  end subroutine hookedSolvePreDerivativeHook
