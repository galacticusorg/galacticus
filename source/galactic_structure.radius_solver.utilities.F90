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

  !!{
  Contains a module that provides utility functions for galactic structure solver.
  !!}

module Galactic_Structure_Radius_Solver_Utilities
  !!{
  Provides utility functions for galactic structure solver.
  !!}
  use :: Galacticus_Nodes          , only : treeNode
  use :: Galactic_Structure_Options, only : enumerationComponentTypeType
  private
  public :: radiusSolverPlausibilities, radiusSolverTasks, radiusSolver, solverGet, solverSet

  abstract interface
     double precision function solverGet(node)
       import treeNode
       type(treeNode), intent(inout) :: node
     end function solverGet
  end interface

  abstract interface
     subroutine solverSet(node,value)
       import treeNode
       type            (treeNode), intent(inout) :: node
       double precision          , intent(in   ) :: value
     end subroutine solverSet
  end interface

  abstract interface
     subroutine radiusSolver(node,component,specificAngularMomentum,radiusGet,radiusSet,velocityGet,velocitySet)
       import treeNode, enumerationComponentTypeType, solverGet, solverSet
       type            (treeNode                    ), intent(inout)          :: node
       type            (enumerationComponentTypeType), intent(in   )          :: component
       double precision                              , intent(in   )          :: specificAngularMomentum
       procedure       (solverGet                   ), intent(in   ), pointer :: radiusGet              , velocityGet
       procedure       (solverSet                   ), intent(in   ), pointer :: radiusSet              , velocitySet
     end subroutine radiusSolver
  end interface

contains

  subroutine radiusSolverPlausibilities(node)
    !!{
    Determine if the \mono{node} is physically plausible and solvable for galactic structure calculations.
    !!}
    implicit none
    type(treeNode), intent(inout) :: node
    
    !![
    <eventHookStatic name="radiusSolverPlausibility">
      <callWith>node</callWith>
    </eventHookStatic>
    !!]
    return
  end subroutine radiusSolverPlausibilities

  subroutine radiusSolverTasks(node,specificAngularMomentumRequired,radiusSolve)
    !!{
    !!}
    implicit none
    type            (treeNode                    ), intent(inout)          :: node
    logical                                       , intent(in   )          :: specificAngularMomentumRequired
    procedure       (radiusSolver                ), intent(in   ), pointer :: radiusSolve
    procedure       (solverGet                   )               , pointer :: radiusGet                      , velocityGet
    procedure       (solverSet                   )               , pointer :: radiusSet                      , velocitySet
    type            (enumerationComponentTypeType)                         :: component
    double precision                                                       :: specificAngularMomentum
    logical                                                                :: componentActive
    !$GLC attributes initialized :: radiusGet, velocityGet, radiusSet, velocitySet
    
    !![
    <eventHookStatic name="radiusSolverTask">
      <callWith>node,componentActive,component,specificAngularMomentumRequired,specificAngularMomentum,radiusGet,radiusSet,velocityGet,velocitySet</callWith>
      <onReturn>if (componentActive) call radiusSolve(node,component,specificAngularMomentum,radiusGet,radiusSet,velocityGet,velocitySet)</onReturn>
    </eventHookStatic>
    !!]
    return
  end subroutine radiusSolverTasks
  
end module Galactic_Structure_Radius_Solver_Utilities
