!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements calculations of timesteps for merger tree evolution.

module Merger_Tree_Timesteps
  !% Implements calculations of timesteps for merger tree evolution.
  implicit none
  private
  public :: Time_Step_Get

contains

  double precision function Time_Step_Get(thisNode,evolveToTime,End_Of_Timestep_Task,report,lockNode,lockType)
    !% Computes a suitable timestep over which to evolve a node in a tree.
    use Merger_Trees_Evolve_Timesteps_Template
    use ISO_Varying_String
    !# <include directive="timeStepsTask" type="moduleUse">
    include 'merger_trees.evolve.timesteps.moduleUse.inc'
    !# </include>
    implicit none
    type            (treeNode                     ), intent(inout)          , pointer :: thisNode
    double precision                               , intent(in   )                    :: evolveToTime
    procedure       (End_Of_Timestep_Task_Template), intent(  out)          , pointer :: End_Of_Timestep_Task
    logical                                        , intent(in   )                    :: report
    type            (treeNode                     ), intent(inout), optional, pointer :: lockNode
    type            (varying_string               ), intent(inout), optional          :: lockType
    procedure       (End_Of_Timestep_Task_Template)                         , pointer :: End_Of_Timestep_Task_Internal
    class           (nodeComponentBasic           )                         , pointer :: thisBasicComponent

    ! Call the function to get the timestep.
    thisBasicComponent => thisNode%basic()
    Time_Step_Get=evolveToTime-thisBasicComponent%time()
    End_Of_Timestep_Task_Internal => null()
    !# <include directive="timeStepsTask" type="functionCall" functionType="void">
    !#  <functionArgs>thisNode,Time_Step_Get,End_Of_Timestep_Task_Internal,report,lockNode,lockType</functionArgs>
    include 'merger_trees.evolve.timesteps.inc'
    !# </include>
    End_Of_Timestep_Task => End_Of_Timestep_Task_Internal
    return
  end function Time_Step_Get

end module Merger_Tree_Timesteps
