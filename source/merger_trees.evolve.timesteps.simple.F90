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






!% Contains a module which implements a simple time-stepping criterion for merger tree evolution.

module Merger_Tree_Timesteps_Simple
  !% Implements a simple time-stepping criterion for merger tree evolution.
  private
  public :: Merger_Tree_Timestep_Simple

  ! Variable inidicating if module is initialized.
  logical          :: timestepSimpleInitialized=.false.

  ! Parameters controlling the size of timesteps.
  double precision :: timestepSimpleRelative,timestepSimpleAbsolute

contains

  !# <timeStepsTask>
  !#  <unitName>Merger_Tree_Timestep_Simple</unitName>
  !# </timeStepsTask>
  subroutine Merger_Tree_Timestep_Simple(thisNode,timeStep,End_Of_Timestep_Task)
    !% Determine a suitable timestep for {\tt thisNode} using the simple method. This simply selects the smaller of {\tt
    !% timestepSimpleAbsolute} and {\tt timestepSimpleRelative}$H^{-1}(t)$.
    use Tree_Nodes
    use Tree_Node_Methods
    use Input_Parameters
    use Cosmology_Functions
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    procedure(),      intent(inout), pointer :: End_Of_Timestep_Task
    double precision, intent(inout)          :: timeStep
    double precision                         :: time,expansionFactor,expansionTimescale,ourTimeStep

    !$omp critical (timestepSimpleInitialize)
    if (.not.timestepSimpleInitialized) then
       !@ <inputParameter>
       !@   <name>timestepSimpleRelative</name>
       !@   <defaultValue>0.1</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The maximum allowed relative change in time for a single step in the evolution of a node.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('timestepSimpleRelative',timestepSimpleRelative,defaultValue=0.1d0)
       !@ <inputParameter>
       !@   <name>timestepSimpleAbsolute</name>
       !@   <defaultValue>1</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The maximum allowed absolute change in time (in Gyr) for a single step in the evolution of a node.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('timestepSimpleAbsolute',timestepSimpleAbsolute,defaultValue=1.0d0)
       timestepSimpleInitialized=.true.
    end if
    !$omp end critical (timestepSimpleInitialize)

    ! Get current cosmic time.
    time=Tree_Node_Time(thisNode)

    ! Find current expansion timescale.
    expansionFactor=Expansion_Factor(time)
    expansionTimescale=1.0d0/Expansion_Rate(expansionFactor)

    ! Determine suitable timestep.
    ourTimeStep=min(timestepSimpleRelative*expansionTimescale,timestepSimpleAbsolute)

    ! Set return value if our timestep is smaller than current one.
    if (ourTimeStep < timeStep) timeStep=ourTimeStep

    return
  end subroutine Merger_Tree_Timestep_Simple

end module Merger_Tree_Timesteps_Simple
