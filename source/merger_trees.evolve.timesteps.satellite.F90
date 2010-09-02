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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


!% Contains a module which implements a time-stepping criterion for merger tree evolution which stops evolution when a merger is
!% about to happen.

module Merger_Tree_Timesteps_Satellite
  private
  public :: Merger_Tree_Timestep_Satellite

  ! Flag indicating whether this module is initialized.
  logical :: mergerTimestepsInitialized=.false.

  ! Flag indicating if this module is limiting timesteps.
  logical :: limitTimesteps

contains

  !# <timeStepsTask>
  !#  <unitName>Merger_Tree_Timestep_Satellite</unitName>
  !# </timeStepsTask>
  subroutine Merger_Tree_Timestep_Satellite(thisNode,timeStep,End_Of_Timestep_Task)
    !% Determines the timestep to go to the time at which the node merges.
    use Tree_Nodes
    use Merger_Trees_Evolve_Timesteps_Template
    implicit none
    type(treeNode),                           intent(inout), pointer :: thisNode
    procedure(End_Of_Timestep_Task_Template), intent(inout), pointer :: End_Of_Timestep_Task
    double precision,                         intent(inout)          :: timeStep
    double precision                                                 :: timeUntilMerging

    ! Initialize the module.
    !$omp critical (Merger_Tree_Timestep_Satellite_Initialize)
    if (.not.mergerTimestepsInitialized) then
       ! Check that the merge time property exists.
       limitTimesteps=associated(Tree_Node_Satellite_Merge_Time)

       ! Flag that the module is initialized.
       mergerTimestepsInitialized=.true.
    end if 
    !$omp end critical (Merger_Tree_Timestep_Satellite_Initialize)

    ! Exit if we are not limiting timesteps.
    if (.not.limitTimesteps) return

    ! Get the time until this node merges.
    timeUntilMerging=Tree_Node_Satellite_Merge_Time(thisNode)

    ! If time is negative, implies this is not a satellite, so return.
    if (timeUntilMerging < 0.0d0) return

    ! Set return value if our timestep is smaller than current one.
    if (timeUntilMerging <= timeStep) then
       timeStep=timeUntilMerging
       End_Of_Timestep_Task => Satellite_Merger_Process
    end if
    return
  end subroutine Merger_Tree_Timestep_Satellite

  subroutine Satellite_Merger_Process(thisTree,thisNode)
    !% Process a satellite node which has undergone a merger with its host node.
    use Merger_Trees
    use Tree_Nodes
    !# <include directive="satelliteMergerTask" type="moduleUse">
    include 'merger_trees.evolve.timesteps.satellite.moduleUse.inc'
    !# </include>
    implicit none
    type(mergerTree), intent(in)             :: thisTree
    type(treeNode),   intent(inout), pointer :: thisNode

    ! Allow arbitrary routines to process the merger.
    !# <include directive="satelliteMergerTask" type="code" action="subroutine">
    !#  <subroutineArgs>thisNode</subroutineArgs>
    include 'merger_trees.evolve.timesteps.satellite.inc'
    !# </include>

    ! Finally remove the satellite node from the host and destroy it.
    call thisNode%removeFromHost()
    call thisNode%destroy
    return
  end subroutine Satellite_Merger_Process

end module Merger_Tree_Timesteps_Satellite
