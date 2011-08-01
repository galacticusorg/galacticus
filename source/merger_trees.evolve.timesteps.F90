!! Copyright 2009, 2010, 2011 Andrew Benson <abenson@caltech.edu>
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


!% Contains a module which implements calculations of timesteps for merger tree evolution.

module Merger_Tree_Timesteps
  !% Implements calculations of timesteps for merger tree evolution.
  implicit none
  private
  public :: Time_Step_Get

contains

  double precision function Time_Step_Get(thisNode,evolveToTime,End_Of_Timestep_Task)
    !% Computes a suitable timestep over which to evolve a node in a tree.
    use Tree_Nodes
    use Merger_Trees
    use Input_Parameters
    use Galacticus_Error
    use Merger_Trees_Evolve_Timesteps_Template
    !# <include directive="timeStepsTask" type="moduleUse">
    include 'merger_trees.evolve.timesteps.moduleUse.inc'
    !# </include>
    implicit none
    type(treeNode),                           intent(inout), pointer :: thisNode
    double precision,                         intent(in)             :: evolveToTime
    procedure(),                              intent(out),   pointer :: End_Of_Timestep_Task
    procedure(End_Of_Timestep_Task_Template),                pointer :: End_Of_Timestep_Task_Internal

    ! Call the function to get the timestep.
    Time_Step_Get=evolveToTime-Tree_Node_Time(thisNode)
    End_Of_Timestep_Task_Internal => null()
    !# <include directive="timeStepsTask" type="code" action="subroutine">
    !#  <subroutineArgs>thisNode,Time_Step_Get,End_Of_Timestep_Task_Internal</subroutineArgs>
    include 'merger_trees.evolve.timesteps.inc'
    !# </include>
    End_Of_Timestep_Task => End_Of_Timestep_Task_Internal
    return
  end function Time_Step_Get

end module Merger_Tree_Timesteps
