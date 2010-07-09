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
