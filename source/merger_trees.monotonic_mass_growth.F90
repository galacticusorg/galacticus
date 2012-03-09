!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@caltech.edu>
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


!% Contains a module which enforces monotonic mass growth along merger tree branches.

module Merger_Trees_Monotonic_Mass_Growth
  !% Enforces monotonic mass growth along merger tree branches.
  implicit none
  private
  public :: Merger_Tree_Monotonic_Mass_Growth
  
  ! Flag indicating if module is initialized.
  logical          :: monotonicGrowthModuleInitialized=.false.
 
  ! Flag indicating if pruning is required.
  logical          :: mergerTreeEnforceMonotonicGrowth
contains

  !# <mergerTreePreEvolveTask>
  !#   <unitName>Merger_Tree_Monotonic_Mass_Growth</unitName>
  !# </mergerTreePreEvolveTask>
  subroutine Merger_Tree_Monotonic_Mass_Growth(thisTree)
    !% Enforce monotonic mass growth along branches of {\tt thisTree}.
    use Merger_Trees
    use Tree_Nodes
    use Input_Parameters
    implicit none
    type(mergerTree), intent(in) :: thisTree
    type(treeNode),   pointer    :: thisNode,progenitorNode
    logical                      :: didModifyTree
    double precision             :: progenitorMass

    ! Check if module is initialized.
    !$omp critical (Merger_Tree_Monotonic_Mass_Growth_Initialize)
    if (.not.monotonicGrowthModuleInitialized) then
       ! Get parameter specifying if monotonic growth should be enforced.
       !@ <inputParameter>
       !@   <name>mergerTreeEnforceMonotonicGrowth</name>
       !@   <defaultValue>false</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Specifies whether or not enforce monotonic mass growth along the branches of merger trees.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeEnforceMonotonicGrowth',mergerTreeEnforceMonotonicGrowth,defaultValue=.false.)
       ! Flag that module is initialized.
       monotonicGrowthModuleInitialized=.true.
    end if
    !$omp end critical (Merger_Tree_Monotonic_Mass_Growth_Initialize)

    ! Enforce monotonic growth if necessary.
    if (mergerTreeEnforceMonotonicGrowth) then
       didModifyTree=.true.
       do while (didModifyTree)
          didModifyTree=.false.
          ! Get root node of the tree.       
          thisNode => thisTree%baseNode
          ! Walk the tree.
          do while (associated(thisNode))
             ! Find nodes that have children.
             if (associated(thisNode%childNode)) then
                ! Find the mass of all progenitor nodes.
                progenitorMass =  0.0d0
                progenitorNode => thisNode%childNode
                do while (associated(progenitorNode))
                   progenitorMass =  progenitorMass+Tree_Node_Mass(progenitorNode)
                   progenitorNode => progenitorNode%siblingNode
                end do
                ! Find nodes which are less massive than the sum of their progenitors.
                if (Tree_Node_Mass(thisNode) < progenitorMass) then
                   call Tree_Node_Mass_Set(thisNode,progenitorMass)
                   didModifyTree=.true.
                end if
             end if
             ! Walk to the next node in the tree.             
             call thisNode%walkTree(thisNode)
          end do
        end do
     end if
    return
  end subroutine Merger_Tree_Monotonic_Mass_Growth
  
end module Merger_Trees_Monotonic_Mass_Growth
