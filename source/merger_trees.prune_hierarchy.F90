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


!% Contains a module which prunes hierarchy below a given depth in merger trees.

module Merger_Trees_Prune_Hierarchy
  !% Prunes hierarchy below a given depth in merger trees.
  private
  public :: Merger_Tree_Prune_Hierarchy
  
  ! Flag indicating if module is initialized.
  logical :: pruneHierarchyModuleInitialized=.false.
 
  ! Depth below which hierarchy should be pruned.
  integer :: mergerTreePruneHierarchyAtDepth
  
contains

  !# <mergerTreePreEvolveTask>
  !#   <unitName>Merger_Tree_Prune_Hierarchy</unitName>
  !# </mergerTreePreEvolveTask>
  subroutine Merger_Tree_Prune_Hierarchy(thisTree)
    !% Prune hierarchy from {\tt thisTree}.
    use Merger_Trees
    use Tree_Nodes
    use Input_Parameters
    implicit none
    type(mergerTree), intent(in) :: thisTree
    type(treeNode),   pointer    :: thisNode,nextNode,destroyNode,previousNode,workNode
    logical                      :: didPruning
    integer                      :: hierarchyDepth

    ! Check if module is initialized.
    !$omp critical (Merger_Tree_Prune_Hierarchy_Initialize)
    if (.not.pruneHierarchyModuleInitialized) then
       ! Get parameter specifying if pruning is required.
       !@ <inputParameter>
       !@   <name>mergerTreePruneHierarchyAtDepth</name>
       !@   <defaultValue>0</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The depth in the hierarchy at which to prune merger trees. (Zero indicates to not prune.)
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreePruneHierarchyAtDepth',mergerTreePruneHierarchyAtDepth,defaultValue=0)
       ! Flag that module is initialized.
       pruneHierarchyModuleInitialized=.true.
    end if
    !$omp end critical (Merger_Tree_Prune_Hierarchy_Initialize)

    ! Prune tree if necessary.
    if (mergerTreePruneHierarchyAtDepth > 0) then
       didPruning=.true.
       do while (didPruning)
          didPruning=.false.
          ! Get root node of the tree.       
          thisNode => thisTree%baseNode
          ! Walk the tree, pruning hierarchy.
          do while (associated(thisNode))
             ! Record the parent node to which we will return.
             previousNode => thisNode%parentNode

             ! Find the depth in the hierarchy.
             hierarchyDepth=0
             workNode => thisNode
             do while (associated(workNode))
                ! Increment hierarchy depth if this node is not the main progenitor.
                if (.not.workNode%isPrimaryProgenitor().and.associated(workNode%parentNode)) hierarchyDepth=hierarchyDepth+1
                workNode => workNode%parentNode
             end do

             ! Prune if this node is sufficiently deep in the hierarchy.
             if (hierarchyDepth >= mergerTreePruneHierarchyAtDepth) then
                didPruning=.true.
                ! Decouple from other nodes.
                if (thisNode%isPrimaryProgenitorOf(previousNode)) then
                   previousNode%childNode => thisNode%siblingNode
                else
                   nextNode => previousNode%childNode
                   do while (.not.associated(nextNode%siblingNode,thisNode))
                      nextNode => nextNode%siblingNode
                   end do
                   nextNode%siblingNode => thisNode%siblingNode
                end if
                ! Should call the "destroyBranch" method on "thisTree" here, but seems to cause a compiler crash under gFortran
                ! v4.4. So, instead, do the branch destroy manually.
                nextNode => thisNode
                do while (associated(nextNode))
                   destroyNode => nextNode
                   call destroyNode%walkBranch(thisNode,nextNode)
                   call destroyNode%destroy
                end do
                ! Return to parent node.
                thisNode => previousNode
             end if
             call thisNode%walkTree(thisNode)
          end do
       end do
    end if

    return
  end subroutine Merger_Tree_Prune_Hierarchy
  
end module Merger_Trees_Prune_Hierarchy
