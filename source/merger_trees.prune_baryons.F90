!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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

!% Contains a module which prunes branches which can contain no baryons from merger trees.

module Merger_Trees_Prune_Baryons
  !% Prunes branches which can contain no baryons from merger trees.
  implicit none
  private
  public :: Merger_Tree_Prune_Baryons

  ! Flag indicating if module is initialized.
  logical :: pruneBaryonsModuleInitialized=.false.

  ! Flag indicating if pruning is required.
  logical :: mergerTreePruneBaryons

contains

  !# <mergerTreePreEvolveTask>
  !#   <unitName>Merger_Tree_Prune_Baryons</unitName>
  !# </mergerTreePreEvolveTask>
  subroutine Merger_Tree_Prune_Baryons(tree)
    !% Prune branches from {\normalfont \ttfamily tree}.
    use Galacticus_Nodes
    use Merger_Trees_Pruning_Utilities
    use Input_Parameters
    use Accretion_Halos
    use Virial_Density_Contrast
    implicit none
    type   (mergerTree                ), intent(in   ), target :: tree
    type   (treeNode                  ), pointer               :: nodeNext              , nodePrevious, &
         &                                                        node
    class  (nodeComponentBasic        ), pointer               :: basic
    type   (mergerTree                ), pointer               :: treeCurrent
    class  (accretionHaloClass        ), pointer               :: accretionHalo_
    class  (virialDensityContrastClass), pointer               :: virialDensityContrast_
    logical                                                    :: didPruning
    double precision                                           :: densityContrast
    
    ! Check if module is initialized.
    if (.not.pruneBaryonsModuleInitialized) then
       !$omp critical (Merger_Tree_Prune_Baryons_Initialize)
       if (.not.pruneBaryonsModuleInitialized) then
          ! Get parameter specifying if pruning is required.
          !# <inputParameter>
          !#   <name>mergerTreePruneBaryons</name>
          !#   <cardinality>1</cardinality>
          !#   <defaultValue>.false.</defaultValue>
          !#   <description>Specifies whether or not to prune merger trees prior to evolution.</description>
          !#   <source>globalParameters</source>
          !#   <type>boolean</type>
          !# </inputParameter>
          ! Flag that module is initialized.
          pruneBaryonsModuleInitialized=.true.
       end if
       !$omp end critical (Merger_Tree_Prune_Baryons_Initialize)
    end if
    ! Prune tree if necessary.
    if (mergerTreePruneBaryons) then
       !# <workaround type="unknown">
       !#  <description>
       !#    When using the percolation virial density contrast type we run into problems with initialization if we do not
       !#    explicitly cause the percolation tables to be initialized here. The problem is a segmentation fault - I have not
       !#    understood why this happens.
       !#  </description>
       ! Compute a density contrast.
       node                   => tree                  %baseNode
       basic                  => node                  %basic          (                         )
       virialDensityContrast_ => virialDensityContrast                 (                         )
       densityContrast        =  virialDensityContrast_%densityContrast(basic%mass(),basic%time())
       !# </workaround>
       ! Get required objects.
       accretionHalo_ => accretionHalo()
       ! Iterate over trees.
       treeCurrent => tree
       do while (associated(treeCurrent))
          didPruning=.true.
          do while (didPruning)
             didPruning=.false.
             ! Get root node of the tree.
             node  => treeCurrent%baseNode
             basic => node       %basic   ()
             ! Check if the tree has any baryons.
             if (.not.accretionHalo_%branchHasBaryons(node)) then
                ! Entire tree can be pruned. Destroy all but this base node. (Leaving just
                ! the base node makes the tree inert - i.e. it can not do anything.)
                node => node%firstChild
                do while (associated(node))
                   nodeNext => node%sibling
                   call node%destroyBranch()
                   deallocate(node)
                   node => nodeNext
                end do
             else
                ! Walk the tree, pruning branches.
                do while (associated(node))
                   basic => node%basic()
                   ! Record the parent node to which we will return.
                   nodePrevious => node%parent
                   if (.not.accretionHalo_%branchHasBaryons(node).and.node%uniqueID() /= nodePrevious%uniqueID()) then
                      didPruning=.true.
                      ! Decouple from other nodes.
                      call Merger_Tree_Prune_Unlink_Parent(node,nodePrevious,.not.accretionHalo_%branchHasBaryons(nodePrevious),.true.)
                      ! Clean the branch.
                      call Merger_Tree_Prune_Clean_Branch(node)
                      ! Destroy the branch.
                      call node%destroyBranch()
                      deallocate(node)
                      ! Return to parent node.
                      node => nodePrevious
                   end if
                   node => node%walkTree()
                end do
             end if
          end do
          ! Move to the next tree.
          treeCurrent => treeCurrent%nextTree
       end do
       ! Uniqueify nodes.
       call Merger_Tree_Prune_Uniqueify_IDs(tree)
    end if
    return
  end subroutine Merger_Tree_Prune_Baryons

end module Merger_Trees_Prune_Baryons
