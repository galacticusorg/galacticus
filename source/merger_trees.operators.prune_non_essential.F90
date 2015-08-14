!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015 Andrew Benson <abenson@obs.carnegiescience.edu>
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

  !% Contains a module which implements a merger tree operator which prunes all branches which do not contain an ``essential''
  !% node.

  use Kind_Numbers
  
  !# <mergerTreeOperator name="mergerTreeOperatorPruneNonEssential">
  !#  <description>
  !#   Contains a module which implements a merger tree operator which prunes branches that do not directly influence an
  !#   ``essential'' node. Any branch which does not connect to the branch into which the node identified by ID {\normalfont
  !#   \ttfamily [essentialNodeI]} descends by time {\normalfont \ttfamily essetialNodeTime]} will be pruned.
  !#  </description>
  !# </mergerTreeOperator>
  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorPruneNonEssential
     !% A merger tree operator class which prunes branches which do not contain an ``essential'' node.
     private
     integer         (kind=kind_int8) :: essentialNodeID
     double precision                 :: essentialNodeTime
   contains
     final     ::            pruneNonEssentialDestructor
     procedure :: operate => pruneNonEssentialOperate
  end type mergerTreeOperatorPruneNonEssential
  
  interface mergerTreeOperatorPruneNonEssential
     !% Constructors for the prune-non-essential merger tree operator class.
     module procedure pruneNonEssentialConstructorParameters
     module procedure pruneNonEssentialConstructorInternal
  end interface mergerTreeOperatorPruneNonEssential

contains

  function pruneNonEssentialConstructorParameters(parameters)
    !% Constructor for the prune-non-essential merger tree operator class which takes a parameter set as input.
    implicit none
    type   (mergerTreeOperatorPruneNonEssential)                :: pruneNonEssentialConstructorParameters
    type   (inputParameters                    ), intent(in   ) :: parameters
    !# <inputParameterList label="allowedParameterNames" />

    call parameters%checkParameters(allowedParameterNames)
    !# <inputParameter>
    !#   <name>essentialNodeID</name>
    !#   <source>parameters</source>
    !#   <variable>pruneNonEssentialConstructorParameters%essentialNodeID</variable>
    !#   <description>ID of the essential node to avoid pruning.</description>
    !#   <type>integer</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>essentialNodeTime</name>
    !#   <source>parameters</source>
    !#   <variable>pruneNonEssentialConstructorParameters%essentialNodeTime</variable>
    !#   <description>Time of the essential node to avoid pruning.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    return
  end function pruneNonEssentialConstructorParameters

  function pruneNonEssentialConstructorInternal(essentialNodeID,essentialNodeTime)
    !% Internal constructor for the prune-non-essential merger tree operator class.
    use Galacticus_Error
    implicit none
    type            (mergerTreeOperatorPruneNonEssential)                :: pruneNonEssentialConstructorInternal
    integer                                              , intent(in   ) :: essentialNodeID
    double precision                                     , intent(in   ) :: essentialNodeTime

    pruneNonEssentialConstructorInternal%essentialNodeID  =essentialNodeID
    pruneNonEssentialConstructorInternal%essentialNodeTime=essentialNodeTime
   return
  end function pruneNonEssentialConstructorInternal

  elemental subroutine pruneNonEssentialDestructor(self)
    !% Destructor for the merger tree operator function class.
    implicit none
    type(mergerTreeOperatorPruneNonEssential), intent(inout) :: self

    ! Nothing to do.
    return
  end subroutine pruneNonEssentialDestructor

  subroutine pruneNonEssentialOperate(self,tree)
    !% Perform a prune-non-essential operation on a merger tree.
    use Merger_Trees_Pruning_Utilities
    implicit none
    class  (mergerTreeOperatorPruneNonEssential), intent(inout)          :: self
    type   (mergerTree                         ), intent(inout), target  :: tree
    type   (mergerTree                         )               , pointer :: treeCurrent
    type   (treeNode                           )               , pointer :: node          , nodeEssential, &
         &                                                                  nodePrevious
    class  (nodeComponentBasic                 )               , pointer :: basic
    
    ! Iterate over trees.
    treeCurrent => tree
    do while (associated(treeCurrent))
       ! Find the essential node.
       nodeEssential => treeCurrent%getNode(self%essentialNodeID)
       if (associated(nodeEssential)) then
          ! Trace the essential node to the required time.
          basic => nodeEssential%basic()
          do while (basic%time() < self%essentialNodeTime .and. associated(nodeEssential%parent))
             nodeEssential => nodeEssential%parent
             basic         => nodeEssential%basic ()
          end do
          ! Get root node of the tree.
          node => treeCurrent%baseNode
          ! Walk the tree, pruning branches.
          do while (associated(node))
             ! Record the parent node to which we will return.
             nodePrevious => node%parent
             if     (                                                &
                  &   associated(nodePrevious)                       &
                  &  .and.                                           &
                  &   .not.                                          &
                  &    (                                             &
                  &      node         %isProgenitorOf(nodeEssential) &
                  &     .or.                                         &
                  &      nodeEssential%isProgenitorOf(node         ) &
                  &    )                                             &
                  & ) then
                ! Decouple from other nodes.                
                call Merger_Tree_Prune_Unlink_Parent(                                                    &
                     &                               node                                              , &
                     &                               nodePrevious                                      , &
                     &                               .not.                                               &
                     &                                    (                                              &
                     &                                      nodePrevious %isProgenitorOf(nodeEssential)  &
                     &                                     .or.                                          &
                     &                                      nodeEssential%isProgenitorOf(nodePrevious )  &
                     &                                    )                                            , &
                     &                               .true.                                              &
                     &                              )
                ! Clean the branch.
                call Merger_Tree_Prune_Clean_Branch(node)
                ! Destroy the branch.
                call treeCurrent%destroyBranch(node)
                ! Return to parent node.
                node => nodePrevious
             end if
             node => node%walkTree()
          end do
       else
          ! Entire tree can be pruned. Destroy all but this base node. (Leaving just
          ! the base node makes the tree inert - i.e. it can not do anything.)
          node => treeCurrent%baseNode%firstChild
          do while (associated(node))
             nodePrevious => node%sibling
             call Merger_Tree_Prune_Clean_Branch(node)
            call treeCurrent%destroyBranch(node)
             node => nodePrevious
          end do
       end if
       ! Move to the next tree.
       treeCurrent => treeCurrent%nextTree
    end do
    ! Uniqueify nodes.
    call Merger_Tree_Prune_Uniqueify_IDs(tree)
    return
  end subroutine pruneNonEssentialOperate
