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

!!{RST
Implements a merger tree operator which validates the structural invariants of merger trees.
!!}

  !![
  <mergerTreeOperator name="mergerTreeOperatorValidateStructure" docformat="rst">
   <description>
   Provides a merger tree operator which validates the structural invariants of merger trees, reporting a fatal error if any are
   violated. The invariants checked (for every node of every tree in the forest) are:

   * the host tree pointer is set, and points to the tree being walked;
   * the node time is positive;
   * a node with no parent is the base node of its tree;
   * a node with a parent appears among that parent's children or satellites;
   * a child node does not exist at a time later than its parent;
   * every child (and satellite) of a node points back to that node as its parent;
   * the parent chain of every node terminates at the base node of the tree (in a bounded number of steps - which also detects
     cycles);
   * a node with a merge target appears on that target's list of mergees, and every mergee of a node points back to that node as
     its merge target;
   * a node whose primary progenitor is a clone of itself (i.e. has the same node index) has a distinct unique ID from that
     clone;
   * every event attached to a node which names a paired node has a matching (same ID) event attached to that paired node.

   This is intended primarily for use in testing and debugging (for example, validating the construction of merger trees read
   from file), but is cheap enough to enable in production models if desired.
   </description>
  </mergerTreeOperator>
  !!]
  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorValidateStructure
     !!{RST
     A merger tree operator class which validates the structural invariants of merger trees.
     !!}
     private
   contains
     procedure :: operatePreEvolution => validateStructureOperatePreEvolution
  end type mergerTreeOperatorValidateStructure

  interface mergerTreeOperatorValidateStructure
     !!{RST
     Constructors for the :galacticus-class:`mergerTreeOperatorValidateStructure` merger tree operator class.
     !!}
     module procedure validateStructureConstructorParameters
  end interface mergerTreeOperatorValidateStructure

contains

  function validateStructureConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`mergerTreeOperatorValidateStructure` merger tree operator class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(mergerTreeOperatorValidateStructure)                :: self
    type(inputParameters                    ), intent(inout) :: parameters

    self=mergerTreeOperatorValidateStructure()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function validateStructureConstructorParameters

  subroutine validateStructureOperatePreEvolution(self,tree)
    !!{RST
    Validate the structural invariants of a merger tree forest, reporting a fatal error listing all violations found.
    !!}
    use            :: Display            , only : displayMessage          , verbosityLevelInfo
    use            :: Error              , only : Error_Report
    use            :: Galacticus_Nodes   , only : mergerTree              , nodeComponentBasic, nodeEvent, treeNode
    use, intrinsic :: ISO_C_Binding      , only : c_size_t
    use            :: ISO_Varying_String , only : assignment(=)           , operator(//)      , var_str  , varying_string
    use            :: Merger_Tree_Walkers, only : mergerTreeWalkerAllNodes
    use            :: String_Handling    , only : operator(//)
    implicit none
    class  (mergerTreeOperatorValidateStructure), intent(inout), target  :: self
    type   (mergerTree                         ), intent(inout), target  :: tree
    type   (mergerTree                         )               , pointer :: treeCurrent
    type   (treeNode                           )               , pointer :: node          , nodeChild, &
         &                                                                  nodeWork
    class  (nodeComponentBasic                 )               , pointer :: basic         , basicParent
    class  (nodeEvent                          )               , pointer :: event         , eventPaired
    type   (mergerTreeWalkerAllNodes           )                         :: treeWalker
    integer(c_size_t                           )                         :: countNodes    , countSteps
    integer                                                              :: countViolations
    integer                                     , parameter              :: countViolationsReportMaximum=20
    logical                                                              :: found         , isChild  , &
         &                                                                  isSatellite
    type   (varying_string                     )                         :: message
    !$GLC attributes unused :: self

    countViolations=0
    message        ='tree structure validation failed:'
    ! Count nodes in the forest - used to bound all chain traversals below (so that corrupted, cyclic chains are detected rather
    ! than walked indefinitely).
    countNodes=0_c_size_t
    treeWalker=mergerTreeWalkerAllNodes(tree,spanForest=.true.)
    do while (treeWalker%next(nodeWork))
       countNodes=countNodes+1_c_size_t
    end do
    ! Iterate over trees in the forest.
    treeCurrent => tree
    do while (associated(treeCurrent))
       ! Validate that the tree has a base node.
       if (.not.associated(treeCurrent%nodeBase)) then
          call violate(var_str('tree [')//treeCurrent%index//'] has no base node')
          treeCurrent => treeCurrent%nextTree
          cycle
       end if
       ! Iterate over all nodes in this tree.
       treeWalker=mergerTreeWalkerAllNodes(treeCurrent,spanForest=.false.)
       do while (treeWalker%next(node))
          basic => node%basic()
          ! Validate that the host tree pointer is set correctly.
          if (.not.associated(node%hostTree,treeCurrent)) &
               & call violate(var_str('node [')//node%index()//'] has an incorrect host tree pointer')
          ! Validate that the node time is positive.
          if (basic%time() <= 0.0d0) &
               & call violate(var_str('node [')//node%index()//'] has a non-positive time')
          ! Validate parentage.
          if (associated(node%parent)) then
             ! Validate that this node appears among its parent's children or satellites.
             isChild    =chainContains(node%parent%firstChild    ,node)
             isSatellite=chainContains(node%parent%firstSatellite,node)
             if (.not.(isChild.or.isSatellite)) &
                  & call violate(var_str('node [')//node%index()//'] does not appear among the children or satellites of its parent ['//node%parent%index()//']')
             ! Validate that a child node does not exist later than its parent.
             if (isChild) then
                basicParent => node%parent%basic()
                if (basic%time() > basicParent%time()) &
                     & call violate(var_str('node [')//node%index()//'] exists later than its parent ['//node%parent%index()//']')
             end if
             ! Validate that the parent chain terminates at the base node of the tree.
             nodeWork   => node
             countSteps =  0_c_size_t
             do while (associated(nodeWork%parent).and.countSteps <= countNodes)
                nodeWork   => nodeWork%parent
                countSteps =  countSteps+1_c_size_t
             end do
             if (countSteps > countNodes) then
                call violate(var_str('the parent chain of node [')//node%index()//'] does not terminate - a cycle exists')
             else if (.not.associated(nodeWork,treeCurrent%nodeBase)) then
                call violate(var_str('the parent chain of node [')//node%index()//'] does not reach the base node of its tree')
             end if
          else
             ! Validate that a parent-less node is the base node of its tree.
             if (.not.associated(node,treeCurrent%nodeBase)) &
                  & call violate(var_str('node [')//node%index()//'] has no parent but is not the base node of its tree')
          end if
          ! Validate that children and satellites point back to this node as their parent.
          nodeChild  => node%firstChild
          countSteps =  0_c_size_t
          do while (associated(nodeChild).and.countSteps <= countNodes)
             if (.not.associated(nodeChild%parent,node)) &
                  & call violate(var_str('child [')//nodeChild%index()//'] of node ['//node%index()//'] does not point back to it as parent')
             nodeChild  => nodeChild%sibling
             countSteps =  countSteps+1_c_size_t
          end do
          if (countSteps > countNodes) call violate(var_str('the child chain of node [')//node%index()//'] does not terminate - a cycle exists')
          nodeChild  => node%firstSatellite
          countSteps =  0_c_size_t
          do while (associated(nodeChild).and.countSteps <= countNodes)
             if (.not.associated(nodeChild%parent,node)) &
                  & call violate(var_str('satellite [')//nodeChild%index()//'] of node ['//node%index()//'] does not point back to it as parent')
             nodeChild  => nodeChild%sibling
             countSteps =  countSteps+1_c_size_t
          end do
          if (countSteps > countNodes) call violate(var_str('the satellite chain of node [')//node%index()//'] does not terminate - a cycle exists')
          ! Validate merge target/mergee consistency.
          if (associated(node%mergeTarget)) then
             nodeChild  => node%mergeTarget%firstMergee
             countSteps =  0_c_size_t
             found      =  .false.
             do while (associated(nodeChild).and.countSteps <= countNodes)
                if (associated(nodeChild,node)) then
                   found=.true.
                   exit
                end if
                nodeChild  => nodeChild%siblingMergee
                countSteps =  countSteps+1_c_size_t
             end do
             if (.not.found) &
                  & call violate(var_str('node [')//node%index()//'] does not appear on the mergee list of its merge target ['//node%mergeTarget%index()//']')
          end if
          nodeChild  => node%firstMergee
          countSteps =  0_c_size_t
          do while (associated(nodeChild).and.countSteps <= countNodes)
             if (.not.associated(nodeChild%mergeTarget,node)) &
                  & call violate(var_str('mergee [')//nodeChild%index()//'] of node ['//node%index()//'] does not point back to it as merge target')
             nodeChild  => nodeChild%siblingMergee
             countSteps =  countSteps+1_c_size_t
          end do
          if (countSteps > countNodes) call violate(var_str('the mergee chain of node [')//node%index()//'] does not terminate - a cycle exists')
          ! Validate that a cloned primary progenitor has a distinct unique ID.
          if (associated(node%firstChild)) then
             if (node%firstChild%index() == node%index() .and. node%firstChild%uniqueID() == node%uniqueID()) &
                  & call violate(var_str('node [')//node%index()//'] and its clone share a unique ID')
          end if
          ! Validate that paired events are attached to both nodes of the pair.
          event => node%event
          do while (associated(event))
             if (associated(event%node)) then
                eventPaired => event%node%event
                found       =  .false.
                do while (associated(eventPaired))
                   if (eventPaired%ID == event%ID) then
                      found=.true.
                      exit
                   end if
                   eventPaired => eventPaired%next
                end do
                if (.not.found) &
                     & call violate(var_str('event [')//event%ID//'] on node ['//node%index()//'] has no paired event on node ['//event%node%index()//']')
             end if
             event => event%next
          end do
       end do
       ! Move to the next tree in the forest.
       treeCurrent => treeCurrent%nextTree
    end do
    ! Report.
    if (countViolations > 0) then
       if (countViolations > countViolationsReportMaximum) &
            & message=message//char(10)//'   [plus '//(countViolations-countViolationsReportMaximum)//' further violations]'
       call Error_Report(message//{introspection:location})
    else
       message='tree ['
       message=message//tree%index//'] structure validated: '//countNodes//' nodes'
       call displayMessage(message,verbosityLevelInfo)
    end if
    return

  contains

    subroutine violate(description)
      !!{RST
      Record a single violation of the tree structure invariants.
      !!}
      implicit none
      type(varying_string), intent(in   ) :: description

      countViolations=countViolations+1
      if (countViolations <= countViolationsReportMaximum) message=message//char(10)//'   '//description
      return
    end subroutine violate

    function chainContains(chainFirst,nodeSought) result(contained)
      !!{RST
      Return true if ``nodeSought`` appears in the sibling chain beginning at ``chainFirst``, walking at most the number of
      nodes in the forest.
      !!}
      implicit none
      logical                                        :: contained
      type   (treeNode), pointer     , intent(in   ) :: chainFirst
      type   (treeNode), pointer     , intent(in   ) :: nodeSought
      type   (treeNode), pointer                     :: nodeChain
      integer(c_size_t)                              :: countStepsChain

      contained       =  .false.
      nodeChain       => chainFirst
      countStepsChain =  0_c_size_t
      do while (associated(nodeChain).and.countStepsChain <= countNodes)
         if (associated(nodeChain,nodeSought)) then
            contained=.true.
            return
         end if
         nodeChain       => nodeChain%sibling
         countStepsChain =  countStepsChain+1_c_size_t
      end do
      return
    end function chainContains

  end subroutine validateStructureOperatePreEvolution
