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

  !!{
  Implements a merger tree operator which restructures the tree onto a fixed grid of timesteps.
  !!}

  use :: Output_Times, only : outputTimesClass

  !![
  <mergerTreeOperator name="mergerTreeOperatorRegridTimes">
   <description>
    A merger tree operator class which will interpolate the merger tree structure onto a new array of timesteps. The timestep
    array is specified via an \refClass{outputTimesClass} object. Along each branch of the tree, new halos are inserted at
    times corresponding to the times in the resulting array. The masses of these nodes are linearly interpolated between the
    existing nodes on the branch. Once these new nodes have been added, all other nodes are removed from the tree\footnote{The
    base node of the tree is never removed, even if it does not lie on one of the times in the constructed array.} The
    processing is useful to construct representations of trees as they would be if only sparse time sampling were available. As
    such, it is useful for exploring how the number of snapshots in merger trees extracted from N-body simulations\index{merger
    tree!N-body} affects the properties of galaxies that form in them.
   </description>
  </mergerTreeOperator>
  !!]
  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorRegridTimes
     !!{
     A merger tree operator class which restructures the tree onto a fixed grid of timesteps.
     !!}
     private
     class           (outputTimesClass), pointer                   :: outputTimes_    => null()
     logical                                                       :: dumpTrees                , removeUngridded, &
          &                                                           propagateLabels
     double precision                                              :: snapTolerance
     double precision                  , allocatable, dimension(:) :: timeGrid
   contains
     final     ::                             regridTimesDestructor
     procedure :: operatePreInitialization => regridTimesOperatePreInitialization
  end type mergerTreeOperatorRegridTimes

  interface mergerTreeOperatorRegridTimes
     !!{
     Constructors for the regrid times merger tree operator class.
     !!}
     module procedure regridTimesConstructorParameters
     module procedure regridTimesConstructorInternal
  end interface mergerTreeOperatorRegridTimes

contains

  function regridTimesConstructorParameters(parameters) result(self)
    !!{
    Constructor for the regrid times merger tree operator class which takes a parameter set as input.
    !!}
    implicit none
    type            (mergerTreeOperatorRegridTimes)                :: self
    type            (inputParameters              ), intent(inout) :: parameters
    class           (outputTimesClass             ), pointer       :: outputTimes_
    logical                                                        :: dumpTrees      , removeUngridded, &
         &                                                            propagateLabels
    double precision                                               :: snapTolerance

    !![
    <inputParameter>
      <name>dumpTrees</name>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
      <description>Specifies whether or not to dump merger trees as they are regridded.</description>
    </inputParameter>
    <inputParameter>
      <name>removeUngridded</name>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
      <description>If true, remove nodes not at gridded times. Otherwise, leave them in place.</description>
    </inputParameter>
    <inputParameter>
      <name>propagateLabels</name>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
      <description>If true, any labels attached to progenitor nodes are propagated to newly inserted nodes.</description>
    </inputParameter>
    <inputParameter>
      <name>snapTolerance</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>The fractional tolerance used in deciding if a node should be snapped to a time on the grid.</description>
    </inputParameter>
    <objectBuilder class="outputTimes" name="outputTimes_" source="parameters"/>
    !!]
    self=mergerTreeOperatorRegridTimes(snapTolerance,dumpTrees,removeUngridded,propagateLabels,outputTimes_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="outputTimes_"/>
    !!]
    return
  end function regridTimesConstructorParameters

  function regridTimesConstructorInternal(snapTolerance,dumpTrees,removeUngridded,propagateLabels,outputTimes_) result(self)
    !!{
    Internal constructor for the regrid times merger tree operator class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (mergerTreeOperatorRegridTimes)                        :: self
    double precision                               , intent(in   )         :: snapTolerance
    logical                                        , intent(in   )         :: dumpTrees      , removeUngridded, &
         &                                                                    propagateLabels
    class           (outputTimesClass             ), intent(in   ), target :: outputTimes_
    integer         (c_size_t                     )                        :: i
    !![
    <constructorAssign variables="dumpTrees, snapTolerance, removeUngridded, propagateLabels, *outputTimes_"/>
    !!]

    ! Validate arguments.
    if (self%removeUngridded .and. self%outputTimes_%count() < 2_c_size_t) call Error_Report('2 or more output times are required'//{introspection:location})
    allocate(self%timeGrid(self%outputTimes_%count()))
    do i=1,self%outputTimes_%count()
       self%timeGrid(i)=self%outputTimes_%time(i)
    end do
    return
  end function regridTimesConstructorInternal

  subroutine regridTimesDestructor(self)
    !!{
    Destructor for the merger tree operator function class.
    !!}
    implicit none
    type(mergerTreeOperatorRegridTimes), intent(inout) :: self

    !![
    <objectDestructor name="self%outputTimes_"/>
    !!]
    return
  end subroutine regridTimesDestructor

  subroutine regridTimesOperatePreInitialization(self,tree)
    !!{
    Perform a regrid times operation on a merger tree.
    !!}
    use            :: Display                , only : displayIndent           , displayMessage               , displayUnindent       , verbosityLevelWorking, &
         &                                            displayMagenta          , displayReset
    use            :: Error                  , only : Error_Report            , Warn
    use            :: Galacticus_Nodes       , only : mergerTree              , nodeComponentBasic           , nodeComponentSatellite, nodeEvent            , &
          &                                           treeNode                , treeNodeList
    use, intrinsic :: ISO_C_Binding          , only : c_size_t
    use            :: ISO_Varying_String     , only : var_str
    use            :: Kind_Numbers           , only : kind_int8
    use            :: Merger_Tree_Walkers    , only : mergerTreeWalkerAllNodes, mergerTreeWalkerIsolatedNodes
    use            :: Merger_Trees_Dump      , only : Merger_Tree_Dump
    use            :: Nodes_Labels           , only : nodeLabelIsPresent      , nodeLabelSet                 , nodeLabelCount
    use            :: Numerical_Comparison   , only : Values_Agree
    use            :: Numerical_Interpolation, only : interpolator
    use            :: String_Handling        , only : operator(//)
    implicit none
    class           (mergerTreeOperatorRegridTimes), intent(inout), target                :: self
    type            (mergerTree                   ), intent(inout), target                :: tree
    type            (treeNode                     )                             , pointer :: nodeChild                       , mergee     , &
         &                                                                                   nodeSibling                     , node
    type            (treeNodeList                 ), allocatable  , dimension(:)          :: newNodes
    integer         (kind_int8                    ), allocatable  , dimension(:)          :: highlightNodes
    class           (nodeComponentBasic           )                             , pointer :: basicChild                      , basicParent, &
         &                                                                                   basic
    class           (nodeComponentSatellite       )                             , pointer :: mergeeSatellite
    type            (mergerTree                   )                             , pointer :: currentTree
    class           (nodeEvent                    )                             , pointer :: event                           , pairedEvent
    type            (mergerTreeWalkerAllNodes     )                                       :: treeWalkerAllNodes
    type            (mergerTreeWalkerIsolatedNodes)                                       :: treeWalkerIsolatedNodes
    logical                                        , save                                 :: mergeTargetWarningIssued=.false.
    logical                                                                               :: snap
    type            (interpolator                 )                                       :: interpolator_
    integer         (c_size_t                     )                                       :: iNow                            , iParent    , &
         &                                                                                   iTime                           , countNodes
    integer                                                                               :: allocErr                        , i
    double precision                                                                      :: massNow                         , massParent , &
         &                                                                                   timeNow                         , timeParent
    integer         (kind=kind_int8               )                                       :: firstNewNode                    , nodeIndex

    ! Build an interpolator.
    interpolator_=interpolator(self%timeGrid)
    ! Iterate over trees.
    currentTree => tree
    do while (associated(currentTree))
       call displayIndent(var_str('Regridding tree ')//currentTree%index,verbosityLevelWorking)
       ! Dump the unprocessed tree if required.
       if (self%dumpTrees) call Merger_Tree_Dump(                              &
            &                                    currentTree                 , &
            &                                    backgroundColor    ='white' , &
            &                                    nodeColor          ='black' , &
            &                                    highlightColor     ='black' , &
            &                                    edgeColor          ='black' , &
            &                                    nodeStyle          ='solid' , &
            &                                    highlightStyle     ='filled', &
            &                                    edgeStyle          ='solid' , &
            &                                    labelNodes         =.false. , &
            &                                    scaleNodesByLogMass=.true.  , &
            &                                    edgeLengthsToTimes =.true.    &
            &                                   )
       ! Iterate through to tree to:
       !  a) Find the current maximum node index in the tree, and;
       !  b) Snap halos to snapshot times if requested, and;
       !  c) Count the number of nodes in the tree.
       nodeIndex =0_kind_int8
       countNodes=0_c_size_t
       treeWalkerAllNodes=mergerTreeWalkerAllNodes(currentTree)
       do while (treeWalkerAllNodes%next(node))
          nodeIndex=max(nodeIndex,node%index())
          ! Count node.
          countNodes=countNodes+1_c_size_t
          ! Check for merge targets being set - these are not supported under the regridding transformation, so issue a warning.
          if (associated(node%mergeTarget).and..not.mergeTargetWarningIssued) then
             !$omp critical (mergeTargetWarning)
             if (.not.mergeTargetWarningIssued) then
                call Warn(                                                                                          &
                     &                                                                         displayMagenta(  )// &
                     &    'WARNING:'                                                         //displayReset  (  )// &
                     &    ' nodes in this tree have merge targets set'                       //char          (10)// &
                     &    '         this is not supported by the regridding operator'        //char          (10)// &
                     &    '         your tree may crash or deadlock'                         //char          (10)// &
                     &    '         to avoid this problem do not preset merge targets, e.g. '//char          (10)// &
                     &    '           <mergerTreeReadPresetMergerNodes value="false"/>'                             &
                     &   )
                mergeTargetWarningIssued=.true.
             end if
             !$omp end critical (mergeTargetWarning)
          end if
          ! Check if this node can be snapped to a grid time.
          if (self%snapTolerance > 0.0d0) then
             ! Get the basic component.
             basic   => node %basic()
             ! Get the time for this node.
             timeNow =  basic%time ()
             ! Find the closest time in the new time grid.
             iNow    =  interpolator_%locate(timeNow,closest=.true.)
             ! Test how close the node is to this time.
             snap    =  Values_Agree(timeNow,self%timeGrid(iNow),relTol=self%snapTolerance)
             ! Only snap if this node is closer to the output time than its parent or child.
             if (snap) then
                if (associated(node%parent    )) then
                   basicParent => node%parent    %basic()
                   snap        =   abs(timeNow-self%timeGrid(iNow)) < abs(basicParent%time()-self%timeGrid(iNow))      &
                        &         .or.                                                                                 &
                        &          .not.Values_Agree(basicParent%time(),self%timeGrid(iNow),relTol=self%snapTolerance)
                end if
                if (associated(node%firstChild)) then
                   basicChild  => node%firstChild%basic()
                   snap        =   abs(timeNow-self%timeGrid(iNow)) < abs(basicChild %time()-self%timeGrid(iNow))      &
                        &         .or.                                                                                 &
                        &          .not.Values_Agree(basicChild %time(),self%timeGrid(iNow),relTol=self%snapTolerance)
                end if
             end if
             ! If this node is to be snapped, do so now.
             if (snap) then
                ! Adjust the time of the node.
                call basic%timeSet(self%timeGrid(iNow))
                ! Check for mergees with their merge times set to match the time of this node.
                mergee => node%firstMergee
                do while (associated(mergee))
                   mergeeSatellite => mergee         %satellite    ()
                   ! Get the merge time for this mergee.
                   timeNow         =  mergeeSatellite%timeOfMerging()
                   ! Find the closest time in the new time grid.
                   iNow    =  interpolator_%locate(timeNow,closest=.true.)
                   if (Values_Agree(timeNow,self%timeGrid(iNow),relTol=self%snapTolerance)) &
                        & call mergeeSatellite%timeOfMergingSet(self%timeGrid(iNow))
                   mergee => mergee%siblingMergee
                end do
                ! Check for events with their event times set to match the time of this node.
                event => node%event
                do while (associated(event))
                   ! Get the merge time for this event.
                   timeNow=event%time
                   ! Find the closest time in the new time grid.
                   iNow   =interpolator_%locate(timeNow,closest=.true.)
                   if (Values_Agree(timeNow,self%timeGrid(iNow),relTol=self%snapTolerance)) then
                      event%time=self%timeGrid(iNow)
                      if (associated(event%node)) then
                         pairedEvent => event%node%event
                         do while (associated(pairedEvent))
                            if (pairedEvent%ID == event%ID) then
                               pairedEvent%time=self%timeGrid(iNow)
                               exit
                            end if
                            pairedEvent => pairedEvent%next
                         end do
                      end if
                   end if
                   event => event%next
                end do
             end if
          end if
       end do
       firstNewNode=nodeIndex+1
       call displayMessage(var_str('Tree contains ')//countNodes//' nodes prior to regridding',verbosityLevelWorking)
       ! Walk the tree, locating branches which intersect grid times.
       treeWalkerIsolatedNodes=mergerTreeWalkerIsolatedNodes(currentTree)
       do while (treeWalkerIsolatedNodes%next(node))
          ! Skip this node if it is the root node.
          if (associated(node%parent)) then
             basic       => node       %basic()
             basicParent => node%parent%basic()
             ! Get the time of this node and its parent.
             timeNow   =basic      %time()
             timeParent=basicParent%time()
             ! Locate these times in the list of grid times.
             if (size(self%timeGrid) > 1) then
                iNow   =interpolator_%locate(timeNow   )
                iParent=interpolator_%locate(timeParent)
             else
                iNow   =0
                iParent=0
                if (timeNow    >= self%timeGrid(1)) iNow   =1
                if (timeParent >  self%timeGrid(1)) iParent=1
             end if
             ! For nodes existing precisely at a grid time, ignore this grid point. (These are,
             ! typically, nodes which have been created at these points.)
             if (timeNow    == self%timeGrid(iNow   )) iNow   =iNow   +1
             if (timeParent == self%timeGrid(iParent)) iParent=iParent-1
             ! If the branch from node to parent spans one or more grid times, insert new nodes
             ! at those points.
             if (iParent > iNow) then
                ! Get masses of these halos.
                massNow   =basic      %mass()
                massParent=basicParent%mass()
                if (node%isPrimaryProgenitor()) then
                   ! Remove the mass in any non-primary progenitors - we don't want to include
                   ! their mass in the estimated mass growth rate of this node.
                   nodeChild => node%parent%firstChild%sibling
                   do while (associated(nodeChild))
                      basicChild => nodeChild%basic()
                      massParent =  massParent-basicChild%mass()
                      nodeChild  => nodeChild%sibling
                   end do
                   ! Do not let the parent mass decrease along the branch.
                   massParent=max(massParent,massNow)
                else
                   ! Halo is not the primary progenitor of its parent. Assume that its mass does
                   ! not grow further.
                   massParent=massNow
                end if
                ! Create new nodes.
                allocate(newNodes(iParent-iNow),stat=allocErr)
                if (allocErr/=0) call Error_Report('unable to allocate new nodes'//{introspection:location})
                do iTime=iNow+1,iParent
                   nodeIndex=nodeIndex+1_kind_int8
                   newNodes(iTime-iNow)%node => treeNode(hostTree=currentTree)
                   call newNodes(iTime-iNow)%node%indexSet(nodeIndex)
                end do
                ! Assign node properties and build links.
                do iTime=iNow+1,iParent
                   ! Assign a time and a mass
                   basic => newNodes(iTime-iNow)%node%basic(autoCreate=.true.)
                   call basic%timeSet(                              self%timeGrid(iTime)                              )
                   call basic%massSet(massNow+(massParent-massNow)*(self%timeGrid(iTime)-timeNow)/(timeParent-timeNow))
                   ! Link to child node.
                   if (iTime > iNow+1 ) newNodes(iTime-iNow)%node%firstChild => newNodes(iTime-iNow-1)%node
                   ! Link to parent node.
                   if (iTime < iParent) newNodes(iTime-iNow)%node%parent     => newNodes(iTime-iNow+1)%node
                   ! Set subsampling rate.
                   call newNodes(iTime-iNow)%node%subsamplingWeightSet(node%subsamplingWeight())
                   ! Propagate labels.
                   if (self%propagateLabels .and. nodeLabelCount() > 0) then
                      do i=1,nodeLabelCount()
                         if (nodeLabelIsPresent(i,node)) call nodeLabelSet(i,newNodes(iTime-iNow)%node)
                      end do
                   end if
                end do
                ! Link final node to the parent.
                newNodes(iParent-iNow)%node%parent  => node%parent
                ! Link final node sibling to current node sibling.
                newNodes(iParent-iNow)%node%sibling => node%sibling
                ! Link the parent to the final node.
                if (node%isPrimaryProgenitor()) then
                   ! Node is the main progenitor of its parent, so simply replace it with the
                   ! final node in our list.
                   node%parent%firstChild  => newNodes(iParent-iNow)%node
                else
                   ! Node is not the main progenitor of its parent, so find the child node that
                   ! has it as a sibling.
                   nodeChild => node%parent%firstChild
                   do while (.not.associated(nodeChild%sibling,node))
                      nodeChild => nodeChild%sibling
                   end do
                   nodeChild%sibling => newNodes(iParent-iNow)%node
                end if
                ! Link the child of the first node to the node being processed.
                newNodes(1)%node%firstChild  => node
                ! Nullify any sibling of the node being processed.
                node%sibling => null()
                ! Link the parent of the node being processed to the first node of the list.
                node%parent  => newNodes(1)%node
                ! Erase the node list.
                deallocate(newNodes)
             end if
          end if
       end do
       ! Dump the intermediate tree if required.
       if (self%dumpTrees) then
          allocate(highlightNodes(nodeIndex-firstNewNode+2))
          highlightNodes(1)=currentTree%nodeBase%index()
          do nodeIndex=1,nodeIndex-firstNewNode+1
             highlightNodes(nodeIndex+1)=firstNewNode+nodeIndex-1
          end do
          call Merger_Tree_Dump(                                    &
               &                currentTree                       , &
               &                highlightNodes     =highlightNodes, &
               &                backgroundColor    ='white'       , &
               &                nodeColor          ='black'       , &
               &                highlightColor     ='black'       , &
               &                edgeColor          ='black'       , &
               &                nodeStyle          ='solid'       , &
               &                highlightStyle     ='filled'      , &
               &                edgeStyle          ='dotted'      , &
               &                labelNodes         =.false.       , &
               &                scaleNodesByLogMass=.true.        , &
               &                edgeLengthsToTimes =.true.          &
               &               )
          deallocate(highlightNodes)
       end if
       ! Walk the tree removing nodes not at grid times.
       if (self%removeUngridded) then
          countNodes              =0_c_size_t
          treeWalkerIsolatedNodes=mergerTreeWalkerIsolatedNodes(currentTree)
          do while (treeWalkerIsolatedNodes%next(node))
             basic => node%basic()
             ! Get the time for this node.
             timeNow=basic%time()
             ! Find the closest time in the new time grid.
             iNow   =interpolator_%locate(timeNow,closest=.true.)
             ! If this node does not lie precisely on the grid then remove it.
             if (associated(node%parent) .and. timeNow /= self%timeGrid(iNow)) then
                if (node%isPrimaryProgenitor()) then
                   ! Handle primary progenitor nodes.
                   if (associated(node%firstChild)) then
                      ! Handle primary progenitors with children
                      nodeChild => node%firstChild
                      ! Assign all children a parent that is the parent of the current node.
                      do while (associated(nodeChild))
                         nodeChild%parent => node%parent
                         if (.not.associated(nodeChild%sibling)) then
                            nodeChild%sibling => node     %sibling
                            nodeChild         => null()
                         else
                            nodeChild         => nodeChild%sibling
                         end if
                      end do
                      ! Assign the current node's parent a child that is the child of the current node.
                      node%parent%firstChild => node%firstChild
                   else
                      ! Handle primary nodes with no children - simply make the parent's main
                      ! progenitor the sibling of the current node.
                      node%parent%firstChild => node%sibling
                   end if
                else
                   ! Handle non-primary nodes.
                   if (associated(node%firstChild)) then
                      ! Handle non-primary nodes with children.
                      ! Assign all children a parent that is the parent of the current node.
                      nodeChild => node%firstChild
                      do while (associated(nodeChild))
                         nodeChild%parent => node%parent
                         if (.not.associated(nodeChild%sibling)) then
                            nodeChild%sibling => node     %sibling
                            nodeChild         => null()
                         else
                            nodeChild         => nodeChild%sibling
                         end if
                      end do
                      ! Find which sibling points to the current node and link in the children of the current node.
                      nodeSibling => node%parent%firstChild
                      do while (.not.associated(nodeSibling%sibling,node))
                         nodeSibling => nodeSibling%sibling
                      end do
                      nodeSibling%sibling => node%firstChild
                   else
                      ! Handle non-primary nodes with no children - just snip it out of the sibling list.
                      nodeSibling => node%parent%firstChild
                      do while (.not.associated(nodeSibling%sibling,node))
                         nodeSibling => nodeSibling%sibling
                      end do
                      nodeSibling%sibling => node%sibling
                   end if
                end if
                ! Destroy the node.
                call node%destroy()
                deallocate(node)
                call treeWalkerIsolatedNodes%previous(node)
             else
                countNodes=countNodes+1_c_size_t
             end if
          end do
       else
          countNodes             =0_c_size_t
          treeWalkerIsolatedNodes=mergerTreeWalkerIsolatedNodes(currentTree)
          do while (treeWalkerIsolatedNodes%next(node))
             countNodes=countNodes+1_c_size_t
          end do
       end if
       call displayMessage(var_str('Tree contains ')//countNodes//' nodes after regridding',verbosityLevelWorking)
       ! Dump the processed tree if required.
       if (self%dumpTrees) call Merger_Tree_Dump(                              &
            &                                    currentTree                 , &
            &                                    backgroundColor    ='white' , &
            &                                    nodeColor          ='black' , &
            &                                    highlightColor     ='black' , &
            &                                    edgeColor          ='black' , &
            &                                    nodeStyle          ='solid' , &
            &                                    highlightStyle     ='filled', &
            &                                    edgeStyle          ='solid' , &
            &                                    labelNodes         =.false. , &
            &                                    scaleNodesByLogMass=.true.  , &
            &                                    edgeLengthsToTimes =.true.    &
            &                                   )
       call displayUnindent('Done',verbosityLevelWorking)
       ! Move to the next tree.
       currentTree => currentTree%nextTree
    end do
    return
  end subroutine regridTimesOperatePreInitialization
