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

!+    Contributions to this file made by: Andrew Benson, Claude.

!!{RST
Contains a program which tests the time-limited isolated nodes branch merger tree walker.
!!}

program Tests_Merger_Tree_Walker_Branch_Time_Limited
  !!{RST
  Tests that the :galacticus-class:`mergerTreeWalkerIsolatedNodesBranch` merger tree walker, when given an earliest time, visits
  every node in the branch that exists at or after that time, including the progenitors of non-primary children whose primary
  sibling precedes the earliest time.
  !!}
  use :: Display                   , only : displayVerbositySet         , verbosityLevelStandard
  use :: Events_Hooks              , only : eventsHooksInitialize
  use :: Functions_Global_Utilities, only : Functions_Global_Set
  use :: Galacticus_Nodes          , only : mergerTree                  , nodeClassHierarchyFinalize      , nodeClassHierarchyInitialize       , nodeComponentBasic, &
       &                                    treeNode                    , treeNodeList
  use :: Input_Parameters          , only : inputParameters
  use :: Kind_Numbers              , only : kind_int8
  use :: Merger_Tree_Walkers       , only : mergerTreeWalkerIsolatedNodesBranch
  use :: Node_Components           , only : Node_Components_Initialize  , Node_Components_Thread_Initialize, Node_Components_Thread_Uninitialize, Node_Components_Uninitialize
  use :: Unit_Tests                , only : Assert                      , Unit_Tests_Begin_Group          , Unit_Tests_End_Group               , Unit_Tests_Finish
  implicit none
  type            (mergerTree                         ), pointer      :: tree
  type            (treeNodeList                       ), dimension(6) :: nodes
  type            (treeNode                           ), pointer      :: node
  class           (nodeComponentBasic                 ), pointer      :: basic
  type            (mergerTreeWalkerIsolatedNodesBranch)               :: treeWalker
  type            (inputParameters                    )               :: parameters
  ! Times (in Gyr) of the five nodes. Node 1 is the branch head. Its primary progenitor (node 2) is early, while its non-primary
  ! progenitor (node 3) is later, and has a single progenitor of its own (node 4) which precedes the earliest time. Node 5 is a
  ! second, early, non-primary progenitor of node 1, and node 6 a third, later, one (so that the walk must skip an early sibling to
  ! reach a later one).
  double precision                                     , dimension(6) :: times       =[10.0d0,2.0d0,6.0d0,3.0d0,1.0d0,8.0d0]
  double precision                                     , parameter    :: timeEarliest=4.0d0
  integer         (kind_int8                          ), dimension(6) :: visitedLimited  , visitedAll
  integer         (kind_int8                          )               :: i
  integer                                                             :: countLimited    , countAll

  ! Initialize.
  call displayVerbositySet              (verbosityLevelStandard)
  call eventsHooksInitialize            (                      )
  call Functions_Global_Set             (                      )
  parameters=inputParameters()
  call nodeClassHierarchyInitialize     (parameters            )
  call Node_Components_Initialize       (parameters            )
  call Node_Components_Thread_Initialize(parameters            )
  call Unit_Tests_Begin_Group("Time-limited isolated nodes branch walker")
  ! Build the tree:
  !
  !   1 (t=10) -> children: 2 (t=2, primary), 3 (t=6), 5 (t=1), 6 (t=8)
  !   3 (t=6)  -> child   : 4 (t=3)
  allocate(tree)
  do i=1,6
     nodes(i)%node => treeNode()
     call nodes(i)%node%indexSet(i)
     basic => nodes(i)%node%basic(autoCreate=.true.)
     call basic%timeSet(times(i))
     call basic%massSet(1.0d10/dble(i))
  end do
  tree%nodeBase                => nodes(1)%node
  tree%event                   => null()
  nodes(1)%node%firstChild     => nodes(2)%node
  nodes(2)%node%parent         => nodes(1)%node
  nodes(2)%node%sibling        => nodes(3)%node
  nodes(3)%node%parent         => nodes(1)%node
  nodes(3)%node%sibling        => nodes(5)%node
  nodes(5)%node%parent         => nodes(1)%node
  nodes(5)%node%sibling        => nodes(6)%node
  nodes(6)%node%parent         => nodes(1)%node
  nodes(3)%node%firstChild     => nodes(4)%node
  nodes(4)%node%parent         => nodes(3)%node
  ! Walk the branch with no time limit. All nodes must be visited, depth-first.
  countAll  =0
  visitedAll=0_kind_int8
  treeWalker=mergerTreeWalkerIsolatedNodesBranch(nodes(1)%node)
  do while (treeWalker%next(node))
     countAll=countAll+1
     if (countAll <= size(visitedAll)) visitedAll(countAll)=node%index()
  end do
  call Assert('all nodes visited when not time limited'             ,countAll               ,6                                      )
  call Assert('depth-first order when not time limited'             ,visitedAll             ,[2_kind_int8,4_kind_int8,3_kind_int8,5_kind_int8,6_kind_int8,1_kind_int8])
  ! Walk the branch with a time limit. Exactly the nodes existing at or after the earliest time (nodes 3, 6 and 1) must be visited.
  countLimited  =0
  visitedLimited=0_kind_int8
  treeWalker    =mergerTreeWalkerIsolatedNodesBranch(nodes(1)%node,timeEarliest=timeEarliest)
  do while (treeWalker%next(node))
     countLimited=countLimited+1
     if (countLimited <= size(visitedLimited)) visitedLimited(countLimited)=node%index()
  end do
  call Assert('number of nodes visited when time limited'           ,countLimited           ,3                                      )
  call Assert('non-primary child after earliest time visited'       ,visitedLimited(1)      ,3_kind_int8                            )
  call Assert('later sibling reached past an early sibling'         ,visitedLimited(2)      ,6_kind_int8                            )
  call Assert('branch head visited last when time limited'          ,visitedLimited(3)      ,1_kind_int8                            )
  ! Clean up.
  call tree%destroy()
  deallocate(tree)
  call Unit_Tests_End_Group               ()
  call Unit_Tests_Finish                  ()
  call Node_Components_Thread_Uninitialize()
  call Node_Components_Uninitialize       ()
  call nodeClassHierarchyFinalize         ()
end program Tests_Merger_Tree_Walker_Branch_Time_Limited
