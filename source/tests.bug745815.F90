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

program Tests_Bug745815
  !!{
  Tests for regression of Bug \#745815 (http://bugs.launchpad.net/galacticus/+bug/745815): Skipping of a node during a tree
  walk.
  !!}
  use :: Display            , only : displayVerbositySet     , verbosityLevelStandard
  use :: Galacticus_Nodes   , only : mergerTree              , treeNode              , treeNodeList
  use :: ISO_Varying_String , only : assignment(=)           , varying_string
  use :: Input_Parameters   , only : inputParameters
  use :: Kind_Numbers       , only : kind_int8
  use :: Merger_Tree_Walkers, only : mergerTreeWalkerAllNodes
  use :: Unit_Tests         , only : Assert                  , Unit_Tests_Begin_Group, Unit_Tests_End_Group, Unit_Tests_Finish
  implicit none
  type   (varying_string          )          :: parameterFile
  type   (treeNodeList            )          :: nodes        (5)
  logical                                    :: nodeFound    (5)
  type   (mergerTree              )          :: tree
  type   (treeNode                ), pointer :: node
  integer(kind=kind_int8          )          :: i
  type   (inputParameters         )          :: parameters
  type   (mergerTreeWalkerAllNodes)          :: treeWalker

  ! Set verbosity level.
  call displayVerbositySet(verbosityLevelStandard)
  ! Begin unit tests.
  call Unit_Tests_Begin_Group("Bug #745815: Node skip during tree-walk")

  ! Open the parameter file.
  parameterFile='testSuite/parameters/bug745815.xml'
  parameters=inputParameters(parameterFile)

  ! Create nodes.
  do i=1,5
     nodes(i)%node => treeNode()
  end do

  ! Set indices of nodes.
  call nodes(1)%node%indexSet(100017990003559_kind_int8)
  call nodes(2)%node%indexSet(100017990003560_kind_int8)
  call nodes(3)%node%indexSet(100017990003561_kind_int8)
  call nodes(4)%node%indexSet(100017990003562_kind_int8)
  call nodes(5)%node%indexSet(100017990003571_kind_int8)

  ! Set child nodes.
  nodes(1)%node%firstChild => nodes(2)%node
  nodes(2)%node%firstChild => nodes(3)%node
  nodes(3)%node%firstChild => nodes(4)%node

  ! Set parent nodes.
  nodes(2)%node%parent => nodes(1)%node
  nodes(3)%node%parent => nodes(2)%node
  nodes(4)%node%parent => nodes(3)%node
  nodes(5)%node%parent => nodes(3)%node

  ! Set satellite nodes.
  nodes(3)%node%firstSatellite => nodes(5)%node

  ! Set the base node of our tree.
  tree%nodeBase => nodes(1)%node

  ! Walk the tree, with satellites.
  nodeFound=.false.
  treeWalker=mergerTreeWalkerAllNodes(tree,spanForest=.false.)
  do while (treeWalker%next(node))
     do i=1,5
        if (nodes(i)%node%index() == node%index()) nodeFound(i)=.true.
     end do
  end do
  call Assert('All nodes walked to',all(nodeFound),.true.)

  ! Destroy nodes.
  do i=1,5
     call nodes(i)%node%destroy()
  end do

  ! End unit tests.
  call Unit_Tests_End_Group()
  call Unit_Tests_Finish()

end program Tests_Bug745815
