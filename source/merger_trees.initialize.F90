!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

!% Contains a module which implements initialization of merger tree structures.

module Merger_Trees_Initialize
  !% Implements initialization of merger tree structures.
  implicit none
  private
  public :: Merger_Tree_Initialize

contains

  subroutine Merger_Tree_Initialize(tree,endTime)
    !% Walk through all nodes of a tree and call any routines that requested to perform initialization tasks.
    use :: Galacticus_Nodes   , only : mergerTree              , nodeComponentBasic, treeNode
    use :: Merger_Tree_Walkers, only : mergerTreeWalkerAllNodes
    !# <include directive="mergerTreeInitializeTask" type="moduleUse">
    include 'merger_trees.initialize.tasks.modules.inc'
    !# </include>
    implicit none
    type            (mergerTree              ), intent(inout) :: tree
    double precision                          , intent(in   ) :: endTime
    type            (treeNode                ), pointer       :: node
    class           (nodeComponentBasic      ), pointer       :: basic
    type            (mergerTreeWalkerAllNodes)                :: treeWalker

    if (tree%initializedUntil < endTime) then
       treeWalker=mergerTreeWalkerAllNodes(tree,spanForest=.false.)
       do while (treeWalker%next(node))
          ! Initialize only nodes that exist before the end time.
          basic => node%basic()
          if (basic%time() > tree%initializedUntil .and. basic%time() <= endTime) then
             ! Call subroutines to perform any necessary initialization of this node.
             !# <include directive="mergerTreeInitializeTask" type="functionCall" functionType="void">
             !#  <functionArgs>node</functionArgs>
             include 'merger_trees.initialize.tasks.inc'
             !# </include>
          end if
       end do
       tree%initializedUntil=endTime
    end if
    return
  end subroutine Merger_Tree_Initialize

end module Merger_Trees_Initialize
