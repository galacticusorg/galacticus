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
Implements a depth-first merger tree walker over all isolated nodes in a given branch.
!!}

  !![
  <mergerTreeWalker name="mergerTreeWalkerIsolatedNodesBranch" docformat="rst">
   <description>
   Provides a merger tree walker which iterates depth-first over all isolated nodes in a given branch. If an earliest time is
   given to the internal constructor, nodes which precede that time are skipped, so that the walk visits the branch head and
   every node in the branch which exists at or after that time.
   </description>
  </mergerTreeWalker>
  !!]
  type, extends(mergerTreeWalkerClass) :: mergerTreeWalkerIsolatedNodesBranch
     !!{RST
     A merger tree walker which iterates depth-first over all isolated nodes in a given branch.
     !!}
     private
     type            (treeNode), pointer :: branchHead   => null(), node        => null()
     logical                             :: nodesRemain_          , timeLimited
     double precision                    :: timeEarliest
   contains
     procedure :: next        => isolatedNodesBranchNext
     procedure :: nodesRemain => isolatedNodesBranchNodesRemain
  end type mergerTreeWalkerIsolatedNodesBranch

  interface mergerTreeWalkerIsolatedNodesBranch
     !!{RST
     Constructors for the :galacticus-class:`mergerTreeWalkerIsolatedNodesBranch` merger tree walker class.
     !!}
     module procedure isolatedNodesBranchParameters
     module procedure isolatedNodesBranchInternal
  end interface mergerTreeWalkerIsolatedNodesBranch

contains

  function isolatedNodesBranchParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`mergerTreeWalkerIsolatedNodesBranch` merger tree walker class which takes a parameter set as input.
    !!}
    use :: Error           , only : Error_Report
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(mergerTreeWalkerIsolatedNodesBranch)                :: self
    type(inputParameters                    ), intent(inout) :: parameters
    !$GLC attributes unused :: self, parameters

    call Error_Report('this class can not be built from parameters'//{introspection:location})
    return
  end function isolatedNodesBranchParameters

  function isolatedNodesBranchInternal(branchHead,timeEarliest) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`mergerTreeWalkerIsolatedNodesBranch` merger tree walker class. If
    ``timeEarliest`` is present, nodes which precede that time are skipped (the branch head is always visited).
    !!}
    implicit none
    type            (mergerTreeWalkerIsolatedNodesBranch)                          :: self
    type            (treeNode                           ), intent(in   ), target   :: branchHead
    double precision                                     , intent(in   ), optional :: timeEarliest

    self%branchHead   => branchHead
    self%node         => null()
    self%nodesRemain_ = .true.
    self%timeLimited  = present(timeEarliest)
    if (self%timeLimited) self%timeEarliest=timeEarliest
    return
  end function isolatedNodesBranchInternal

  logical function isolatedNodesBranchNext(self,node)
    !!{RST
    This function will update the given ``node`` to the next node which should be visited in a tree branch to perform a depth-first walk. Once the entire branch has been walked, a ``null()`` pointer will be set, and a value of ``false`` returned indicating that there are no more nodes to walk. Each node will be visited once and once only if the branch is walked in this way (or, for a time-limited walk, the branch head and each node existing at or after the earliest time).
    !!}
    use :: Galacticus_Nodes, only : treeNode
    implicit none
    class(mergerTreeWalkerIsolatedNodesBranch), intent(inout)          :: self
    type (treeNode                           ), intent(inout), pointer :: node
    type (treeNode                           )               , pointer :: nodeNext

    ! If the node is currently pointing to the head node of the branch, then the tree walk is complete.
    if (associated(self%node,self%branchHead)) then
       node                    => null()
       self%nodesRemain_       =  .false.
       isolatedNodesBranchNext =  .false.
       return
    end if
    ! If the node is currently null, set to the head node of the branch, and descend to children. Otherwise, move to the next
    ! sibling (and descend to its children), or, if there is none, to the parent. If the walk is time limited, nodes which
    ! precede the earliest time are skipped, so that every node existing at or after that time is visited. (A node's children
    ! need not all exist at the same time, so stopping the descent at an early first child would miss later siblings and their
    ! progenitors.)
    if (.not.associated(self%node)) then
       self%node => self%branchHead
       nodeNext  => isolatedNodesBranchNextInTime(self,self%node%firstChild)
       do while (associated(nodeNext))
          self%node => nodeNext
          nodeNext  => isolatedNodesBranchNextInTime(self,self%node%firstChild)
       end do
    else
       nodeNext => isolatedNodesBranchNextInTime(self,self%node%sibling)
       if (associated(nodeNext)) then
          do while (associated(nodeNext))
             self%node => nodeNext
             nodeNext  => isolatedNodesBranchNextInTime(self,self%node%firstChild)
          end do
       else
          self%node => self%node%parent
       end if
    end if
    node                    => self%node
    isolatedNodesBranchNext =  .true.
    return
  end function isolatedNodesBranchNext

  function isolatedNodesBranchNextInTime(self,node) result(nodeNext)
    !!{RST
    Return a pointer to the first node in the sibling list starting at ``node`` which exists at or after the earliest time of a
    time-limited walk (or simply to ``node`` if the walk is not time limited). A null pointer is returned if there is no such node.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    type (treeNode                           ), pointer                :: nodeNext
    class(mergerTreeWalkerIsolatedNodesBranch), intent(in   )          :: self
    type (treeNode                           ), intent(in   ), pointer :: node
    class(nodeComponentBasic                 )               , pointer :: basic

    nodeNext => node
    if (.not.self%timeLimited) return
    do while (associated(nodeNext))
       basic => nodeNext%basic()
       if (basic%time() >= self%timeEarliest) return
       nodeNext => nodeNext%sibling
    end do
    return
  end function isolatedNodesBranchNextInTime

  logical function isolatedNodesBranchNodesRemain(self)
    !!{RST
    Returns true if nodes remain to be visited in the branch.
    !!}
    implicit none
    class(mergerTreeWalkerIsolatedNodesBranch), intent(inout) :: self

    isolatedNodesBranchNodesRemain=self%nodesRemain_
    return
  end function isolatedNodesBranchNodesRemain
