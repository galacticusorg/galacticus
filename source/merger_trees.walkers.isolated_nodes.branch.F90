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
Implements a depth-first merger tree walker over all isolated nodes in a given branch.
!!}

  !![
  <mergerTreeWalker name="mergerTreeWalkerIsolatedNodesBranch">
   <description>Provides a merger tree walker which iterates depth-first over all isolated nodes in a given branch.</description>
  </mergerTreeWalker>
  !!]
  type, extends(mergerTreeWalkerClass) :: mergerTreeWalkerIsolatedNodesBranch
     !!{
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
     !!{
     Constructors for the \refClass{mergerTreeWalkerIsolatedNodesBranch} merger tree walker class.
     !!}
     module procedure isolatedNodesBranchParameters
     module procedure isolatedNodesBranchInternal
  end interface mergerTreeWalkerIsolatedNodesBranch

contains

  function isolatedNodesBranchParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerTreeWalkerIsolatedNodesBranch} merger tree walker class which takes a parameter set as input.
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
    !!{
    Internal constructor for the \refClass{mergerTreeWalkerIsolatedNodesBranch} merger tree walker class.
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
    !!{
    This function will update the given {\normalfont \ttfamily node} to the next node which should be visited in a tree branch
    to perform a depth-first walk. Once the entire branch has been walked, a {\normalfont \ttfamily null()} pointer will be
    set, and a value of {\normalfont \ttfamily false} returned indicating that there are no more nodes to walk. Each node will
    be visited once and once only if the branch is walked in this way.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class(mergerTreeWalkerIsolatedNodesBranch), intent(inout)          :: self
    type (treeNode                           ), intent(inout), pointer :: node
    class(nodeComponentBasic                 )               , pointer :: basic

    ! If the node is currently pointing to the head node of the branch, then the tree walk is complete.
    if (associated(self%node,self%branchHead)) then
       node                    => null()
       self%nodesRemain_       =  .false.
       isolatedNodesBranchNext =  .false.
       return
    end if
    ! If the node is currently null, set to the head node of the branch, and descend to children.
    if (.not.associated(self%node)) then
       self%node => self%branchHead
       do while (associated(self%node%firstChild))
          if (self%timeLimited) then
             basic => self%node%firstChild%basic()
             if (basic%time() < self%timeEarliest) exit
          end if
          self%node => self%node%firstChild
       end do
    else
       if (associated(self%node%sibling)) then
          self%node => self%node%sibling
          do while (associated(self%node%firstChild))
          if (self%timeLimited) then
             basic => self%node%firstChild%basic()
             if (basic%time() < self%timeEarliest) exit
          end if
             self%node => self%node%firstChild
          end do
       else
          self%node => self%node%parent
       end if
    end if
    node                    => self%node
    isolatedNodesBranchNext =  .true.
    return
  end function isolatedNodesBranchNext

  logical function isolatedNodesBranchNodesRemain(self)
    !!{
    Returns true if nodes remain to be visited in the branch.
    !!}
    implicit none
    class(mergerTreeWalkerIsolatedNodesBranch), intent(inout) :: self

    isolatedNodesBranchNodesRemain=self%nodesRemain_
    return
  end function isolatedNodesBranchNodesRemain
