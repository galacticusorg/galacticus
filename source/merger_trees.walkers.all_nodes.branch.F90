!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
Implements a depth-first merger tree walker over all all nodes in a given branch.
!!}

  !![
  <mergerTreeWalker name="mergerTreeWalkerAllNodesBranch">
   <description>Provides a merger tree walker which iterates depth-first over all all nodes in a given branch.</description>
  </mergerTreeWalker>
  !!]
  type, extends(mergerTreeWalkerClass) :: mergerTreeWalkerAllNodesBranch
     !!{
     A merger tree walker which iterates depth-first over all all nodes in a given branch.
     !!}
     private
     type   (treeNode), pointer :: branchHead   => null(), node => null()
     logical                    :: nodesRemain_
   contains
     !![
     <methods>
       <method description="Descend through the hierarchy to the deepest node along the current branch." method="descend" />
     </methods>
     !!]
     procedure :: next        => allNodesBranchNext
     procedure :: nodesRemain => allNodesBranchNodesRemain
     procedure :: descend     => allNodesBranchDescend
 end type mergerTreeWalkerAllNodesBranch

  interface mergerTreeWalkerAllNodesBranch
     !!{
     Constructors for the {\normalfont \ttfamily allNodesBranch} merger tree walker class.
     !!}
     module procedure allNodesBranchParameters
     module procedure allNodesBranchInternal
  end interface mergerTreeWalkerAllNodesBranch

contains

  function allNodesBranchParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily allNodesBranch} merger tree walker class which takes a parameter set as input.
    !!}
    use :: Error           , only : Error_Report
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(mergerTreeWalkerAllNodesBranch)                :: self
    type(inputParameters               ), intent(inout) :: parameters
    !$GLC attributes unused :: self, parameters

    call Error_Report('this class can not be built from parameters'//{introspection:location})
    return
  end function allNodesBranchParameters

  function allNodesBranchInternal(branchHead) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily allNodesBranch} merger tree walker class.
    !!}
    implicit none
    type(mergerTreeWalkerAllNodesBranch)                          :: self
    type(treeNode                      ), intent(in   ), target   :: branchHead

    self%branchHead   => branchHead
    self%node         => null()
    self%nodesRemain_ = .true.
    return
  end function allNodesBranchInternal

  logical function allNodesBranchNext(self,node)
    !!{
    This function will update the given {\normalfont \ttfamily node} to the next node which should be visited in a tree branch
    to perform a depth-first walk. Once the entire branch has been walked, a {\normalfont \ttfamily null()} pointer will be
    set, and a value of {\normalfont \ttfamily false} returned indicating that there are no more nodes to walk. Each node will
    be visited once and once only if the branch is walked in this way.
    !!}
    implicit none
    class(mergerTreeWalkerAllNodesBranch), intent(inout)          :: self
    type (treeNode                      ), intent(inout), pointer :: node

    ! If the node is currently pointing to the head node of the branch, then the tree walk is complete.
    if (associated(self%node,self%branchHead)) then
       node               => null()
       self%nodesRemain_  =  .false.
       allNodesBranchNext =  .false.
       return
    end if
    ! If the node is currently null, set to the head node of the branch, and descend to children.
    if (.not.associated(self%node)) then
       ! This is the base of the branch.
       self%node => self%branchHead
       ! Descend through satellites and children.
       call self%descend()
    else
       if (associated(self%node%sibling)) then
          self%node => self%node%sibling
          call self%descend()
       else
          ! About to move back up the tree. Check if the node we're moving up from is a satellite.
          if (self%node%isSatellite()) then
             ! It is a satellite. Therefore, the parent may have children that have yet to be
             ! visited. Check if the parent has children.
             if (associated(self%node%parent%firstChild)) then
                ! Parent does have children, so move to the first one.
                self%node => self%node%parent%firstChild
                ! Descend through satellites and children.
                call self%descend()
             else
                ! Parent has no children, so move to the parent.
                self%node => self%node%parent
             end if
          else
             ! It is not a satellite, so all satellites and children have been processed.
             self%node => self%node%parent
          end if
       end if
    end if
    node              => self%node
    allNodesBranchNext =  .true.
    return
  end function allNodesBranchNext

  logical function allNodesBranchNodesRemain(self)
    !!{
    Returns true if nodes remain to be visited in the branch.
    !!}
    implicit none
    class(mergerTreeWalkerAllNodesBranch), intent(inout) :: self

    allNodesBranchNodesRemain=self%nodesRemain_
    return
  end function allNodesBranchNodesRemain

  subroutine allNodesBranchDescend(self)
    !!{
    Descend to the deepest progenitor (satellites and children) of the current branch.
    !!}
    implicit none
    class(mergerTreeWalkerAllNodesBranch), intent(inout) :: self

    ! Descend through satellites and children.
    do while (associated(self%node%firstSatellite).or.associated(self%node%firstChild))
       if (associated(self%node%firstSatellite)) then
          self%node => self%node%firstSatellite
       else
          self%node => self%node%firstChild
       end if
    end do
    return
  end subroutine allNodesBranchDescend

