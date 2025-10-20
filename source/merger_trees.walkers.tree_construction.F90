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
  Implements a tree walker for trees under construction.
  !!}
  use :: Galacticus_Nodes, only : mergerTree, treeNode

  !![
  <mergerTreeWalker name="mergerTreeWalkerTreeConstruction">
   <description>Provides a merger tree walker for trees under construction.</description>
  </mergerTreeWalker>
  !!]
  type, extends(mergerTreeWalkerClass) :: mergerTreeWalkerTreeConstruction
     !!{
     A merger tree walker for trees under construction.
     !!}
     private
     type   (mergerTree), pointer :: tree         => null()
     type   (treeNode  ), pointer :: node         => null()
     logical                      :: nodesRemain_
   contains
     !![
     <methods>
       <method description="Set the walker to the given node." method="setNode"/>
     </methods>
     !!]
     procedure :: next        => treeConstructionNext
     procedure :: nodesRemain => treeConstructionNodesRemain
     procedure :: setNode     => treeConstructionSetNode
  end type mergerTreeWalkerTreeConstruction

  interface mergerTreeWalkerTreeConstruction
     !!{
     Constructors for the \refClass{mergerTreeWalkerTreeConstruction} merger tree walker class.
     !!}
     module procedure treeConstructionParameters
     module procedure treeConstructionInternal
  end interface mergerTreeWalkerTreeConstruction

contains

  function treeConstructionParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerTreeWalkerTreeConstruction} merger tree walker class which takes a parameter set as input.
    !!}
    use :: Error           , only : Error_Report
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(mergerTreeWalkerTreeConstruction)                :: self
    type(inputParameters                 ), intent(inout) :: parameters
    !$GLC attributes unused :: self, parameters

    call Error_Report('this class can not be built from parameters'//{introspection:location})
    return
  end function treeConstructionParameters

  function treeConstructionInternal(tree) result(self)
    !!{
    Internal constructor for the \refClass{mergerTreeWalkerTreeConstruction} merger tree walker class.
    !!}
    implicit none
    type   (mergerTreeWalkerTreeConstruction)                        :: self
    type   (mergerTree                      ), intent(in   ), target :: tree

    self%tree        => tree
    self%node        => null()
    self%nodesRemain_ = .true.
    return
  end function treeConstructionInternal

  recursive logical function treeConstructionNext(self,node)
    !!{
    This function will update the given {\normalfont \ttfamily node} to the next node which should be visited in a tree to
    perform a walk suitable for trees under construction.
    !!}
    implicit none
    class(mergerTreeWalkerTreeConstruction), intent(inout)          :: self
    type (treeNode                        ), intent(inout), pointer :: node

    if (.not.self%nodesRemain_) then
       node                 => null()
       treeConstructionNext =  .false.
       return
    end if
    treeConstructionNext=.true.
    if (.not.associated(self%node)) then
       ! Initial node. This should be the base node, unless that already has children, in which case descend following normal
       ! rules.
       self%node => self%tree%nodeBase
       if (.not.associated(self%node%firstChild)) then
          node => self%node
          return
       end if
    end if
    if (associated(self%node%firstChild)) then
       ! Move to the primary child if one exists.
       do while (associated(self%node%firstChild))
          self%node => self%node%firstChild
       end do
    else
       if (associated(self%node%sibling)) then
          self%node => self%node%sibling
          do while (associated(self%node%firstChild))
             self%node => self%node%firstChild
          end do
       else
          do while (associated(self%node))
             if (associated(self%node%parent)) then
                self%node => self%node%parent
                if (associated(self%node%sibling)) then
                   self%node => self%node%sibling
                   do while (associated(self%node%firstChild))
                      self%node => self%node%firstChild
                   end do
                   exit
                end if
             else
                self%node            => null()
                self%nodesRemain_    = .false.
                treeConstructionNext = .false.
             end if
          end do
       end if
    end if
    node => self%node
    return
  end function treeConstructionNext

  logical function treeConstructionNodesRemain(self)
    !!{
    Returns true if more nodes remain to be walked to.
    !!}
    implicit none
    class(mergerTreeWalkerTreeConstruction), intent(inout) :: self

    treeConstructionNodesRemain=self%nodesRemain_
    return
  end function treeConstructionNodesRemain

  subroutine treeConstructionSetNode(self,node)
    !!{
    Set the current node for the walker.
    !!}
    implicit none
    class(mergerTreeWalkerTreeConstruction), intent(inout)         :: self
    type (treeNode                        ), intent(in   ), target :: node

    self%nodesRemain_ =  .true.
    self%node         => node
    self%tree         => node%hostTree
    return
  end subroutine treeConstructionSetNode

