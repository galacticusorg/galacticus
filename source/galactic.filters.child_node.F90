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
Implements a galactic filter which applies another filter to a child node of the given node.
!!}
  
  !![
  <galacticFilter name="galacticFilterChildNode">
   <description>
     Applies a filter to a child node of the given node. The {\normalfont \ttfamily [childRank]} parameter specifies which child
     to use---a rank of 1 means the first child, a rank of $N$ means the sibling of the $N-1$ rank child. If a child of the
     specified rank does not exist this filter fails to pass.
   </description>
  </galacticFilter>
  !!]
  type, extends(galacticFilterClass) :: galacticFilterChildNode
     !!{
     A galactic filter which applies another filter to a child node of the given node.
     !!}
     private
     class  (galacticFilterClass), pointer :: galacticFilter_ => null()
     integer                               :: childRank
   contains
     final     ::           childNodeDestructor
     procedure :: passes => childNodePasses
  end type galacticFilterChildNode

  interface galacticFilterChildNode
     !!{
     Constructors for the \refClass{galacticFilterChildNode} galactic filter class.
     !!}
     module procedure childNodeConstructorParameters
     module procedure childNodeConstructorInternal
  end interface galacticFilterChildNode

contains

  function childNodeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{galacticFilterChildNode} galactic filter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (galacticFilterChildNode)                :: self
    type   (inputParameters        ), intent(inout) :: parameters
    class  (galacticFilterClass    ), pointer       :: galacticFilter_
    integer                                         :: childRank

    !![
    <inputParameter>
      <name>childRank</name>
      <source>parameters</source>
      <description>The rank of the child to use---a rank of 1 means the first child, a rank of $N$ means the sibling of the $N-1$ rank child.</description>
    </inputParameter>
    <objectBuilder class="galacticFilter" name="galacticFilter_" source="parameters"/>
    !!]
    self=galacticFilterChildNode(childRank,galacticFilter_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="galacticFilter_"/>
    !!]
    return
  end function childNodeConstructorParameters
  
  function childNodeConstructorInternal(childRank,galacticFilter_) result(self)
    !!{
    Internal constructor for the \refClass{galacticFilterChildNode} galactic filter class.
    !!}
    implicit none
    type   (galacticFilterChildNode)                        :: self
    integer                         , intent(in   )         :: childRank
    class  (galacticFilterClass    ), intent(in   ), target :: galacticFilter_
    !![
    <constructorAssign variables="childRank, *galacticFilter_"/>
    !!]

    return
  end function childNodeConstructorInternal
  
  subroutine childNodeDestructor(self)
    !!{
    Destructor for the \refClass{galacticFilterChildNode} galactic filter class.
    !!}
    implicit none
    type(galacticFilterChildNode), intent(inout) :: self

    !![
    <objectDestructor name="self%galacticFilter_"/>
    !!]
    return
  end subroutine childNodeDestructor
  
  logical function childNodePasses(self,node)
    !!{
    Implement a filter on child node properties.
    !!}
    implicit none
    class  (galacticFilterChildNode), intent(inout)         :: self
    type   (treeNode               ), intent(inout), target :: node
    type   (treeNode               ), pointer               :: nodeChild
    integer                                                 :: i

    ! Find the child of the requested rank.
    i         =  1
    nodeChild => node%firstChild
    do while (i < self%childRank .and. associated(nodeChild))
       i=i+1
       nodeChild => nodeChild%sibling
    end do
    ! If the child exists, apply a filter to it, otherwise return false.
    childNodePasses=.false.
    if (associated(nodeChild)) childNodePasses=self%galacticFilter_%passes(nodeChild)
    return
  end function childNodePasses
