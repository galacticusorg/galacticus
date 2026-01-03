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
Implements a galactic filter which applies another filter to a parent node of the given node.
!!}
  
  !![
  <galacticFilter name="galacticFilterParentNode">
   <description>
     Applies a filter to a parent node of the given node. If a parent of the specified rank does not exist this filter fails to
     pass.
   </description>
  </galacticFilter>
  !!]
  type, extends(galacticFilterClass) :: galacticFilterParentNode
     !!{
     A galactic filter which applies another filter to a parent node of the given node.
     !!}
     private
     class(galacticFilterClass), pointer :: galacticFilter_ => null()
   contains
     final     ::           parentNodeDestructor
     procedure :: passes => parentNodePasses
  end type galacticFilterParentNode

  interface galacticFilterParentNode
     !!{
     Constructors for the \refClass{galacticFilterParentNode} galactic filter class.
     !!}
     module procedure parentNodeConstructorParameters
     module procedure parentNodeConstructorInternal
  end interface galacticFilterParentNode

contains

  function parentNodeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{galacticFilterParentNode} galactic filter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (galacticFilterParentNode)                :: self
    type   (inputParameters         ), intent(inout) :: parameters
    class  (galacticFilterClass     ), pointer       :: galacticFilter_

    !![
    <objectBuilder class="galacticFilter" name="galacticFilter_" source="parameters"/>
    !!]
    self=galacticFilterParentNode(galacticFilter_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="galacticFilter_"/>
    !!]
    return
  end function parentNodeConstructorParameters
  
  function parentNodeConstructorInternal(galacticFilter_) result(self)
    !!{
    Internal constructor for the \refClass{galacticFilterParentNode} galactic filter class.
    !!}
    implicit none
    type   (galacticFilterParentNode)                        :: self
    class  (galacticFilterClass     ), intent(in   ), target :: galacticFilter_
    !![
    <constructorAssign variables="*galacticFilter_"/>
    !!]

    return
  end function parentNodeConstructorInternal
  
  subroutine parentNodeDestructor(self)
    !!{
    Destructor for the \refClass{galacticFilterParentNode} galactic filter class.
    !!}
    implicit none
    type(galacticFilterParentNode), intent(inout) :: self

    !![
    <objectDestructor name="self%galacticFilter_"/>
    !!]
    return
  end subroutine parentNodeDestructor
  
  logical function parentNodePasses(self,node)
    !!{
    Implement a filter on parent node properties.
    !!}
    implicit none
    class  (galacticFilterParentNode), intent(inout)         :: self
    type   (treeNode                ), intent(inout), target :: node
    type   (treeNode                ), pointer               :: nodeParent

    nodeParent => node%parent
    ! If the parent exists, apply a filter to it, otherwise return false.
    parentNodePasses=.false.
    if (associated(nodeParent)) parentNodePasses=self%galacticFilter_%passes(nodeParent)
    return
  end function parentNodePasses
