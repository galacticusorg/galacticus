!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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
Contains a module which implements a galactic filter which applies another filter to all descendant nodes of the given node and passes if any descendant passes.
!!}
  
  use :: Cosmology_Functions, only : cosmologyFunctionsClass

  !![
  <galacticFilter name="galacticFilterAnyDescendantNode">
   <description>Applies a filter to all descendant nodes of the given node and passes if any descendant passes.</description>
  </galacticFilter>
  !!]
  type, extends(galacticFilterClass) :: galacticFilterAnyDescendantNode
     !!{
     A galactic filter which applies another filter to all descendant nodes of the given node and passes if any descendant passes.
     !!}
     private
     class  (galacticFilterClass), pointer :: galacticFilter_ => null()
     logical                               :: allowSelf
   contains
     final     ::           anyDescendantNodeDestructor
     procedure :: passes => anyDescendantNodePasses
  end type galacticFilterAnyDescendantNode

  interface galacticFilterAnyDescendantNode
     !!{
     Constructors for the ``anyDescendantNode'' galactic filter class.
     !!}
     module procedure anyDescendantNodeConstructorParameters
     module procedure anyDescendantNodeConstructorInternal
  end interface galacticFilterAnyDescendantNode

contains

  function anyDescendantNodeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``anyDescendantNode'' galactic filter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (galacticFilterAnyDescendantNode)                :: self
    type   (inputParameters                ), intent(inout) :: parameters
    class  (galacticFilterClass            ), pointer       :: galacticFilter_
    logical                                                 :: allowSelf
         
    !![
    <inputParameter>
      <name>allowSelf</name>
      <source>parameters</source>
      <description>If true, the node itself is considered as a descendant, otherwise the node itself is excluded from the descendant node search.</description>
    </inputParameter>
    <objectBuilder class="galacticFilter" name="galacticFilter_" source="parameters"/>
    !!]
    self=galacticFilterAnyDescendantNode(allowSelf,galacticFilter_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="galacticFilter_"    />
    !!]
    return
  end function anyDescendantNodeConstructorParameters
  
  function anyDescendantNodeConstructorInternal(allowSelf,galacticFilter_) result(self)
    !!{
    Internal constructor for the ``anyDescendantNode'' galactic filter class.
    !!}
    implicit none
    type   (galacticFilterAnyDescendantNode)                        :: self
    logical                                 , intent(in   )         :: allowSelf
    class  (galacticFilterClass            ), intent(in   ), target :: galacticFilter_
    !![
    <constructorAssign variables="allowSelf, *galacticFilter_"/>
    !!]

    return
  end function anyDescendantNodeConstructorInternal
  
  subroutine anyDescendantNodeDestructor(self)
    !!{
    Destructor for  the ``anyDescendantNode'' galactic filter class.
    !!}
    implicit none
    type(galacticFilterAnyDescendantNode), intent(inout) :: self

    !![
    <objectDestructor name="self%galacticFilter_"/>
    !!]
    return
  end subroutine anyDescendantNodeDestructor
  
  logical function anyDescendantNodePasses(self,node)
    !!{
    Implement a filter on descendant node properties.
    !!}
    implicit none
    class(galacticFilterAnyDescendantNode), intent(inout)         :: self
    type (treeNode                       ), intent(inout), target :: node
    type (treeNode                       ), pointer               :: nodeDescendant

    anyDescendantNodePasses=.false.
    if (self%allowSelf) then
       nodeDescendant => node
    else
       nodeDescendant => node%parent
    end if
    do while (associated(nodeDescendant))
       anyDescendantNodePasses =  self          %galacticFilter_%passes(nodeDescendant)
       nodeDescendant          => nodeDescendant%parent
       if (anyDescendantNodePasses) exit
    end do
    return
  end function anyDescendantNodePasses
