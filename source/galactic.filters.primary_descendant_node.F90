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
  Implements a galactic filter which applies another filter to all descendant nodes of the given node and
  passes if any primary descendant passes.  
  !!}
  
  !![
  <galacticFilter name="galacticFilterPrimaryDescendantNode">
   <description>Applies a filter to all descendant nodes of the given node and passes if any primary descendant passes.</description>
  </galacticFilter>
  !!]
  type, extends(galacticFilterClass) :: galacticFilterPrimaryDescendantNode
     !!{
     A galactic filter which applies another filter to all descendant nodes of the given node and passes if any primary descendant passes.
     !!}
     private
     class  (galacticFilterClass), pointer :: galacticFilter_ => null()
     logical                               :: allowSelf
   contains
     final     ::           primaryDescendantNodeDestructor
     procedure :: passes => primaryDescendantNodePasses
  end type galacticFilterPrimaryDescendantNode

  interface galacticFilterPrimaryDescendantNode
     !!{
     Constructors for the {\normalfont \ttfamily primaryDescendantNode} galactic filter class.
     !!}
     module procedure primaryDescendantNodeConstructorParameters
     module procedure primaryDescendantNodeConstructorInternal
  end interface galacticFilterPrimaryDescendantNode

contains

  function primaryDescendantNodeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily primaryDescendantNode} galactic filter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (galacticFilterPrimaryDescendantNode)                :: self
    type   (inputParameters                    ), intent(inout) :: parameters
    class  (galacticFilterClass                ), pointer       :: galacticFilter_
    logical                                                     :: allowSelf
         
    !![
    <inputParameter>
      <name>allowSelf</name>
      <source>parameters</source>
      <description>If true, the node itself is considered as a descendant, otherwise the node itself is excluded from the descendant node search.</description>
    </inputParameter>
    <objectBuilder class="galacticFilter" name="galacticFilter_" source="parameters"/>
    !!]
    self=galacticFilterPrimaryDescendantNode(allowSelf,galacticFilter_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="galacticFilter_"    />
    !!]
    return
  end function primaryDescendantNodeConstructorParameters
  
  function primaryDescendantNodeConstructorInternal(allowSelf,galacticFilter_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily primaryDescendantNode} galactic filter class.
    !!}
    implicit none
    type   (galacticFilterPrimaryDescendantNode)                        :: self
    logical                                     , intent(in   )         :: allowSelf
    class  (galacticFilterClass                ), intent(in   ), target :: galacticFilter_
    !![
    <constructorAssign variables="allowSelf, *galacticFilter_"/>
    !!]

    return
  end function primaryDescendantNodeConstructorInternal
  
  subroutine primaryDescendantNodeDestructor(self)
    !!{
    Destructor for  the {\normalfont \ttfamily primaryDescendantNode} galactic filter class.
    !!}
    implicit none
    type(galacticFilterPrimaryDescendantNode), intent(inout) :: self

    !![
    <objectDestructor name="self%galacticFilter_"/>
    !!]
    return
  end subroutine primaryDescendantNodeDestructor
  
  logical function primaryDescendantNodePasses(self,node)
    !!{
    Implement a filter on descendant node properties.
    !!}
    implicit none
    class(galacticFilterPrimaryDescendantNode), intent(inout)         :: self
    type (treeNode                           ), intent(inout), target :: node
    type (treeNode                           ), pointer               :: nodeDescendant

    primaryDescendantNodePasses=.false.
    if (self%allowSelf) then
       nodeDescendant    => node
    else
       if (node%isPrimaryProgenitor()) then
          nodeDescendant => node%parent
       else
          nodeDescendant => null()
       end if
    end if
    do while (associated(nodeDescendant))
       primaryDescendantNodePasses=self%galacticFilter_%passes(nodeDescendant)
       if (nodeDescendant%isPrimaryProgenitor()) then
          nodeDescendant => nodeDescendant%parent
       else
          nodeDescendant => null()
       end if
       if (primaryDescendantNodePasses) exit
    end do
    return
  end function primaryDescendantNodePasses
