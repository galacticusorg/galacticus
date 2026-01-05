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
Implements a galactic filter which applies another filter to the isolated host node of the given node.
!!}
  
  !![
  <galacticFilter name="galacticFilterIsolatedHostNode">
   <description>Applies a filter to the isolated host node of the given node.</description>
  </galacticFilter>
  !!]
  type, extends(galacticFilterClass) :: galacticFilterIsolatedHostNode
     !!{
     A galactic filter which applies another filter to the isolated host node of the given node.
     !!}
     private
     class(galacticFilterClass), pointer :: galacticFilter_ => null()
   contains
     final     ::           isolatedHostNodeDestructor
     procedure :: passes => isolatedHostNodePasses
  end type galacticFilterIsolatedHostNode

  interface galacticFilterIsolatedHostNode
     !!{
     Constructors for the \refClass{galacticFilterIsolatedHostNode} galactic filter class.
     !!}
     module procedure isolatedHostNodeConstructorParameters
     module procedure isolatedHostNodeConstructorInternal
  end interface galacticFilterIsolatedHostNode

contains

  function isolatedHostNodeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{galacticFilterIsolatedHostNode} galactic filter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (galacticFilterIsolatedHostNode)                :: self
    type   (inputParameters               ), intent(inout) :: parameters
    class  (galacticFilterClass           ), pointer       :: galacticFilter_

    !![
    <objectBuilder class="galacticFilter" name="galacticFilter_" source="parameters"/>
    !!]
    self=galacticFilterIsolatedHostNode(galacticFilter_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="galacticFilter_"/>
    !!]
    return
  end function isolatedHostNodeConstructorParameters
  
  function isolatedHostNodeConstructorInternal(galacticFilter_) result(self)
    !!{
    Internal constructor for the \refClass{galacticFilterIsolatedHostNode} galactic filter class.
    !!}
    implicit none
    type (galacticFilterIsolatedHostNode)                        :: self
    class(galacticFilterClass           ), intent(in   ), target :: galacticFilter_
    !![
    <constructorAssign variables="*galacticFilter_"/>
    !!]

    return
  end function isolatedHostNodeConstructorInternal
  
  subroutine isolatedHostNodeDestructor(self)
    !!{
    Destructor for the \refClass{galacticFilterIsolatedHostNode} galactic filter class.
    !!}
    implicit none
    type(galacticFilterIsolatedHostNode), intent(inout) :: self

    !![
    <objectDestructor name="self%galacticFilter_"/>
    !!]
    return
  end subroutine isolatedHostNodeDestructor
  
  logical function isolatedHostNodePasses(self,node) result(passes)
    !!{
    Implement a filter on parent node properties.
    !!}
    implicit none
    class(galacticFilterIsolatedHostNode), intent(inout)         :: self
    type (treeNode                      ), intent(inout), target :: node
    type (treeNode                      ), pointer               :: nodeHost

    nodeHost => node
    do while (nodeHost%isSatellite())
       nodeHost => nodeHost%parent
    end do
    passes=self%galacticFilter_%passes(nodeHost)
    return
  end function isolatedHostNodePasses
