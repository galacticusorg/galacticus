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
Implements a galactic filter which applies another filter to satellite nodes of the given node and returns true if \emph{any}
satellite node passes.
!!}
  
  !![
  <galacticFilter name="galacticFilterAnySatelliteNode">
   <description>
    Applies a filter to satellite nodes of the given node and returns true if \emph{any} satellite node passes.
   </description>
  </galacticFilter>
  !!]
  type, extends(galacticFilterClass) :: galacticFilterAnySatelliteNode
     !!{
     A galactic filter which applies another filter to satellite nodes of the given node and returns true if \emph{any}
     satellite node passes.
     !!}
     private
     class(galacticFilterClass), pointer :: galacticFilter_ => null()
   contains
     final     ::           anySatelliteNodeDestructor
     procedure :: passes => anySatelliteNodePasses
  end type galacticFilterAnySatelliteNode

  interface galacticFilterAnySatelliteNode
     !!{
     Constructors for the \refClass{galacticFilterAnySatelliteNode} galactic filter class.
     !!}
     module procedure anySatelliteNodeConstructorParameters
     module procedure anySatelliteNodeConstructorInternal
  end interface galacticFilterAnySatelliteNode

contains

  function anySatelliteNodeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{galacticFilterAnySatelliteNode} galactic filter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (galacticFilterAnySatelliteNode)                :: self
    type (inputParameters               ), intent(inout) :: parameters
    class(galacticFilterClass           ), pointer       :: galacticFilter_
         
    !![
    <objectBuilder class="galacticFilter" name="galacticFilter_" source="parameters"/>
    !!]
    self=galacticFilterAnySatelliteNode(galacticFilter_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="galacticFilter_"/>
    !!]
    return
  end function anySatelliteNodeConstructorParameters
  
  function anySatelliteNodeConstructorInternal(galacticFilter_) result(self)
    !!{
    Internal constructor for the \refClass{galacticFilterAnySatelliteNode} galactic filter class.
    !!}
    implicit none
    type (galacticFilterAnySatelliteNode)                        :: self
    class(galacticFilterClass           ), intent(in   ), target :: galacticFilter_
    !![
    <constructorAssign variables="*galacticFilter_"/>
    !!]

    return
  end function anySatelliteNodeConstructorInternal
  
  subroutine anySatelliteNodeDestructor(self)
    !!{
    Destructor for the \refClass{galacticFilterAnySatelliteNode} galactic filter class.
    !!}
    implicit none
    type(galacticFilterAnySatelliteNode), intent(inout) :: self

    !![
    <objectDestructor name="self%galacticFilter_"/>
    !!]
    return
  end subroutine anySatelliteNodeDestructor
  
  logical function anySatelliteNodePasses(self,node)
    !!{
    Implement a filter on satellite node properties which passes if any satellite is passed.
    !!}
    implicit none
    class(galacticFilterAnySatelliteNode), intent(inout)         :: self
    type (treeNode                      ), intent(inout), target :: node
    type (treeNode                      ), pointer               :: nodeSatellite

    anySatelliteNodePasses =  .false.
    nodeSatellite          => node%firstSatellite
    do while (associated(nodeSatellite))
       anySatelliteNodePasses=self%galacticFilter_%passes(nodeSatellite)
       if (anySatelliteNodePasses) return
       nodeSatellite => nodeSatellite%sibling
    end do
    return
  end function anySatelliteNodePasses
