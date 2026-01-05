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
  Implements a node operator class that relocates subhalos to the position of their host.
  !!}

  !![
  <nodeOperator name="nodeOperatorPositionToHost">
   <description>A node operator class that relocates subhalos to the position of their host.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorPositionToHost
     !!{
     A node operator class that relocates subhalos to the position of their host.
     !!}
     private
   contains
     procedure :: nodePromote => positionToHostNodePromote
     procedure :: nodesMerge  => positionToHostNodesMerge
     procedure :: autoHook    => positionToHostAutoHook
  end type nodeOperatorPositionToHost
  
  interface nodeOperatorPositionToHost
     !!{
     Constructors for the \refClass{nodeOperatorPositionToHost} node operator class.
     !!}
     module procedure positionToHostConstructorParameters
  end interface nodeOperatorPositionToHost
  
contains
  
  function positionToHostConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorPositionToHost} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorPositionToHost)                :: self
    type (inputParameters           ), intent(inout) :: parameters

    self=nodeOperatorPositionToHost()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function positionToHostConstructorParameters

  subroutine positionToHostAutoHook(self)
    !!{
    Attach to various event hooks.
    !!}
    use :: Events_Hooks, only : openMPThreadBindingAtLevel, satelliteHostChangeEvent
    implicit none
    class(nodeOperatorPositionToHost), intent(inout) :: self

    call satelliteHostChangeEvent%attach(self,positionToHostSatelliteHostChange,openMPThreadBindingAtLevel,label='positionToHost')
    return
  end subroutine positionToHostAutoHook

  subroutine positionToHostDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorPositionToHost} node operator class.
    !!}
    use :: Events_Hooks, only : openMPThreadBindingAtLevel, satelliteHostChangeEvent
    implicit none
    type(nodeOperatorPositionToHost), intent(inout) :: self

    if (satelliteHostChangeEvent%isAttached(self,positionToHostSatelliteHostChange)) call satelliteHostChangeEvent%detach(self,positionToHostSatelliteHostChange)
    return
  end subroutine positionToHostDestructor

  subroutine positionToHostSatelliteHostChange(self,node)
    !!{
    Update the maximum host mass of this node in response to a change in host.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(*       ), intent(inout)         :: self
    type (treeNode), intent(inout), target :: node
    
    select type (self)
    class is (nodeOperatorPositionToHost)
       call self%nodePromote(node)
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine positionToHostSatelliteHostChange
  
  subroutine positionToHostNodePromote(self,node)
    !!{
    Relocate any subhalos to the position of the new node on node promotion.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentPosition
    implicit none
    class(nodeOperatorPositionToHost), intent(inout) :: self
    type (treeNode                  ), intent(inout) :: node
    class(nodeComponentPosition     ), pointer       :: positionParent, positionSatellite
    type (treeNode                  ), pointer       :: nodeSatellite
    
    positionParent => node%parent        %position()
    nodeSatellite  => node%firstSatellite
    do while (associated(nodeSatellite))
       positionSatellite => nodeSatellite%position()
       call positionSatellite%positionSet(positionParent%position())
       call positionSatellite%velocitySet(positionParent%velocity())
       nodeSatellite     => nodeSatellite%sibling
    end do
    return
  end subroutine positionToHostNodePromote

  subroutine positionToHostNodesMerge(self,node)
    !!{
    Relocate the position of the node to that of its host on a node merger.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentPosition
    implicit none
    class(nodeOperatorPositionToHost), intent(inout) :: self
    type (treeNode                  ), intent(inout) :: node
    type (treeNode                  ), pointer       :: nodeHost
    class(nodeComponentPosition     ), pointer       :: positionHost, position
    
    nodeHost     => node    %parent
    position     => node    %position()
    positionHost => nodeHost%position()
    call position%positionSet(positionHost%position())
    call position%velocitySet(positionHost%velocity())
    return
  end subroutine positionToHostNodesMerge
