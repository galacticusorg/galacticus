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
Implements a node operator class that records the index of the node in which a node was last a satellite.
!!}

  !![
  <nodeOperator name="nodeOperatorIndexLastHost">
   <description>A node operator class that records the index of the node in which a node was last a satellite.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorIndexLastHost
     !!{
     A node operator class that records the index of the node in which a node was last a satellite.
     !!}
     private
     integer :: indexLastHostID
   contains
     final     ::                       indexLastHostDestructor
     procedure :: nodeTreeInitialize => indexLastHostNodeTreeInitialize
     procedure :: nodesMerge         => indexLastHostNodesMerge
     procedure :: nodePromote        => indexLastHostNodePromote
     procedure :: autoHook           => indexLastHostAutoHook
  end type nodeOperatorIndexLastHost

  interface nodeOperatorIndexLastHost
     !!{
     Constructors for the \refClass{nodeOperatorIndexLastHost} node operator class.
     !!}
     module procedure indexLastHostConstructorParameters
     module procedure indexLastHostConstructorInternal
  end interface nodeOperatorIndexLastHost

contains

  function indexLastHostConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorIndexLastHost} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodeOperatorIndexLastHost)                :: self
    type(inputParameters           ), intent(inout) :: parameters
    
    self=nodeOperatorIndexLastHost()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function indexLastHostConstructorParameters

  function indexLastHostConstructorInternal() result(self)
    !!{
    Constructor for the \refClass{nodeOperatorIndexLastHost} node operator class which takes a parameter set as input.
    !!}
    implicit none
    type(nodeOperatorIndexLastHost) :: self
    
    !![
    <addMetaProperty component="basic" name="nodeIndexLastHost" type="longInteger" id="self%indexLastHostID" isCreator="yes"/>
    !!]
    return
  end function indexLastHostConstructorInternal

  subroutine indexLastHostAutoHook(self)
    !!{
    Attach to various event hooks.
    !!}
    use :: Events_Hooks, only : openMPThreadBindingAtLevel, satelliteHostChangeEvent
    implicit none
    class(nodeOperatorIndexLastHost), intent(inout) :: self

    call satelliteHostChangeEvent%attach(self,indexLastHostSatelliteHostChange,openMPThreadBindingAtLevel,label='indexLastHost')
    return
  end subroutine indexLastHostAutoHook

  subroutine indexLastHostDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorIndexLastHost} node operator class.
    !!}
    use :: Events_Hooks, only : satelliteHostChangeEvent
    implicit none
    type(nodeOperatorIndexLastHost), intent(inout) :: self

    if (satelliteHostChangeEvent%isAttached(self,indexLastHostSatelliteHostChange)) call satelliteHostChangeEvent%detach(self,indexLastHostSatelliteHostChange)
    return
  end subroutine indexLastHostDestructor

  subroutine indexLastHostNodeTreeInitialize(self,node)
    !!{
    Initialize host node indices.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class  (nodeOperatorIndexLastHost), intent(inout), target  :: self
    type   (treeNode                 ), intent(inout), target  :: node
    class  (nodeComponentBasic       )               , pointer :: basic
    integer(kind_int8                )                         :: indexHost

    basic     => node%basic()
    indexHost =  -1_kind_int8
    if (node%isSatellite()) indexHost=node%parent%index()
    call basic%longIntegerRank0MetaPropertySet(self%indexLastHostID,indexHost)
    return
  end subroutine indexLastHostNodeTreeInitialize

  subroutine indexLastHostSatelliteHostChange(self,node)
    !!{
    Update the host halo index of this node in response to a change in host.
    !!}
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(*                 ), intent(inout)         :: self
    type (treeNode          ), intent(inout), target :: node
    class(nodeComponentBasic), pointer               :: basic

    select type (self)
    class is (nodeOperatorIndexLastHost)
       basic => node%basic()
       call basic%longIntegerRank0MetaPropertySet(self%indexLastHostID,node%parent%index())
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine indexLastHostSatelliteHostChange

  subroutine indexLastHostNodePromote(self,node)
    !!{
    Update the maximum host mass of this node as a result of node promotion.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(nodeOperatorIndexLastHost), intent(inout) :: self
    type (treeNode                 ), intent(inout) :: node
    type (treeNode                 ), pointer       :: nodeSatellite
    class(nodeComponentBasic       ), pointer       :: basic
    
    nodeSatellite => node%firstSatellite
    do while (associated(nodeSatellite))
       basic => nodeSatellite%basic()
       call basic%longIntegerRank0MetaPropertySet(self%indexLastHostID,nodeSatellite%parent%index())
       nodeSatellite => nodeSatellite%sibling
    end do
    return
  end subroutine indexLastHostNodePromote
  
  subroutine indexLastHostNodesMerge(self,node)
    !!{
    Update the maximum host mass of this node as a result of node merging.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(nodeOperatorIndexLastHost), intent(inout) :: self
    type (treeNode                 ), intent(inout) :: node
    class(nodeComponentBasic       ), pointer       :: basic
    
    basic => node%basic()
    call basic%longIntegerRank0MetaPropertySet(self%indexLastHostID,node%parent%index())
    return
  end subroutine indexLastHostNodesMerge
