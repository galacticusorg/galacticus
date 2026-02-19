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
  Implements a node operator class that tracks the maximum host halo mass which a node has occupied.
  !!}
  
  !![
  <nodeOperator name="nodeOperatorMassHostMaximum">
    <description>
      A node operator class that tracks the maximum host halo mass which a node has occupied. Intended to be paired with the
      \refClass{nodePropertyExtractorMassHostMaximum} class to extract these masses for output.
    </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorMassHostMaximum
     !!{
     A node operator class that tracks the maximum host halo mass which a node has occupied.
     !!}
     private
     integer :: massHostMaximumID
   contains
     !![
     <methods>
       <method method="update" description="Update the maximum host mass of this node."/>
     </methods>
     !!]
     final     ::                              massHostMaximumDestructor
     procedure :: nodeInitialize            => massHostMaximumNodeInitialize
     procedure :: nodesMerge                => massHostMaximumNodesMerge
     procedure :: nodePromote               => massHostMaximumNodePromote
     procedure :: update                    => massHostMaximumUpdate
     procedure :: differentialEvolutionPost => massHostMaximumDifferentialEvolutionPost
     procedure :: autoHook                  => massHostMaximumAutoHook
  end type nodeOperatorMassHostMaximum
  
  interface nodeOperatorMassHostMaximum
     !!{
     Constructors for the \refClass{nodeOperatorMassHostMaximum} node operator class.
     !!}
     module procedure massHostMaximumConstructorParameters
     module procedure massHostMaximumConstructorInternal
  end interface nodeOperatorMassHostMaximum
  
contains

  function massHostMaximumConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorMassHostMaximum} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodeOperatorMassHostMaximum)                :: self
    type(inputParameters            ), intent(inout) :: parameters
    
    self=nodeOperatorMassHostMaximum()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function massHostMaximumConstructorParameters

  function massHostMaximumConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorMassHostMaximum} node operator class.
    !!}
    implicit none
    type(nodeOperatorMassHostMaximum) :: self
    
    !![
    <addMetaProperty component="basic" name="massHostMaximum" id="self%massHostMaximumID" isEvolvable="no" isCreator="yes"/>
    !!]
    return
  end function massHostMaximumConstructorInternal

  subroutine massHostMaximumAutoHook(self)
    !!{
    Attach to various event hooks.
    !!}
    use :: Events_Hooks, only : openMPThreadBindingAtLevel, satelliteHostChangeEvent
    implicit none
    class(nodeOperatorMassHostMaximum), intent(inout) :: self

    call satelliteHostChangeEvent%attach(self,massHostMaximumSatelliteHostChange,openMPThreadBindingAtLevel,label='massHostMaximum')
    return
  end subroutine massHostMaximumAutoHook

  subroutine massHostMaximumDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorMassHostMaximum} node operator class.
    !!}
    use :: Events_Hooks, only : satelliteHostChangeEvent
    implicit none
    type(nodeOperatorMassHostMaximum), intent(inout) :: self

    if (satelliteHostChangeEvent%isAttached(self,massHostMaximumSatelliteHostChange)) call satelliteHostChangeEvent%detach(self,massHostMaximumSatelliteHostChange)
    return
  end subroutine massHostMaximumDestructor

  subroutine massHostMaximumUpdate(self,node,isMerging)
    !!{
    Update the maximum host mass of this node.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class  (nodeOperatorMassHostMaximum), intent(inout) :: self
    type   (treeNode                   ), intent(inout) :: node
    logical                             , intent(in   ) :: isMerging
    class  (nodeComponentBasic         ), pointer       :: basic    , basicHost
    
    if (.not.(node%isSatellite() .or. isMerging)) return
    basic     => node       %basic()
    basicHost => node%parent%basic()
    call basic%floatRank0MetaPropertySet(                                                                 &
         &                                                                       self%massHostMaximumID , &
         &                               max(                                                             &
         &                                   basic    %floatRank0MetaPropertyGet(self%massHostMaximumID), &
         &                                   basicHost%mass                     (                      )  &
         &                                  )                                                             &
         &                              )
    return
  end subroutine massHostMaximumUpdate

  subroutine massHostMaximumSatelliteHostChange(self,node)
    !!{
    Update the maximum host mass of this node in response to a change in host.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(*       ), intent(inout)         :: self
    type (treeNode), intent(inout), target :: node
    
    select type (self)
    class is (nodeOperatorMassHostMaximum)
      call self%update(node,isMerging=.false.)
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine massHostMaximumSatelliteHostChange
  
  subroutine massHostMaximumNodeInitialize(self,node)
    !!{
    Initialize the maximum host mass of this node.
    !!}
    implicit none
    class(nodeOperatorMassHostMaximum), intent(inout), target :: self
    type (treeNode                   ), intent(inout), target :: node
    
    call self%update(node,isMerging=.false.)
    return
  end subroutine massHostMaximumNodeInitialize

  subroutine massHostMaximumNodePromote(self,node)
    !!{
    Update the maximum host mass of this node as a result of node promotion.
    !!}
     use :: Galacticus_Nodes, only : nodeComponentBasic
   implicit none
    class(nodeOperatorMassHostMaximum), intent(inout) :: self
    type (treeNode                   ), intent(inout) :: node
    type (treeNode                   ), pointer       :: nodeSatellite
    class(nodeComponentBasic         ), pointer       :: basicSatellite, basicHost

    call self%update(node,isMerging=.false.)
    nodeSatellite => node%firstSatellite
    do while (associated(nodeSatellite))
       basicSatellite => nodeSatellite       %basic()
       basicHost      => node         %parent%basic()
       call basicSatellite%floatRank0MetaPropertySet(                                                                      &
            &                                                                                     self%massHostMaximumID , &
            &                                        max(                                                                  &
            &                                            basicSatellite%floatRank0MetaPropertyGet(self%massHostMaximumID), &
            &                                            basicHost     %mass                     (                      )  &
            &                                           )                                                                  &
            &                                       )  
       nodeSatellite => nodeSatellite%sibling
    end do
    return
  end subroutine massHostMaximumNodePromote
  
  subroutine massHostMaximumNodesMerge(self,node)
    !!{
    Update the maximum host mass of this node as a result of node merging.
    !!}
    implicit none
    class(nodeOperatorMassHostMaximum), intent(inout) :: self
    type (treeNode                   ), intent(inout) :: node
    
    call self%update(node,isMerging=.true.)
    return
  end subroutine massHostMaximumNodesMerge

  subroutine massHostMaximumDifferentialEvolutionPost(self,node)
    !!{
    Update the maximum host mass of this node, and of any satellite nodes.
    !!}
    implicit none
    class(nodeOperatorMassHostMaximum), intent(inout) :: self
    type (treeNode                   ), intent(inout) :: node
    type (treeNode                   ), pointer       :: nodeSatellite
    
    call self%update(node,isMerging=.false.)
    nodeSatellite => node%firstSatellite
    do while (associated(nodeSatellite))
       call self%update(nodeSatellite,isMerging=.false.)
       nodeSatellite => nodeSatellite%sibling
    end do
    return
  end subroutine massHostMaximumDifferentialEvolutionPost
