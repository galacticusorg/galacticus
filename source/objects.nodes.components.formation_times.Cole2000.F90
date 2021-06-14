!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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

!% Contains a module of halo formation time methods.

module Node_Component_Formation_Times_Cole2000
  !% Implement tracking of halo formation times.
  implicit none
  private
  public :: Node_Component_Formation_Times_Cole2000_Initialize         , Node_Component_Formation_Times_Cole2000_Rate_Compute     , &
       &    Node_Component_Formation_Times_Cole2000_Tree_Initialize    , Node_Component_Formation_Times_Cole2000_Thread_Initialize, &
       &    Node_Component_Formation_Times_Cole2000_Thread_Uninitialize

  !# <component>
  !#  <class>formationTime</class>
  !#  <name>Cole2000</name>
  !#  <isDefault>false</isDefault>
  !#  <properties>
  !#   <property>
  !#     <name>formationTime</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="false" isGettable="true" isEvolvable="false" />
  !#     <getFunction bindsTo="component">FormationTimeCole2000FormationTime</getFunction>
  !#   </property>
  !#  </properties>
  !#  <functions>objects.nodes.components.formation_times.Cole2000.bound_functions.inc</functions>
  !# </component>

  ! Factor by which mass must increase to trigger a new formation event.
  double precision :: haloReformationMassFactor

  ! Switch indicating whether or not halo reformation should only be checked for at node promotion events.
  logical          :: haloReformationOnPromotionOnly

contains

  !# <nodeComponentInitializationTask>
  !#  <unitName>Node_Component_Formation_Times_Cole2000_Initialize</unitName>
  !# </nodeComponentInitializationTask>
  subroutine Node_Component_Formation_Times_Cole2000_Initialize(parameters_)
    !% Initializes the tree node formation time tracking module.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters_

    !# <inputParameter>
    !#   <name>haloReformationMassFactor</name>
    !#   <defaultValue>2.0d0</defaultValue>
    !#   <description>Factor by which halo mass must have increased to trigger a new formation event.</description>
    !#   <source>parameters_</source>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>haloReformationOnPromotionOnly</name>
    !#   <defaultValue>.false.</defaultValue>
    !#   <description>Specifies whether halo reformation should occur only at node promotion events, or at the precise time that
    !#      the halo mass has increased sufficiently in mass.</description>
    !#   <source>parameters_</source>
    !# </inputParameter>
    return
  end subroutine Node_Component_Formation_Times_Cole2000_Initialize

  !# <nodeComponentThreadInitializationTask>
  !#  <unitName>Node_Component_Formation_Times_Cole2000_Thread_Initialize</unitName>
  !# </nodeComponentThreadInitializationTask>
  subroutine Node_Component_Formation_Times_Cole2000_Thread_Initialize(parameters_)
    !% Initializes the tree node hot halo methods module.
    use :: Events_Hooks    , only : nodePromotionEvent           , openMPThreadBindingAtLevel, dependencyRegEx, dependencyDirectionAfter
    use :: Galacticus_Nodes, only : defaultFormationTimeComponent
    use :: Input_Parameters, only : inputParameter               , inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters_
    type(dependencyRegEx), dimension(1)  :: dependencies
    !$GLC attributes unused :: parameters_
    
    ! Check if this implementation is selected.
    if (defaultFormationTimeComponent%cole2000IsActive()) then
       dependencies(1)=dependencyRegEx(dependencyDirectionAfter,'^nodeComponentHotHalo')
       call nodePromotionEvent%attach(defaultFormationTimeComponent,nodePromotion,openMPThreadBindingAtLevel,label='nodeComponentFormationTimeCole2000',dependencies=dependencies)
    end if
    return
  end subroutine Node_Component_Formation_Times_Cole2000_Thread_Initialize

  !# <nodeComponentThreadUninitializationTask>
  !#  <unitName>Node_Component_Formation_Times_Cole2000_Thread_Uninitialize</unitName>
  !# </nodeComponentThreadUninitializationTask>
  subroutine Node_Component_Formation_Times_Cole2000_Thread_Uninitialize()
    !% Uninitializes the tree node hot halo methods module.
    use :: Events_Hooks    , only : nodePromotionEvent
    use :: Galacticus_Nodes, only : defaultFormationTimeComponent
    implicit none

    if (defaultFormationTimeComponent%cole2000IsActive()) &
         & call nodePromotionEvent%detach(defaultFormationTimeComponent,nodePromotion)
    return
  end subroutine Node_Component_Formation_Times_Cole2000_Thread_Uninitialize

  !# <rateComputeTask>
  !#  <unitName>Node_Component_Formation_Times_Cole2000_Rate_Compute</unitName>
  !# </rateComputeTask>
  subroutine Node_Component_Formation_Times_Cole2000_Rate_Compute(node,interrupt,interruptProcedure,propertyType)
    !% Check for need to update the formation time of a node in the {\normalfont \ttfamily Cole2000} formation time component.
    use :: Galacticus_Nodes, only : defaultFormationTimeComponent, interruptTask, nodeComponentBasic, nodeComponentFormationTime, &
          &                         propertyTypeInactive         , treeNode
    implicit none
    type     (treeNode                   ), intent(inout)          :: node
    logical                               , intent(inout)          :: interrupt
    procedure(interruptTask              ), intent(inout), pointer :: interruptProcedure
    integer                               , intent(in   )          :: propertyType
    class    (nodeComponentFormationTime )               , pointer :: formationTime
    class    (nodeComponentBasic         )               , pointer :: basicFormation    , basic

    ! Return immediately if inactive variables are requested.
    if (propertyType == propertyTypeInactive) return
    ! Return immediately if this class is not in use.
    if (.not.defaultFormationTimeComponent%cole2000IsActive()) return
    ! Get the hot halo component.
    formationTime => node%formationTime()
    if (defaultFormationTimeComponent%cole2000IsActive()) then
       ! Check if the halo has grown sufficiently in mass to trigger a new formation event.
       if (.not.haloReformationOnPromotionOnly) then
          basic      => node              %basic()
          basicFormation => node%formationNode%basic()
          if (basic%mass() > haloReformationMassFactor*basicFormation%mass()) then
             interrupt=.true.
             interruptProcedure => Node_Component_Formation_Times_Cole2000_Create
             return
          end if
       end if
    end if
    return
  end subroutine Node_Component_Formation_Times_Cole2000_Rate_Compute

  subroutine nodePromotion(self,node)
    !% Check if the node has undergone a formation event.
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentFormationTime, treeNode
    implicit none
    class(*                         ), intent(inout) :: self
    type (treeNode                  ), intent(inout) :: node
    class(nodeComponentFormationTime), pointer       :: formationTime
    class(nodeComponentBasic        ), pointer       :: basicFormation, basicParent
    !$GLC attributes unused :: self
    
    if (haloReformationOnPromotionOnly) then
       formationTime  => node              %formationTime()
       basicParent    => node%parent       %basic        ()
       basicFormation => node%formationNode%basic        ()
       if (basicParent%mass() > haloReformationMassFactor*basicFormation%mass()) &
            & call Node_Component_Formation_Times_Cole2000_Create(node)
    end if
    return
  end subroutine nodePromotion

  subroutine Node_Component_Formation_Times_Cole2000_Create(node)
    !% Creates a halo formation time component for {\normalfont \ttfamily node}. This function is also used to ``reform'' the halo, since it
    !% simply resets the formation time and mass to the current values.
    use :: Events_Halo_Formation, only : Event_Halo_Formation
    use :: Galacticus_Nodes     , only : nodeComponentFormationTime, treeNode
    implicit none
    type (treeNode                  ), intent(inout), target :: node
    class(nodeComponentFormationTime), pointer               :: formationTime

    ! Trigger a halo formation event.
    call Event_Halo_Formation(node)

    ! Create the component.
    formationTime => node%formationTime(autoCreate=.true.)
    ! Make a copy of the formation node, and decouple it from the tree, using the parentNode pointer to point to the node of which
    ! it is the formation node.
    if (associated(node%formationNode)) call node%formationNode%destroy()
    allocate(node%formationNode)
    call node%copyNodeTo(node%formationNode,skipFormationNode=.true.)
    node%formationNode%parent         => node
    node%formationNode%firstChild     => null()
    node%formationNode%sibling        => null()
    node%formationNode%firstSatellite => null()
    node%formationNode%firstMergee    => null()
    node%formationNode%mergeTarget    => null()
    node%formationNode%siblingMergee  => null()
    node%formationNode%formationNode  => null()
    return
  end subroutine Node_Component_Formation_Times_Cole2000_Create

  !# <mergerTreeInitializeTask>
  !#  <unitName>Node_Component_Formation_Times_Cole2000_Tree_Initialize</unitName>
  !# </mergerTreeInitializeTask>
  subroutine Node_Component_Formation_Times_Cole2000_Tree_Initialize(node)
    !% Initialize the formation node pointer for any childless node.
    use :: Galacticus_Nodes, only : defaultFormationTimeComponent, treeNode
    implicit none
    type(treeNode), intent(inout), pointer :: node

    ! If this method is selected and the node has no child then initialize it.
    if (defaultFormationTimeComponent%cole2000IsActive().and..not.associated(node%firstChild)) &
         & call Node_Component_Formation_Times_Cole2000_Create(node)

    return
  end subroutine Node_Component_Formation_Times_Cole2000_Tree_Initialize

end module Node_Component_Formation_Times_Cole2000
