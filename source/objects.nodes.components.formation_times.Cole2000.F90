!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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
  use Galacticus_Nodes
  implicit none
  private
  public :: Node_Component_Formation_Times_Cole2000_Initialize    , Node_Component_Formation_Time_Cole2000_Rate_Compute  , &
       &    Node_Component_Formation_Time_Cole2000_Tree_Initialize, Node_Component_Formation_Time_Cole2000_Node_Promotion

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

  ! Record of whether this module has been initialized.
  logical          :: moduleInitialized             =.false.

  ! Factor by which mass must increase to trigger a new formation event.
  double precision :: haloReformationMassFactor

  ! Switch indicating whether or not halo reformation should only be checked for at node promotion events.
  logical          :: haloReformationOnPromotionOnly

contains

  subroutine Node_Component_Formation_Times_Cole2000_Initialize()
    !% Initializes the tree node formation time tracking module.
    use Input_Parameters
    implicit none

    ! Initialize the module if necessary.
    !$omp critical (Node_Component_Formation_Times_Cole2000_Initialize)
    if (.not.moduleInitialized) then
       !@ <inputParameter>
       !@   <name>haloReformationMassFactor</name>
       !@   <defaultValue>2.0</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Factor by which halo mass must have increased to trigger a new formation event.
       !@   </description>
       !@   <type>double</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('haloReformationMassFactor',haloReformationMassFactor,defaultValue=2.0d0)
       !@ <inputParameter>
       !@   <name>haloReformationOnPromotionOnly</name>
       !@   <defaultValue>false</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Specifies whether halo reformation should occur only at node promotion events, or at the precise time that
       !@     the halo mass has increased sufficiently in mass.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('haloReformationOnPromotionOnly',haloReformationOnPromotionOnly,defaultValue=.false.)

       ! Record that the module is now initialized.
       moduleInitialized=.true.
    end if
    !$omp end critical (Node_Component_Formation_Times_Cole2000_Initialize)
    return
  end subroutine Node_Component_Formation_Times_Cole2000_Initialize

  !# <rateComputeTask>
  !#  <unitName>Node_Component_Formation_Time_Cole2000_Rate_Compute</unitName>
  !# </rateComputeTask>
  subroutine Node_Component_Formation_Time_Cole2000_Rate_Compute(node,odeConverged,interrupt,interruptProcedure)
    !% Check for need to update the formation time of a node in the {\normalfont \ttfamily Cole2000} formation time component.
    implicit none
    type     (treeNode                   ), intent(inout), pointer :: node
    logical                               , intent(in   )          :: odeConverged
    logical                               , intent(inout)          :: interrupt
    procedure(interruptTask              ), intent(inout), pointer :: interruptProcedure
    class    (nodeComponentFormationTime )               , pointer :: formationTime
    class    (nodeComponentBasic         )               , pointer :: basicFormation    , basic
    !GCC$ attributes unused :: odeConverged
    
    ! Get the hot halo component.
    formationTime => node%formationTime()
    if (defaultFormationTimeComponent%cole2000IsActive()) then
       ! Check if the halo has grown sufficiently in mass to trigger a new formation event.
       if (.not.haloReformationOnPromotionOnly) then
          basic      => node              %basic()
          basicFormation => node%formationNode%basic()
          if (basic%mass() > haloReformationMassFactor*basicFormation%mass()) then
             interrupt=.true.
             interruptProcedure => Node_Component_Formation_Time_Cole2000_Create
             return
          end if
       end if
    end if
    return
  end subroutine Node_Component_Formation_Time_Cole2000_Rate_Compute

  !# <nodePromotionTask>
  !#  <unitName>Node_Component_Formation_Time_Cole2000_Node_Promotion</unitName>
  !#  <after>re:Node_Component_Hot_Halo_.*_Promote</after>
  !# </nodePromotionTask>
  subroutine Node_Component_Formation_Time_Cole2000_Node_Promotion(node)
    implicit none
    type (treeNode                  ), intent(inout), pointer :: node
    class(nodeComponentFormationTime)               , pointer :: formationTime
    class(nodeComponentBasic        )               , pointer :: basicFormation, basicParent

    ! Get the formation time component.
    formationTime => node%formationTime()
    ! Ensure that it is of specified class.
    select type (formationTime)
    class is (nodeComponentFormationTimeCole2000)
       if (haloReformationOnPromotionOnly) then
          basicParent    => node%parent       %basic()
          basicFormation => node%formationNode%basic()
          if (basicParent%mass() > haloReformationMassFactor*basicFormation%mass()) &
               & call Node_Component_Formation_Time_Cole2000_Create(node)
       end if
    end select
    return
  end subroutine Node_Component_Formation_Time_Cole2000_Node_Promotion

  subroutine Node_Component_Formation_Time_Cole2000_Create(node)
    !% Creates a halo formation time component for {\normalfont \ttfamily node}. This function is also used to ``reform'' the halo, since it
    !% simply resets the formation time and mass to the current values.
    use Events_Halo_Formation
    implicit none
    type (treeNode                  ), intent(inout), pointer :: node
    class(nodeComponentFormationTime)               , pointer :: formationTime

    ! Ensure that this module has been initialized.
    call Node_Component_Formation_Times_Cole2000_Initialize()

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
  end subroutine Node_Component_Formation_Time_Cole2000_Create

  !# <mergerTreeInitializeTask>
  !#  <unitName>Node_Component_Formation_Time_Cole2000_Tree_Initialize</unitName>
  !# </mergerTreeInitializeTask>
  subroutine Node_Component_Formation_Time_Cole2000_Tree_Initialize(node)
    !% Initialize the formation node pointer for any childless node.
    implicit none
    type(treeNode), intent(inout), pointer :: node

    ! If this method is selected and the node has no child then initialize it.
    if (defaultFormationTimeComponent%cole2000IsActive().and..not.associated(node%firstChild)) &
         & call Node_Component_Formation_Time_Cole2000_Create(node)

    return
  end subroutine Node_Component_Formation_Time_Cole2000_Tree_Initialize

end module Node_Component_Formation_Times_Cole2000
