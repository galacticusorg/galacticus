!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements an extension to the very simple hot halo node component by including an outflowed reservoir
!% with delayed reincorporation.

module Node_Component_Hot_Halo_VS_Delayed
  !% Implements an extension to the very simple hot halo node component by including an outflowed reservoir
  !% with delayed reincorporation.
  use Galacticus_Nodes
  implicit none
  private
  public :: Node_Component_Hot_Halo_VS_Delayed_Node_Merger     , Node_Component_Hot_Halo_VS_Delayed_Rate_Compute   , &
       &    Node_Component_Hot_Halo_VS_Delayed_Scale_Set       , Node_Component_Hot_Halo_VS_Delayed_Tree_Initialize, &
       &    Node_Component_Hot_Halo_VS_Delayed_Satellite_Merger, Node_Component_Hot_Halo_VS_Delayed_Promote        , &
       &    Node_Component_Hot_Halo_VS_Delayed_Post_Evolve

  !# <component>
  !#  <class>hotHalo</class>
  !#  <name>verySimpleDelayed</name>
  !#  <extends>
  !#   <class>hotHalo</class>
  !#   <name>verySimple</name>
  !#  </extends>
  !#  <isDefault>no</isDefault>
  !#  <properties>
  !#   <property>
  !#     <name>outflowedMass</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#     <output unitsInSI="massSolar" comment="Mass of gas in the hot halo."/>
  !#   </property>
  !#  </properties>
  !# </component>

  ! Record of whether this module has been initialized.
  logical :: moduleInitialized=.false.

contains

  subroutine Node_Component_Hot_Halo_VS_Delayed_Initialize()
    !% Initializes the very simple hot halo component module.
    implicit none
    type(nodeComponentHotHaloVerySimpleDelayed) :: hotHaloComponent

    ! Initialize the module if necessary.
    !$omp critical (Node_Component_Hot_Halo_VS_Delayed_Initialize)
    if (.not.moduleInitialized) then
       ! Bind outflowing material pipes to the functions that will handle input of outflowing material to the hot halo.
       call hotHaloComponent%outflowingMassRateFunction(Node_Component_Hot_Halo_VS_Delayed_Outflowing_Mass_Rate)
       ! Record that the module is now initialized.
       moduleInitialized=.true.
    end if
    !$omp end critical (Node_Component_Hot_Halo_VS_Delayed_Initialize)
    return
  end subroutine Node_Component_Hot_Halo_VS_Delayed_Initialize

  subroutine Node_Component_Hot_Halo_VS_Delayed_Outflowing_Mass_Rate(self,rate,interrupt,interruptProcedure)
    !% Accept outflowing gas from a galaxy and deposit it into very simple hot halo.
    implicit none
    class           (nodeComponentHotHalo), intent(inout)                    :: self
    double precision                      , intent(in   )                    :: rate
    logical                               , intent(inout), optional          :: interrupt
    procedure       (                    ), intent(inout), optional, pointer :: interruptProcedure

    ! Funnel the outflow gas into the outflowed reservoir.
    call self%outflowedMassRate(rate)
    return
  end subroutine Node_Component_Hot_Halo_VS_Delayed_Outflowing_Mass_Rate
  
  !# <rateComputeTask>
  !#  <unitName>Node_Component_Hot_Halo_VS_Delayed_Rate_Compute</unitName>
  !# </rateComputeTask>
  subroutine Node_Component_Hot_Halo_VS_Delayed_Rate_Compute(node,interrupt,interruptProcedure)
    !% Compute the very simple hot halo component mass rate of change.
    use Hot_Halo_Outflows_Reincorporations
    implicit none
    type            (treeNode                          ), intent(inout), pointer :: node
    logical                                             , intent(inout)          :: interrupt
    procedure       (                                  ), intent(inout), pointer :: interruptProcedure
    class           (nodeComponentHotHalo              )               , pointer :: hotHalo
    class           (hotHaloOutflowReincorporationClass)               , pointer :: hotHaloOutflowReincorporation_
    double precision                                                             :: outflowReturnRate

    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    select type (hotHalo)
    class is (nodeComponentHotHaloVerySimpleDelayed)
       ! Move outflowed material back to the hot reservoir.
       hotHaloOutflowReincorporation_ => hotHaloOutflowReincorporation      (    )
       outflowReturnRate              =  hotHaloOutflowReincorporation_%rate(node)
       call hotHalo%outflowedMassRate(-outflowReturnRate)
       call hotHalo%massRate         (+outflowReturnRate)
    end select
    return
  end subroutine Node_Component_Hot_Halo_VS_Delayed_Rate_Compute

  !# <scaleSetTask>
  !#  <unitName>Node_Component_Hot_Halo_VS_Delayed_Scale_Set</unitName>
  !# </scaleSetTask>
  subroutine Node_Component_Hot_Halo_VS_Delayed_Scale_Set(node)
    !% Set scales for properties of {\tt node}.
    implicit none
    type            (treeNode            ), intent(inout), pointer :: node
    double precision                      , parameter              :: scaleMassRelative=1.0d-2
    class           (nodeComponentHotHalo)               , pointer :: hotHalo
    class           (nodeComponentBasic  )               , pointer :: basic
    double precision                                               :: massVirial

    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    ! Ensure that it is of our class.
    select type (hotHalo)
    class is (nodeComponentHotHaloVerySimpleDelayed)
       ! The the basic component.
       basic      => node %basic()
       ! Get virial properties.
       massVirial =  basic%mass ()
       ! Set the scale.
       call hotHalo%outflowedMassScale(massVirial*scaleMassRelative)
    end select
    return
  end subroutine Node_Component_Hot_Halo_VS_Delayed_Scale_Set

  !# <mergerTreeInitializeTask>
  !#  <unitName>Node_Component_Hot_Halo_VS_Delayed_Tree_Initialize</unitName>
  !#  <after>Node_Component_Hot_Halo_Very_Simple_Tree_Initialize</after>
  !# </mergerTreeInitializeTask>
  subroutine Node_Component_Hot_Halo_VS_Delayed_Tree_Initialize(node)
    !% Initialize the contents of the very simple hot halo component.
    use Cosmology_Parameters
    implicit none
    type (treeNode                ), intent(inout), pointer :: node
    class(nodeComponentHotHalo    )               , pointer :: hotHalo

    ! If the very simple hot halo is not active, then return immediately.
    if (.not.defaultHotHaloComponent%verySimpleDelayedIsActive()) return

    ! Ensure that this module has been initialized.
    call Node_Component_Hot_Halo_VS_Delayed_Initialize()
    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    ! Ensure that it is of our class.
    select type (hotHalo)
    type is (nodeComponentHotHaloVerySimpleDelayed)
       call hotHalo%outflowedMassSet(0.0d0)
    end select
    return
  end subroutine Node_Component_Hot_Halo_VS_Delayed_Tree_Initialize

  !# <satelliteMergerTask>
  !#  <unitName>Node_Component_Hot_Halo_VS_Delayed_Satellite_Merger</unitName>
  !# </satelliteMergerTask>
  subroutine Node_Component_Hot_Halo_VS_Delayed_Satellite_Merger(node)
    !% Remove any hot halo associated with {\tt node} before it merges with its host halo.
    implicit none
    type (treeNode            ), intent(inout), pointer :: node
    type (treeNode            )               , pointer :: hostNode
    class(nodeComponentHotHalo)               , pointer :: hostHotHalo, hotHalo

    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    ! Ensure that it is of unspecified class.
    select type (hotHalo)
    class is (nodeComponentHotHaloVerySimpleDelayed)
       ! Find the node to merge with.
       hostNode    => node    %mergesWith()
       hostHotHalo => hostNode%hotHalo   ()
       ! Move the hot halo to the host.
       call hostHotHalo%outflowedMassSet(                             &
            &                            +hostHotHalo%outflowedMass() &
            &                            +    hotHalo%outflowedMass() &
            &                           )
       call     hotHalo%outflowedMassSet(                             &
            &                            +0.0d0                       &
            &                           )
    end select
    return
  end subroutine Node_Component_Hot_Halo_VS_Delayed_Satellite_Merger

  !# <nodePromotionTask>
  !#  <unitName>Node_Component_Hot_Halo_VS_Delayed_Promote</unitName>
  !# </nodePromotionTask>
  subroutine Node_Component_Hot_Halo_VS_Delayed_Promote(node)
    !% Ensure that {\tt node} is ready for promotion to its parent. In this case, we simply update the hot halo mass of {\tt
    !% node} to account for any hot halo already in the parent.
    implicit none
    type (treeNode            ), intent(inout), pointer :: node
    type (treeNode            )               , pointer :: parentNode
    class(nodeComponentHotHalo)               , pointer :: parentHotHalo, hotHalo

    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    ! Ensure that it is of specified class.
    select type (hotHalo)
    class is (nodeComponentHotHaloVerySimpleDelayed)
       ! Get the parent node of this node.
       parentNode    => node      %parent
       parentHotHalo => parentNode%hotHalo()
       ! If the parent node has a hot halo component, then add it to that of this node, and perform other changes needed prior to
       ! promotion.
       select type (parentHotHalo)
       class is (nodeComponentHotHaloVerySimpleDelayed)
          call hotHalo%outflowedMassSet(                               &
               &                        +      hotHalo%outflowedMass() &
               &                        +parentHotHalo%outflowedMass() &
               &                       )
       end select
    end select
    return
  end subroutine Node_Component_Hot_Halo_VS_Delayed_Promote

  !# <postEvolveTask>
  !#  <unitName>Node_Component_Hot_Halo_VS_Delayed_Post_Evolve</unitName>
  !# </postEvolveTask>
  subroutine Node_Component_Hot_Halo_VS_Delayed_Post_Evolve(node)
    !% Do processing of the node required after evolution.
    implicit none
    type (treeNode            ), intent(inout), pointer :: node
    type (treeNode            )               , pointer :: parentNode
    class(nodeComponentHotHalo)               , pointer :: parentHotHalo, hotHalo

    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    select type (hotHalo)
    class is (nodeComponentHotHaloVerySimpleDelayed)
       ! Check if this node is a satellite.
       if (node%isSatellite()) then
          ! Transfer any outflowed gas to the hot halo of the parent node.
          parentNode => node%parent
          do while (parentNode%isSatellite())
             parentNode => parentNode%parent
          end do
          parentHotHalo => parentNode%hotHalo()
          call parentHotHalo%outflowedMassSet(parentHotHalo%outflowedMass()+hotHalo%outflowedMass())
          call       hotHalo%outflowedMassSet(                                                0.0d0)
       end if
    end select
    return
  end subroutine Node_Component_Hot_Halo_VS_Delayed_Post_Evolve

  !# <nodeMergerTask>
  !#  <unitName>Node_Component_Hot_Halo_VS_Delayed_Node_Merger</unitName>
  !# </nodeMergerTask>
  subroutine Node_Component_Hot_Halo_VS_Delayed_Node_Merger(node)
    !% Starve {\tt node} by transferring its hot halo to its parent.
    implicit none
    type (treeNode            ), intent(inout), pointer :: node
    type (treeNode            )               , pointer :: parentNode
    class(nodeComponentHotHalo)               , pointer :: parentHotHalo, hotHalo

    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    select type (hotHalo)
    class is (nodeComponentHotHaloVerySimpleDelayed)
       ! Find the parent node and its hot halo component.
       parentNode    => node      %parent
       parentHotHalo => parentNode%hotHalo()
       ! Move the hot halo to the parent. We leave the hot halo in place even if it is starved, since outflows will accumulate
       ! to this hot halo (and will be moved to the parent at the end of the evolution timestep).
       call parentHotHalo%outflowedMassSet(parentHotHalo%outflowedMass()+hotHalo%outflowedMass())
       call       hotHalo%outflowedMassSet(                                                0.0d0)
    end select
    return
  end subroutine Node_Component_Hot_Halo_VS_Delayed_Node_Merger

end module Node_Component_Hot_Halo_VS_Delayed
