!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
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

!% Contains a module which implements a very simple hot halo node component.

module Node_Component_Hot_Halo_Very_Simple
  !% Implements a very simple hot halo node component.
  use Galacticus_Nodes
  implicit none
  private
  public :: Node_Component_Hot_Halo_Very_Simple_Reset            , Node_Component_Hot_Halo_Very_Simple_Rate_Compute   , &
       &    Node_Component_Hot_Halo_Very_Simple_Scale_Set        , Node_Component_Hot_Halo_Very_Simple_Tree_Initialize, &
       &    Node_Component_Hot_Halo_Very_Simple_Satellite_Merging, Node_Component_Hot_Halo_Very_Simple_Promote        , &
       &    Node_Component_Hot_Halo_Very_Simple_Post_Evolve      , Node_Component_Hot_Halo_Very_Simple_Node_Merger

  !# <component>
  !#  <class>hotHalo</class>
  !#  <name>verySimple</name>
  !#  <isDefault>no</isDefault>
  !#  <properties>
  !#   <property>
  !#     <name>mass</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#     <output unitsInSI="massSolar" comment="Mass of gas in the hot halo."/>
  !#   </property>
  !#   <property>
  !#     <name>unaccretedMass</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" />
  !#     <output unitsInSI="massSolar" comment="Mass of gas that failed to accrete into the hot halo."/>
  !#   </property>
  !#   <property>
  !#     <name>outflowingMass</name>
  !#     <attributes isSettable="false" isGettable="false" isEvolvable="true" isDeferred="rate" isVirtual="true" />
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#   </property>
  !#   <property>
  !#     <name>hotHaloCoolingMass</name>
  !#     <attributes isSettable="false" isGettable="false" isEvolvable="true" isDeferred="rate" bindsTo="top" isVirtual="true" />
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#   </property>
  !#   <property>
  !#     <name>outerRadius</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="false" isGettable="true" isEvolvable="false" isDeferred="get" isVirtual="true" />
  !#   </property>
  !#  </properties>
  !# </component>

  ! Quantities stored to avoid repeated computation.
  logical          :: gotCoolingRate   =.false.
  double precision :: coolingRate
  !$omp threadprivate(gotCoolingRate,coolingRate)
  ! Record of whether this module has been initialized.
  logical          :: moduleInitialized=.false.

contains

  subroutine Node_Component_Hot_Halo_Very_Simple_Initialize()
    !% Initializes the very simple hot halo component module.
    implicit none
    type(nodeComponentHotHaloVerySimple) :: hotHaloComponent

    ! Initialize the module if necessary.
    !$omp critical (Node_Component_Hot_Halo_Very_Simple_Initialize)
    if (.not.moduleInitialized) then

       ! Bind outflowing material pipes to the functions that will handle input of outflowing material to the hot halo.
       call hotHaloComponent%outflowingMassRateFunction(Node_Component_Hot_Halo_Very_Simple_Outflowing_Mass_Rate)

       ! Bind outer radius function.
       call hotHaloComponent%       outerRadiusFunction(Node_Component_Hot_Halo_Very_Simple_Outer_Radius        )

       ! Record that the module is now initialized.
       moduleInitialized=.true.

    end if
    !$omp end critical (Node_Component_Hot_Halo_Very_Simple_Initialize)
    return
  end subroutine Node_Component_Hot_Halo_Very_Simple_Initialize

  !# <calculationResetTask>
  !# <unitName>Node_Component_Hot_Halo_Very_Simple_Reset</unitName>
  !# </calculationResetTask>
  subroutine Node_Component_Hot_Halo_Very_Simple_Reset(thisNode)
    !% Remove memory of stored computed values as we're about to begin computing derivatives anew.
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    gotCoolingRate=.false.
    return
  end subroutine Node_Component_Hot_Halo_Very_Simple_Reset

  subroutine Node_Component_Hot_Halo_Very_Simple_Push_To_Cooling_Pipes(thisNode,massRate,interrupt,interruptProcedure)
    !% Push mass through the cooling pipes at the given rate.
    implicit none
    type            (treeNode            ), intent(inout), pointer :: thisNode
    double precision                      , intent(in   )          :: massRate
    logical                               , intent(inout)          :: interrupt
    procedure       (                    ), intent(inout), pointer :: interruptProcedure
    class           (nodeComponentHotHalo)               , pointer :: thisHotHaloComponent

    ! Get the hot halo component.
    thisHotHaloComponent => thisNode%hotHalo()
    select type (thisHotHaloComponent)
    class is (nodeComponentHotHaloVerySimple)
       ! Ignore zero rates.
       if (massRate /= 0.0d0 .and. thisHotHaloComponent%mass() > 0.0d0) then
          ! Remove mass from the hot component.
          call    thisHotHaloComponent%massRate       (-massRate                             )
          ! Pipe the mass rate to whatever component claimed it.
          if (thisHotHaloComponent%hotHaloCoolingMassRateIsAttached()) then
             call thisHotHaloComponent%hotHaloCoolingMassRate(+massRate,interrupt,interruptProcedure)
             if (interrupt) return
          end if
       end if
    end select
    return
  end subroutine Node_Component_Hot_Halo_Very_Simple_Push_To_Cooling_Pipes

  subroutine Node_Component_Hot_Halo_Very_Simple_Outflowing_Mass_Rate(self,rate,interrupt,interruptProcedure)
    !% Accept outflowing gas from a galaxy and deposit it into very simple hot halo.
    implicit none
    class           (nodeComponentHotHalo), intent(inout)                    :: self
    double precision                      , intent(in   )                    :: rate
    logical                               , intent(inout), optional          :: interrupt
    procedure       (                    ), intent(inout), optional, pointer :: interruptProcedure

    ! Funnel the outflow gas into the hot halo.
    call self%massRate(rate)
    return
  end subroutine Node_Component_Hot_Halo_Very_Simple_Outflowing_Mass_Rate
  
  double precision function Node_Component_Hot_Halo_Very_Simple_Outer_Radius(self)
    !% Return the outer radius of the hot halo. Assumes a simple model in which this always equals the virial radius.
    use Dark_Matter_Halo_Scales
    implicit none
    class(nodeComponentHotHaloVerySimple), intent(inout) :: self
    class(darkMatterHaloScaleClass      ), pointer       :: darkMatterHaloScale_

    darkMatterHaloScale_ => darkMatterHaloScale()
    Node_Component_Hot_Halo_Very_Simple_Outer_Radius=darkMatterHaloScale_%virialRadius(self%hostNode)
    return
  end function Node_Component_Hot_Halo_Very_Simple_Outer_Radius

  !# <rateComputeTask>
  !#  <unitName>Node_Component_Hot_Halo_Very_Simple_Rate_Compute</unitName>
  !# </rateComputeTask>
  subroutine Node_Component_Hot_Halo_Very_Simple_Rate_Compute(thisNode,interrupt,interruptProcedure)
    !% Compute the very simple hot halo component mass rate of change.
    use Accretion_Halos
    implicit none
    type            (treeNode            ), intent(inout), pointer :: thisNode
    logical                               , intent(inout)          :: interrupt
    procedure       (                    ), intent(inout), pointer :: interruptProcedure
    class           (nodeComponentHotHalo)               , pointer :: thisHotHaloComponent
    class           (accretionHaloClass  )               , pointer :: accretionHalo_
    double precision                                               :: massAccretionRate   , failedMassAccretionRate
    
    ! Get the hot halo component.
    thisHotHaloComponent => thisNode%hotHalo()
    select type (thisHotHaloComponent)
    class is (nodeComponentHotHaloVerySimple)
       ! Get required objects.
       accretionHalo_ => accretionHalo()
       ! Find the rate of gas mass accretion onto the halo.
       massAccretionRate      =accretionHalo_%accretionRate      (thisNode,accretionModeTotal)
       failedMassAccretionRate=accretionHalo_%failedAccretionRate(thisNode,accretionModeTotal)
       ! Apply accretion rates.
       if (      massAccretionRate > 0.0d0 .or. thisHotHaloComponent%mass() > 0.0d0) &
            & call thisHotHaloComponent%          massRate(      massAccretionRate,interrupt,interruptProcedure)
       if (failedMassAccretionRate > 0.0d0 .or. thisHotHaloComponent%mass() > 0.0d0) &
            & call thisHotHaloComponent%unaccretedMassRate(failedMassAccretionRate,interrupt,interruptProcedure)
       ! Next compute the cooling rate in this halo.
       call Node_Component_Hot_Halo_Very_Simple_Cooling_Rate         (thisNode                                         )
       ! Pipe the cooling rate to which ever component claimed it.
       call Node_Component_Hot_Halo_Very_Simple_Push_To_Cooling_Pipes(thisNode,coolingRate,interrupt,interruptProcedure)
    end select
    return
  end subroutine Node_Component_Hot_Halo_Very_Simple_Rate_Compute

  !# <scaleSetTask>
  !#  <unitName>Node_Component_Hot_Halo_Very_Simple_Scale_Set</unitName>
  !# </scaleSetTask>
  subroutine Node_Component_Hot_Halo_Very_Simple_Scale_Set(thisNode)
    !% Set scales for properties of {\normalfont \ttfamily thisNode}.
    implicit none
    type            (treeNode            ), intent(inout), pointer :: thisNode
    double precision                      , parameter              :: scaleMassRelative   =1.0d-2
    class           (nodeComponentHotHalo)               , pointer :: thisHotHaloComponent
    class           (nodeComponentBasic  )               , pointer :: thisBasicComponent
    double precision                                               :: massVirial

    ! Get the hot halo component.
    thisHotHaloComponent => thisNode%hotHalo()
    ! Ensure that it is of the very simple class.
    select type (thisHotHaloComponent)
    class is (nodeComponentHotHaloVerySimple)
       ! The the basic component.
       thisBasicComponent => thisNode%basic()
       ! Get virial properties.
       massVirial=thisBasicComponent%mass()
       ! Set the scale.
       call thisHotHaloComponent%          massScale(massVirial*scaleMassRelative)
       call thisHotHaloComponent%unaccretedMassScale(massVirial*scaleMassRelative)
    end select
    return
  end subroutine Node_Component_Hot_Halo_Very_Simple_Scale_Set

  !# <mergerTreeInitializeTask>
  !#  <unitName>Node_Component_Hot_Halo_Very_Simple_Tree_Initialize</unitName>
  !#  <after>darkMatterProfile</after>
  !# </mergerTreeInitializeTask>
  subroutine Node_Component_Hot_Halo_Very_Simple_Tree_Initialize(thisNode)
    !% Initialize the contents of the very simple hot halo component.
    use Cosmology_Parameters
    use Accretion_Halos
    implicit none
    type            (treeNode            ), intent(inout), pointer :: thisNode
    type            (treeNode            )               , pointer :: childNode
    class           (nodeComponentHotHalo)               , pointer :: currentHotHaloComponent, thisHotHaloComponent
    class           (nodeComponentBasic  )               , pointer :: childBasicComponent    , currentBasicComponent
    class           (accretionHaloClass  )               , pointer :: accretionHalo_
    double precision                                               :: hotHaloMass            , failedHotHaloMass

    ! If the very simple hot halo is not active, then return immediately.
    if (associated(thisNode%firstChild).or..not.defaultHotHaloComponent%verySimpleIsActive()) return

    ! Ensure that this module has been initialized.
    call Node_Component_Hot_Halo_Very_Simple_Initialize()

    ! Get the hot halo component.
    currentHotHaloComponent => thisNode%hotHalo()
    ! Ensure that it is of unspecified class.
    select type (currentHotHaloComponent)
    type is (nodeComponentHotHalo)
       ! Get required objects.
       accretionHalo_ => accretionHalo()
       ! Get the mass of hot gas accreted and the mass that failed to accrete.
       hotHaloMass      =accretionHalo_%accretedMass      (thisNode,accretionModeTotal)
       failedHotHaloMass=accretionHalo_%failedAccretedMass(thisNode,accretionModeTotal)
       ! If either is non-zero, then create a hot halo component and add these masses to it.
       if (hotHaloMass > 0.0d0 .or. failedHotHaloMass > 0.0d0) then
          call Node_Component_Hot_Halo_Very_Simple_Create(thisNode)
          thisHotHaloComponent => thisNode%hotHalo()
          call thisHotHaloComponent%           massSet(      hotHaloMass )
          call thisHotHaloComponent% unaccretedMassSet(failedHotHaloMass )
       end if
    end select
    return
  end subroutine Node_Component_Hot_Halo_Very_Simple_Tree_Initialize
  
  !# <satelliteMergerTask>
  !#  <unitName>Node_Component_Hot_Halo_Very_Simple_Satellite_Merging</unitName>
  !# </satelliteMergerTask>
  subroutine Node_Component_Hot_Halo_Very_Simple_Satellite_Merging(thisNode)
    !% Remove any hot halo associated with {\normalfont \ttfamily thisNode} before it merges with its host halo.
    implicit none
    type (treeNode            ), intent(inout), pointer :: thisNode
    type (treeNode            )               , pointer :: hostNode
    class(nodeComponentHotHalo)               , pointer :: hostHotHaloComponent, thisHotHaloComponent

    ! Get the hot halo component.
    thisHotHaloComponent => thisNode%hotHalo()
    ! Ensure that it is of unspecified class.
    select type (thisHotHaloComponent)
    class is (nodeComponentHotHaloVerySimple)

       ! Find the node to merge with.
       hostNode             => thisNode%mergesWith()
       hostHotHaloComponent => hostNode%hotHalo   ()

       ! Move the hot halo to the host.
       call hostHotHaloComponent%                    massSet(                                                 &
            &                                                 hostHotHaloComponent%mass                    () &
            &                                                +thisHotHaloComponent%mass                    () &
            &                                               )
       call thisHotHaloComponent%                    massSet(                                                 &
            &                                                 0.0d0                                           &
            &                                               )
    end select
    return
  end subroutine Node_Component_Hot_Halo_Very_Simple_Satellite_Merging

  !# <nodePromotionTask>
  !#  <unitName>Node_Component_Hot_Halo_Very_Simple_Promote</unitName>
  !# </nodePromotionTask>
  subroutine Node_Component_Hot_Halo_Very_Simple_Promote(thisNode)
    !% Ensure that {\normalfont \ttfamily thisNode} is ready for promotion to its parent. In this case, we simply update the hot halo mass of {\tt
    !% thisNode} to account for any hot halo already in the parent.
    implicit none
    type (treeNode            ), intent(inout), pointer :: thisNode
    type (treeNode            )               , pointer :: parentNode
    class(nodeComponentHotHalo)               , pointer :: parentHotHaloComponent, thisHotHaloComponent

    ! Get the hot halo component.
    thisHotHaloComponent => thisNode%hotHalo()
    ! Ensure that it is of specified class.
    select type (thisHotHaloComponent)
    class is (nodeComponentHotHaloVerySimple)
       ! Get the parent node of this node.
       parentNode             => thisNode  %parent
       parentHotHaloComponent => parentNode%hotHalo()
       ! If the parent node has a hot halo component, then add it to that of this node, and perform other changes needed prior to
       ! promotion.
       select type (parentHotHaloComponent)
       class is (nodeComponentHotHaloVerySimple)
          call thisHotHaloComponent%unaccretedMassSet(                                         &
               &                                      +thisHotHaloComponent  %unaccretedMass() &
               &                                      +parentHotHaloComponent%unaccretedMass() &
               &                                     )
          call thisHotHaloComponent%          massSet(                                         &
               &                                      +thisHotHaloComponent  %          mass() &
               &                                      +parentHotHaloComponent%          mass() &
               &                                     )
       end select
    end select
    return
  end subroutine Node_Component_Hot_Halo_Very_Simple_Promote

  !# <postEvolveTask>
  !#  <unitName>Node_Component_Hot_Halo_Very_Simple_Post_Evolve</unitName>
  !# </postEvolveTask>
  subroutine Node_Component_Hot_Halo_Very_Simple_Post_Evolve(thisNode)
    !% Do processing of the node required after evolution.
    implicit none
    type (treeNode            ), intent(inout), pointer :: thisNode
    type (treeNode            )               , pointer :: parentNode
    class(nodeComponentHotHalo)               , pointer :: parentHotHaloComponent, thisHotHaloComponent

    ! Get the hot halo component.
    thisHotHaloComponent => thisNode%hotHalo()
    select type (thisHotHaloComponent)
    class is (nodeComponentHotHaloVerySimple)
       ! Check if this node is a satellite.
       if (thisNode%isSatellite()) then
          ! Transfer any outflowed gas to the hot halo of the parent node.
          parentNode => thisNode%parent
          do while (parentNode%isSatellite())
             parentNode => parentNode%parent
          end do
          parentHotHaloComponent => parentNode%hotHalo()
          call parentHotHaloComponent%massSet(parentHotHaloComponent%mass()+thisHotHaloComponent%mass())
          call   thisHotHaloComponent%massSet(                                                    0.0d0)
       end if
    end select
    return
  end subroutine Node_Component_Hot_Halo_Very_Simple_Post_Evolve

  !# <nodeMergerTask>
  !#  <unitName>Node_Component_Hot_Halo_Very_Simple_Node_Merger</unitName>
  !# </nodeMergerTask>
  subroutine Node_Component_Hot_Halo_Very_Simple_Node_Merger(thisNode)
    !% Starve {\normalfont \ttfamily thisNode} by transferring its hot halo to its parent.
    implicit none
    type (treeNode            ), intent(inout), pointer :: thisNode
    type (treeNode            )               , pointer :: parentNode
    class(nodeComponentHotHalo)               , pointer :: parentHotHaloComponent, thisHotHaloComponent

    ! Get the hot halo component.
    thisHotHaloComponent => thisNode%hotHalo()
    select type (thisHotHaloComponent)
    class is (nodeComponentHotHaloVerySimple)
       ! Find the parent node and its hot halo component.
       parentNode => thisNode%parent
       parentHotHaloComponent => parentNode%hotHalo()
       ! Any gas that failed to be accreted by this halo is always transferred to the parent.
       call parentHotHaloComponent%unaccretedMassSet(parentHotHaloComponent%unaccretedMass()+thisHotHaloComponent%unaccretedMass())
       call   thisHotHaloComponent%unaccretedMassSet(                                                                        0.0d0)
       ! Move the hot halo to the parent. We leave the hot halo in place even if it is starved, since outflows will accumulate
       ! to this hot halo (and will be moved to the parent at the end of the evolution timestep).
       call parentHotHaloComponent%          massSet(parentHotHaloComponent%          mass()+thisHotHaloComponent%          mass())
       call   thisHotHaloComponent%          massSet(                                                                        0.0d0)
    end select
    return
  end subroutine Node_Component_Hot_Halo_Very_Simple_Node_Merger

  subroutine Node_Component_Hot_Halo_Very_Simple_Cooling_Rate(thisNode)
    !% Get and store the cooling rate for {\normalfont \ttfamily thisNode}.
    use Cooling_Rates
    implicit none
    type (treeNode            ), intent(inout), pointer :: thisNode
    class(nodeComponentHotHalo)               , pointer :: thisHotHaloComponent

    if (.not.gotCoolingRate) then
       ! Get the hot halo component.
       thisHotHaloComponent => thisNode%hotHalo()
       if (thisHotHaloComponent%mass() > 0.0d0) then
          ! Get the cooling time.
          coolingRate=Cooling_Rate(thisNode)
       else
          coolingRate=0.0d0
       end if
       ! Flag that cooling rate has now been computed.
       gotCoolingRate=.true.
    end if
    return
  end subroutine Node_Component_Hot_Halo_Very_Simple_Cooling_Rate

  subroutine Node_Component_Hot_Halo_Very_Simple_Create(thisNode)
    !% Creates a very simple hot halo component for {\normalfont \ttfamily thisNode}.
    implicit none
    type (treeNode            ), intent(inout), pointer :: thisNode
    class(nodeComponentHotHalo)               , pointer :: thisHotHaloComponent

    ! Ensure that this module has been initialized.
    call Node_Component_Hot_Halo_Very_Simple_Initialize()

    ! Create the component.
    thisHotHaloComponent => thisNode%hotHalo(autoCreate=.true.)
    return
  end subroutine Node_Component_Hot_Halo_Very_Simple_Create

end module Node_Component_Hot_Halo_Very_Simple
