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

!% Contains a module which implements an extension to the very simple hot halo node component by including an outflowed reservoir
!% with delayed reincorporation.

module Node_Component_Hot_Halo_VS_Delayed
  !% Implements an extension to the very simple hot halo node component by including an outflowed reservoir
  !% with delayed reincorporation.
  use :: Hot_Halo_Outflows_Reincorporations, only : hotHaloOutflowReincorporationClass
  implicit none
  private
  public :: Node_Component_Hot_Halo_VS_Delayed_Node_Merger        , Node_Component_Hot_Halo_VS_Delayed_Rate_Compute     , &
       &    Node_Component_Hot_Halo_VS_Delayed_Scale_Set          , Node_Component_Hot_Halo_VS_Delayed_Tree_Initialize  , &
       &    Node_Component_Hot_Halo_VS_Delayed_Thread_Uninitialize, Node_Component_Hot_Halo_VS_Delayed_Thread_Initialize, &
       &    Node_Component_Hot_Halo_VS_Delayed_Initialize

  !# <component>
  !#  <class>hotHalo</class>
  !#  <name>verySimpleDelayed</name>
  !#  <extends>
  !#   <class>hotHalo</class>
  !#   <name>verySimple</name>
  !#  </extends>
  !#  <isDefault>false</isDefault>
  !#  <properties>
  !#   <property>
  !#     <name>outflowedMass</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#     <output unitsInSI="massSolar" comment="Mass of gas in the hot halo."/>
  !#   </property>
  !#   <property>
  !#     <name>outflowedAbundances</name>
  !#     <type>abundances</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#     <output unitsInSI="massSolar" comment="Mass of metals in the outflowed phase of the hot halo."/>
  !#   </property>
  !#  </properties>
  !# </component>

  ! Options controlling the numerical implementation.
  double precision                                               :: hotHaloVerySimpleDelayedMassScaleRelative

  ! Objects used by this component.
  class           (hotHaloOutflowReincorporationClass ), pointer :: hotHaloOutflowReincorporation_
  !$omp threadprivate(hotHaloOutflowReincorporation_)
  
contains

  !# <nodeComponentInitializationTask>
  !#  <unitName>Node_Component_Hot_Halo_VS_Delayed_Initialize</unitName>
  !# </nodeComponentInitializationTask>
  subroutine Node_Component_Hot_Halo_VS_Delayed_Initialize(parameters_)
    !% Initializes the very simple hot halo component module.
    use :: Galacticus_Nodes, only : defaultHotHaloComponent, nodeComponentHotHaloVerySimpleDelayed
    use :: Input_Parameters, only : inputParameter         , inputParameters
    implicit none
    type(inputParameters                      ), intent(inout) :: parameters_
    type(nodeComponentHotHaloVerySimpleDelayed)                :: hotHalo

    !$omp critical (Node_Component_Hot_Halo_Very_Simple_Delayed_Initialize)
    if (defaultHotHaloComponent%verySimpleDelayedIsActive()) then
       ! Bind outflowing material pipes to the functions that will handle input of outflowing material to the hot halo.
       call hotHalo%      outflowingMassRateFunction(Node_Component_Hot_Halo_VS_Delayed_Outflowing_Mass_Rate      )
       call hotHalo%outflowingAbundancesRateFunction(Node_Component_Hot_Halo_VS_Delayed_Outflowing_Abundances_Rate)
       ! Read parameters controlling the physical implementation.
       !# <inputParameter>
       !#   <name>hotHaloVerySimpleDelayedMassScaleRelative</name>
       !#   <defaultValue>1.0d-2</defaultValue>
       !#   <description>The mass scale, relative to the total mass of the node, below which calculations in the delayed very simple hot halo component are allowed to become inaccurate.</description>
       !#   <source>parameters_</source>
       !# </inputParameter>
    end if
    !$omp end critical (Node_Component_Hot_Halo_Very_Simple_Delayed_Initialize)
    return
  end subroutine Node_Component_Hot_Halo_VS_Delayed_Initialize

  !# <nodeComponentThreadInitializationTask>
  !#  <unitName>Node_Component_Hot_Halo_VS_Delayed_Thread_Initialize</unitName>
  !# </nodeComponentThreadInitializationTask>
  subroutine Node_Component_Hot_Halo_VS_Delayed_Thread_Initialize(parameters_)
    !% Initializes the tree node very simple disk profile module.
    use :: Events_Hooks    , only : nodePromotionEvent     , satelliteMergerEvent    , postEvolveEvent, openMPThreadBindingAtLevel, &
         &                          dependencyRegEx        , dependencyDirectionAfter
    use :: Galacticus_Nodes, only : defaultHotHaloComponent
    use :: Input_Parameters, only : inputParameter         , inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters_
    type(dependencyRegEx), dimension(1)  :: dependencies
    !$GLC attributes unused :: parameters_

    if (defaultHotHaloComponent%verySimpleDelayedIsActive()) then
       !# <objectBuilder class="hotHaloOutflowReincorporation" name="hotHaloOutflowReincorporation_" source="parameters_"/>
       dependencies(1)=dependencyRegEx(dependencyDirectionAfter,'^remnantStructure:')
       call nodePromotionEvent  %attach(defaultHotHaloComponent,nodePromotion  ,openMPThreadBindingAtLevel,label='nodeComponentHotHaloVerySimpleDelayed'                          )
       call satelliteMergerEvent%attach(defaultHotHaloComponent,satelliteMerger,openMPThreadBindingAtLevel,label='nodeComponentHotHaloVerySimpleDelayed',dependencies=dependencies)
       call postEvolveEvent     %attach(defaultHotHaloComponent,postEvolve     ,openMPThreadBindingAtLevel,label='nodeComponentHotHaloVerySimpleDelayed'                          )
    end if
    return
  end subroutine Node_Component_Hot_Halo_VS_Delayed_Thread_Initialize

  !# <nodeComponentThreadUninitializationTask>
  !#  <unitName>Node_Component_Hot_Halo_VS_Delayed_Thread_Uninitialize</unitName>
  !# </nodeComponentThreadUninitializationTask>
  subroutine Node_Component_Hot_Halo_VS_Delayed_Thread_Uninitialize()
    !% Uninitializes the tree node very simple disk profile module.
    use :: Events_Hooks    , only : nodePromotionEvent     , satelliteMergerEvent, postEvolveEvent
    use :: Galacticus_Nodes, only : defaultHotHaloComponent
    implicit none

    if (defaultHotHaloComponent%verySimpleDelayedIsActive()) then
       !# <objectDestructor name="hotHaloOutflowReincorporation_"/>
       call nodePromotionEvent  %detach(defaultHotHaloComponent,nodePromotion  )
       call satelliteMergerEvent%detach(defaultHotHaloComponent,satelliteMerger)
       call postEvolveEvent     %detach(defaultHotHaloComponent,postEvolve     )
    end if
    return
  end subroutine Node_Component_Hot_Halo_VS_Delayed_Thread_Uninitialize

  subroutine Node_Component_Hot_Halo_VS_Delayed_Outflowing_Mass_Rate(self,rate,interrupt,interruptProcedure)
    !% Accept outflowing gas from a galaxy and deposit it into very simple hot halo.
    use :: Galacticus_Nodes, only : nodeComponentHotHalo
    implicit none
    class           (nodeComponentHotHalo), intent(inout)                    :: self
    double precision                      , intent(in   )                    :: rate
    logical                               , intent(inout), optional          :: interrupt
    procedure       (                    ), intent(inout), optional, pointer :: interruptProcedure
    !$GLC attributes unused :: interrupt, interruptProcedure

    ! Funnel the outflowing gas into the outflowed reservoir.
    call self%outflowedMassRate(rate)
    return
  end subroutine Node_Component_Hot_Halo_VS_Delayed_Outflowing_Mass_Rate

  subroutine Node_Component_Hot_Halo_VS_Delayed_Outflowing_Abundances_Rate(self,rate,interrupt,interruptProcedure)
    !% Accept outflowing gas abundances from a galaxy and deposit them into very simple hot halo.
    use :: Abundances_Structure, only : abundances
    use :: Galacticus_Nodes    , only : nodeComponentHotHalo
    implicit none
    class    (nodeComponentHotHalo), intent(inout)                    :: self
    type     (abundances          ), intent(in   )                    :: rate
    logical                        , intent(inout), optional          :: interrupt
    procedure(                    ), intent(inout), optional, pointer :: interruptProcedure
    !$GLC attributes unused :: interrupt, interruptProcedure

    ! Funnel the outflowing gas abundances into the outflowed reservoir.
    call self%outflowedAbundancesRate(rate)
    return
  end subroutine Node_Component_Hot_Halo_VS_Delayed_Outflowing_Abundances_Rate

  !# <rateComputeTask>
  !#  <unitName>Node_Component_Hot_Halo_VS_Delayed_Rate_Compute</unitName>
  !# </rateComputeTask>
  subroutine Node_Component_Hot_Halo_VS_Delayed_Rate_Compute(node,interrupt,interruptProcedure,propertyType)
    !% Compute the very simple hot halo component mass rate of change.
    use :: Abundances_Structure, only : abundances                   , operator(*)                          , zeroAbundances
    use :: Galacticus_Nodes    , only : nodeComponentHotHalo         , nodeComponentHotHaloVerySimpleDelayed, propertyTypeInactive, treeNode, &
         &                              defaultHotHaloComponent
    implicit none
    type            (treeNode             ), intent(inout)          :: node
    logical                                , intent(inout)          :: interrupt
    procedure       (                     ), intent(inout), pointer :: interruptProcedure
    integer                                , intent(in   )          :: propertyType
    class           (nodeComponentHotHalo )               , pointer :: hotHalo
    type            (abundances           ), save                   :: abundancesReturnRate
    !$omp threadprivate(abundancesReturnRate)
    double precision                                                :: outflowReturnRate
    !$GLC attributes unused :: interrupt, interruptProcedure

    ! Return immediately if inactive variables are requested.
    if (propertyType == propertyTypeInactive) return
    ! Return immediately if this class is not in use.
    if (.not.defaultHotHaloComponent%verySimpleDelayedIsActive()) return
    ! Don't reincorporate gas for satellites - we don't want it to be able to re-infall back onto the satellite.
    if (node%isSatellite()) return
    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    select type (hotHalo)
    class is (nodeComponentHotHaloVerySimpleDelayed)
       ! Move outflowed material back to the hot reservoir.
       outflowReturnRate=hotHaloOutflowReincorporation_%rate(node)
       if (hotHalo%outflowedMass() > 0.0d0) then
          abundancesReturnRate           =  +outflowReturnRate             &
               &                            *hotHalo%outflowedAbundances() &
               &                            /hotHalo%outflowedMass      ()
       else
          abundancesReturnRate           =  zeroAbundances
       end if
       call hotHalo%outflowedMassRate      (-   outflowReturnRate)
       call hotHalo%outflowedAbundancesRate(-abundancesReturnRate)
       call hotHalo%massRate               (+   outflowReturnRate)
       call hotHalo%abundancesRate         (+abundancesReturnRate)
    end select
    return
  end subroutine Node_Component_Hot_Halo_VS_Delayed_Rate_Compute

  !# <scaleSetTask>
  !#  <unitName>Node_Component_Hot_Halo_VS_Delayed_Scale_Set</unitName>
  !# </scaleSetTask>
  subroutine Node_Component_Hot_Halo_VS_Delayed_Scale_Set(node)
    !% Set scales for properties of {\normalfont \ttfamily node}.
    use :: Abundances_Structure, only : unitAbundances
    use :: Galacticus_Nodes    , only : nodeComponentBasic     , nodeComponentHotHalo, nodeComponentHotHaloVerySimpleDelayed, treeNode, &
         &                              defaultHotHaloComponent
    implicit none
    type            (treeNode            ), intent(inout), pointer :: node
    class           (nodeComponentHotHalo)               , pointer :: hotHalo
    class           (nodeComponentBasic  )               , pointer :: basic
    double precision                                               :: massVirial

    ! Check if we are the default method.
    if (.not.defaultHotHaloComponent%verySimpleDelayedIsActive()) return
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
       call hotHalo%outflowedMassScale      (               massVirial*hotHaloVerySimpleDelayedMassScaleRelative)
       call hotHalo%outflowedAbundancesScale(unitAbundances*massVirial*hotHaloVerySimpleDelayedMassScaleRelative)
    end select
    return
  end subroutine Node_Component_Hot_Halo_VS_Delayed_Scale_Set

  !# <mergerTreeInitializeTask>
  !#  <unitName>Node_Component_Hot_Halo_VS_Delayed_Tree_Initialize</unitName>
  !#  <after>darkMatterProfile</after>
  !#  <after>Node_Component_Hot_Halo_Very_Simple_Tree_Initialize</after>
  !# </mergerTreeInitializeTask>
  subroutine Node_Component_Hot_Halo_VS_Delayed_Tree_Initialize(node)
    !% Initialize the contents of the very simple hot halo component.
    use :: Abundances_Structure, only : zeroAbundances
    use :: Galacticus_Nodes    , only : defaultHotHaloComponent, nodeComponentHotHalo, nodeComponentHotHaloVerySimpleDelayed, treeNode
    implicit none
    type (treeNode                ), intent(inout), pointer :: node
    class(nodeComponentHotHalo    )               , pointer :: hotHalo

    ! If the very simple hot halo is not active, then return immediately.
    if (.not.defaultHotHaloComponent%verySimpleDelayedIsActive()) return

    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    ! Ensure that it is of our class.
    select type (hotHalo)
    type is (nodeComponentHotHaloVerySimpleDelayed)
       call hotHalo%outflowedMassSet      (         0.0d0)
       call hotHalo%outflowedAbundancesSet(zeroAbundances)
    end select
    return
  end subroutine Node_Component_Hot_Halo_VS_Delayed_Tree_Initialize

  subroutine satelliteMerger(self,node)
    !% Remove any hot halo associated with {\normalfont \ttfamily node} before it merges with its host halo.
    use :: Abundances_Structure, only : zeroAbundances
    use :: Galacticus_Nodes    , only : nodeComponentHotHalo, nodeComponentHotHaloVerySimpleDelayed, treeNode
    implicit none
    class(*                   ), intent(inout) :: self
    type (treeNode            ), intent(inout) :: node
    type (treeNode            ), pointer       :: nodeHost
    class(nodeComponentHotHalo), pointer       :: hotHaloHost, hotHalo
    !$GLC attributes unused :: self

    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    ! Ensure that it is of unspecified class.
    select type (hotHalo)
    class is (nodeComponentHotHaloVerySimpleDelayed)
       ! Find the node to merge with.
       nodeHost    => node    %mergesWith(                 )
       hotHaloHost => nodeHost%hotHalo   (autoCreate=.true.)
       ! Move the hot halo to the host.
       call hotHaloHost%outflowedMassSet      (                                   &
            &                                  +hotHaloHost%outflowedMass      () &
            &                                  +    hotHalo%outflowedMass      () &
            &                                 )
       call hotHaloHost%outflowedAbundancesSet(                                   &
            &                                  +hotHaloHost%outflowedAbundances() &
            &                                  +    hotHalo%outflowedAbundances() &
            &                                 )
       call     hotHalo%outflowedMassSet      (                                   &
            &                                  +0.0d0                             &
            &                                 )
       call     hotHalo%outflowedAbundancesSet(                                   &
            &                                  +zeroAbundances                    &
            &                                 )
    end select
    return
  end subroutine satelliteMerger

  subroutine nodePromotion(self,node)
    !% Ensure that {\normalfont \ttfamily node} is ready for promotion to its parent. In this case, we simply update the hot halo mass of {\normalfont \ttfamily
    !% node} to account for any hot halo already in the parent.
    use :: Galacticus_Nodes, only : nodeComponentHotHalo, nodeComponentHotHaloVerySimpleDelayed, treeNode
    implicit none
    class(*                   ), intent(inout)          :: self
    type (treeNode            ), intent(inout), target  :: node
    type (treeNode            )               , pointer :: nodeParent
    class(nodeComponentHotHalo)               , pointer :: hotHaloParent, hotHalo
    !$GLC attributes unused :: self
    
    hotHalo       => node      %hotHalo(                 )
    nodeParent    => node      %parent
    hotHaloParent => nodeParent%hotHalo(autoCreate=.true.)
    ! If the parent node has a hot halo component, then add it to that of this node, and perform other changes needed prior to
    ! promotion.
    select type (hotHaloParent)
    class is (nodeComponentHotHaloVerySimpleDelayed)
       call hotHalo%outflowedMassSet      (                                     &
            &                              +      hotHalo%outflowedMass      () &
            &                              +hotHaloParent%outflowedMass      () &
            &                             )
       call hotHalo%outflowedAbundancesSet(                                     &
            &                              +      hotHalo%outflowedAbundances() &
            &                              +hotHaloParent%outflowedAbundances() &
            &                             )
    end select
    return
  end subroutine nodePromotion

  subroutine postEvolve(self,node)
    !% Do processing of the node required after evolution.
    use :: Abundances_Structure, only : zeroAbundances
    use :: Galacticus_Nodes    , only : nodeComponentHotHalo, nodeComponentHotHaloVerySimpleDelayed, treeNode
    implicit none
    class(*                   ), intent(inout) :: self
    type (treeNode            ), intent(inout) :: node
    type (treeNode            ), pointer       :: nodeParent
    class(nodeComponentHotHalo), pointer       :: hotHaloParent, hotHalo
    !$GLC attributes unused :: self

    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    select type (hotHalo)
    class is (nodeComponentHotHaloVerySimpleDelayed)
       ! Check if this node is a satellite.
       if (node%isSatellite()) then
          ! Transfer any outflowed gas to the hot halo of the parent node.
          nodeParent => node%parent
          do while (nodeParent%isSatellite())
             nodeParent => nodeParent%parent
          end do
          hotHaloParent => nodeParent%hotHalo(autoCreate=.true.)
          call hotHaloParent%outflowedMassSet      (hotHaloParent%outflowedMass      ()+hotHalo%outflowedMass      ())
          call hotHaloParent%outflowedAbundancesSet(hotHaloParent%outflowedAbundances()+hotHalo%outflowedAbundances())
          call       hotHalo%outflowedMassSet      (                                                            0.0d0)
          call       hotHalo%outflowedAbundancesSet(                                                   zeroAbundances)
       end if
    end select
    return
  end subroutine postEvolve

  !# <nodeMergerTask>
  !#  <unitName>Node_Component_Hot_Halo_VS_Delayed_Node_Merger</unitName>
  !# </nodeMergerTask>
  subroutine Node_Component_Hot_Halo_VS_Delayed_Node_Merger(node)
    !% Starve {\normalfont \ttfamily node} by transferring its hot halo to its parent.
    use :: Abundances_Structure, only : zeroAbundances
    use :: Galacticus_Nodes    , only : nodeComponentHotHalo, nodeComponentHotHaloVerySimpleDelayed, treeNode, defaultHotHaloComponent
    implicit none
    type (treeNode            ), intent(inout) :: node
    type (treeNode            ), pointer       :: nodeParent
    class(nodeComponentHotHalo), pointer       :: hotHaloParent, hotHalo

    ! Return immediately if this class is not in use.
    if (.not.defaultHotHaloComponent%verySimpleDelayedIsActive()) return
    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    select type (hotHalo)
    class is (nodeComponentHotHaloVerySimpleDelayed)
       ! Find the parent node and its hot halo component.
       nodeParent    => node      %parent
       hotHaloParent => nodeParent%hotHalo(autoCreate=.true.)
       ! Move the hot halo to the parent. We leave the hot halo in place even if it is starved, since outflows will accumulate
       ! to this hot halo (and will be moved to the parent at the end of the evolution timestep).
       call hotHaloParent%outflowedMassSet      (hotHaloParent%outflowedMass      ()+hotHalo%outflowedMass      ())
       call hotHaloParent%outflowedAbundancesSet(hotHaloParent%outflowedAbundances()+hotHalo%outflowedAbundances())
       call       hotHalo%outflowedMassSet      (                                                            0.0d0)
       call       hotHalo%outflowedAbundancesSet(                                                   zeroAbundances)
    end select
    return
  end subroutine Node_Component_Hot_Halo_VS_Delayed_Node_Merger

end module Node_Component_Hot_Halo_VS_Delayed
