!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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
Contains a module which implements a very simple hot halo node component.
!!}

module Node_Component_Hot_Halo_Very_Simple
  !!{
  Implements a very simple hot halo node component.
  !!}
  use :: Accretion_Halos        , only : accretionHaloClass
  use :: Cooling_Rates          , only : coolingRateClass
  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass
  implicit none
  private
  public :: Node_Component_Hot_Halo_Very_Simple_Reset            , Node_Component_Hot_Halo_Very_Simple_Rate_Compute       , &
       &    Node_Component_Hot_Halo_Very_Simple_Scale_Set        , Node_Component_Hot_Halo_Very_Simple_Tree_Initialize    , &
       &    Node_Component_Hot_Halo_Very_Simple_Thread_Initialize, Node_Component_Hot_Halo_Very_Simple_Thread_Uninitialize, &
       &    Node_Component_Hot_Halo_Very_Simple_State_Store      , Node_Component_Hot_Halo_Very_Simple_State_Restore      , &
       &    Node_Component_Hot_Halo_Very_Simple_Node_Merger

  !![
  <component>
   <class>hotHalo</class>
   <name>verySimple</name>
   <isDefault>false</isDefault>
   <properties>
    <property>
      <name>mass</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" />
      <output unitsInSI="massSolar" comment="Mass of gas in the hot halo."/>
    </property>
    <property>
      <name>abundances</name>
      <type>abundances</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" />
      <output unitsInSI="massSolar" comment="Mass of metals in the hot halo."/>
    </property>
    <property>
      <name>unaccretedMass</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" />
      <output unitsInSI="massSolar" comment="Mass of gas that failed to accrete into the hot halo."/>
    </property>
    <property>
      <name>outflowingMass</name>
      <attributes isSettable="false" isGettable="false" isEvolvable="true" isDeferred="rate" isVirtual="true" />
      <type>double</type>
      <rank>0</rank>
    </property>
    <property>
      <name>outflowingAbundances</name>
      <attributes isSettable="false" isGettable="false" isEvolvable="true" isDeferred="rate" isVirtual="true" />
      <type>abundances</type>
      <rank>0</rank>
    </property>
    <property>
      <name>hotHaloCoolingMass</name>
      <attributes isSettable="false" isGettable="false" isEvolvable="true" isDeferred="rate" bindsTo="top" isVirtual="true" />
      <type>double</type>
      <rank>0</rank>
    </property>
    <property>
      <name>hotHaloCoolingAbundances</name>
      <attributes isSettable="false" isGettable="false" isEvolvable="true" isDeferred="rate" bindsTo="top" isVirtual="true" />
      <type>abundances</type>
      <rank>0</rank>
    </property>
    <property>
      <name>outerRadius</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="false" isGettable="true" isEvolvable="false" isDeferred="get" isVirtual="true" />
    </property>
   </properties>
  </component>
  !!]

  ! Objects used by this component.
  class(darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_
  class(coolingRateClass        ), pointer :: coolingRate_
  class(accretionHaloClass      ), pointer :: accretionHalo_
  !$omp threadprivate(darkMatterHaloScale_,coolingRate_,accretionHalo_)

  ! Quantities stored to avoid repeated computation.
  logical          :: gotCoolingRate   =.false.
  double precision :: rateCooling
  !$omp threadprivate(gotCoolingRate,rateCooling)

contains

  subroutine Node_Component_Hot_Halo_Very_Simple_Initialize()
    !!{
    Initializes the very simple hot halo component module.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentHotHaloVerySimple
    implicit none
    type(nodeComponentHotHaloVerySimple) :: hotHalo

    ! Bind outflowing material pipes to the functions that will handle input of outflowing material to the hot halo.
    call hotHalo%      outflowingMassRateFunction(Node_Component_Hot_Halo_Very_Simple_Outflowing_Mass_Rate      )
    call hotHalo%outflowingAbundancesRateFunction(Node_Component_Hot_Halo_Very_Simple_Outflowing_Abundances_Rate)
    ! Bind outer radius function.
    call hotHalo%             outerRadiusFunction(Node_Component_Hot_Halo_Very_Simple_Outer_Radius              )
    return
  end subroutine Node_Component_Hot_Halo_Very_Simple_Initialize

  !![
  <nodeComponentThreadInitializationTask>
   <unitName>Node_Component_Hot_Halo_Very_Simple_Thread_Initialize</unitName>
  </nodeComponentThreadInitializationTask>
  !!]
  subroutine Node_Component_Hot_Halo_Very_Simple_Thread_Initialize(parameters_)
    !!{
    Initializes the tree node very simple disk profile module.
    !!}
    use :: Events_Hooks    , only : nodePromotionEvent     , satelliteMergerEvent    , postEvolveEvent, openMPThreadBindingAtLevel, &
         &                          dependencyRegEx        , dependencyDirectionAfter
    use :: Galacticus_Nodes, only : defaultHotHaloComponent
    use :: Input_Parameters, only : inputParameter         , inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters_
    type(dependencyRegEx), dimension(1)  :: dependencies

    if (defaultHotHaloComponent%verySimpleIsActive()) then
       !![
       <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters_"/>
       <objectBuilder class="coolingRate"         name="coolingRate_"         source="parameters_"/>
       <objectBuilder class="accretionHalo"       name="accretionHalo_"       source="parameters_"/>
       !!]
       dependencies(1)=dependencyRegEx(dependencyDirectionAfter,'^remnantStructure:')
       call nodePromotionEvent  %attach(defaultHotHaloComponent,nodePromotion  ,openMPThreadBindingAtLevel,label='nodeComponentHotHaloVerySimple'                          )
       call satelliteMergerEvent%attach(defaultHotHaloComponent,satelliteMerger,openMPThreadBindingAtLevel,label='nodeComponentHotHaloVerySimple',dependencies=dependencies)
       call postEvolveEvent     %attach(defaultHotHaloComponent,postEvolve     ,openMPThreadBindingAtLevel,label='nodeComponentHotHaloVerySimple'                          )
    end if
    return
  end subroutine Node_Component_Hot_Halo_Very_Simple_Thread_Initialize

  !![
  <nodeComponentThreadUninitializationTask>
   <unitName>Node_Component_Hot_Halo_Very_Simple_Thread_Uninitialize</unitName>
  </nodeComponentThreadUninitializationTask>
  !!]
  subroutine Node_Component_Hot_Halo_Very_Simple_Thread_Uninitialize()
    !!{
    Uninitializes the tree node very simple disk profile module.
    !!}
    use :: Events_Hooks    , only : nodePromotionEvent     , satelliteMergerEvent, postEvolveEvent
    use :: Galacticus_Nodes, only : defaultHotHaloComponent
    implicit none

    if (defaultHotHaloComponent%verySimpleIsActive()) then
       !![
       <objectDestructor name="darkMatterHaloScale_"/>
       <objectDestructor name="coolingRate_"        />
       <objectDestructor name="accretionHalo_"      />
       !!]
       if (nodePromotionEvent  %isAttached(defaultHotHaloComponent,nodePromotion  )) call nodePromotionEvent  %detach(defaultHotHaloComponent,nodePromotion  )
       if (satelliteMergerEvent%isAttached(defaultHotHaloComponent,satelliteMerger)) call satelliteMergerEvent%detach(defaultHotHaloComponent,satelliteMerger)
       if (postEvolveEvent     %isAttached(defaultHotHaloComponent,postEvolve     )) call postEvolveEvent     %detach(defaultHotHaloComponent,postEvolve     )
    end if
    return
  end subroutine Node_Component_Hot_Halo_Very_Simple_Thread_Uninitialize

  !![
  <calculationResetTask>
  <unitName>Node_Component_Hot_Halo_Very_Simple_Reset</unitName>
  </calculationResetTask>
  !!]
  subroutine Node_Component_Hot_Halo_Very_Simple_Reset(node)
    !!{
    Remove memory of stored computed values as we're about to begin computing derivatives anew.
    !!}
    use :: Galacticus_Nodes, only : treeNode
    implicit none
    type(treeNode), intent(inout) :: node
    !$GLC attributes unused :: node

    gotCoolingRate=.false.
    return
  end subroutine Node_Component_Hot_Halo_Very_Simple_Reset

  subroutine Node_Component_Hot_Halo_Very_Simple_Push_To_Cooling_Pipes(node,massRate,interrupt,interruptProcedure)
    !!{
    Push mass through the cooling pipes at the given rate.
    !!}
    use :: Abundances_Structure, only : abundances          , operator(*)
    use :: Galacticus_Nodes    , only : nodeComponentHotHalo, nodeComponentHotHaloVerySimple, treeNode
    implicit none
    type            (treeNode            ), intent(inout)          :: node
    double precision                      , intent(in   )          :: massRate
    logical                               , intent(inout)          :: interrupt
    procedure       (                    ), intent(inout), pointer :: interruptProcedure
    class           (nodeComponentHotHalo)               , pointer :: hotHalo
    type            (abundances          )                         :: abundancesCoolingRate

    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    select type (hotHalo)
    class is (nodeComponentHotHaloVerySimple)
       ! Ignore zero rates.
       if (massRate /= 0.0d0 .and. hotHalo%mass() > 0.0d0) then
          ! Remove mass from the hot component.
          abundancesCoolingRate=+massRate             &
               &                *hotHalo%abundances() &
               &                /hotHalo%mass      ()
          call    hotHalo%massRate                    (-massRate                                          )
          call    hotHalo%abundancesRate              (-abundancesCoolingRate                             )
          ! Pipe the mass rate to whatever component claimed it.
          if (hotHalo%hotHaloCoolingMassRateIsAttached      ()) &
               & call hotHalo%hotHaloCoolingMassRate      (+massRate             ,interrupt,interruptProcedure)
          if (hotHalo%hotHaloCoolingAbundancesRateIsAttached()) &
               & call hotHalo%hotHaloCoolingAbundancesRate(+abundancesCoolingRate,interrupt,interruptProcedure)
       end if
    end select
    return
  end subroutine Node_Component_Hot_Halo_Very_Simple_Push_To_Cooling_Pipes

  subroutine Node_Component_Hot_Halo_Very_Simple_Outflowing_Mass_Rate(self,rate,interrupt,interruptProcedure)
    !!{
    Accept outflowing gas from a galaxy and deposit it into very simple hot halo.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentHotHalo
    implicit none
    class           (nodeComponentHotHalo), intent(inout)                    :: self
    double precision                      , intent(in   )                    :: rate
    logical                               , intent(inout), optional          :: interrupt
    procedure       (                    ), intent(inout), optional, pointer :: interruptProcedure
    !$GLC attributes unused :: interrupt, interruptProcedure

    ! Funnel the outflow gas into the hot halo.
    call self%massRate(rate)
    return
  end subroutine Node_Component_Hot_Halo_Very_Simple_Outflowing_Mass_Rate

  subroutine Node_Component_Hot_Halo_Very_Simple_Outflowing_Abundances_Rate(self,rate,interrupt,interruptProcedure)
    !!{
    Accept outflowing gas abundances from a galaxy and deposit them into very simple hot halo.
    !!}
    use :: Abundances_Structure, only : abundances
    use :: Galacticus_Nodes    , only : nodeComponentHotHalo
    implicit none
    class    (nodeComponentHotHalo), intent(inout)                    :: self
    type     (abundances          ), intent(in   )                    :: rate
    logical                        , intent(inout), optional          :: interrupt
    procedure(                    ), intent(inout), optional, pointer :: interruptProcedure
    !$GLC attributes unused :: interrupt, interruptProcedure

    ! Funnel the outflow gas abundances into the hot halo.
    call self%abundancesRate(rate)
    return
  end subroutine Node_Component_Hot_Halo_Very_Simple_Outflowing_Abundances_Rate

  double precision function Node_Component_Hot_Halo_Very_Simple_Outer_Radius(self)
    !!{
    Return the outer radius of the hot halo. Assumes a simple model in which this always equals the virial radius.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentHotHaloVerySimple
    implicit none
    class(nodeComponentHotHaloVerySimple), intent(inout) :: self

    Node_Component_Hot_Halo_Very_Simple_Outer_Radius=darkMatterHaloScale_%radiusVirial(self%hostNode)
    return
  end function Node_Component_Hot_Halo_Very_Simple_Outer_Radius

  !![
  <rateComputeTask>
   <unitName>Node_Component_Hot_Halo_Very_Simple_Rate_Compute</unitName>
  </rateComputeTask>
  !!]
  subroutine Node_Component_Hot_Halo_Very_Simple_Rate_Compute(node,interrupt,interruptProcedure,propertyType)
    !!{
    Compute the very simple hot halo component mass rate of change.
    !!}
    use :: Accretion_Halos , only : accretionModeTotal
    use :: Galacticus_Nodes, only : nodeComponentHotHalo   , nodeComponentHotHaloVerySimple, propertyTypeInactive, treeNode, &
         &                          defaultHotHaloComponent
    implicit none
    type            (treeNode            ), intent(inout)          :: node
    logical                               , intent(inout)          :: interrupt
    procedure       (                    ), intent(inout), pointer :: interruptProcedure
    integer                               , intent(in   )          :: propertyType
    class           (nodeComponentHotHalo)               , pointer :: hotHalo
    double precision                                               :: massAccretionRate   , failedMassAccretionRate

    ! Return immediately if inactive variables are requested.
    if (propertyType == propertyTypeInactive) return
    ! Return immediately if this class is not in use.
    if (.not.defaultHotHaloComponent%verySimpleIsActive()) return
    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    select type (hotHalo)
    class is (nodeComponentHotHaloVerySimple)
       ! Find the rate of gas mass accretion onto the halo.
       massAccretionRate      =accretionHalo_%accretionRate      (node,accretionModeTotal)
       failedMassAccretionRate=accretionHalo_%failedAccretionRate(node,accretionModeTotal)
       ! Apply accretion rates.
       if (      massAccretionRate > 0.0d0 .or. hotHalo%mass() > 0.0d0) &
            & call hotHalo%          massRate(      massAccretionRate,interrupt,interruptProcedure)
       if (failedMassAccretionRate > 0.0d0 .or. hotHalo%mass() > 0.0d0) &
            & call hotHalo%unaccretedMassRate(failedMassAccretionRate,interrupt,interruptProcedure)

       ! Next compute the cooling rate in this halo.
       call Node_Component_Hot_Halo_Very_Simple_Cooling_Rate         (node                                         )
       ! Pipe the cooling rate to which ever component claimed it.
       call Node_Component_Hot_Halo_Very_Simple_Push_To_Cooling_Pipes(node,rateCooling,interrupt,interruptProcedure)
    end select
    return
  end subroutine Node_Component_Hot_Halo_Very_Simple_Rate_Compute

  !![
  <scaleSetTask>
   <unitName>Node_Component_Hot_Halo_Very_Simple_Scale_Set</unitName>
  </scaleSetTask>
  !!]
  subroutine Node_Component_Hot_Halo_Very_Simple_Scale_Set(node)
    !!{
    Set scales for properties of {\normalfont \ttfamily node}.
    !!}
    use :: Abundances_Structure, only : unitAbundances
    use :: Galacticus_Nodes    , only : nodeComponentBasic     , nodeComponentHotHalo, nodeComponentHotHaloVerySimple, treeNode, &
         &                              defaultHotHaloComponent
    implicit none
    type            (treeNode            ), intent(inout), pointer :: node
    double precision                      , parameter              :: scaleMassRelative=1.0d-2
    class           (nodeComponentHotHalo)               , pointer :: hotHalo
    class           (nodeComponentBasic  )               , pointer :: basic
    double precision                                               :: massVirial

    ! Check if we are the default method.
    if (.not.defaultHotHaloComponent%verySimpleIsActive()) return
    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    ! Ensure that it is of the very simple class.
    select type (hotHalo)
    class is (nodeComponentHotHaloVerySimple)
       ! The the basic component.
       basic => node%basic()
       ! Get virial properties.
       massVirial=basic%mass()
       ! Set the scale.
       call hotHalo%          massScale(               massVirial*scaleMassRelative)
       call hotHalo%unaccretedMassScale(               massVirial*scaleMassRelative)
       call hotHalo%    abundancesScale(unitAbundances*massVirial*scaleMassRelative)
    end select
    return
  end subroutine Node_Component_Hot_Halo_Very_Simple_Scale_Set

  !![
  <mergerTreeInitializeTask>
   <unitName>Node_Component_Hot_Halo_Very_Simple_Tree_Initialize</unitName>
   <after>darkMatterProfile</after>
  </mergerTreeInitializeTask>
  !!]
  subroutine Node_Component_Hot_Halo_Very_Simple_Tree_Initialize(node)
    !!{
    Initialize the contents of the very simple hot halo component.
    !!}
    use :: Abundances_Structure, only : zeroAbundances
    use :: Accretion_Halos     , only : accretionModeTotal
    use :: Galacticus_Nodes    , only : defaultHotHaloComponent, nodeComponentHotHalo, nodeEvent, nodeEventSubhaloPromotion, &
          &                             treeNode
    implicit none
    type            (treeNode            ), intent(inout), pointer :: node
    class           (nodeComponentHotHalo)               , pointer :: hotHaloCurrent, hotHalo
    class           (nodeEvent           )               , pointer :: event
    double precision                                               :: massHotHalo   , massFailedHotHalo

    ! If the very simple hot halo is not active, then return immediately.
    if (associated(node%firstChild).or..not.defaultHotHaloComponent%verySimpleIsActive()) return
    ! Search for a subhalo promotion events associated with this node.
    event => node%event
    do while (associated(event))
       ! Check if this event:
       !  a) is a subhalo promotion event;
       !  b) has no associated task (which means this is the node being promoted to, not the node being promoted itself).
       ! Do not assign any mass to such nodes, as they should receive gas from the node which is promoted to them.
       select type (event)
       type is (nodeEventSubhaloPromotion)
          if (.not.associated(event%task)) return
       end select
       event => event%next
    end do

    ! Ensure that this module has been initialized.
    call Node_Component_Hot_Halo_Very_Simple_Initialize()

    ! Get the hot halo component.
    hotHaloCurrent => node%hotHalo()
    ! Ensure that it is of unspecified class.
    select type (hotHaloCurrent)
    type is (nodeComponentHotHalo)
       ! Get the mass of hot gas accreted and the mass that failed to accrete.
       massHotHalo      =accretionHalo_%accretedMass      (node,accretionModeTotal)
       massFailedHotHalo=accretionHalo_%failedAccretedMass(node,accretionModeTotal)
       ! If either is non-zero, then create a hot halo component and add these masses to it.
       if (massHotHalo > 0.0d0 .or. massFailedHotHalo > 0.0d0) then
          call Node_Component_Hot_Halo_Very_Simple_Create(node)
          hotHalo => node%hotHalo()
          call hotHalo%          massSet(      massHotHalo)
          call hotHalo%unaccretedMassSet(massFailedHotHalo)
          call hotHalo%    abundancesSet(   zeroAbundances)
       end if
    end select
    return
  end subroutine Node_Component_Hot_Halo_Very_Simple_Tree_Initialize

  subroutine satelliteMerger(self,node)
    !!{
    Remove any hot halo associated with {\normalfont \ttfamily node} before it merges with its host halo.
    !!}
    use :: Abundances_Structure, only : abundances          , zeroAbundances
    use :: Galacticus_Nodes    , only : nodeComponentHotHalo, nodeComponentHotHaloVerySimple, treeNode
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
    class is (nodeComponentHotHaloVerySimple)

       ! Find the node to merge with.
       nodeHost    => node%mergesWith(                 )
       hotHaloHost => nodeHost%hotHalo   (autoCreate=.true.)

       ! Move the hot halo to the host.
       call hotHaloHost%                    massSet(                                   &
            &                                                 hotHaloHost%mass      () &
            &                                                +hotHalo    %mass      () &
            &                                               )
       call hotHaloHost%              abundancesSet(                                   &
            &                                                 hotHaloHost%abundances() &
            &                                                +hotHalo    %abundances() &
            &                                               )
       call hotHalo%                    massSet(                                       &
            &                                                 0.0d0                    &
            &                                               )
       call hotHalo%              abundancesSet(                                       &
            &                                                 zeroAbundances           &
            &                                               )
    end select
    return
  end subroutine satelliteMerger

  subroutine nodePromotion(self,node)
    !!{
    Ensure that {\normalfont \ttfamily node} is ready for promotion to its parent. In this case, we simply update the hot halo mass of {\normalfont \ttfamily
    node} to account for any hot halo already in the parent.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentHotHalo, nodeComponentHotHaloVerySimple, treeNode
    implicit none
    class(*                   ), intent(inout)          :: self
    type (treeNode            ), intent(inout), target  :: node
    type (treeNode            )               , pointer :: nodeParent
    class(nodeComponentHotHalo)               , pointer :: hotHaloParent, hotHalo
    !$GLC attributes unused :: self

    hotHalo => node%hotHalo()
    ! Ensure that it is of specified class.
    select type (hotHalo)
    class is (nodeComponentHotHaloVerySimple)
       nodeParent    => node      %parent
       hotHaloParent => nodeParent%hotHalo(autoCreate=.true.)
       ! If the parent node has a hot halo component, then add it to that of this node, and perform other changes needed prior to
       ! promotion.
       select type (hotHaloParent)
          class is (nodeComponentHotHaloVerySimple)
          call hotHalo%unaccretedMassSet(                                &
               &                         +hotHalo      %unaccretedMass() &
               &                         +hotHaloParent%unaccretedMass() &
               &                        )
          call hotHalo%          massSet(                                &
               &                         +hotHalo      %          mass() &
               &                         +hotHaloParent%          mass() &
               &                        )
          call hotHalo%    abundancesSet(                                &
               &                         +hotHalo      %    abundances() &
               &                         +hotHaloParent%    abundances() &
               &                        )
       end select
    end select
    return
  end subroutine nodePromotion

  subroutine postEvolve(self,node)
    !!{
    Do processing of the node required after evolution.
    !!}
    use :: Abundances_Structure, only : abundances          , zeroAbundances
    use :: Galacticus_Nodes    , only : nodeComponentHotHalo, nodeComponentHotHaloVerySimple, treeNode
    implicit none
    class(*                   ), intent(inout) :: self
    type (treeNode            ), intent(inout) :: node
    type (treeNode            ), pointer       :: nodeParent
    class(nodeComponentHotHalo), pointer       :: hotHaloParent, hotHalo
    !$GLC attributes unused :: self

    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    select type (hotHalo)
    class is (nodeComponentHotHaloVerySimple)
       ! Check if this node is a satellite.
       if (node%isSatellite()) then
          ! Transfer any outflowed gas to the hot halo of the parent node.
          nodeParent => node%parent
          do while (nodeParent%isSatellite())
             nodeParent => nodeParent%parent
          end do
          hotHaloParent => nodeParent%hotHalo(autoCreate=.true.)
          call hotHaloParent%      massSet(hotHaloParent%mass      ()+hotHalo%mass      ())
          call hotHaloParent%abundancesSet(hotHaloParent%abundances()+hotHalo%abundances())
          call hotHalo      %      massSet(                                          0.0d0)
          call hotHalo      %abundancesSet(                                 zeroAbundances)
       end if
    end select
    return
  end subroutine postEvolve

  !![
  <nodeMergerTask>
   <unitName>Node_Component_Hot_Halo_Very_Simple_Node_Merger</unitName>
  </nodeMergerTask>
  !!]
  subroutine Node_Component_Hot_Halo_Very_Simple_Node_Merger(node)
    !!{
    Starve {\normalfont \ttfamily node} by transferring its hot halo to its parent.
    !!}
    use :: Abundances_Structure, only : abundances        , operator(*)         , zeroAbundances
    use :: Accretion_Halos     , only : accretionModeHot  , accretionModeTotal
    use :: Galacticus_Nodes    , only : nodeComponentBasic, nodeComponentHotHalo, nodeComponentHotHaloVerySimple, treeNode, defaultHotHaloComponent
    implicit none
    type            (treeNode            ), intent(inout) :: node
    type            (treeNode            ), pointer       :: nodeParent
    class           (nodeComponentHotHalo), pointer       :: hotHaloParent       , hotHalo
    class           (nodeComponentBasic  ), pointer       :: basic               , basicParent
    type            (abundances          ), save          :: massMetalsAccreted  , fractionMetalsAccreted, &
         &                                                   massMetalsReaccreted
    !$omp threadprivate(massMetalsAccreted,fractionMetalsAccreted,massMetalsReaccreted)
    double precision                                      :: massAccreted        , massUnaccreted        , &
         &                                                   fractionAccreted    , massReaccreted

    ! Return immediately if this class is not in use.
    if (.not.defaultHotHaloComponent%verySimpleIsActive()) return
    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    select type (hotHalo)
    class is (nodeComponentHotHaloVerySimple)
       ! Find the parent node and its hot halo component.
       nodeParent    => node      %parent
       hotHaloParent => nodeParent%hotHalo(autoCreate=.true.)
       ! Get the basic components.
       basic       => node      %basic()
       basicParent => nodeParent%basic()
       ! Any gas that failed to be accreted by this halo is always transferred to the parent.
       call hotHaloParent%unaccretedMassSet(hotHaloParent%unaccretedMass()+hotHalo%unaccretedMass())
       call       hotHalo%unaccretedMassSet(                                                  0.0d0)
       ! Move the hot halo to the parent. We leave the hot halo in place even if it is starved, since outflows will accumulate
       ! to this hot halo (and will be moved to the parent at the end of the evolution timestep).
       call hotHaloParent%          massSet(hotHaloParent%          mass()+hotHalo%          mass())
       call hotHaloParent%    abundancesSet(hotHaloParent%    abundances()+hotHalo%    abundances())
       call       hotHalo%          massSet(                                                  0.0d0)
       call       hotHalo%    abundancesSet(                                         zeroAbundances)
       ! Finally, since the parent node is undergoing mass growth through this merger we potentially return some of the unaccreted
       ! gas to the hot phase.
       !! First, find the masses of hot and failed mass the node would have if it formed instantaneously.
       massAccreted  =accretionHalo_%accretedMass      (nodeParent,accretionModeTotal)
       massUnaccreted=accretionHalo_%failedAccretedMass(nodeParent,accretionModeTotal)
       !! Find the fraction of mass that would be successfully accreted.
       fractionAccreted=+  massAccreted   &
            &           /(                &
            &             +massAccreted   &
            &             +massUnaccreted &
            &            )
       !! Find the change in the unaccreted mass.
       massReaccreted=+hotHaloParent   %unaccretedMass() &
            &         *fractionAccreted                  &
            &         *basic           %          mass() &
            &         /basicParent     %          mass()
       !! Reaccrete the gas,
       call hotHaloParent%unaccretedMassSet(hotHaloParent%unaccretedMass()-massReaccreted)
       call hotHaloParent%          massSet(hotHaloParent%          mass()+massReaccreted)
       ! Compute the reaccreted metals.
       !! First, find the metal mass the node would have if it formed instantaneously.
       massMetalsAccreted=accretionHalo_%accretedMassMetals(nodeParent,accretionModeHot)
       !! Find the mass fraction of metals that would be successfully accreted.
       fractionMetalsAccreted=+  massMetalsAccreted &
            &                 /(                    &
            &                   +massAccreted       &
            &                   +massUnaccreted     &
            &                  )
       !! Find the change in the unaccreted mass.
       massMetalsReaccreted=+hotHaloParent   %unaccretedMass() &
            &               *fractionMetalsAccreted            &
            &               *basic           %          mass() &
            &               /basicParent     %          mass()
       !! Reaccrete the metals.
       call hotHaloParent%abundancesSet(hotHaloParent%abundances()+massMetalsReaccreted)
    end select
    return
  end subroutine Node_Component_Hot_Halo_Very_Simple_Node_Merger

  subroutine Node_Component_Hot_Halo_Very_Simple_Cooling_Rate(node)
    !!{
    Get and store the cooling rate for {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentHotHalo, treeNode
    implicit none
    type (treeNode            ), intent(inout) :: node
    class(nodeComponentHotHalo), pointer       :: hotHalo

    if (.not.gotCoolingRate) then
       ! Get the hot halo component.
       hotHalo => node%hotHalo()
       if (hotHalo%mass() > 0.0d0) then
          ! Get the cooling time.
          rateCooling=coolingRate_%rate(node)
       else
          rateCooling=0.0d0
       end if
       ! Flag that cooling rate has now been computed.
       gotCoolingRate=.true.
    end if
    return
  end subroutine Node_Component_Hot_Halo_Very_Simple_Cooling_Rate

  subroutine Node_Component_Hot_Halo_Very_Simple_Create(node)
    !!{
    Creates a very simple hot halo component for {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentHotHalo, treeNode
    implicit none
    type (treeNode            ), intent(inout), pointer :: node
    class(nodeComponentHotHalo)               , pointer :: hotHalo

    ! Ensure that this module has been initialized.
    call Node_Component_Hot_Halo_Very_Simple_Initialize()

    ! Create the component.
    hotHalo => node%hotHalo(autoCreate=.true.)
    return
  end subroutine Node_Component_Hot_Halo_Very_Simple_Create

  !![
  <galacticusStateStoreTask>
   <unitName>Node_Component_Hot_Halo_Very_Simple_State_Store</unitName>
  </galacticusStateStoreTask>
  !!]
  subroutine Node_Component_Hot_Halo_Very_Simple_State_Store(stateFile,gslStateFile,stateOperationID)
    !!{
    Store object state,
    !!}
    use            :: Display      , only : displayMessage, verbosityLevelInfo
    use, intrinsic :: ISO_C_Binding, only : c_ptr         , c_size_t
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call displayMessage('Storing state for: componentHotHalo -> verySimple',verbosity=verbosityLevelInfo)
    !![
    <stateStore variables="darkMatterHaloScale_ coolingRate_ accretionHalo_"/>
    !!]
    return
  end subroutine Node_Component_Hot_Halo_Very_Simple_State_Store

  !![
  <galacticusStateRetrieveTask>
   <unitName>Node_Component_Hot_Halo_Very_Simple_State_Restore</unitName>
  </galacticusStateRetrieveTask>
  !!]
  subroutine Node_Component_Hot_Halo_Very_Simple_State_Restore(stateFile,gslStateFile,stateOperationID)
    !!{
    Retrieve object state.
    !!}
    use            :: Display      , only : displayMessage, verbosityLevelInfo
    use, intrinsic :: ISO_C_Binding, only : c_ptr         , c_size_t
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call displayMessage('Retrieving state for: componentHotHalo -> verySimple',verbosity=verbosityLevelInfo)
    !![
    <stateRestore variables="darkMatterHaloScale_ coolingRate_ accretionHalo_"/>
    !!]
    return
  end subroutine Node_Component_Hot_Halo_Very_Simple_State_Restore

end module Node_Component_Hot_Halo_Very_Simple
