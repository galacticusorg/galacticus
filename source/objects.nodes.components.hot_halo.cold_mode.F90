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
Contains a module which implements an extension to the standard hot halo node component which
supports a cold mode reservoir.
!!}

module Node_Component_Hot_Halo_Cold_Mode
  !!{
  Implements an extension to the standard hot halo node component which supports a cold mode
  reservoir.
  !!}
  use :: Accretion_Halos                      , only : accretionHaloClass
  use :: Cosmology_Parameters                 , only : cosmologyParametersClass
  use :: Dark_Matter_Halo_Scales              , only : darkMatterHaloScaleClass
  use :: Hot_Halo_Cold_Mode_Mass_Distributions, only : hotHaloColdModeMassDistributionClass
  implicit none
  private
  public :: Node_Component_Hot_Halo_Cold_Mode_Initialize       , Node_Component_Hot_Halo_Cold_Mode_Node_Merger        , &
       &    Node_Component_Hot_Halo_Cold_Mode_Scale_Set        , Node_Component_Hot_Halo_Cold_Mode_Tree_Initialize    , &
       &    Node_Component_Hot_Halo_Cold_Mode_Thread_Initialize, Node_Component_Hot_Halo_Cold_Mode_Thread_Uninitialize, &
       &    Node_Component_Hot_Halo_Cold_Mode_State_Store      , Node_Component_Hot_Halo_Cold_Mode_State_Restore

  !![
  <component>
   <class>hotHalo</class>
   <name>coldMode</name>
   <extends>
     <class>hotHalo</class>
     <name>standard</name>
   </extends>
   <isDefault>false</isDefault>
   <properties>
    <property>
      <name>massCold</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" />
      <output unitsInSI="massSolar" comment="Mass of cold-mode gas in the hot halo."/>
    </property>
    <property>
      <name>abundancesCold</name>
      <type>abundances</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" />
      <output unitsInSI="massSolar" comment="Mass of metals in the cold-mode of the hot halo."/>
    </property>
    <property>
      <name>angularMomentumCold</name>
      <type>double</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" createIfNeeded="true" />
      <output unitsInSI="massSolar*megaParsec*kilo" comment="Angular momentum of cold-mode gas in the hot halo."/>
    </property>
    <property>
      <name>massTotal</name>
      <attributes isSettable="false" isGettable="true" isEvolvable="false" isVirtual="true" />
      <type>double</type>
      <rank>0</rank>
      <getFunction>Node_Component_Hot_Halo_Cold_Mode_Mass_Total</getFunction>
    </property>
   </properties>
   <bindings>
     <binding method="massDistribution" isDeferred="true" >
      <interface>
       <type>class(massDistributionClass), pointer</type>
       <rank>0</rank>
       <module>Galactic_Structure_Options, only : enumerationWeightByType, enumerationComponentTypeType, enumerationMassTypeType</module>
       <module>Mass_Distributions        , only : massDistributionClass                                                         </module>
       <self pass="true" intent="inout" />
       <argument>type   (enumerationComponentTypeType), intent(in   ), optional :: componentType</argument>
       <argument>type   (enumerationMassTypeType     ), intent(in   ), optional :: massType     </argument>
       <argument>type   (enumerationWeightByType     ), intent(in   ), optional :: weightBy     </argument>
       <argument>integer                              , intent(in   ), optional :: weightIndex  </argument>
      </interface>
     </binding>
     <binding method="massBaryonic" function="Node_Component_Hot_Halo_Cole_Mode_Mass_Baryonic"/>
   </bindings>
   <functions>objects.nodes.components.hot_halo.cold_mode.bound_functions.inc</functions>
  </component>
  !!]

  ! Objects used by this component.
  class(accretionHaloClass                  ), pointer :: accretionHalo_
  class(cosmologyParametersClass            ), pointer :: cosmologyParameters_
  class(darkMatterHaloScaleClass            ), pointer :: darkMatterHaloScale_
  class(hotHaloColdModeMassDistributionClass), pointer :: hotHaloColdModeMassDistribution_
  !$omp threadprivate(accretionHalo_,cosmologyParameters_,darkMatterHaloScale_,hotHaloColdModeMassDistribution_)

  ! Internal count of abundances.
  integer :: abundancesCount

  ! A threadprivate object used to track to which thread events are attached.
  integer :: thread
  !$omp threadprivate(thread)

  ! Procedure pointer to mass distribution function.
  procedure(Node_Component_Hot_Halo_Cold_Mode_Mass_Distribution), pointer :: Node_Component_Hot_Halo_Cold_Mode_Mass_Distribution_
  
contains

  !![
  <nodeComponentInitializationTask function="Node_Component_Hot_Halo_Cold_Mode_Initialize"/>
  !!]
  subroutine Node_Component_Hot_Halo_Cold_Mode_Initialize(parameters)
    !!{
    Initializes the tree node hot halo methods module.
    !!}
    use :: Abundances_Structure, only : Abundances_Property_Count
    use :: Galacticus_Nodes    , only : defaultHotHaloComponent  , nodeComponentHotHaloColdMode
    use :: Input_Parameters    , only : inputParameter           , inputParameters
    implicit none
    type(inputParameters             ), intent(inout) :: parameters
    type(nodeComponentHotHaloColdMode)                :: hotHalo
    type(inputParameters             )                :: subParameters

    ! Initialize the module if necessary.
    !$omp critical (Node_Component_Hot_Halo_Cold_Mode_Initialize)
    if (defaultHotHaloComponent%coldModeIsActive()) then
       ! Get numbers of abundance properties.
       abundancesCount=Abundances_Property_Count()
       ! Find our parameters.
       subParameters=parameters%subParameters('componentHotHalo')
       ! Bind the mass distribution function.
       Node_Component_Hot_Halo_Cold_Mode_Mass_Distribution_ => Node_Component_Hot_Halo_Cold_Mode_Mass_Distribution
       call hotHalo%massDistributionFunction(Node_Component_Hot_Halo_Cold_Mode_Mass_Distribution_)
    end if
    !$omp end critical (Node_Component_Hot_Halo_Cold_Mode_Initialize)
    return
  end subroutine Node_Component_Hot_Halo_Cold_Mode_Initialize

  !![
  <nodeComponentThreadInitializationTask function="Node_Component_Hot_Halo_Cold_Mode_Thread_Initialize"/>
  !!]
  subroutine Node_Component_Hot_Halo_Cold_Mode_Thread_Initialize(parameters)
    !!{
    Initializes the tree node hot halo cold mode methods module.
    !!}
    use :: Events_Hooks                         , only : nodePromotionEvent      , satelliteMergerEvent, openMPThreadBindingAtLevel, dependencyRegEx, &
         &                                               dependencyDirectionAfter, haloFormationEvent
    use :: Galacticus_Nodes                     , only : defaultHotHaloComponent
    use :: Hot_Halo_Cold_Mode_Density_Core_Radii, only : hotHaloColdModeCoreRadii
    use :: Input_Parameters                     , only : inputParameter          , inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters
    type(dependencyRegEx), dimension(1)  :: dependencies
    type(inputParameters)                :: subParameters

    if (defaultHotHaloComponent%coldModeIsActive()) then
       ! Find our parameters.
       subParameters=parameters%subParameters('componentHotHalo')
       !![
       <objectBuilder class="cosmologyParameters"             name="cosmologyParameters_"             source="subParameters"/>
       <objectBuilder class="darkMatterHaloScale"             name="darkMatterHaloScale_"             source="subParameters"/>
       <objectBuilder class="accretionHalo"                   name="accretionHalo_"                   source="subParameters"/>
       <objectBuilder class="hotHaloColdModeMassDistribution" name="hotHaloColdModeMassDistribution_" source="subParameters"/>
       !!]
       dependencies(1)=dependencyRegEx(dependencyDirectionAfter,'^remnantStructure:')
       call nodePromotionEvent  %attach(thread,nodePromotion  ,openMPThreadBindingAtLevel,label='nodeComponentHotHaloColdMode'                          )
       call satelliteMergerEvent%attach(thread,satelliteMerger,openMPThreadBindingAtLevel,label='nodeComponentHotHaloColdMode',dependencies=dependencies)
       call haloFormationEvent  %attach(thread,haloFormation  ,openMPThreadBindingAtLevel,label='nodeComponentHotHaloColdMode'                          )
    end if
    return
  end subroutine Node_Component_Hot_Halo_Cold_Mode_Thread_Initialize

  !![
  <nodeComponentThreadUninitializationTask function="Node_Component_Hot_Halo_Cold_Mode_Thread_Uninitialize"/>
  !!]
  subroutine Node_Component_Hot_Halo_Cold_Mode_Thread_Uninitialize()
    !!{
    Uninitializes the tree node hot halo cold mode methods module.
    !!}
    use :: Events_Hooks    , only : nodePromotionEvent     , satelliteMergerEvent, haloFormationEvent
    use :: Galacticus_Nodes, only : defaultHotHaloComponent
    implicit none

    if (defaultHotHaloComponent%coldModeIsActive()) then
       !![
       <objectDestructor name="cosmologyParameters_"            />
       <objectDestructor name="darkMatterHaloScale_"            />
       <objectDestructor name="accretionHalo_"                  />
       <objectDestructor name="hotHaloColdModeMassDistribution_"/>
       !!]
       if (nodePromotionEvent  %isAttached(thread,nodePromotion  )) call nodePromotionEvent  %detach(thread,nodePromotion  )
       if (satelliteMergerEvent%isAttached(thread,satelliteMerger)) call satelliteMergerEvent%detach(thread,satelliteMerger)
       if (haloFormationEvent  %isAttached(thread,haloFormation  )) call haloFormationEvent  %detach(thread,haloFormation  )
    end if
    return
  end subroutine Node_Component_Hot_Halo_Cold_Mode_Thread_Uninitialize

  !![
  <scaleSetTask function="Node_Component_Hot_Halo_Cold_Mode_Scale_Set"/>
  !!]
  subroutine Node_Component_Hot_Halo_Cold_Mode_Scale_Set(node)
    !!{
    Set scales for properties of \mono{node}.
    !!}
    use :: Abundances_Structure, only : unitAbundances
    use :: Galacticus_Nodes    , only : nodeComponentBasic     , nodeComponentHotHalo, nodeComponentHotHaloColdMode, treeNode, &
         &                              defaultHotHaloComponent
    implicit none
    type            (treeNode            ), intent(inout), pointer :: node
    class           (nodeComponentHotHalo)               , pointer :: hotHalo
    class           (nodeComponentBasic  )               , pointer :: basic
    double precision                      , parameter              :: scaleMassRelative   =1.0d-3
    double precision                      , parameter              :: scaleRadiusRelative =1.0d+0
    double precision                                               :: massVirial                 , radiusVirial, &
         &                                                            velocityVirial

    ! Check if we are the default method.
    if (.not.defaultHotHaloComponent%coldModeIsActive()) return
    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    ! Ensure that it is of the cold mode class.
    select type (hotHalo)
    class is (nodeComponentHotHaloColdMode)
       ! The the basic component.
       basic => node%basic()
       ! Get virial properties.
       massVirial    =basic%mass()
       radiusVirial  =darkMatterHaloScale_%radiusVirial  (node)
       velocityVirial=darkMatterHaloScale_%velocityVirial(node)
       call    hotHalo%           massColdScale(               massVirial                            *scaleMassRelative)
       call    hotHalo%     abundancesColdScale(unitAbundances*massVirial                            *scaleMassRelative)
       call    hotHalo%angularMomentumColdScale(               massVirial*radiusVirial*velocityVirial*scaleMassRelative)
    end select
    return
  end subroutine Node_Component_Hot_Halo_Cold_Mode_Scale_Set

  !![
  <mergerTreeInitializeTask function="Node_Component_Hot_Halo_Cold_Mode_Tree_Initialize" after="Node_Component_Hot_Halo_Standard_Tree_Initialize"/>
  !!]
  subroutine Node_Component_Hot_Halo_Cold_Mode_Tree_Initialize(node)
    !!{
    Initialize the contents of the hot halo component for any sub-resolution accretion (i.e. the gas that would have been
    accreted if the merger tree had infinite resolution).
    !!}
    use :: Accretion_Halos , only : accretionModeCold
    use :: Galacticus_Nodes, only : defaultHotHaloComponent  , nodeComponentBasic, nodeComponentHotHalo, nodeEvent, &
          &                         nodeEventSubhaloPromotion, treeNode          , nodeComponentSpin
    implicit none
    type            (treeNode            ), intent(inout), pointer :: node
    class           (nodeComponentHotHalo)               , pointer :: hotHalo
    class           (nodeComponentBasic  )               , pointer :: basic
    class           (nodeComponentSpin   )               , pointer :: spin
    class           (nodeEvent           )               , pointer :: event
    double precision                                               :: angularMomentum, coldModeMass

    ! If the node has a child or the standard hot halo is not active, then return immediately.
    if (associated(node%firstChild).or..not.defaultHotHaloComponent%coldModeIsActive()) return
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
    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    ! Get the mass of cold mode gas accreted.
    coldModeMass=accretionHalo_%accretedMass(node,accretionModeCold)
    ! If non-zero, then create a hot halo component and add to it.
    if (coldModeMass > 0.0d0) then
       ! Ensure that it is of unspecified class.
       hotHalo => node%hotHalo(autoCreate=.true.)
       basic   => node%basic  (                 )
       spin    => node%spin   (                 )
       call hotHalo%massColdSet(coldModeMass)
       ! Also add the appropriate angular momentum.
       angularMomentum=+      coldModeMass      &
            &          *spin %angularMomentum() &
            &          /basic%mass           ()
       call hotHalo%angularMomentumColdSet(angularMomentum)
       ! Add the appropriate abundances.
       call hotHalo%abundancesColdSet(accretionHalo_%accretedMassMetals(node,accretionModeCold))
    end if
    return
  end subroutine Node_Component_Hot_Halo_Cold_Mode_Tree_Initialize

  !![
  <nodeMergerTask function="Node_Component_Hot_Halo_Cold_Mode_Node_Merger" before="Node_Component_Hot_Halo_Standard_Node_Merger"/>
  !!]
  subroutine Node_Component_Hot_Halo_Cold_Mode_Node_Merger(node)
    !!{
    Starve \mono{node} by transferring its hot halo to its parent.
    !!}
    use :: Abundances_Structure                 , only : abundances                     , operator(*)            , zeroAbundances
    use :: Accretion_Halos                      , only : accretionModeCold              , accretionModeTotal
    use :: Galactic_Structure_Options           , only : componentTypeAll               , massTypeBaryonic
    use :: Galacticus_Nodes                     , only : nodeComponentBasic             , nodeComponentHotHalo   , nodeComponentHotHaloColdMode, nodeComponentSpin, &
          &                                              treeNode                       , defaultHotHaloComponent
    use :: Node_Component_Hot_Halo_Standard_Data, only : fractionBaryonLimitInNodeMerger, starveSatellites
    use :: Mass_Distributions                   , only : massDistributionClass
    implicit none
    type            (treeNode             ), intent(inout) :: node
    type            (treeNode             ), pointer       :: nodeParent
    class           (nodeComponentHotHalo ), pointer       :: hotHaloParent          , hotHalo
    class           (nodeComponentSpin    ), pointer       :: spinParent
    class           (nodeComponentBasic   ), pointer       :: basicParent            , basic
    class           (massDistributionClass), pointer       :: massDistribution_
    double precision                                       :: baryonicMassCurrent    , baryonicMassMaximum   , &
         &                                                    fractionRemove         , massAccretedCold      , &
         &                                                    massAccreted           , massUnaccreted        , &
         &                                                    angularMomentumAccreted, massReaccreted        , &
         &                                                    fractionAccreted
    type            (abundances           ), save          :: massMetalsAccreted     , fractionMetalsAccreted, &
         &                                                    massMetalsReaccreted
    !$omp threadprivate(massMetalsAccreted,fractionMetalsAccreted,massMetalsReaccreted)

    ! Return immediately if this class is not in use.
    if (.not.defaultHotHaloComponent%coldModeIsActive()) return
    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    ! Ensure that it is of cold mode class.
    select type (hotHalo)
    class is (nodeComponentHotHaloColdMode)
       ! Find the parent node and its hot halo and angular momentum components.
       nodeParent    => node      %parent
       hotHaloParent => nodeParent%hotHalo(autoCreate=.true.)
       spinParent    => nodeParent%spin   (                 )
       basicParent   => nodeParent%basic  (                 )
       basic         => node      %basic  (                 )
       ! Since the parent node is undergoing mass growth through this merger we potentially return some of the unaccreted gas to
       ! the hot phase.
       !! First, find the masses of hot and failed mass the node would have if it formed instantaneously.
       massAccretedCold=accretionHalo_%      accretedMass(nodeParent,accretionModeCold )
       massAccreted    =accretionHalo_%      accretedMass(nodeParent,accretionModeTotal)
       massUnaccreted  =accretionHalo_%failedAccretedMass(nodeParent,accretionModeTotal)
       !! Find the fraction of mass that would be successfully accreted.
       fractionAccreted=+  massAccretedCold &
            &           /(                  &
            &             +massAccreted     &
            &             +massUnaccreted   &
            &            )
       !! Find the change in the unaccreted mass.
       massReaccreted=+hotHaloParent   %unaccretedMass() &
            &         *fractionAccreted                  &
            &         *basic           %          mass() &
            &         /basicParent     %          mass()
       !! Reaccrete the gas.
       call hotHaloParent%unaccretedMassSet(hotHaloParent%unaccretedMass()-massReaccreted)
       call hotHaloParent%      massColdSet(hotHaloParent%      massCold()+massReaccreted)
       ! Compute the reaccreted angular momentum.
       angularMomentumAccreted=+            massReaccreted    &
            &                  *spinParent %angularMomentum() &
            &                  /basicParent%mass           ()
       call hotHaloParent%angularMomentumColdSet(hotHaloParent%angularMomentumCold()+angularMomentumAccreted)
       ! Compute the reaccreted metals.
       !! First, find the metal mass the node would have if it formed instantaneously.
       massMetalsAccreted=accretionHalo_%accretedMassMetals(nodeParent,accretionModeCold)
       !! Find the mass fraction of metals that would be successfully accreted.
       fractionMetalsAccreted=+  massMetalsAccreted &
            &                 /(                    &
            &                   +massAccreted       &
            &                   +massUnaccreted     &
            &                  )
       !! Find the change in the unaccreted mass.
       massMetalsReaccreted=+hotHaloParent         %unaccretedMass() &
            &               *fractionMetalsAccreted                  &
            &               *basic                 %          mass() &
            &               /basicParent           %          mass()
       !! Reaccrete the metals.
       call hotHaloParent%abundancesColdSet(hotHaloParent%abundancesCold()+massMetalsReaccreted)
       ! Determine if starvation is to be applied.
       if (starveSatellites) then
          ! Move the hot halo to the parent. We leave the hot halo in place even if it is starved, since outflows will accumulate to
          ! this hot halo (and will be moved to the parent at the end of the evolution timestep).
          call hotHaloParent%           massColdSet(                                      &
               &                                     hotHaloParent %massCold           () &
               &                                    +hotHalo       %massCold           () &
               &                                   )
          call hotHaloParent%angularMomentumColdSet(                                      &
               &                                     hotHaloParent %angularMomentumCold() &
               &                                    +hotHalo       %massCold           () &
               &                                    *spinParent    %angularMomentum    () &
               &                                    /basicParent   %mass               () &
               &                                   )
          call hotHalo      %           massColdSet(                                      &
               &                                     0.0d0                                &
               &                                   )
          call hotHalo      %angularMomentumColdSet(                                      &
               &                                     0.0d0                                &
               &                                   )
          call hotHaloParent%     abundancesColdSet(                                      &
               &                                     hotHaloParent %abundancesCold     () &
               &                                    +hotHalo       %abundancesCold     () &
               &                                   )
          call hotHalo      %     abundancesColdSet(                                      &
               &                                     zeroAbundances                       &
               &                                   )
          ! Check if the baryon fraction in the parent hot halo exceeds the universal value. If it does, mitigate this by moving
          ! some of the mass to the failed accretion reservoir.
          if (fractionBaryonLimitInNodeMerger) then
             massDistribution_   =>  nodeParent          %massDistribution(massType=massTypeBaryonic)
             baryonicMassMaximum =  +basicParent         %mass            (                         ) &
                  &                 *cosmologyParameters_%omegaBaryon     (                         ) &
                  &                 /cosmologyParameters_%omegaMatter     (                         )
             baryonicMassCurrent =  +massDistribution_   %massTotal       (                         )
             !![
	     <objectDestructor name="massDistribution_"/>
	     !!]
             if (baryonicMassCurrent > baryonicMassMaximum .and. hotHaloParent%mass()+hotHaloParent%massCold() > 0.0d0) then
                fractionRemove=min((baryonicMassCurrent-baryonicMassMaximum)/hotHaloParent%massTotal(),1.0d0)
                call hotHaloParent%     unaccretedMassSet(                                                            &
                     &                                     hotHaloParent%unaccretedMass     ()                        &
                     &                                    +hotHaloParent%massCold           ()*       fractionRemove  &
                     &                                   )
                call hotHaloParent%           massColdSet( hotHaloParent%massCold           ()*(1.0d0-fractionRemove))
                call hotHaloParent%angularMomentumColdSet( hotHaloParent%angularMomentumCold()*(1.0d0-fractionRemove))
                call hotHaloParent%     abundancesColdSet( hotHaloParent%abundancesCold     ()*(1.0d0-fractionRemove))
             end if
          end if
       end if
    end select
    return
  end subroutine Node_Component_Hot_Halo_Cold_Mode_Node_Merger

  subroutine satelliteMerger(self,node)
    !!{
    Remove any cold mode gas associated with \mono{node} before it merges with its host halo.
    !!}
    use :: Abundances_Structure                 , only : abundances          , zeroAbundances
    use :: Galacticus_Nodes                     , only : nodeComponentHotHalo, nodeComponentHotHaloColdMode, nodeComponentSpin, nodeComponentBasic, &
         &                                               treeNode
    use :: Node_Component_Hot_Halo_Standard_Data, only : starveSatellites
    implicit none
    class(*                   ), intent(inout) :: self
    type (treeNode            ), intent(inout) :: node
    type (treeNode            ), pointer       :: nodeHost
    class(nodeComponentBasic  ), pointer       :: basicHost
    class(nodeComponentHotHalo), pointer       :: hotHaloHost, hotHalo
    class(nodeComponentSpin   ), pointer       :: spinHost
    !$GLC attributes unused :: self

    ! Return immediately if satellites are starved, as in that case there is no hot halo to transfer.
    if (starveSatellites) return
    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    ! Ensure that it is of unspecified class.
    select type (hotHalo)
    class is (nodeComponentHotHaloColdMode)
       ! Find the node with which to merge.
       nodeHost    => node    %mergesWith(                 )
       hotHaloHost => nodeHost%hotHalo   (autoCreate=.true.)
       basicHost   => nodeHost%basic     (                 )
       spinHost    => nodeHost%spin      (                 )
       ! Move the cold mode to the host.
       call hotHaloHost%               massSet(                                 &
            &                                   hotHaloHost  %mass           () &
            &                                  +hotHalo      %massCold       () &
            &                                 )
       call hotHaloHost%    angularMomentumSet(                                 &
            &                                   hotHaloHost  %angularMomentum() &
            &                                  +hotHalo      %massCold       () &
            &                                  *spinHost     %angularMomentum() &
            &                                  /basicHost    %mass           () &
            &                                 )
       call hotHalo    %           massColdSet(                                 &
            &                                   0.0d0                           &
            &                                 )
       call hotHalo    %angularMomentumColdSet(                                 &
            &                                   0.0d0                           &
            &                                 )
       call hotHaloHost%         abundancesSet(                                 &
            &                                   hotHaloHost  %abundances     () &
            &                                  +hotHalo      %abundancesCold () &
            &                                 )
       call hotHalo    %     abundancesColdSet(                                 &
            &                                  zeroAbundances                   &
            &                                 )
    end select
    return
  end subroutine satelliteMerger

  subroutine nodePromotion(self,node)
    !!{
    Ensure that \mono{node} is ready for promotion to its parent. In this case, we simply
    update the cold mode mass of \mono{node} to account for any cold mode gas already in the
    parent.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentHotHalo, nodeComponentHotHaloColdMode, treeNode
    implicit none
    class(*                   ), intent(inout) :: self
    type (treeNode            ), intent(inout) :: node
    type (treeNode            ), pointer       :: nodeParent
    class(nodeComponentHotHalo), pointer       :: hotHaloParent, hotHalo
    !$GLC attributes unused :: self

    hotHalo       => node      %hotHalo(autoCreate=.true.)
    nodeParent    => node      %parent
    hotHaloParent => nodeParent%hotHalo(autoCreate=.true.)
    ! If the parent node has a hot halo component, then add its cold mode to that of this node,
    ! and perform other changes needed prior to promotion.
    select type (hotHaloParent)
    class is (nodeComponentHotHaloColdMode)
       call hotHalo%           massColdSet(                                      &
            &                                hotHalo      %massCold           () &
            &                               +hotHaloParent%massCold           () &
            &                              )
       call hotHalo%angularMomentumColdSet(                                      &
            &                                hotHalo      %angularMomentumCold() &
            &                               +hotHaloParent%angularMomentumCold() &
            &                              )
       call hotHalo%     abundancesColdSet(                                      &
            &                                hotHalo      %abundancesCold     () &
            &                               +hotHaloParent%abundancesCold     () &
            &                              )
    end select
    return
  end subroutine nodePromotion

  subroutine haloFormation(self,node)
    !!{
    Updates the hot halo gas distribution at a formation event, if requested.
    !!}
    use :: Abundances_Structure                 , only : abundances              , zeroAbundances
    use :: Galacticus_Nodes                     , only : nodeComponentHotHalo    , nodeComponentHotHaloColdMode, treeNode
    use :: Node_Component_Hot_Halo_Standard_Data, only : outflowReturnOnFormation
    implicit none
    class(*                   ), intent(inout) :: self
    type (treeNode            ), intent(inout) :: node
    class(nodeComponentHotHalo), pointer       :: hotHalo
    !$GLC attributes unused :: self

    ! Return immediately if return of outflowed gas on formation events is not requested.
    if (.not.outflowReturnOnFormation) return
    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    ! Ensure that it is of unspecified class.
    select type (hotHalo)
    class is (nodeComponentHotHaloColdMode)
       ! Transfer mass, angular momentum and abundances.
       call hotHalo%                    massSet(                                    &
            &                                    hotHalo%         mass           () &
            &                                   +hotHalo%outflowedMass           () &
            &                                  )
       call hotHalo%         angularMomentumSet(                                    &
            &                                    hotHalo%         angularMomentum() &
            &                                   +hotHalo%outflowedAngularMomentum() &
            &                                  )
       call hotHalo%              abundancesSet(                                    &
            &                                    hotHalo%         abundances     () &
            &                                   +hotHalo%outflowedAbundances     () &
            &                                  )
       call hotHalo%           outflowedMassSet(                                    &
            &                                    0.0d0                              &
            &                                  )
       call hotHalo%outflowedAngularMomentumSet(                                    &
            &                                    0.0d0                              &
            &                                  )
       call hotHalo%     outflowedAbundancesSet(                                    &
            &                                    zeroAbundances                     &
            &                                  )
    end select
    return
  end subroutine haloFormation

  !![
  <stateStoreTask function="Node_Component_Hot_Halo_Cold_Mode_State_Store"/>
  !!]
  subroutine Node_Component_Hot_Halo_Cold_Mode_State_Store(stateFile,gslStateFile,stateOperationID)
    !!{
    Store object state,
    !!}
    use            :: Display      , only : displayMessage, verbosityLevelInfo
    use, intrinsic :: ISO_C_Binding, only : c_ptr         , c_size_t
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call displayMessage('Storing state for: componentHotHalo -> coldMode',verbosity=verbosityLevelInfo)
    !![
    <stateStore variables="accretionHalo_ cosmologyParameters_ hotHaloColdModeMassDistribution_"/>
    !!]
    return
  end subroutine Node_Component_Hot_Halo_Cold_Mode_State_Store

  !![
  <stateRetrieveTask function="Node_Component_Hot_Halo_Cold_Mode_State_Restore"/>
  !!]
  subroutine Node_Component_Hot_Halo_Cold_Mode_State_Restore(stateFile,gslStateFile,stateOperationID)
    !!{
    Retrieve object state.
    !!}
    use            :: Display      , only : displayMessage, verbosityLevelInfo
    use, intrinsic :: ISO_C_Binding, only : c_ptr         , c_size_t
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile

    call displayMessage('Retrieving state for: componentHotHalo -> coldMode',verbosity=verbosityLevelInfo)
    !![
    <stateRestore variables="accretionHalo_ cosmologyParameters_ hotHaloColdModeMassDistribution_"/>
    !!]
    return
  end subroutine Node_Component_Hot_Halo_Cold_Mode_State_Restore

  function Node_Component_Hot_Halo_Cold_Mode_Mass_Distribution(self,componentType,massType,weightBy,weightIndex) result(massDistribution_)
    !!{
    Return the mass distribution associated with the hot halo.
    !!}
    use :: Galacticus_Nodes          , only : nodeComponentHotHaloStandard, nodeComponentHotHaloColdMode
    use :: Galactic_Structure_Options, only : enumerationWeightByType     , enumerationComponentTypeType, enumerationMassTypeType  , componentTypeColdHalo, &
         &                                    massTypeGaseous
    use :: Mass_Distributions        , only : massDistributionClass       , kinematicsDistributionLocal , massDistributionComposite, massDistributionList , &
         &                                    massDistributionMatches_
    implicit none
    class  (massDistributionClass       ), pointer                 :: massDistributionHotMode   , massDistributionColdMode, &
         &                                                            massDistribution_
    type   (kinematicsDistributionLocal ), pointer                 :: kinematicsDistribution_
    type   (massDistributionComposite   ), pointer                 :: massDistributionTotal
    type   (massDistributionList        ), pointer                 :: massDistributionComponents
    class  (nodeComponentHotHaloStandard), intent(inout)           :: self
    type   (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type   (enumerationMassTypeType     ), intent(in   ), optional :: massType
    type   (enumerationWeightByType     ), intent(in   ), optional :: weightBy
    integer                              , intent(in   ), optional :: weightIndex

    select type (self)
    class is (nodeComponentHotHaloColdMode)
       if (massDistributionMatches_(componentTypeColdHalo,massTypeGaseous,componentType,massType)) then
          massDistributionColdMode => hotHaloColdModeMassDistribution_                             %get             (self%hostNode         ,weightBy,weightIndex)
       else
          massDistributionColdMode => null()
       end if
       massDistributionHotMode     => self                            %nodeComponentHotHaloStandard%massDistribution(componentType,massType,weightBy,weightIndex)
       if (associated(massDistributionColdMode)) then
          allocate(kinematicsDistribution_)
          !![
	  <referenceConstruct object="kinematicsDistribution_" constructor="kinematicsDistributionLocal(alpha=1.0d0/sqrt(2.0d0))"/>
          !!]
          call massDistributionColdMode%setKinematicsDistribution(kinematicsDistribution_)
          !![
	  <objectDestructor name="kinematicsDistribution_"/>
          !!]
       end if
       if (.not.associated(massDistributionColdMode)) then
          if (.not.associated(massDistributionHotMode)) then
             massDistribution_ => null()
          else
             massDistribution_ => massDistributionHotMode
          end if
       else
          if (.not.associated(massDistributionHotMode)) then
             massDistribution_ => massDistributionColdMode
          else          
             allocate(massDistributionTotal          )
             allocate(massDistributionComponents     )
             allocate(massDistributionComponents%next)
             massDistributionComponents     %massDistribution_ => massDistributionHotMode
             massDistributionComponents%next%massDistribution_ => massDistributionColdMode
             !![
	     <referenceConstruct object="massDistributionTotal" constructor="massDistributionComposite(massDistributionComponents)"/>
	     <objectDestructor name="massDistributionHotMode" />
	     <objectDestructor name="massDistributionColdMode"/>
             !!]
             nullify(massDistributionComponents)
          end if
       end if
    class default
       call Error_Report('unexpected class'//{introspection:location})
    end select
    return
  end function Node_Component_Hot_Halo_Cold_Mode_Mass_Distribution

end module Node_Component_Hot_Halo_Cold_Mode
