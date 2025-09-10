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
  Implements a node operator class that implements accretion of gas into the \gls{cgm}.
  !!}

  use :: Accretion_Halos, only : accretionHaloClass

  !![
  <nodeOperator name="nodeOperatorCGMAccretion">
   <description>
    A node operator class that implements accretion of gas into the \gls{cgm}.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorCGMAccretion
     !!{
     A node operator class that implements accretion of gas into the \gls{cgm}.
     !!}
     private
     class  (accretionHaloClass), pointer :: accretionHalo_       => null()
     logical                              :: allowNegativeCGMMass          , angularMomentumAlwaysGrows
     integer                              :: countChemicals
   contains
     final     ::                          cgmAccretionDestructor
     procedure :: autoHook              => cgmAccretionAutoHook
     procedure :: nodeInitialize        => cgmAccretionNodeInitialize
     procedure :: nodesMerge            => cgmAccretionNodesMerge
     procedure :: nodePromote           => cgmAccretionNodePromote
     procedure :: differentialEvolution => cgmAccretionDifferentialEvolution
  end type nodeOperatorCGMAccretion
  
  interface nodeOperatorCGMAccretion
     !!{
     Constructors for the \refClass{nodeOperatorCGMAccretion} node operator class.
     !!}
     module procedure cgmAccretionConstructorParameters
     module procedure cgmAccretionConstructorInternal
  end interface nodeOperatorCGMAccretion
  
contains
  
  function cgmAccretionConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorCGMAccretion} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (nodeOperatorCGMAccretion)                :: self
    type   (inputParameters         ), intent(inout) :: parameters
    class  (accretionHaloClass      ), pointer       :: accretionHalo_
    logical                                          :: allowNegativeCGMMass, angularMomentumAlwaysGrows
    
    !![
    <inputParameter>
      <name>allowNegativeCGMMass</name>
      <defaultValue>.true.</defaultValue>
      <description>
	If true, allow negative mass in the \gls{cgm}. If false, rates that would drive the \gls{cgm} mass to be negative are
	truncated to zero.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>angularMomentumAlwaysGrows</name>
      <defaultValue>.false.</defaultValue>
      <description>
	Specifies whether or not negative rates of accretion of angular momentum into the hot halo will be treated as positive for
        the purposes of computing the hot halo angular momentum.
      </description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="accretionHalo" name="accretionHalo_" source="parameters"/>
    !!]
    self=nodeOperatorCGMAccretion(allowNegativeCGMMass,angularMomentumAlwaysGrows,accretionHalo_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="accretionHalo_"/>
    !!]
    return
  end function cgmAccretionConstructorParameters

  function cgmAccretionConstructorInternal(allowNegativeCGMMass,angularMomentumAlwaysGrows,accretionHalo_) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorCGMAccretion} node operator class.
    !!}
    use :: Chemical_Abundances_Structure, only : Chemicals_Property_Count
    implicit none
    type   (nodeOperatorCGMAccretion)                        :: self
    class  (accretionHaloClass      ), intent(in   ), target :: accretionHalo_
    logical                          , intent(in   )         :: allowNegativeCGMMass, angularMomentumAlwaysGrows
    !![
    <constructorAssign variables="allowNegativeCGMMass, angularMomentumAlwaysGrows, *accretionHalo_"/>
    !!]

    self%countChemicals=Chemicals_Property_Count()
    return
  end function cgmAccretionConstructorInternal

  subroutine cgmAccretionAutoHook(self)
    !!{
    Attach to various event hooks.
    !!}
    use :: Events_Hooks, only : satelliteMergerEvent, openMPThreadBindingAtLevel, dependencyDirectionAfter, dependencyRegEx
    implicit none
    class(nodeOperatorCGMAccretion), intent(inout) :: self
    type (dependencyRegEx         ), dimension(1)  :: dependenciesSatelliteMerger 

    dependenciesSatelliteMerger(1)=dependencyRegEx(dependencyDirectionAfter,'^remnantStructure:')
    call satelliteMergerEvent%attach(self,satelliteMerger,openMPThreadBindingAtLevel,label='cgmAccretion',dependencies=dependenciesSatelliteMerger)
    return
  end subroutine cgmAccretionAutoHook
  
  subroutine cgmAccretionDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorCGMAccretion} node operator class.
    !!}
    use :: Events_Hooks, only : satelliteMergerEvent
    implicit none
    type(nodeOperatorCGMAccretion), intent(inout) :: self

    if (satelliteMergerEvent%isAttached(self,satelliteMerger)) call satelliteMergerEvent%detach(self,satelliteMerger)
    !![
    <objectDestructor name="self%accretionHalo_"/>
    !!]
    return
  end subroutine cgmAccretionDestructor

  subroutine cgmAccretionNodeInitialize(self,node)
    !!{
    Initialize the \gls{cgm} content of a node.
    !!}
    use :: Accretion_Halos     , only : accretionModeTotal       , accretionModeHot
    use :: Galacticus_Nodes    , only : nodeComponentHotHalo     , nodeComponentBasic, nodeComponentSpin, nodeEvent, &
         &                              nodeEventSubhaloPromotion
    implicit none
    class           (nodeOperatorCGMAccretion), intent(inout), target  :: self
    type            (treeNode                ), intent(inout), target  :: node
    class           (nodeComponentBasic      )               , pointer :: basic
    class           (nodeComponentHotHalo    )               , pointer :: hotHalo
    class           (nodeComponentSpin       )               , pointer :: spin
    class           (nodeEvent               )               , pointer :: event
    double precision                                                   :: massHotHalo    , massFailedHotHalo, &
         &                                                                angularMomentum

    ! If this is not a branch tip, then return immediately.
    if (associated(node%firstChild)) return
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
    ! Get the mass of hot gas accreted and the mass that failed to accrete.
    massHotHalo      =self%accretionHalo_%accretedMass      (node,accretionModeHot  )
    massFailedHotHalo=self%accretionHalo_%failedAccretedMass(node,accretionModeTotal)
    ! If either is non-zero, then create a hot halo component and add these masses to it.
    if (massHotHalo > 0.0d0 .or. massFailedHotHalo > 0.0d0) then
       hotHalo => node%hotHalo(autoCreate=.true.)
       basic   => node%basic  (                 )
       spin    => node%spin   (                 )
       ! Compute the initial angular momentum.
       angularMomentum=+spin %angularMomentum() &
            &          *      massHotHalo       &
            &          /basic%mass           ()
       ! Set the other initial properties.
       call        hotHalo%                massSet(                          massHotHalo                                 )
       call        hotHalo%      unaccretedMassSet(                          massFailedHotHalo                           )
       if (hotHalo%     angularMomentumIsSettable()) &
            & call hotHalo%     angularMomentumSet(                          angularMomentum                             )
       call        hotHalo%          abundancesSet(self%accretionHalo_%      accretedMassMetals   (node,accretionModeHot))
       if (hotHalo%unaccretedAbundancesIsSettable()) &
            & call hotHalo%unaccretedAbundancesSet(self%accretionHalo_%failedAccretedMassMetals   (node,accretionModeHot))
       if (hotHalo%           chemicalsIsSettable()) &
            & call hotHalo%           chemicalsSet(self%accretionHalo_%      accretedMassChemicals(node,accretionModeHot))
    end if
    return
  end subroutine cgmAccretionNodeInitialize
  
  subroutine cgmAccretionNodePromote(self,node)
    !!{
    Update the \gls{cgm} content of a node as a result of promotion.
    !!}
    use :: Abundances_Structure         , only : zeroAbundances
    use :: Chemical_Abundances_Structure, only : zeroChemicalAbundances
    use :: Galacticus_Nodes             , only : nodeComponentHotHalo, nodeComponentHotHaloStandard
    implicit none
    class(nodeOperatorCGMAccretion), intent(inout) :: self
    type (treeNode                ), intent(inout) :: node
    type (treeNode                ), pointer       :: nodeParent
    class(nodeComponentHotHalo    ), pointer       :: hotHaloParent, hotHalo

    nodeParent    => node      %parent
    hotHaloParent => nodeParent%hotHalo()
    ! If the parent node has a hot halo component, then add it to that of this node, and perform other changes needed prior to
    ! promotion.
    select type (hotHaloParent)
    type is (nodeComponentHotHalo)
       ! The parent has no hot halo component - nothing to do.
    class default
       hotHalo => node%hotHalo(autoCreate=.true.)
       ! If mass is non-positive, set mass and all related quantities to zero.
       if (hotHalo%         mass() <= 0.0d0) then
          call        hotHalo%         massSet           (                 0.0d0)
          if (hotHalo%         angularMomentumIsSettable())                       &
               & call hotHalo%         angularMomentumSet(                 0.0d0)
          call        hotHalo%         abundancesSet     (        zeroAbundances)
          if (hotHalo%               chemicalsIsSettable())                       &
               & call hotHalo%         chemicalsSet      (zeroChemicalAbundances)
       end if
       ! Transfer CGM from the parent node.
       call        hotHalo%          unaccretedMassSet(                                           &
            &                                          +hotHalo      %unaccretedMass           () &
            &                                          +hotHaloParent%unaccretedMass           () &
            &                                         )
       call        hotHalo%                    massSet(                                           &
            &                                          +hotHalo      %          mass           () &
            &                                          +hotHaloParent%          mass           () &
            &                                         )
       if (hotHalo%         angularMomentumIsSettable())                                          &
            & call hotHalo%         angularMomentumSet(                                           &
            &                                           hotHalo      %          angularMomentum() &
            &                                          +hotHaloParent%          angularMomentum() &
            &                                         )
       if (hotHalo%    unaccretedAbundancesIsSettable())                                          &
            & call hotHalo%    unaccretedAbundancesSet(                                           &
            &                                           hotHalo      %unaccretedAbundances     () &
            &                                          +hotHaloParent%unaccretedAbundances     () &
            &                                         )
       call        hotHalo% abundancesSet             (                                           &
            &                                          +hotHalo      %          abundances     () &
            &                                          +hotHaloParent%          abundances     () &
            &                                         )
       if (hotHalo%               chemicalsIsSettable())                                          &
            & call hotHalo%               chemicalsSet(                                           &
            &                                           hotHalo      %          chemicals      () &
            &                                          +hotHaloParent%          chemicals      () &
            &                                         )
    end select
    return
  end subroutine cgmAccretionNodePromote

  subroutine cgmAccretionNodesMerge(self,node)
    !!{
    Update the \gls{cgm} content of a node as a result of a merger.
    !!}
    use :: Abundances_Structure         , only : abundances        , zeroAbundances        , operator(*)
    use :: Chemical_Abundances_Structure, only : chemicalAbundances, zeroChemicalAbundances, operator(*)      , operator(>)
    use :: Accretion_Halos              , only : accretionModeHot  , accretionModeTotal
    use :: Galacticus_Nodes             , only : nodeComponentBasic, nodeComponentHotHalo  , nodeComponentSpin
    implicit none
    class           (nodeOperatorCGMAccretion), intent(inout) :: self
    type            (treeNode                ), intent(inout) :: node
    type            (treeNode                ), pointer       :: nodeParent
    class           (nodeComponentHotHalo    ), pointer       :: hotHaloParent          , hotHalo
    class           (nodeComponentBasic      ), pointer       :: basicParent            , basic
    class           (nodeComponentSpin       ), pointer       :: spinParent
    type            (abundances              ), save          :: massMetalsAccreted     , fractionMetalsAccreted   , &
         &                                                       massMetalsReaccreted
    type            (chemicalAbundances      ), save          :: massChemicalsAccreted  , fractionChemicalsAccreted, &
         &                                                       massChemicalsReaccreted
    !$omp threadprivate(massMetalsAccreted,fractionMetalsAccreted,massMetalsReaccreted,massChemicalsAccreted,fractionChemicalsAccreted,massChemicalsReaccreted)
    logical                                                   :: massTotalNonZero
    double precision                                          :: massAccreted           , massUnaccreted           , &
         &                                                       fractionAccreted       , massReaccreted           , &
         &                                                       angularMomentumAccreted, massAccretedHot

    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    select type (hotHalo)
    type is (nodeComponentHotHalo)
       ! No hot halo exists - nothing to do.
    class default
       ! Find the parent node and its hot halo component.
       nodeParent    => node      %parent
       hotHaloParent => nodeParent%hotHalo(autoCreate=.true.)
       ! Get the basic components.
       basic       => node      %basic()
       basicParent => nodeParent%basic()
       spinParent  => nodeParent%spin ()
       ! Any gas that failed to be accreted by this halo is always transferred to the parent.
       call    hotHaloParent%      unaccretedMassSet(hotHaloParent%unaccretedMass      ()+hotHalo%unaccretedMass      ())
       call    hotHalo      %      unaccretedMassSet(0.0d0                                                              )
       if (hotHaloParent%unaccretedAbundancesIsSettable()) then
          call hotHaloParent%unaccretedAbundancesSet(hotHaloParent%unaccretedAbundances()+hotHalo%unaccretedAbundances())
          call hotHalo      %unaccretedAbundancesSet(zeroAbundances                                                     )
       end if
       ! Since the parent node is undergoing mass growth through this merger we potentially return some of the unaccreted
       ! gas to the hot phase.
       !! First, find the masses of hot and failed mass the node would have if it formed instantaneously.
       massAccretedHot =self%accretionHalo_%      accretedMass(nodeParent,accretionModeHot  )
       massAccreted    =self%accretionHalo_%      accretedMass(nodeParent,accretionModeTotal)
       massUnaccreted  =self%accretionHalo_%failedAccretedMass(nodeParent,accretionModeTotal)
       massTotalNonZero=+massAccreted+massUnaccreted > 0.0d0
       !! Find the fraction of mass that would be successfully accreted.
       if (massAccretedHot > 0.0d0) then
          if (.not.massTotalNonZero) call Error_Report('mass of hot-mode gas accreted is non-zero, but total mass is zero'//{introspection:location})
          fractionAccreted=+  massAccretedHot &
               &           /(                 &
               &             +massAccreted    &
               &             +massUnaccreted  &
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
          if (hotHaloParent%unaccretedAbundancesIsSettable()) then
             !! First, find the metal mass the node would have if it formed instantaneously.
             massMetalsAccreted=self%accretionHalo_%accretedMassMetals(nodeParent,accretionModeHot)
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
             call hotHaloParent%unaccretedAbundancesSet(hotHaloParent%unaccretedAbundances()-massMetalsReaccreted)
             call hotHaloParent%          abundancesSet(hotHaloParent%          abundances()+massMetalsReaccreted)
          end if
          ! Compute the reaccreted angular momentum.
          if (hotHaloParent%angularMomentumIsSettable()) then
             angularMomentumAccreted=+            massReaccreted    &
                  &                  *spinParent %angularMomentum() &
                  &                  /basicParent%mass           ()
             call hotHaloParent%angularMomentumSet(hotHaloParent%angularMomentum()+angularMomentumAccreted)
          end if
       end if
       ! Compute the reaccreted chemicals.
       !! First, find the chemicals mass that would be successfully accreted.
       massChemicalsAccreted=self%accretionHalo_%accretedMassChemicals(nodeParent,accretionModeHot)
       !! Find the mass fraction of chemicals that would be successfully accreted.
       if (massChemicalsAccreted > zeroChemicalAbundances) then
          if (.not.massTotalNonZero) call Error_Report('mass of hot-mode chemicals accreted is non-zero, but total mass is zero'//{introspection:location})
          fractionChemicalsAccreted=+  massChemicalsAccreted &
               &                    /(                       &
               &                      +massAccreted          &
               &                      +massUnaccreted        &
               &                     )
          !! Find the change in the unaccreted mass.
          massChemicalsReaccreted=+hotHaloParent   %unaccretedMass() &
               &                  *fractionChemicalsAccreted         &
               &                  *basic           %          mass() &
               &                  /basicParent     %          mass()
          !! Reaccrete the chemicals.
          call hotHaloParent%chemicalsSet(hotHaloParent%chemicals()+massChemicalsReaccreted)
       end if
    end select
    return
  end subroutine cgmAccretionNodesMerge
  
  subroutine satelliteMerger(self,node)
    !!{
    Remove any hot halo associated with {\normalfont \ttfamily node} before it merges with its host halo.
    !!}
    use :: Abundances_Structure         , only : zeroAbundances
    use :: Chemical_Abundances_Structure, only : zeroChemicalAbundances
    use :: Error                        , only : Error_Report
    use :: Abundances_Structure         , only : zeroAbundances
    use :: Galacticus_Nodes             , only : nodeComponentHotHalo  , nodeComponentSpin, nodeComponentBasic
    implicit none
    class(*                   ), intent(inout)         :: self
    type (treeNode            ), intent(inout), target :: node
    type (treeNode            ), pointer               :: nodeHost
    class(nodeComponentBasic  ), pointer               :: basicHost
    class(nodeComponentHotHalo), pointer               :: hotHaloHost, hotHalo
    class(nodeComponentSpin   ), pointer               :: spinHost

    select type (self)
    class is (nodeOperatorCGMAccretion)
       ! Get the hot halo component.
       hotHalo => node%hotHalo()
       select type (hotHalo)
       type is (nodeComponentHotHalo)
          ! No hot halo exists - nothing to do.
       class default
          ! Find the node to merge with.
          nodeHost    => node    %mergesWith(                 )
          basicHost   => nodeHost%basic     (                 )
          hotHaloHost => nodeHost%hotHalo   (autoCreate=.true.)
          spinHost    => nodeHost%spin      (                 )
          ! Move the hot halo to the host.
          call        hotHaloHost%                    massSet(                                         &
               &                                              +hotHaloHost%          mass           () &
               &                                              +hotHalo    %          mass           () &
               &                                             )
          call        hotHaloHost%          unaccretedMassSet(                                         &
               &                                               hotHaloHost%unaccretedMass           () &
               &                                              +hotHalo    %unaccretedMass           () &
               &                                             )
          if (hotHaloHost%         angularMomentumIsSettable())                                        &
               & call hotHaloHost%         angularMomentumSet(                                         &
               &                                               hotHaloHost%          angularMomentum() &
               &                                              +hotHalo    %          mass           () &
               &                                              *spinHost   %          angularMomentum() &
               &                                              /basicHost  %          mass           () &
               &                                             )
          call        hotHaloHost%              abundancesSet(                                         &
               &                                              +hotHaloHost%          abundances     () &
               &                                              +hotHalo    %          abundances     () &
               &                                             )
          if (hotHaloHost%    unaccretedAbundancesIsSettable())                                        &
               & call hotHaloHost%    unaccretedAbundancesSet(                                         &
               &                                               hotHaloHost%unaccretedAbundances     () &
               &                                              +hotHalo    %unaccretedAbundances     () &
               &                                             )
          if (hotHaloHost%               chemicalsIsSettable())                                        &
               & call hotHaloHost%               chemicalsSet(                                         &
               &                                               hotHaloHost%          chemicals      () &
               &                                              +hotHalo    %          chemicals      () &
               &                                             )
          call        hotHalo    %                    massSet(                                         &
               &                                               0.0d0                                   &
               &                                             )
          if (hotHalo   %          angularMomentumIsSettable())                                        &
               & call hotHalo    %         angularMomentumSet(                                         &
               &                                               0.0d0                                   &
               &                                             )
          call        hotHalo    %              abundancesSet(                                         &
               &                                               zeroAbundances                          &
               &                                             )
          if (hotHalo        %unaccretedAbundancesIsSettable())                                        &
               & call hotHalo    %    unaccretedAbundancesSet(                                         &
               &                                               zeroAbundances                          &
               &                                             )
          if (hotHalo    %               chemicalsIsSettable())                                        &
               & call hotHalo    %               chemicalsSet(                                         &
               &                                               zeroChemicalAbundances                  &
               &                                             )
       end select
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine satelliteMerger
  
  subroutine cgmAccretionDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Perform accretion of gas into the \gls{cgm}.
    !!}
    use :: Abundances_Structure         , only : abundances
    use :: Chemical_Abundances_Structure, only : chemicalAbundances
    use :: Accretion_Halos              , only : accretionModeTotal  , accretionModeHot
    use :: Galacticus_Nodes             , only : nodeComponentHotHalo, nodeComponentSpin, nodeComponentBasic, propertyInactive
    implicit none
    class           (nodeOperatorCGMAccretion), intent(inout), target  :: self
    type            (treeNode                ), intent(inout), target  :: node
    logical                                   , intent(inout)          :: interrupt
    procedure       (interruptTask           ), intent(inout), pointer :: functionInterrupt
    integer                                   , intent(in   )          :: propertyType
    class           (nodeComponentBasic      )               , pointer :: basic
    class           (nodeComponentHotHalo    )               , pointer :: hotHalo
    class           (nodeComponentSpin       )               , pointer :: spin
    type            (abundances              ), save                   :: rateMetalsAccretion         , rateMetalsAccretionFailed
    type            (chemicalAbundances      ), save                   :: rateChemicalsAccretion
    !$omp threadprivate(rateMetalsAccretion,rateMetalsAccretionFailed,rateChemicalsAccretion)
    double precision                                                   :: rateMassAccretion           , rateMassAccretionFailed  , &
         &                                                                rateAngularMomentumAccretion
    
    ! Return immediately if inactive variables are requested.
    if (propertyInactive(propertyType)) return
    ! Get the components needed.
    basic   => node%basic  ()
    hotHalo => node%hotHalo()
    ! Find the rate of gas mass accretion onto the halo. We take all of the unaccreted gas, including any which is nominally cold
    ! mode, since we do not care what mode it should be in (since it is not actually accreted). It will be assigned to the
    ! relevant mode later when reaccreted. Negative accretion is allowed here as (for example) we may need to record the amount of
    ! negative accretion during periods where a halo is declining in mass, such that this loss of mass must be accounted for when
    ! the halo starts growing again before we allow gas to actually accrete into the hot component.
    rateMassAccretion               =   self%accretionHalo_%      accretionRate         (node,accretionModeHot  )
    rateMassAccretionFailed         =   self%accretionHalo_%failedAccretionRate         (node,accretionModeTotal)
    ! Get the rates of abundances accretion onto the halo.
    rateMetalsAccretion             =   self%accretionHalo_%      accretionRateMetals   (node,accretionModeHot  )
    rateMetalsAccretionFailed       =   self%accretionHalo_%failedAccretionRateMetals   (node,accretionModeHot  )
    ! Get the rates of chemicals accretion onto the halo.
    if (self%countChemicals > 0) &
         & rateChemicalsAccretion   =   self%accretionHalo_%      accretionRateChemicals(node,accretionModeHot  )
    ! Get the rate of angular momentum accretion onto the halo.
    if (basic%accretionRate() /= 0.0d0 .and. hotHalo%angularMomentumIsSettable()) then
       spin                         =>  node %spin                     ()
       rateAngularMomentumAccretion =  +spin %angularMomentumGrowthRate() &
            &                          *      rateMassAccretion           &
            &                          /basic%accretionRate            ()
       if (self%angularMomentumAlwaysGrows) &
            & rateAngularMomentumAccretion=abs(rateAngularMomentumAccretion)
    else
       rateAngularMomentumAccretion =  +0.0d0
    end if
    ! Apply accretion rates.
    if (      rateMassAccretion > 0.0d0 .or. hotHalo%mass() > 0.0d0 .or. self%allowNegativeCGMMass) &
         & call hotHalo%                massRate(      rateMassAccretion     ,interrupt,functionInterrupt)
    if (rateMassAccretionFailed > 0.0d0 .or. hotHalo%mass() > 0.0d0 .or. self%allowNegativeCGMMass) &
         & call hotHalo%      unaccretedMassRate(rateMassAccretionFailed     ,interrupt,functionInterrupt)
    call        hotHalo%          abundancesRate(rateMetalsAccretion         ,interrupt,functionInterrupt)
    call        hotHalo%unaccretedAbundancesRate(rateMetalsAccretionFailed   ,interrupt,functionInterrupt)
    if (self%countChemicals > 0) &
         & call hotHalo%           chemicalsRate(rateChemicalsAccretion      ,interrupt,functionInterrupt)
    call        hotHalo%     angularMomentumRate(rateAngularMomentumAccretion,interrupt,functionInterrupt)
    return
  end subroutine cgmAccretionDifferentialEvolution
  
