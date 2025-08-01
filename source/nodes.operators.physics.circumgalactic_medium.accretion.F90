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
     class(accretionHaloClass), pointer :: accretionHalo_ => null()
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
    type (nodeOperatorCGMAccretion)                :: self
    type (inputParameters         ), intent(inout) :: parameters
    class(accretionHaloClass      ), pointer       :: accretionHalo_

    !![
    <objectBuilder class="accretionHalo" name="accretionHalo_" source="parameters"/>
    !!]
    self=nodeOperatorCGMAccretion(accretionHalo_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="accretionHalo_"/>
    !!]
    return
  end function cgmAccretionConstructorParameters

  function cgmAccretionConstructorInternal(accretionHalo_) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorCGMAccretion} node operator class.
    !!}
    implicit none
    type (nodeOperatorCGMAccretion)                        :: self
    class(accretionHaloClass      ), intent(in   ), target :: accretionHalo_
    !![
    <constructorAssign variables="*accretionHalo_"/>
    !!]

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
    use :: Abundances_Structure, only : zeroAbundances
    use :: Accretion_Halos     , only : accretionModeTotal
    use :: Galacticus_Nodes    , only : nodeComponentHotHalo, nodeEvent, nodeEventSubhaloPromotion
    implicit none
    class           (nodeOperatorCGMAccretion), intent(inout), target  :: self
    type            (treeNode                ), intent(inout), target  :: node
    class           (nodeComponentHotHalo    )               , pointer :: hotHalo
    class           (nodeEvent               )               , pointer :: event
    double precision                                                   :: massHotHalo, massFailedHotHalo

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
    massHotHalo      =self%accretionHalo_%accretedMass      (node,accretionModeTotal)
    massFailedHotHalo=self%accretionHalo_%failedAccretedMass(node,accretionModeTotal)
    ! If either is non-zero, then create a hot halo component and add these masses to it.
    if (massHotHalo > 0.0d0 .or. massFailedHotHalo > 0.0d0) then
       hotHalo => node%hotHalo(autoCreate=.true.)
       call hotHalo%          massSet(      massHotHalo)
       call hotHalo%unaccretedMassSet(massFailedHotHalo)
       call hotHalo%    abundancesSet(   zeroAbundances)
    end if
    return
  end subroutine cgmAccretionNodeInitialize
  
  subroutine cgmAccretionNodePromote(self,node)
    !!{
    Update the \gls{cgm} content of a node as a result of promotion.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentHotHalo
    implicit none
    class(nodeOperatorCGMAccretion), intent(inout) :: self
    type (treeNode                ), intent(inout) :: node
    type (treeNode                ), pointer       :: nodeParent
    class(nodeComponentHotHalo    ), pointer       :: hotHaloParent, hotHalo

    nodeParent    => node      %parent
    hotHalo       => node      %hotHalo()
    hotHaloParent => nodeParent%hotHalo()
    ! If the parent node has a hot halo component, then add it to that of this node, and perform other changes needed prior to
    ! promotion.
    select type (hotHaloParent)
    type is (nodeComponentHotHalo)
       ! The parent has no hot halo component - nothing to do.
    class default
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
    return
  end subroutine cgmAccretionNodePromote

  subroutine cgmAccretionNodesMerge(self,node)
    !!{
    Update the \gls{cgm} content of a node as a result of a merger.
    !!}
    use :: Abundances_Structure, only : abundances        , operator(*)         , zeroAbundances
    use :: Accretion_Halos     , only : accretionModeHot  , accretionModeTotal
    use :: Galacticus_Nodes    , only : nodeComponentBasic, nodeComponentHotHalo
    implicit none
    class           (nodeOperatorCGMAccretion), intent(inout) :: self
    type            (treeNode                ), intent(inout) :: node
    type            (treeNode                ), pointer       :: nodeParent
    class           (nodeComponentHotHalo    ), pointer       :: hotHaloParent       , hotHalo
    class           (nodeComponentBasic      ), pointer       :: basic               , basicParent
    type            (abundances              ), save          :: massMetalsAccreted  , fractionMetalsAccreted, &
         &                                                       massMetalsReaccreted
    !$omp threadprivate(massMetalsAccreted,fractionMetalsAccreted,massMetalsReaccreted)
    double precision                                          :: massAccreted        , massUnaccreted        , &
         &                                                       fractionAccreted    , massReaccreted

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
       massAccreted  =self%accretionHalo_%accretedMass      (nodeParent,accretionModeTotal)
       massUnaccreted=self%accretionHalo_%failedAccretedMass(nodeParent,accretionModeTotal)
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
       call hotHaloParent%abundancesSet(hotHaloParent%abundances()+massMetalsReaccreted)
    end select
    return
  end subroutine cgmAccretionNodesMerge
  
  subroutine satelliteMerger(self,node)
    !!{
    Remove any hot halo associated with {\normalfont \ttfamily node} before it merges with its host halo.
    !!}
    use :: Error               , only : Error_Report
    use :: Abundances_Structure, only : zeroAbundances
    use :: Galacticus_Nodes    , only : nodeComponentHotHalo
    implicit none
    class(*                   ), intent(inout)         :: self
    type (treeNode            ), intent(inout), target :: node
    type (treeNode            ), pointer               :: nodeHost
    class(nodeComponentHotHalo), pointer               :: hotHaloHost, hotHalo

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
          hotHaloHost => nodeHost%hotHalo   (autoCreate=.true.)
          ! Move the hot halo to the host.
          call hotHaloHost%      massSet(                          &
               &                         +hotHaloHost%mass      () &
               &                         +hotHalo    %mass      () &
               &                        )
          call hotHaloHost%abundancesSet(                          &
               &                         +hotHaloHost%abundances() &
               &                         +hotHalo    %abundances() &
               &                        )
          call hotHalo    %      massSet(                          &
               &                         +0.0d0                    &
               &                        )
          call hotHalo    %abundancesSet(                          &
               &                          zeroAbundances           &
               &                        )
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
    use :: Accretion_Halos , only : accretionModeTotal
    use :: Galacticus_Nodes, only : nodeComponentHotHalo, propertyInactive
    implicit none
    class           (nodeOperatorCGMAccretion), intent(inout), target  :: self
    type            (treeNode                ), intent(inout), target  :: node
    logical                                   , intent(inout)          :: interrupt
    procedure       (interruptTask           ), intent(inout), pointer :: functionInterrupt
    integer                                   , intent(in   )          :: propertyType
    class           (nodeComponentHotHalo    )               , pointer :: hotHalo
    double precision                                                   :: massAccretionRate, failedMassAccretionRate
    
    ! Return immediately if inactive variables are requested.
    if (propertyInactive(propertyType)) return
    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    ! Find the rate of gas mass accretion onto the halo.
    massAccretionRate      =self%accretionHalo_%accretionRate      (node,accretionModeTotal)
    failedMassAccretionRate=self%accretionHalo_%failedAccretionRate(node,accretionModeTotal)
    ! Apply accretion rates.
    if (      massAccretionRate > 0.0d0 .or. hotHalo%mass() > 0.0d0) &
         & call hotHalo%          massRate(      massAccretionRate,interrupt,functionInterrupt)
    if (failedMassAccretionRate > 0.0d0 .or. hotHalo%mass() > 0.0d0) &
         & call hotHalo%unaccretedMassRate(failedMassAccretionRate,interrupt,functionInterrupt)
    return
  end subroutine cgmAccretionDifferentialEvolution
  
