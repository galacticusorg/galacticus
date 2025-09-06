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
  Implements a node operator class that implements reincorporation of outflowing gas into the \gls{cgm}.
  !!}

  use :: Hot_Halo_Outflows_Reincorporations, only : hotHaloOutflowReincorporationClass

  !![
  <nodeOperator name="nodeOperatorCGMOutflowReincorporation">
   <description>
    A node operator class that implements reincorporation of outflowing gas into the \gls{cgm}.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorCGMOutflowReincorporation
     !!{
     A node operator class that implements reincorporation of outflowing gas into the \gls{cgm}.
     !!}
     private
     class  (hotHaloOutflowReincorporationClass), pointer :: hotHaloOutflowReincorporation_ => null()
     logical                                              :: includeSatellites
   contains
     final     ::                          cgmOutflowReincorporationDestructor
     procedure :: autoHook              => cgmOutflowReincorporationAutoHook
     procedure :: nodePromote           => cgmOutflowReincorporationNodePromote
     procedure :: differentialEvolution => cgmOutflowReincorporationDifferentialEvolution
  end type nodeOperatorCGMOutflowReincorporation
  
  interface nodeOperatorCGMOutflowReincorporation
     !!{
     Constructors for the \refClass{nodeOperatorCGMOutflowReincorporation} node operator class.
     !!}
     module procedure cgmOutflowReincorporationConstructorParameters
     module procedure cgmOutflowReincorporationConstructorInternal
  end interface nodeOperatorCGMOutflowReincorporation
  
contains
  
  function cgmOutflowReincorporationConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorCGMOutflowReincorporation} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (nodeOperatorCGMOutflowReincorporation)                :: self
    type   (inputParameters                      ), intent(inout) :: parameters
    class  (hotHaloOutflowReincorporationClass   ), pointer       :: hotHaloOutflowReincorporation_
    logical                                                       :: includeSatellites

    !![
    <inputParameter>
      <name>includeSatellites</name>
      <defaultValue>.true.</defaultValue>
      <description>
	If true, perform reincorporation for all galaxies. Otherwise, reincorporation is performed only for central galaxies.
      </description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="hotHaloOutflowReincorporation" name="hotHaloOutflowReincorporation_" source="parameters"/>
    !!]
    self=nodeOperatorCGMOutflowReincorporation(includeSatellites,hotHaloOutflowReincorporation_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="hotHaloOutflowReincorporation_"/>
    !!]
    return
  end function cgmOutflowReincorporationConstructorParameters

  function cgmOutflowReincorporationConstructorInternal(includeSatellites,hotHaloOutflowReincorporation_) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorCGMOutflowReincorporation} node operator class.
    !!}
    implicit none
    type   (nodeOperatorCGMOutflowReincorporation)                        :: self
    class  (hotHaloOutflowReincorporationClass   ), intent(in   ), target :: hotHaloOutflowReincorporation_
    logical                                       , intent(in   )         :: includeSatellites
    !![
    <constructorAssign variables="includeSatellites, *hotHaloOutflowReincorporation_"/>
    !!]

    return
  end function cgmOutflowReincorporationConstructorInternal

  subroutine cgmOutflowReincorporationAutoHook(self)
    !!{
    Attach to various event hooks.
    !!}
    use :: Events_Hooks, only : satelliteMergerEvent, openMPThreadBindingAtLevel, dependencyDirectionAfter, dependencyRegEx
    implicit none
    class(nodeOperatorCGMOutflowReincorporation), intent(inout) :: self
    type (dependencyRegEx                      ), dimension(1)  :: dependenciesSatelliteMerger 

    dependenciesSatelliteMerger(1)=dependencyRegEx(dependencyDirectionAfter,'^remnantStructure:')
    call satelliteMergerEvent%attach(self,satelliteMerger,openMPThreadBindingAtLevel,label='cgmOutflowReincorporation',dependencies=dependenciesSatelliteMerger)
    return
  end subroutine cgmOutflowReincorporationAutoHook
  
  subroutine cgmOutflowReincorporationDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorCGMOutflowReincorporation} node operator class.
    !!}
    use :: Events_Hooks, only : satelliteMergerEvent
    implicit none
    type(nodeOperatorCGMOutflowReincorporation), intent(inout) :: self

    if (satelliteMergerEvent%isAttached(self,satelliteMerger)) call satelliteMergerEvent%detach(self,satelliteMerger)
    !![
    <objectDestructor name="self%hotHaloOutflowReincorporation_"/>
    !!]
    return
  end subroutine cgmOutflowReincorporationDestructor

  subroutine cgmOutflowReincorporationNodePromote(self,node)
    !!{
    Update the \gls{cgm} outflowed content of a node as a result of promotion.
    !!}
    use :: Galacticus_Nodes             , only : nodeComponentHotHalo
    use :: Abundances_Structure         , only : zeroAbundances
    use :: Chemical_Abundances_Structure, only : zeroChemicalAbundances 
    implicit none
    class(nodeOperatorCGMOutflowReincorporation), intent(inout) :: self
    type (treeNode                             ), intent(inout) :: node
    type (treeNode                             ), pointer       :: nodeParent
    class(nodeComponentHotHalo                 ), pointer       :: hotHaloParent, hotHalo

    nodeParent    => node      %parent
    hotHalo       => node      %hotHalo()
    hotHaloParent => nodeParent%hotHalo()
    ! If the parent node has a hot halo component, then add it to that of this node, and perform other changes needed prior to
    ! promotion.
    select type (hotHaloParent)
    type is (nodeComponentHotHalo)
       ! The parent has no hot halo component - nothing to do.
    class default
       ! If outflowed mass is non-positive, set mass and all related quantities to zero.
       if (hotHalo%outflowedMass() <= 0.0d0) then
          call        hotHalo%outflowedMassSet           (                 0.0d0)
          if (hotHalo%outflowedAngularMomentumIsSettable())                       &
               & call hotHalo%outflowedAngularMomentumSet(                 0.0d0)
          call        hotHalo%outflowedAbundancesSet     (        zeroAbundances)
          if (hotHalo%      outflowedChemicalsIsSettable())                       &
               &  call hotHalo%outflowedChemicalsSet     (zeroChemicalAbundances)
       end if
       ! Transfer CGM from the parent node.
       if (hotHalo%           outflowedMassIsSettable())                                          &
            & call hotHalo%           outflowedMassSet(                                           &
            &                                           hotHalo      % outflowedMass           () &
            &                                          +hotHaloParent% outflowedMass           () &
            &                                         )
       if (hotHalo%outflowedAngularMomentumIsSettable())                                          &
            & call hotHalo%outflowedAngularMomentumSet(                                           &
            &                                           hotHalo      % outflowedAngularMomentum() &
            &                                          +hotHaloParent% outflowedAngularMomentum() &
            &                                         )
       if (hotHalo%     outflowedAbundancesIsSettable())                                          &
            & call hotHalo%     outflowedAbundancesSet(                                           &
            &                                           hotHalo      % outflowedAbundances     () &
            &                                          +hotHaloParent% outflowedAbundances     () &
            &                                         )
       if (hotHalo%      outflowedChemicalsIsSettable())                                          &
            & call hotHalo%      outflowedChemicalsSet(                                           &
            &                                           hotHalo      % outflowedChemicals      () &
            &                                          +hotHaloParent% outflowedChemicals      () &
            &                                         )
    end select
    return
  end subroutine cgmOutflowReincorporationNodePromote

  subroutine satelliteMerger(self,node)
    !!{
    Remove any hot halo associated with {\normalfont \ttfamily node} before it merges with its host halo.
    !!}
    use :: Error                        , only : Error_Report
    use :: Abundances_Structure         , only : zeroAbundances
    use :: Chemical_Abundances_Structure, only : zeroChemicalAbundances 
    use :: Galacticus_Nodes             , only : nodeComponentHotHalo  , nodeComponentSpin, nodeComponentBasic
    implicit none
    class(*                   ), intent(inout)         :: self
    type (treeNode            ), intent(inout), target :: node
    type (treeNode            ), pointer               :: nodeHost
    class(nodeComponentBasic  ), pointer               :: basicHost
    class(nodeComponentHotHalo), pointer               :: hotHaloHost, hotHalo
    class(nodeComponentSpin   ), pointer               :: spinHost

    select type (self)
    class is (nodeOperatorCGMOutflowReincorporation)
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
          call        hotHaloHost%           outflowedMassSet(                                       &
               &                                              +hotHaloHost%outflowedMass          () &
               &                                              +    hotHalo%outflowedMass          () &
               &                                             )
          if (hotHaloHost%outflowedAngularMomentumIsSettable())                                       &
               & call hotHaloHost%outflowedAngularMomentumSet(                                        &
               &                                               hotHaloHost%outflowedAngularMomentum() &
               &                                              +hotHalo    %outflowedMass           () &
               &                                              *spinHost   %         angularMomentum() &
               &                                              /basicHost  %         mass           () &
               &                                             )
          if (hotHaloHost%     outflowedAbundancesIsSettable())                                       &
               & call hotHaloHost%     outflowedAbundancesSet(                                        &
               &                                               hotHaloHost%outflowedAbundances     () &
               &                                              +hotHalo    %outflowedAbundances     () &
               &                                             )
          if (hotHaloHost%      outflowedChemicalsIsSettable())                                       &
               & call hotHaloHost%      outflowedChemicalsSet(                                        &
               &                                               hotHaloHost%outflowedChemicals      () &
               &                                              +hotHalo    %outflowedChemicals      () &
               &                                             )
          call        hotHalo%               outflowedMassSet(                                        &
               &                                               0.0d0                                  &
               &                                             )
          if (hotHalo    %outflowedAngularMomentumIsSettable())                                       &
               & call hotHalo    %outflowedAngularMomentumSet(                                        &
               &                                               0.0d0                                  &
               &                                             )
          if (hotHalo    %     outflowedAbundancesIsSettable())                                       &
               & call hotHalo    %     outflowedAbundancesSet(                                        &
               &                                               zeroAbundances                         &
               &                                             )
          if (hotHalo    %      outflowedChemicalsIsSettable())                                       &
               & call hotHalo    %      outflowedChemicalsSet(                                        &
               &                                               zeroChemicalAbundances                 &
               &                                             )
       end select
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine satelliteMerger
  
  subroutine cgmOutflowReincorporationDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Perform accretion of gas into the \gls{cgm}.
    !!}
    use :: Abundances_Structure         , only : abundances          , operator(*)
    use :: Chemical_Abundances_Structure, only : chemicalAbundances  , operator(*)
    use :: Galacticus_Nodes             , only : nodeComponentHotHalo, propertyInactive
    implicit none
    class           (nodeOperatorCGMOutflowReincorporation), intent(inout), target  :: self
    type            (treeNode                             ), intent(inout), target  :: node
    logical                                                , intent(inout)          :: interrupt
    procedure       (interruptTask                        ), intent(inout), pointer :: functionInterrupt
    integer                                                , intent(in   )          :: propertyType
    class           (nodeComponentHotHalo                 )               , pointer :: hotHalo
    type            (abundances                           ), save                   :: rateAbundancesReturn
    type            (chemicalAbundances                   ), save                   :: rateChemicalsReturn
    !$omp threadprivate(rateAbundancesReturn,rateChemicalsReturn)
    double precision                                                                :: massOutflowed            , rateMassReturn, &
         &                                                                             rateAngularMomentumReturn
    !$GLC attributes unused :: interrupt, functionInterrupt

    ! Return immediately if inactive variables are requested.
    if (propertyInactive(propertyType)) return
    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    ! Don't reincorporate gas for satellites if requested.
    if (.not.self%includeSatellites.and.node%isSatellite()) return
    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    select type (hotHalo)
    type is (nodeComponentHotHalo)
       ! No hot halo exists - nothing to do.
    class default
       ! Move outflowed material back to the hot reservoir.
       massOutflowed =hotHalo                               %outflowedMass(    )
       rateMassReturn=self   %hotHaloOutflowReincorporation_%rate         (node)
       if (massOutflowed > 0.0d0) then
          rateAbundancesReturn     =+        rateMassReturn                                            &
               &                    *hotHalo%outflowedAbundances     ()                                &
               &                    /        massOutflowed
          rateAngularMomentumReturn=+hotHalo%outflowedAngularMomentum()*(rateMassReturn/massOutflowed)
          rateAbundancesReturn     =+hotHalo%outflowedAbundances     ()*(rateMassReturn/massOutflowed)
          rateChemicalsReturn      =+hotHalo%outflowedChemicals      ()*(rateMassReturn/massOutflowed)
          call    hotHalo%           outflowedMassRate(-           rateMassReturn)
          call    hotHalo%                    massRate(+           rateMassReturn)
          call    hotHalo%outflowedAngularMomentumRate(-rateAngularMomentumReturn)
          call    hotHalo%         angularMomentumRate( rateAngularMomentumReturn)
          call    hotHalo%     outflowedAbundancesRate(-     rateAbundancesReturn)
          call    hotHalo%              abundancesRate(      rateAbundancesReturn)
          if (hotHalo%outflowedChemicalsIsSettable()) then
             call hotHalo%      outflowedChemicalsRate(-      rateChemicalsReturn)
             call hotHalo%               chemicalsRate(       rateChemicalsReturn)
          end if
       end if
    end select
    return
  end subroutine cgmOutflowReincorporationDifferentialEvolution
  
