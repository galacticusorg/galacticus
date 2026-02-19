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
     class(hotHaloOutflowReincorporationClass), pointer :: hotHaloOutflowReincorporation_ => null()
   contains
     final     ::                          cgmOutflowReincorporationDestructor
     procedure :: autoHook              => cgmOutflowReincorporationAutoHook
     procedure :: nodesMerge            => cgmOutflowReincorporationNodesMerge
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
    type (nodeOperatorCGMOutflowReincorporation)                :: self
    type (inputParameters                      ), intent(inout) :: parameters
    class(hotHaloOutflowReincorporationClass   ), pointer       :: hotHaloOutflowReincorporation_

    !![
    <objectBuilder class="hotHaloOutflowReincorporation" name="hotHaloOutflowReincorporation_" source="parameters"/>
    !!]
    self=nodeOperatorCGMOutflowReincorporation(hotHaloOutflowReincorporation_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="hotHaloOutflowReincorporation_"/>
    !!]
    return
  end function cgmOutflowReincorporationConstructorParameters

  function cgmOutflowReincorporationConstructorInternal(hotHaloOutflowReincorporation_) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorCGMOutflowReincorporation} node operator class.
    !!}
    implicit none
    type (nodeOperatorCGMOutflowReincorporation)                        :: self
    class(hotHaloOutflowReincorporationClass   ), intent(in   ), target :: hotHaloOutflowReincorporation_
    !![
    <constructorAssign variables="*hotHaloOutflowReincorporation_"/>
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
    use :: Galacticus_Nodes, only : nodeComponentHotHalo
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
       call hotHalo%outflowedMassSet      (                                     &
            &                              +hotHalo      %outflowedMass      () &
            &                              +hotHaloParent%outflowedMass      () &
            &                             )
       call hotHalo%outflowedAbundancesSet(                                     &
            &                              +hotHalo      %outflowedAbundances() &
            &                              +hotHaloParent%outflowedAbundances() &
            &                             )
    end select
    return
  end subroutine cgmOutflowReincorporationNodePromote

  subroutine cgmOutflowReincorporationNodesMerge(self,node)
    !!{
    Update the \gls{cgm} content of a node as a result of a merger.
    !!}
    use :: Abundances_Structure, only : zeroAbundances
    use :: Galacticus_Nodes    , only : nodeComponentHotHalo
    implicit none
    class(nodeOperatorCGMOutflowReincorporation), intent(inout) :: self
    type (treeNode                             ), intent(inout) :: node
    type (treeNode                             ), pointer       :: nodeParent
    class(nodeComponentHotHalo                 ), pointer       :: hotHaloParent, hotHalo

    ! Get the hot halo component.
    hotHalo       => node%hotHalo()
    select type (hotHalo)
    type is (nodeComponentHotHalo)
       ! No hot halo in the merging node - nothing to do.
    class default
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
  end subroutine cgmOutflowReincorporationNodesMerge
  
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
    class is (nodeOperatorCGMOutflowReincorporation)
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
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine satelliteMerger
  
  subroutine cgmOutflowReincorporationDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Perform accretion of gas into the \gls{cgm}.
    !!}
    use :: Abundances_Structure, only : abundances          , operator(*)     , zeroAbundances
    use :: Galacticus_Nodes    , only : nodeComponentHotHalo, propertyInactive
    implicit none
    class           (nodeOperatorCGMOutflowReincorporation), intent(inout), target  :: self
    type            (treeNode                             ), intent(inout), target  :: node
    logical                                                , intent(inout)          :: interrupt
    procedure       (interruptTask                        ), intent(inout), pointer :: functionInterrupt
    integer                                                , intent(in   )          :: propertyType
    class           (nodeComponentHotHalo                 )               , pointer :: hotHalo
    type            (abundances                           ), save                   :: abundancesReturnRate
    !$omp threadprivate(abundancesReturnRate)
    double precision                                                                :: outflowReturnRate
    !$GLC attributes unused :: interrupt, interruptProcedure

    ! Return immediately if inactive variables are requested.
    if (propertyInactive(propertyType)) return
    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    ! Don't reincorporate gas for satellites - we don't want it to be able to re-infall back onto the satellite.
    if (node%isSatellite()) return
    ! Get the hot halo component.
    hotHalo => node%hotHalo()
    select type (hotHalo)
    type is (nodeComponentHotHalo)
       ! No hot halo exists - nothing to do.
    class default
       ! Move outflowed material back to the hot reservoir.
       outflowReturnRate=self%hotHaloOutflowReincorporation_%rate(node)
       if (hotHalo%outflowedMass() > 0.0d0) then
          abundancesReturnRate           =  +        outflowReturnRate     &
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
  end subroutine cgmOutflowReincorporationDifferentialEvolution
  
