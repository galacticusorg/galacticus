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
  Implements a node operator class that tracks the mass and metals arriving in the \gls{cgm} from outflows.
  !!}

  !![
  <nodeOperator name="nodeOperatorTrackOutflowedMass">
   <description>A node operator class that tracks the mass and metals arriving in the \gls{cgm} from outflows.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorTrackOutflowedMass
     !!{
     A node operator class that tracks the mass and metals arriving in the \gls{cgm} from outflows.
     !!}
     private
     integer :: massOutflowedID, massMetalsOutflowedID
   contains
     final     ::                                trackOutflowedMassDestructor
     procedure :: differentialEvolutionScales => trackOutflowedMassDifferentialEvolutionScales
     procedure :: autoHook                    => trackOutflowedMassAutoHook
  end type nodeOperatorTrackOutflowedMass
  
  interface nodeOperatorTrackOutflowedMass
     !!{
     Constructors for the \refClass{nodeOperatorTrackOutflowedMass} node operator class.
     !!}
     module procedure trackOutflowedMassConstructorParameters
     module procedure trackOutflowedMassConstructorInternal
  end interface nodeOperatorTrackOutflowedMass
  
contains

  function trackOutflowedMassConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorTrackOutflowedMass} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodeOperatorTrackOutflowedMass)                :: self
    type(inputParameters               ), intent(inout) :: parameters

    self=nodeOperatorTrackOutflowedMass()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function trackOutflowedMassConstructorParameters

  function trackOutflowedMassConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorTrackOutflowedMass} node operator class.
    !!}
    implicit none
    type(nodeOperatorTrackOutflowedMass) :: self

    !![
    <addMetaProperty component="hotHalo" name="massOutflowed"       id="self%massOutflowedID"       isEvolvable="yes" isCreator="yes"/>
    <addMetaProperty component="hotHalo" name="massMetalsOutflowed" id="self%massMetalsOutflowedID" isEvolvable="yes" isCreator="yes"/>
    !!]
    return
  end function trackOutflowedMassConstructorInternal

  subroutine trackOutflowedMassAutoHook(self)
    !!{
    Attach to various event hooks.
    !!}
    use :: Events_Hooks, only : openMPThreadBindingAtLevel, hotHaloMassEjectionEvent, hotHaloMassInflowEvent, hotHaloMassReincorporationEvent
    implicit none
    class(nodeOperatorTrackOutflowedMass), intent(inout) :: self

    call hotHaloMassEjectionEvent       %attach(self,trackOutflowedMassHotHaloMassRemoval        ,openMPThreadBindingAtLevel,label='trackOutflowedMass')
    call hotHaloMassInflowEvent         %attach(self,trackOutflowedMassHotHaloMassRemoval        ,openMPThreadBindingAtLevel,label='trackOutflowedMass')
    call hotHaloMassReincorporationEvent%attach(self,trackOutflowedMassHotHaloMassReincorporation,openMPThreadBindingAtLevel,label='trackOutflowedMass')
    return
  end subroutine trackOutflowedMassAutoHook

  subroutine trackOutflowedMassDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorTrackOutflowedMass} node operator class.
    !!}
    use :: Events_Hooks, only : hotHaloMassEjectionEvent, hotHaloMassInflowEvent, hotHaloMassReincorporationEvent
    implicit none
    type(nodeOperatorTrackOutflowedMass), intent(inout) :: self

    if (hotHaloMassEjectionEvent       %isAttached(self,trackOutflowedMassHotHaloMassRemoval        )) call hotHaloMassEjectionEvent       %detach(self,trackOutflowedMassHotHaloMassRemoval        )
    if (hotHaloMassInflowEvent         %isAttached(self,trackOutflowedMassHotHaloMassRemoval        )) call hotHaloMassInflowEvent         %detach(self,trackOutflowedMassHotHaloMassRemoval        )
    if (hotHaloMassReincorporationEvent%isAttached(self,trackOutflowedMassHotHaloMassReincorporation)) call hotHaloMassReincorporationEvent%detach(self,trackOutflowedMassHotHaloMassReincorporation)
    return
  end subroutine trackOutflowedMassDestructor

  subroutine trackOutflowedMassDifferentialEvolutionScales(self,node)
    !!{
    Set absolute ODE solver scale for the outflowed mass in the \gls{cgm}.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentHotHalo, nodeComponentBasic
    use :: Numerical_Constants_Astronomical, only : gigaYear            , massSolar
    implicit none
    class           (nodeOperatorTrackOutflowedMass), intent(inout) :: self
    type            (treeNode                      ), intent(inout) :: node
    double precision                                , parameter     :: scaleMassRelative=1.0d-3, scaleMetallicityRelative=1.0d-3
    class           (nodeComponentHotHalo          ), pointer       :: hotHalo
    class           (nodeComponentBasic            ), pointer       :: basic
    double precision                                                :: massVirial
   
    hotHalo => node%hotHalo()
    select type (hotHalo)
    type is (nodeComponentHotHalo)
       ! Hot halo does not exist - nothing to do here.
    class default
       basic      => node %basic()
       massVirial =  basic%mass ()
       call hotHalo%floatRank0MetaPropertyScale(self%massOutflowedID      ,massVirial*scaleMassRelative                         )
       call hotHalo%floatRank0MetaPropertyScale(self%massMetalsOutflowedID,massVirial*scaleMassRelative*scaleMetallicityRelative)
    end select
    return
  end subroutine trackOutflowedMassDifferentialEvolutionScales

  subroutine trackOutflowedMassHotHaloMassRemoval(self,hotHalo,massRate)
    !!{
    Respond to mass removal from the hot halo component.    
    !!}
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : nodeComponentHotHalo
    implicit none
    class           (*                   ), intent(inout) :: self
    class           (nodeComponentHotHalo), intent(inout) :: hotHalo
    double precision                      , intent(in   ) :: massRate

    select type (self)
    class is (nodeOperatorTrackOutflowedMass)
       call hotHalo%floatRank0MetaPropertyRate(self%massOutflowedID      ,-hotHalo%floatRank0MetaPropertyGet(self%massOutflowedID      )*massRate/hotHalo%mass())
       call hotHalo%floatRank0MetaPropertyRate(self%massMetalsOutflowedID,-hotHalo%floatRank0MetaPropertyGet(self%massMetalsOutflowedID)*massRate/hotHalo%mass())
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine trackOutflowedMassHotHaloMassRemoval

  subroutine trackOutflowedMassHotHaloMassReincorporation(self,hotHalo,rateMassReturn,rateAbundancesReturn)
    !!{
    Respond to mass removal from the hot halo component.    
    !!}
    use :: Error               , only : Error_Report
    use :: Galacticus_Nodes    , only : nodeComponentHotHalo
    use :: Abundances_Structure, only : abundances
    implicit none
    class           (*                    ), intent(inout) :: self
    class           (nodeComponentHotHalo ), intent(inout) :: hotHalo
    double precision                       , intent(in   ) :: rateMassReturn
    type            (abundances           ), intent(in   ) :: rateAbundancesReturn

    select type (self)
    class is (nodeOperatorTrackOutflowedMass)
       call hotHalo%floatRank0MetaPropertyRate(self%massOutflowedID      ,rateMassReturn                    )
       call hotHalo%floatRank0MetaPropertyRate(self%massMetalsOutflowedID,rateAbundancesReturn%metallicity())
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine trackOutflowedMassHotHaloMassReincorporation
