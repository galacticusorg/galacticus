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
  Implements a node operator class that triggers merging of satellites based on a merging time.
  !!}

  use :: Dark_Matter_Halos_Mass_Loss_Rates, only : darkMatterHaloMassLossRateClass

  !![
  <nodeOperator name="nodeOperatorSatelliteMassLoss">
   <description>A node operator class that triggers merging of satellites based on a merging time.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorSatelliteMassLoss
     !!{
     A node operator class that triggers merging of satellites based on a merging time.
     !!}
     private
     class  (darkMatterHaloMassLossRateClass), pointer :: darkMatterHaloMassLossRate_ => null()
     logical                                           :: massBoundIsInactive
   contains
     !![
     <methods>
       <method method="massBoundSet" description="Set the time of merging for a satellite node."/>
     </methods>
     !!]
     final     ::                          satelliteMassLossDestructor
     procedure :: autoHook              => satelliteMassLossAutoHook
     procedure :: nodeInitialize        => satelliteMassLossNodeInitialize
     procedure :: nodePromote           => satelliteMassLossNodePromote
     procedure :: nodesMerge            => satelliteMassLossNodesMerge
     procedure :: massBoundSet          => satelliteMassLossMassBoundSet
     procedure :: differentialEvolution => satelliteMassLossDifferentialEvolution
  end type nodeOperatorSatelliteMassLoss
  
  interface nodeOperatorSatelliteMassLoss
     !!{
     Constructors for the \refClass{nodeOperatorSatelliteMassLoss} node operator class.
     !!}
     module procedure satelliteMassLossConstructorParameters
     module procedure satelliteMassLossConstructorInternal
  end interface nodeOperatorSatelliteMassLoss

contains

  function satelliteMassLossConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorSatelliteMassLoss} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (nodeOperatorSatelliteMassLoss  )                :: self
    type   (inputParameters                ), intent(inout) :: parameters
    class  (darkMatterHaloMassLossRateClass), pointer       :: darkMatterHaloMassLossRate_
    logical                                                 :: massBoundIsInactive
    
    !![
    <inputParameter>
      <name>massBoundIsInactive</name>
      <defaultValue>.false.</defaultValue>
      <description>Specifies whether or not the bound mass of the satellite component is inactive (i.e. does not appear in any ODE being solved).</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="darkMatterHaloMassLossRate" name="darkMatterHaloMassLossRate_" source="parameters"/>
    !!]
    self=nodeOperatorSatelliteMassLoss(massBoundIsInactive,darkMatterHaloMassLossRate_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloMassLossRate_"/>
    !!]
    return
  end function satelliteMassLossConstructorParameters

  function satelliteMassLossConstructorInternal(massBoundIsInactive,darkMatterHaloMassLossRate_) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorSatelliteMassLoss} node operator class.
    !!}
    implicit none
    type   (nodeOperatorSatelliteMassLoss  )                        :: self
    class  (darkMatterHaloMassLossRateClass), intent(in   ), target :: darkMatterHaloMassLossRate_
    logical                                 , intent(in   )         :: massBoundIsInactive
    !![
    <constructorAssign variables="massBoundIsInactive, *darkMatterHaloMassLossRate_"/>
    !!]
    
    return
  end function satelliteMassLossConstructorInternal

  subroutine satelliteMassLossAutoHook(self)
    !!{
    Attach to various event hooks.
    !!}
    use :: Events_Hooks, only : subhaloPromotionEvent, openMPThreadBindingAtLevel
    implicit none
    class(nodeOperatorSatelliteMassLoss), intent(inout) :: self
    
    call subhaloPromotionEvent%attach(self,subhaloPromotion,openMPThreadBindingAtLevel,label='satelliteMassLoss')
    return
  end subroutine satelliteMassLossAutoHook

  subroutine satelliteMassLossDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorSatelliteMassLoss} node operator class.
    !!}
    use :: Events_Hooks, only : subhaloPromotionEvent
    implicit none
    type(nodeOperatorSatelliteMassLoss), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloMassLossRate_"/>
    !!]
    if (subhaloPromotionEvent%isAttached(self,subhaloPromotion)) call subhaloPromotionEvent%detach(self,subhaloPromotion)
    return
  end subroutine satelliteMassLossDestructor

  subroutine satelliteMassLossNodeInitialize(self,node)
    !!{
    Initialize bound mass of any initial satellites in a tree.
    !!}
    implicit none
    class(nodeOperatorSatelliteMassLoss), intent(inout), target :: self
    type (treeNode                     ), intent(inout), target :: node

    if     (                                               &
         &                     node%isSatellite        ()  &
         &  .or.                                           &
         &   (                                             &
         &     .not.           node%isPrimaryProgenitor()  &
         &    .and.                                        &
         &          associated(node%parent               ) &
         &   )                                             &
         & ) call self%massBoundSet(node)
    return
  end subroutine satelliteMassLossNodeInitialize

  subroutine subhaloPromotion(self,node,nodePromotion)
    !!{
    Handle cases where a subhalo is promoted to be an isolated host.    
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSatellite, nodeComponentBasic
    implicit none
    class(*                     ), intent(inout)          :: self
    type (treeNode              ), intent(inout), pointer :: node          , nodePromotion
    class(nodeComponentSatellite)               , pointer :: satellite
    class(nodeComponentBasic    )               , pointer :: basicPromotion

    select type (self)
    class is (nodeOperatorSatelliteMassLoss)
       satellite => node%satellite()
       if (satellite%boundMassIsSettable()) then
          basicPromotion => nodePromotion%basic()
          call satellite%boundMassSet(basicPromotion%mass())
       end if
    end select
    return
  end subroutine subhaloPromotion

  subroutine satelliteMassLossNodePromote(self,node)
    !!{
    Act on promotion of a node. Set the bound mass to the basic mass of the parent.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSatellite, nodeComponentBasic
    implicit none
    class(nodeOperatorSatelliteMassLoss), intent(inout) :: self
    type (treeNode                     ), intent(inout) :: node
    class(nodeComponentSatellite       ), pointer       :: satellite
    class(nodeComponentBasic           ), pointer       :: basicParent

    satellite => node%satellite()
    if (satellite%boundMassIsSettable()) then
       basicParent => node%parent%basic()
       call satellite%boundMassSet(basicParent%mass())
    end if
    return
  end subroutine satelliteMassLossNodePromote

  subroutine satelliteMassLossNodesMerge(self,node)
    !!{
    Act on node merging tree.
    !!}
    implicit none
    class(nodeOperatorSatelliteMassLoss), intent(inout) :: self
    type (treeNode                     ), intent(inout) :: node

    call self%massBoundSet(node)
    return
  end subroutine satelliteMassLossNodesMerge

  subroutine satelliteMassLossMassBoundSet(self,node)
    !!{
    Set the time of merging for the given {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentSatellite
    use :: Kepler_Orbits   , only : keplerOrbit       , keplerOrbitMasses        , keplerOrbitRadius            , keplerOrbitTheta, &
         &                          keplerOrbitPhi    , keplerOrbitVelocityRadial, keplerOrbitVelocityTangential
    implicit none
    class(nodeOperatorSatelliteMassLoss), intent(inout) :: self
    type (treeNode                     ), intent(inout) :: node
    class(nodeComponentBasic           ), pointer       :: basic
    class(nodeComponentSatellite       ), pointer       :: satellite

    basic     => node%basic    ()
    satellite => node%satellite()
    if (satellite%boundMassIsSettable()) &
         & call satellite%boundMassSet(basic%mass())
    return
  end subroutine satelliteMassLossMassBoundSet

  subroutine satelliteMassLossDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Set scales for differential evolution.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSatellite, nodeComponentBasic, propertyEvaluate
    implicit none
    class           (nodeOperatorSatelliteMassLoss), intent(inout), target  :: self
    type            (treeNode                     ), intent(inout), target  :: node
    logical                                        , intent(inout)          :: interrupt
    procedure       (interruptTask                ), intent(inout), pointer :: functionInterrupt
    integer                                        , intent(in   )          :: propertyType
    class           (nodeComponentBasic           )               , pointer :: basic
    class           (nodeComponentSatellite       )               , pointer :: satellite
    double precision                                                        :: massRate
    !$GLC attributes unused :: interrupt, functionInterrupt, propertyType

    if (propertyEvaluate(propertyType,self%massBoundIsInactive)) then
       satellite => node%satellite()
       if (node%isSatellite()) then
          massRate =  self %darkMatterHaloMassLossRate_%rate         (node)
       else
          basic    => node                             %basic        (    )
          massRate =  basic                            %accretionRate(    )
       end if
       call satellite%boundMassRate(massRate)
    end if
    return
  end subroutine satelliteMassLossDifferentialEvolution
  
