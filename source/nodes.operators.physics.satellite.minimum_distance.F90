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
  Implements a node operator class that tracks the minimum distance from the center that a satellite has ever reached in its
  current host halo.
  !!}

  !![
  <enumeration>
   <name>relativeTo</name>
   <description>Options for which host halo to compute the minimum distance to.</description>
   <encodeFunction>yes</encodeFunction>
   <validator>yes</validator>
   <visibility>private</visibility>
   <entry label="immediateHost"/>
   <entry label="isolatedHost" />
  </enumeration>
  !!]

  !![
  <nodeOperator name="nodeOperatorSatelliteMinimumDistance">
    <description>
      A node operator class that tracks the minimum distance from the center that a satellite has ever reached in its current host
      halo.
    </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorSatelliteMinimumDistance
     !!{
     A node operator class that tracks the minimum distance from the center that a satellite has ever reached in its current host
     halo.
     !!}
     private
     integer                            :: satelliteDistanceMinimumID, isolatedHostID
     type   (enumerationRelativeToType) :: relativeTo
   contains
     !![
     <methods>
       <method method="distanceRelative" description="Compute the distance to the center of the relevant host halo."/>
     </methods>
     !!]
     final     ::                              satelliteMinimumDistanceDestructor
     procedure :: differentialEvolutionPost => satelliteMinimumDistanceDifferentialEvolutionPost
     procedure :: nodesMerge                => satelliteMinimumDistanceNodesMerge
     procedure :: autoHook                  => satelliteMinimumDistanceAutoHook
     procedure :: distanceRelative          => satelliteMinimumDistanceDistanceRelative
  end type nodeOperatorSatelliteMinimumDistance
  
  interface nodeOperatorSatelliteMinimumDistance
     !!{
     Constructors for the \refClass{nodeOperatorSatelliteMinimumDistance} node operator class.
     !!}
     module procedure satelliteMinimumDistanceConstructorParameters
     module procedure satelliteMinimumDistanceConstructorInternal
  end interface nodeOperatorSatelliteMinimumDistance
  
contains

  function satelliteMinimumDistanceConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorSatelliteMinimumDistance} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodeOperatorSatelliteMinimumDistance)                :: self
    type(inputParameters                     ), intent(inout) :: parameters
    type(varying_string                      )                :: relativeTo

    !![
    <inputParameter>
      <name>relativeTo</name>
      <defaultValue>var_str('immediateHost')</defaultValue>
      <description>
	Specifies to which host halo the minimum distance should be referenced. ``{\normalfont \ttfamily immediateHost}'' computes
	the minimum distance to the first host (so, for a sub-subhalo, this would be the subhalo in which it is
	orbiting). ``{\normalfont \ttfamily isolatedHost}'' computes the minimum distance to the final (isolated halo) host.
      </description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=nodeOperatorSatelliteMinimumDistance(enumerationRelativeToEncode(char(relativeTo),includesPrefix=.false.))
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function satelliteMinimumDistanceConstructorParameters

  function satelliteMinimumDistanceConstructorInternal(relativeTo) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorSatelliteMinimumDistance} node operator class.
    !!}
    implicit none
    type(nodeOperatorSatelliteMinimumDistance)                :: self
    type(enumerationRelativeToType           ), intent(in   ) :: relativeTo
    !![
    <constructorAssign variables="relativeTo"/>
    !!]
    
    !![
    <addMetaProperty component="satellite" name="satelliteDistanceMinimum" id="self%satelliteDistanceMinimumID" isEvolvable="yes"  isCreator="yes"/>
    <addMetaProperty component="satellite" name="isolatedHostID"           id="self%isolatedHostID"             type="longInteger" isCreator="yes"/>
    !!]
    return
  end function satelliteMinimumDistanceConstructorInternal

  subroutine satelliteMinimumDistanceAutoHook(self)
    !!{
    Attach to various event hooks.
    !!}
    use :: Events_Hooks, only : satelliteHostChangeEvent, subhaloPromotionEvent, openMPThreadBindingAtLevel
    implicit none
    class(nodeOperatorSatelliteMinimumDistance), intent(inout) :: self

    call satelliteHostChangeEvent%attach(self,satelliteHostChange ,openMPThreadBindingAtLevel,label='hierarchy')
    call    subhaloPromotionEvent%attach(self,nodeSubhaloPromotion,openMPThreadBindingAtLevel,label='hierarchy')
    return
  end subroutine satelliteMinimumDistanceAutoHook

  subroutine satelliteMinimumDistanceDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorSatelliteMinimumDistance} node operator class.
    !!}
    use :: Events_Hooks, only : satelliteHostChangeEvent, subhaloPromotionEvent
    implicit none
    type(nodeOperatorSatelliteMinimumDistance), intent(inout) :: self

    if (satelliteHostChangeEvent%isAttached(self,satelliteHostChange )) call satelliteHostChangeEvent%detach(self,satelliteHostChange )
    if (   subhaloPromotionEvent%isAttached(self,nodeSubhaloPromotion)) call    subhaloPromotionEvent%detach(self,nodeSubhaloPromotion)
    return
  end subroutine satelliteMinimumDistanceDestructor

  subroutine satelliteMinimumDistanceNodesMerge(self,node)
    !!{
    Update the minimum distance of approach when two nodes merge.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSatellite
    implicit none
    class(nodeOperatorSatelliteMinimumDistance), intent(inout) :: self
    type (treeNode                            ), intent(inout) :: node
    type (treeNode                            ), pointer       :: nodeIsolated
    class(nodeComponentSatellite              ), pointer       :: satellite

    satellite    => node%satellite     ()
    nodeIsolated => node%isolatedParent()
    call satellite%      floatRank0MetaPropertySet(self%satelliteDistanceMinimumID,self        %distanceRelative(node))
    call satellite%longIntegerRank0MetaPropertySet(self%            isolatedHostID,nodeIsolated%uniqueID        (    ))
    return
  end subroutine satelliteMinimumDistanceNodesMerge

  subroutine satelliteMinimumDistanceDifferentialEvolutionPost(self,node)
    !!{
    Update the minimum distance of approach after differential evolution of the node.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSatellite
    implicit none
    class(nodeOperatorSatelliteMinimumDistance), intent(inout) :: self
    type (treeNode                            ), intent(inout) :: node
    type (treeNode                            ), pointer       :: nodeIsolated
    class(nodeComponentSatellite              ), pointer       :: satellite

    if (.not.node%isSatellite()) return
    satellite    => node%satellite     ()
    nodeIsolated => node%isolatedParent()
    if (satellite%longIntegerRank0MetaPropertyGet(self%isolatedHostID) == 0_kind_int8) then
       ! This node has not been seen before - so simply set the minimum distance to the current distance.
       call satellite%      floatRank0MetaPropertySet(self%satelliteDistanceMinimumID,self        %distanceRelative(node))
       call satellite%longIntegerRank0MetaPropertySet(self%            isolatedHostID,nodeIsolated%uniqueID        (    ))
    else
       ! The node has been seen before, so set the minimum distance to the smaller of the current and prior minimum.
       call satellite%floatRank0MetaPropertySet(                                                                          &
            &                                                                           self%satelliteDistanceMinimumID , &
            &                                   min(                                                                      &
            &                                       satellite%floatRank0MetaPropertyGet(self%satelliteDistanceMinimumID), &
            &                                       self     %distanceRelative         (node                           )  &
            &                                   )                                                                         &
            &                                  )
    end if
    return
  end subroutine satelliteMinimumDistanceDifferentialEvolutionPost

  subroutine satelliteHostChange(self,node)
    !!{
    Update the minimum distance of approach when a satellite changes host node.
    !!}
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : nodeComponentSatellite
    implicit none
    class(*                     ), intent(inout)          :: self
    type (treeNode              ), intent(inout), target  :: node
    type (treeNode              )               , pointer :: nodeIsolated
    class(nodeComponentSatellite)               , pointer :: satellite

    select type (self)
    class is (nodeOperatorSatelliteMinimumDistance)
       satellite    => node%satellite     ()
       nodeIsolated => node%isolatedParent()
       ! Only apply if this node has already been operated on. This avoids opererating in cases where the node should be filtered
       ! out.
       if     (                                                                                                                                                &
            &                                                       satellite%longIntegerRank0MetaPropertyGet(self%isolatedHostID) /= 0_c_size_t               &
            &  .and.                                                                                                                                           &
            &   (                                                                                                                                              &
            &      self%relativeTo == relativeToImmediateHost                                                                                                  &
            &    .or.                                                                                                                                          &
            &     (self%relativeTo == relativeToIsolatedHost  .and. satellite%longIntegerRank0MetaPropertyGet(self%isolatedHostID) /= nodeIsolated%uniqueID()) &
            &   )&
            & ) then
          call satellite%      floatRank0MetaPropertySet(self%satelliteDistanceMinimumID,self        %distanceRelative(node))
          call satellite%longIntegerRank0MetaPropertySet(self%            isolatedHostID,nodeIsolated%uniqueID        (    ))
       end if
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine satelliteHostChange
  
  subroutine nodeSubhaloPromotion(self,node,nodePromotion)
    !!{
    Reset the minimum distance of approach when a satellite is promoted to be an isolated halo.
    !!}
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : nodeComponentSatellite
    implicit none
    class(*                     ), intent(inout)          :: self
    type (treeNode              ), intent(inout), pointer :: node        , nodePromotion
    type (treeNode              )               , pointer :: nodeIsolated
    class(nodeComponentSatellite)               , pointer :: satellite
    !$GLC attributes unused :: nodePromotion

    select type (self)
    class is (nodeOperatorSatelliteMinimumDistance)
       satellite    => node%satellite     ()
       nodeIsolated => node%isolatedParent()
       ! Only apply if this node has already been operated on. This avoids opererating in cases where the node should be filtered
       ! out.
       if (satellite%longIntegerRank0MetaPropertyGet(self%isolatedHostID) /= 0_c_size_t) then
          call satellite%      floatRank0MetaPropertySet(self%satelliteDistanceMinimumID,-1.0d0                 )
          call satellite%longIntegerRank0MetaPropertySet(self%            isolatedHostID,nodeIsolated%uniqueID())
       end if
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select    
    return
  end subroutine nodeSubhaloPromotion

  double precision function satelliteMinimumDistanceDistanceRelative(self,node) result(distance)
    !!{
    Compute the current distance to the relevant host halo.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSatellite
    use :: Vectors         , only : Vector_Magnitude
    implicit none
    class           (nodeOperatorSatelliteMinimumDistance), intent(inout)          :: self
    type            (treeNode                            ), intent(inout), target  :: node
    type            (treeNode                            )               , pointer :: nodeWork
    class           (nodeComponentSatellite              )               , pointer :: satellite
    double precision                                      , dimension(3)           :: position

    position = 0.0d0
    nodeWork => node
    do while (associated(nodeWork))
       satellite =>  nodeWork %satellite()
       position  =  +          position    &
            &       +satellite%position ()
       if (self%relativeTo == relativeToImmediateHost) then
          nodeWork => null()
       else if (nodeWork%isSatellite()) then
          nodeWork => nodeWork%parent
       else
          nodeWork => null()
       end if
    end do
    distance=Vector_Magnitude(position)
    return
  end function satelliteMinimumDistanceDistanceRelative
