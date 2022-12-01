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
  Implements a node operator class that tracks the minimum distance from the center that a satellite has ever reached in its
  current host halo.
  !!}

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
     integer :: satelliteDistanceMinimumID
   contains
     final     ::                                satelliteMinimumDistanceDestructor
     procedure :: differentialEvolutionPost   => satelliteMinimumDistanceDifferentialEvolutionPost
     procedure :: nodesMerge                  => satelliteMinimumDistanceNodesMerge
     procedure :: autoHook                    => satelliteMinimumDistanceAutoHook
  end type nodeOperatorSatelliteMinimumDistance
  
  interface nodeOperatorSatelliteMinimumDistance
     !!{
     Constructors for the {\normalfont \ttfamily satelliteMinimumDistance} node operator class.
     !!}
     module procedure satelliteMinimumDistanceConstructorParameters
     module procedure satelliteMinimumDistanceConstructorInternal
  end interface nodeOperatorSatelliteMinimumDistance
  
contains

  function satelliteMinimumDistanceConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily satelliteMinimumDistance} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorSatelliteMinimumDistance)                :: self
    type (inputParameters                     ), intent(inout) :: parameters

    self=nodeOperatorSatelliteMinimumDistance()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function satelliteMinimumDistanceConstructorParameters

  function satelliteMinimumDistanceConstructorInternal() result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily satelliteMinimumDistance} node operator class.
    !!}
    implicit none
    type (nodeOperatorSatelliteMinimumDistance) :: self

    !![
    <addMetaProperty component="satellite" name="satelliteDistanceMinimum" id="self%satelliteDistanceMinimumID" isEvolvable="yes" isCreator="yes"/>
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
    Destructor for the {\normalfont \ttfamily satelliteMinimumDistance} node operator class.
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
    use :: Vectors         , only : Vector_Magnitude
    implicit none
    class(nodeOperatorSatelliteMinimumDistance), intent(inout) :: self
    type (treeNode                            ), intent(inout) :: node
    class(nodeComponentSatellite              ), pointer       :: satellite

    satellite => node%satellite()
    call satellite%floatRank0MetaPropertySet(self%satelliteDistanceMinimumID,Vector_Magnitude(satellite%position()))
    return
  end subroutine satelliteMinimumDistanceNodesMerge

  subroutine satelliteMinimumDistanceDifferentialEvolutionPost(self,node)
    !!{
    Update the minimum distance of approach after differential evolution of the node.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSatellite
    use :: Vectors         , only : Vector_Magnitude
    implicit none
    class(nodeOperatorSatelliteMinimumDistance), intent(inout) :: self
    type (treeNode                            ), intent(inout) :: node
    class(nodeComponentSatellite              ), pointer       :: satellite

    if (.not.node%isSatellite()) return
    satellite => node%satellite()
    call satellite%floatRank0MetaPropertySet(                                                                                 &
         &                                                                           self     %satelliteDistanceMinimumID   , &
         &                                   min(                                                                             &
         &                                       satellite%floatRank0MetaPropertyGet(self     %satelliteDistanceMinimumID  ), &
         &                                       Vector_Magnitude                   (satellite%position                  ())  &
         &                                   )                                                                                &
         &                                  )  
    return
  end subroutine satelliteMinimumDistanceDifferentialEvolutionPost

  subroutine satelliteHostChange(self,node)
    !!{
    Update the minimum distance of approach when a satellite changes host node.
    !!}
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : nodeComponentSatellite, treeNode
    use :: Vectors         , only : Vector_Magnitude
    implicit none
    class(*                     ), intent(inout) :: self
    type (treeNode              ), intent(inout) :: node
    class(nodeComponentSatellite), pointer       :: satellite

    select type (self)
    class is (nodeOperatorSatelliteMinimumDistance)
       satellite => node%satellite()
       call satellite%floatRank0MetaPropertySet(self%satelliteDistanceMinimumID,Vector_Magnitude(satellite%position())) 
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
    use :: Galacticus_Nodes, only : nodeComponentSatellite, treeNode
    implicit none
    class(*                     ), intent(inout)          :: self
    type (treeNode              ), intent(inout), pointer :: node     , nodePromotion
    class(nodeComponentSatellite)               , pointer :: satellite
    !$GLC attributes unused :: nodePromotion

    select type (self)
    class is (nodeOperatorSatelliteMinimumDistance)
       satellite => node%satellite()
       call satellite%floatRank0MetaPropertySet(self%satelliteDistanceMinimumID,-1.0d0)
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine nodeSubhaloPromotion
