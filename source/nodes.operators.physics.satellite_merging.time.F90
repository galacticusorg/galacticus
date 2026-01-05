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
  Implements a node operator class that triggers merging of satellites based on a merging time.
  !!}

  use :: Satellite_Merging_Timescales, only : satelliteMergingTimescalesClass
  use :: Virial_Orbits               , only : virialOrbitClass
  use :: Kind_Numbers                , only : kind_int8

  !![
  <nodeOperator name="nodeOperatorSatelliteMergingTime">
   <description>A node operator class that triggers merging of satellites based on a merging time.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorSatelliteMergingTime
     !!{
     A node operator class that triggers merging of satellites based on a merging time.
     !!}
     private
     class  (virialOrbitClass               ), pointer :: virialOrbit_                => null()
     class  (satelliteMergingTimescalesClass), pointer :: satelliteMergingTimescales_ => null()
     integer(kind_int8                      )          :: uniqueIDProcessed
     logical                                           :: resetOnHaloFormation
   contains
     !![
     <methods>
       <method method="timeMergingSet" description="Set the time of merging for a satellite node."/>
     </methods>
     !!]
     final     ::                       satelliteMergingTimeDestructor
     procedure :: autoHook           => satelliteMergingTimeAutoHook
     procedure :: nodeTreeInitialize => satelliteMergingTimeNodeTreeInitialize
     procedure :: nodesMerge         => satelliteMergingTimeNodesMerge
     procedure :: timeMergingSet     => satelliteMergingTimeTimeMergingSet
  end type nodeOperatorSatelliteMergingTime
  
  interface nodeOperatorSatelliteMergingTime
     !!{
     Constructors for the \refClass{nodeOperatorSatelliteMergingTime} node operator class.
     !!}
     module procedure satelliteMergingTimeConstructorParameters
     module procedure satelliteMergingTimeConstructorInternal
  end interface nodeOperatorSatelliteMergingTime

contains

  function satelliteMergingTimeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorSatelliteMergingTime} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (nodeOperatorSatelliteMergingTime)                :: self
    type   (inputParameters                 ), intent(inout) :: parameters
    class  (virialOrbitClass                ), pointer       :: virialOrbit_
    class  (satelliteMergingTimescalesClass ), pointer       :: satelliteMergingTimescales_
    logical                                                  :: resetOnHaloFormation

    !![
    <inputParameter>
      <name>resetOnHaloFormation</name>
      <defaultValue>.false.</defaultValue>
      <description>Specifies whether satellite virial orbital parameters should be reset on halo formation events.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="virialOrbit"                name="virialOrbit_"                source="parameters"/>
    <objectBuilder class="satelliteMergingTimescales" name="satelliteMergingTimescales_" source="parameters"/>
    !!]
    self=nodeOperatorSatelliteMergingTime(resetOnHaloFormation,virialOrbit_,satelliteMergingTimescales_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="virialOrbit_"               />
    <objectDestructor name="satelliteMergingTimescales_"/>
    !!]
    return
  end function satelliteMergingTimeConstructorParameters

  function satelliteMergingTimeConstructorInternal(resetOnHaloFormation,virialOrbit_,satelliteMergingTimescales_) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorSatelliteMergingTime} node operator class.
    !!}
    implicit none
    type   (nodeOperatorSatelliteMergingTime)                        :: self
    class  (virialOrbitClass                ), intent(in   ), target :: virialOrbit_
    class  (satelliteMergingTimescalesClass ), intent(in   ), target :: satelliteMergingTimescales_
    logical                                  , intent(in   )         :: resetOnHaloFormation
    !![
    <constructorAssign variables="resetOnHaloFormation, *virialOrbit_, *satelliteMergingTimescales_"/>
    !!]
    
    return
  end function satelliteMergingTimeConstructorInternal

  subroutine satelliteMergingTimeAutoHook(self)
    !!{
    Attach to various event hooks.
    !!}
    use :: Events_Hooks, only : satelliteHostChangeEvent, haloFormationEvent, openMPThreadBindingAtLevel
    implicit none
    class(nodeOperatorSatelliteMergingTime), intent(inout) :: self
    
    call satelliteHostChangeEvent%attach(self,satelliteHostChange,openMPThreadBindingAtLevel,label='satelliteMergingTime')
    call       haloFormationEvent%attach(self,haloFormation      ,openMPThreadBindingAtLevel,label='satelliteMergingTime')
    return
  end subroutine satelliteMergingTimeAutoHook

  subroutine satelliteMergingTimeDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorSatelliteMergingTime} node operator class.
    !!}
    use :: Events_Hooks, only : satelliteHostChangeEvent, haloFormationEvent
    implicit none
    type(nodeOperatorSatelliteMergingTime), intent(inout) :: self

    !![
    <objectDestructor name="self%virialOrbit_"               />
    <objectDestructor name="self%satelliteMergingTimescales_"/>
    !!]
    if (satelliteHostChangeEvent%isAttached(self,satelliteHostChange)) call satelliteHostChangeEvent%detach(self,satelliteHostChange)
    if (      haloFormationEvent%isAttached(self,haloFormation      )) call       haloFormationEvent%detach(self,haloFormation      )
    return
  end subroutine satelliteMergingTimeDestructor

  subroutine satelliteMergingTimeNodeTreeInitialize(self,node)
    !!{
    Initialize merging time of any initial satellites in a tree.
    !!}
    implicit none
    class(nodeOperatorSatelliteMergingTime), intent(inout), target  :: self
    type (treeNode                        ), intent(inout), target  :: node

    if     (                                               &
         &                     node%isSatellite        ()  &
         &  .or.                                           &
         &   (                                             &
         &     .not.           node%isPrimaryProgenitor()  &
         &    .and.                                        &
         &          associated(node%parent               ) &
         &   )                                             &
         & ) call self%timeMergingSet(node)
    return
  end subroutine satelliteMergingTimeNodeTreeInitialize
  
  subroutine satelliteMergingTimeNodesMerge(self,node)
    !!{
    Act on node merging tree.
    !!}
    implicit none
    class(nodeOperatorSatelliteMergingTime), intent(inout) :: self
    type (treeNode                        ), intent(inout) :: node

    call self%timeMergingSet(node)
    return
  end subroutine satelliteMergingTimeNodesMerge
  
  subroutine satelliteHostChange(self,node)
    !!{
    Handle cases where a satellite switches host node.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(*       ), intent(inout)         :: self
    type (treeNode), intent(inout), target :: node

    select type (self)
    class is (nodeOperatorSatelliteMergingTime)
       call self%timeMergingSet(node)
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine satelliteHostChange

  subroutine haloFormation(self,node)
    !!{
    Handle cases where a host halo reforms.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(*       ), intent(inout) :: self
    type (treeNode), intent(inout) :: node
    type (treeNode), pointer       :: nodeSatellite

    select type (self)
    class is (nodeOperatorSatelliteMergingTime)
       ! Return immediately if orbits are not to be reset.
       if (.not.self%resetOnHaloFormation) return
       ! Iterate over all satellites.
       nodeSatellite => node%firstSatellite
       do while (associated(nodeSatellite))
          ! Set a new merging time for this satellite.
          call self%timeMergingSet(nodeSatellite)
          nodeSatellite => nodeSatellite%sibling
       end do
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine haloFormation

  subroutine satelliteMergingTimeTimeMergingSet(self,node)
    !!{
    Set the time of merging for the given {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentSatellite
    use :: Kepler_Orbits   , only : keplerOrbit       , keplerOrbitMasses        , keplerOrbitRadius            , keplerOrbitTheta, &
         &                          keplerOrbitPhi    , keplerOrbitVelocityRadial, keplerOrbitVelocityTangential
    implicit none
    class           (nodeOperatorSatelliteMergingTime), intent(inout) :: self
    type            (treeNode                        ), intent(inout) :: node
    class           (nodeComponentBasic              ), pointer       :: basic
    class           (nodeComponentSatellite          ), pointer       :: satellite
    type            (treeNode                        ), pointer       :: nodeHost
    logical                                           , parameter     :: acceptUnboundOrbits=.false.
    double precision                                                  :: timeUntilMerging           , timeHost
    type            (keplerOrbit                     )                :: orbit

    ! Get the satellite component.
    satellite => node%satellite(autoCreate=.true.)
    ! Get an orbit for this satellite.
    if (node%isSatellite()) then
       nodeHost => node %parent
       basic    => node %basic            ()
       timeHost =  basic%time             ()
    else
       nodeHost => node %parent%firstChild
       basic    => node %parent%basic     ()
       timeHost =  basic%time             ()
    end if
    call orbit%reset()
    if (satellite%virialOrbitIsGettable()) &
         & orbit=satellite%virialOrbit()
    if (orbit%isDefined()) then
       call orbit%reset(keep=[keplerOrbitMasses,keplerOrbitRadius,keplerOrbitTheta,keplerOrbitPhi,keplerOrbitVelocityRadial,keplerOrbitVelocityTangential])
    else
       orbit=self%virialOrbit_%orbit(node,nodeHost,acceptUnboundOrbits)
    end if
    if (satellite%virialOrbitIsSettable())      &
         & call satellite%virialOrbitSet(orbit)
    ! Compute and store a time until merging.
    timeUntilMerging=self%satelliteMergingTimescales_%timeUntilMerging(node,orbit)
    if (timeUntilMerging >= 0.0d0) then
       call satellite%timeOfMergingSet(                  &
            &                          +timeUntilMerging &
            &                          +timeHost         &
            &                         )
    end if
    return
  end subroutine satelliteMergingTimeTimeMergingSet
  
