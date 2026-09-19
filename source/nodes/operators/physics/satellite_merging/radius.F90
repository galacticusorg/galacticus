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

  !+    Contributions to this file made by: Claude.

  !!{RST
  Implements an abstract node operator class which triggers merging of satellites once their orbital radius falls below a threshold.
  !!}

  use :: Kepler_Orbits, only : keplerOrbitCount

  !![
  <nodeOperator name="nodeOperatorSatelliteMergingRadius" abstract="yes" docformat="rst">
   <description>
   An abstract node operator class which merges a satellite with its host once the satellite's orbital radius falls below a threshold radius, and which optionally records the Keplerian orbital elements of merged subhalos. This is an abstract class---the threshold radius itself must be provided by a concrete class. ``recordMergedSubhaloProperties`` records the orbital elements of merged subhalos; ``recordFirstLevelOnly`` restricts that recording to first-level subhalos relative to the final host.
   </description>
  </nodeOperator>
  !!]
  type, abstract, extends(nodeOperatorClass) :: nodeOperatorSatelliteMergingRadius
     !!{RST
     An abstract node operator class which triggers merging of satellites once their orbital radius falls below a threshold.
     !!}
     private
     logical :: recordMergedSubhaloProperties                  , recordFirstLevelOnly
     integer :: mergedSubhaloIDs             (keplerOrbitCount), nodeHierarchyLevelMaximumID
   contains
     !![
     <methods docformat="rst">
       <method description="Compute the radius at which the satellite will be merged."                    method="radiusMerge"         />
       <method description="Initialize the meta-properties used to record properties of merged subhalos." method="recordingInitialize" />
     </methods>
     !!]
     procedure                                              :: differentialEvolution => satelliteMergingRadiusDifferentialEvolution
     procedure                                              :: recordingInitialize   => satelliteMergingRadiusRecordingInitialize
     procedure(satelliteMergingRadiusRadiusMerge), deferred :: radiusMerge
  end type nodeOperatorSatelliteMergingRadius

  abstract interface
     double precision function satelliteMergingRadiusRadiusMerge(self,node)
       !!{RST
       Abstract interface for the merging-radius method of radius-triggered satellite merging node operator classes.
       !!}
       import nodeOperatorSatelliteMergingRadius, treeNode
       class(nodeOperatorSatelliteMergingRadius), intent(inout) :: self
       type (treeNode                          ), intent(inout) :: node
     end function satelliteMergingRadiusRadiusMerge
  end interface

  ! Sub-module-scope pointer to self used in callback function.
  class(nodeOperatorSatelliteMergingRadius), pointer :: self_
  !$omp threadprivate(self_)

contains

  subroutine satelliteMergingRadiusRecordingInitialize(self)
    !!{RST
    Initialize the meta-properties used to record the properties of merged subhalos.
    !!}
    use :: Kepler_Orbits, only : keplerOrbitTimeInitial     , keplerOrbitMassSatellite, keplerOrbitMassHost, keplerOrbitRadius, &
         &                       keplerOrbitRadiusPericenter, keplerOrbitTimeCurrent
    implicit none
    class(nodeOperatorSatelliteMergingRadius), intent(inout) :: self

    if (self%recordMergedSubhaloProperties) then
       !![
       <addMetaProperty component="basic" name="mergedSubhaloTimeCurrent"                     id="self%mergedSubhaloIDs(keplerOrbitTimeCurrent     %ID)" rank="1" isCreator="yes"/>
       <addMetaProperty component="basic" name="mergedSubhaloTimeInitial"                     id="self%mergedSubhaloIDs(keplerOrbitTimeInitial     %ID)" rank="1" isCreator="yes"/>
       <addMetaProperty component="basic" name="mergedSubhaloMassSatellite"                   id="self%mergedSubhaloIDs(keplerOrbitMassSatellite   %ID)" rank="1" isCreator="yes"/>
       <addMetaProperty component="basic" name="mergedSubhaloMassHost"                        id="self%mergedSubhaloIDs(keplerOrbitMassHost        %ID)" rank="1" isCreator="yes"/>
       <addMetaProperty component="basic" name="mergedSubhaloRadius"                          id="self%mergedSubhaloIDs(keplerOrbitRadius          %ID)" rank="1" isCreator="yes"/>
       <addMetaProperty component="basic" name="mergedSubhaloRadiusPericenter"                id="self%mergedSubhaloIDs(keplerOrbitRadiusPericenter%ID)" rank="1" isCreator="yes"/>
       <addMetaProperty component="basic" name="nodeHierarchyLevelMaximum"     type="integer" id="self%nodeHierarchyLevelMaximumID"                               isCreator="no" />
       !!]
    end if
    return
  end subroutine satelliteMergingRadiusRecordingInitialize

  subroutine satelliteMergingRadiusDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{RST
    Trigger merging of a satellite halo based on its orbital radius.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSatellite
    use :: Vectors         , only : Vector_Magnitude
    implicit none
    class           (nodeOperatorSatelliteMergingRadius), intent(inout), target  :: self
    type            (treeNode                          ), intent(inout), target  :: node
    logical                                             , intent(inout)          :: interrupt
    procedure       (interruptTask                     ), intent(inout), pointer :: functionInterrupt
    integer                                             , intent(in   )          :: propertyType
    class           (nodeComponentSatellite            )               , pointer :: satellite
    double precision                                    , dimension(3)           :: position
    double precision                                                             :: radius
    !$GLC attributes unused :: propertyType
    
    if (.not.node%isSatellite()) return
    satellite => node     %satellite(        )
    position  =  satellite%position (        )
    radius    =  Vector_Magnitude   (position)
    ! Test for merging.
    if     (                                 &
         &   radius > 0.0d0                  &
         &  .and.                            &
         &   radius < self%radiusMerge(node) &
         & ) then
       ! Merging criterion met - trigger an interrupt.
       interrupt         =  .true.
       functionInterrupt => mergerTrigger
       self_             => self
    end if
    return
  end subroutine satelliteMergingRadiusDifferentialEvolution

  subroutine mergerTrigger(node,timeEnd)
    !!{RST
    Trigger a merger of the satellite by setting the time until merging to zero.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSatellite     , nodeComponentBasic    , treeNode
    use :: Kepler_Orbits   , only : keplerOrbit                , keplerOrbitTimeInitial, keplerOrbitMassSatellite, keplerOrbitMassHost, &
         &                          keplerOrbitRadiusPericenter, keplerOrbitRadius     , keplerOrbitTimeCurrent
    implicit none
    type            (treeNode              ), intent(inout), target      :: node
    double precision                        , intent(in   ), optional    :: timeEnd
    type            (treeNode              )               , pointer     :: nodeHost
    class           (nodeComponentBasic    )               , pointer     :: basic          , basicHost
    class           (nodeComponentSatellite)               , pointer     :: satellite
    double precision                        , dimension(:) , allocatable :: propertyCurrent, propertyNew
    double precision                                                     :: property
    type            (keplerOrbit           )                             :: orbit
    integer                                                              :: i              , ID
    !$GLC attributes unused :: timeEnd

    ! Set the time of merging to the current time.
    basic     => node%basic    ()
    satellite => node%satellite()
    call satellite%timeOfMergingSet(basic%time())
    ! Record properties of the merging subhalo if necessary.
    if (self_%recordMergedSubhaloProperties) then
       ! Find the node to merge with.
       nodeHost  => node    %mergesWith()
       basicHost => nodeHost%basic     ()
       ! Only record if we are recording mergers from all levels of the hierarchy, or if this is a first level subhalo relative to the host.
       if     (                                                                             &
            &   .not.self_%recordFirstLevelOnly                                             &
            &  .or.                                                                         &
            &    basic    %integerRank0MetaPropertyGet(self_%nodeHierarchyLevelMaximumID)   &
            &   ==                                                                          &
            &    basicHost%integerRank0MetaPropertyGet(self_%nodeHierarchyLevelMaximumID)+1 &
            & ) then
          ! Get the virial orbit of the halo about to merge.
          orbit=satellite%virialOrbit()
          ! Append the orbit data.
          do i=1,6
             select case (i)
             case (1)
                ID      =keplerOrbitTimeInitial     %ID
                property=basic%timeLastIsolated()
             case (2)
                ID      =keplerOrbitTimeCurrent     %ID
                property=basic%time            ()
             case (3)
                ID      =keplerOrbitMassSatellite   %ID
                property=orbit%massSatellite   ()
             case (4)
                ID      =keplerOrbitMassHost        %ID
                property=orbit%massHost        ()
             case (5)
                ID      =keplerOrbitRadius          %ID
                property=orbit%radius          ()
             case (6)
                ID      =keplerOrbitRadiusPericenter%ID
                property=orbit%radiusPericenter()
             end select
             propertyCurrent=basicHost%floatRank1MetaPropertyGet(self_%mergedSubhaloIDs(ID))
             allocate(propertyNew(size(propertyCurrent)+1_c_size_t))
             propertyNew(1_c_size_t:size(propertyCurrent))=propertyCurrent(:)
             propertyNew(size(propertyNew))=property
             call basicHost%floatRank1MetaPropertySet(self_%mergedSubhaloIDs(ID),propertyNew)
             deallocate(propertyCurrent)
             deallocate(propertyNew    )
          end do
       end if
    end if
    return
  end subroutine mergerTrigger
