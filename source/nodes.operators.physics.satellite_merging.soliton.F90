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

  !+    Contributions to this file made by: Yu Zhao

  !!{
  Implements a node operator class that triggers merging of satellites based on their orbital radius.
  !!}

  use :: Kepler_Orbits, only : keplerOrbitCount

  !![
  <nodeOperator name="nodeOperatorSatelliteMergingSoliton">
   <description>A node operator class that triggers merging of satellites based on their orbital radius.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorSatelliteMergingSoliton
     !!{
     A node operator that triggers satellite merging in FDM models when the orbital radius falls below the sum of the soliton core radii of the host and satellite.
     !!}
     private
     logical :: recordMergedSubhaloProperties                  , recordFirstLevelOnly
     integer :: mergedSubhaloIDs             (keplerOrbitCount), nodeHierarchyLevelMaximumID
     integer :: radiusCoreID                                   , massCoreID                 , &
                randomOffsetID                                 , massCoreNormalID
   contains
     !![
     <methods>
       <method description="Compute the radius at which the satellite will be merged in FDM models." method="radiusMerge" />
     </methods>
     !!]
     final     ::                          satelliteMergingSolitonDestructor
     procedure :: differentialEvolution => satelliteMergingSolitonDifferentialEvolution
     procedure :: radiusMerge           => satelliteMergingSolitonRadiusMerge
     procedure :: autoHook              => satelliteMergingSolitonAutoHook
  end type nodeOperatorSatelliteMergingSoliton
  
  interface nodeOperatorSatelliteMergingSoliton
     !!{
     Constructors for the \refClass{nodeOperatorSatelliteMergingSoliton} node operator class.
     !!}
     module procedure satelliteMergingSolitonConstructorParameters
     module procedure satelliteMergingSolitonConstructorInternal
  end interface nodeOperatorSatelliteMergingSoliton

  ! Sub-module-scope pointer to self used in callback function.
  class(nodeOperatorSatelliteMergingSoliton), pointer :: self_
  !$omp threadprivate(self)
  
contains

  function satelliteMergingSolitonConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorSatelliteMergingSoliton} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (nodeOperatorSatelliteMergingSoliton)                :: self
    type   (inputParameters                    ), intent(inout) :: parameters
    logical                                                     :: recordMergedSubhaloProperties, recordFirstLevelOnly

    !![
    <inputParameter>
      <name>recordMergedSubhaloProperties</name>
      <defaultValue>.false.</defaultValue>
      <description>If true, record the orbital properties of subhalo that merge.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>recordFirstLevelOnly</name>
      <defaultValue>.false.</defaultValue>
      <description>If true, record only mergers with first-level subhalos relative to the host.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=nodeOperatorSatelliteMergingSoliton(recordMergedSubhaloProperties,recordFirstLevelOnly)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function satelliteMergingSolitonConstructorParameters

  function satelliteMergingSolitonConstructorInternal(recordMergedSubhaloProperties,recordFirstLevelOnly) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorSatelliteMergingSoliton} node operator class.
    !!}
    use :: Kepler_Orbits, only : keplerOrbitTimeInitial     , keplerOrbitMassSatellite, keplerOrbitMassHost, keplerOrbitRadius, &
         &                       keplerOrbitRadiusPericenter, keplerOrbitTimeCurrent
    implicit none
    type   (nodeOperatorSatelliteMergingSoliton)                :: self
    logical                                     , intent(in   ) :: recordMergedSubhaloProperties, recordFirstLevelOnly
    !![
    <constructorAssign variables="recordMergedSubhaloProperties, recordFirstLevelOnly"/>
    !!]
    
    !![
    <addMetaProperty component="darkMatterProfile" name="solitonRandomOffset"   id="self%randomOffsetID"   isEvolvable="no"  isCreator="no"/>
    <addMetaProperty component="darkMatterProfile" name="solitonRadiusCore"     id="self%radiusCoreID"     isEvolvable="no"  isCreator="no"/>
    <addMetaProperty component="darkMatterProfile" name="solitonMassCore"       id="self%massCoreID"       isEvolvable="no"  isCreator="no"/>
    <addMetaProperty component="darkMatterProfile" name="solitonMassCoreNormal" id="self%massCoreNormalID" isEvolvable="yes" isCreator="no"/>
    !!]
    
    if (recordMergedSubhaloProperties) then
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
  end function satelliteMergingSolitonConstructorInternal

  subroutine satelliteMergingSolitonAutoHook(self)
    !!{
    Attach to various event hooks.
    !!}
    use :: Events_Hooks, only : satelliteMergerEvent, openMPThreadBindingAtLevel, dependencyDirectionAfter, dependencyRegEx
    implicit none
    class(nodeOperatorSatelliteMergingSoliton), intent   (inout) :: self
    type (dependencyRegEx                    ), dimension(    1) :: dependenciesSatelliteMerger

    dependenciesSatelliteMerger(1)=dependencyRegEx(dependencyDirectionAfter,'^remnantStructure:')
    call satelliteMergerEvent%attach(self,satelliteMerger,openMPThreadBindingAtLevel,label='satelliteMergingSoliton',dependencies=dependenciesSatelliteMerger)
    return
  end subroutine satelliteMergingSolitonAutoHook

  subroutine satelliteMergingSolitonDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorSatelliteMergingSoliton} node operator class.
    !!}
    use :: Events_Hooks, only : satelliteMergerEvent
    implicit none
    type(nodeOperatorSatelliteMergingSoliton), intent(inout) :: self

    if (satelliteMergerEvent%isAttached(self,satelliteMerger)) call satelliteMergerEvent%detach(self,satelliteMerger)
    return
  end subroutine satelliteMergingSolitonDestructor
  
  subroutine satelliteMergingSolitonDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Trigger merging of a satellite halo based on its orbital radius.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSatellite
    use :: Vectors         , only : Vector_Magnitude
    implicit none
    class           (nodeOperatorSatelliteMergingSoliton), intent   (inout), target  :: self
    type            (treeNode                           ), intent   (inout), target  :: node
    logical                                              , intent   (inout)          :: interrupt
    procedure       (interruptTask                      ), intent   (inout), pointer :: functionInterrupt
    integer                                              , intent   (in   )          :: propertyType
    class           (nodeComponentSatellite             )                  , pointer :: satellite
    double precision                                     , dimension(3    )          :: position
    double precision                                                                 :: radius
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
  end subroutine satelliteMergingSolitonDifferentialEvolution

  subroutine mergerTrigger(node,timeEnd)
    !!{
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
             propertyNew(           size(propertyNew    ))=property
             call basicHost%floatRank1MetaPropertySet(self_%mergedSubhaloIDs(ID),propertyNew)
             deallocate(propertyCurrent)
             deallocate(propertyNew    )
          end do
       end if
    end if
    return
  end subroutine mergerTrigger

  double precision function satelliteMergingSolitonRadiusMerge(self,node) result(radiusMerge)
    !!{
    Compute the merging radius for a node.
    !!}
    use :: Galacticus_Nodes          , only : nodeComponentSatellite, nodeComponentDarkMatterProfile, treeNode
    use :: Galactic_Structure_Options, only : massTypeGalactic
    use :: Mass_Distributions        , only : massDistributionClass
    implicit none
    class           (nodeOperatorSatelliteMergingSoliton), intent(inout) :: self
    class           (nodeComponentDarkMatterProfile     ), pointer       :: darkMatterProfile, darkMatterProfileHost
    type            (treeNode                           ), intent(inout) :: node
    type            (treeNode                           ), pointer       :: nodeHost
    double precision                                                     :: radiusCoreHost   , radiusCoreSatellite

    ! Find the host node profile.
    nodeHost              => node    %mergesWith       ()
    darkMatterProfileHost => nodeHost%darkMatterProfile()
    ! Get the satellite profile.
    darkMatterProfile     => node    %darkMatterProfile()
    ! Compute the merging radius.
    radiusCoreHost     =+darkMatterProfileHost%floatRank0MetaPropertyGet(self%radiusCoreID)
    radiusCoreSatellite=+darkMatterProfile    %floatRank0MetaPropertyGet(self%radiusCoreID)
    radiusMerge        =+radiusCoreHost      &
         &              +radiusCoreSatellite
    return
  end function satelliteMergingSolitonRadiusMerge

  subroutine satelliteMerger(self,node)
    !!{
    Merge the solitonic cores of the satellite and host halos.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSatellite, nodeComponentDarkMatterProfile, treeNode
    implicit none
    class           (*                             ), intent(inout)         :: self
    type            (treeNode                      ), intent(inout), target :: node
    type            (treeNode                      ), pointer               :: nodeHost
    class           (nodeComponentDarkMatterProfile), pointer               :: darkMatterProfile         , darkMatterProfileHost
    double precision                                , parameter             :: fractionMassRetained=0.7d0
    double precision                                                        :: massCoreHost              , massCoreSatellite      , &
            &                                                                  massCoreNormalHost        , massCoreNormalSatellite
    !$GLC attributes unused :: self

    ! Find the host node profile.
    nodeHost              => node    %mergesWith       ()
    darkMatterProfileHost => nodeHost%darkMatterProfile()
    ! Get the satellite profile.
    darkMatterProfile     => node    %darkMatterProfile()
    ! Compute the new core mass.
    massCoreNormalHost     = darkMatterProfileHost%floatRank0MetaPropertyGet(self_%massCoreNormalID)
    massCoreNormalSatellite= darkMatterProfile    %floatRank0MetaPropertyGet(self_%massCoreNormalID)
    massCoreHost           = darkMatterProfileHost%floatRank0MetaPropertyGet(self_%massCoreID      )
    massCoreSatellite      = darkMatterProfile    %floatRank0MetaPropertyGet(self_%massCoreID      )
    call darkMatterProfileHost%floatRank0MetaPropertySet(                                                                    &
         &                                                self_%massCoreNormalID                                           , &
         &                                               +fractionMassRetained*(massCoreNormalHost+massCoreNormalSatellite)  &
         &                                              )
    call darkMatterProfileHost%floatRank0MetaPropertySet(                                                                    &
         &                                                self_%massCoreID                                                 , &
         &                                               +fractionMassRetained*(massCoreHost      +massCoreSatellite      )  &
         &                                              )
    call darkMatterProfileHost%floatRank0MetaPropertySet(                                                                    &
         &                                                self_%randomOffsetID                                             , &
         &                                               +fractionMassRetained*(massCoreHost      +massCoreSatellite      )  &
         &                                               -fractionMassRetained*(massCoreNormalHost+massCoreNormalSatellite)  &
         &                                              )
    return
  end subroutine satelliteMerger


