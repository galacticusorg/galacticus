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
  Implements a node operator class that triggers merging of satellites based on their orbital radius.
  !!}

  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass
  use :: Kepler_Orbits          , only : keplerOrbitCount

  !![
  <nodeOperator name="nodeOperatorSatelliteMergingRadiusTrigger">
   <description>A node operator class that triggers merging of satellites based on their orbital radius.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorSatelliteMergingRadiusTrigger
     !!{
     A node operator class that triggers merging of satellites based on their orbital radius.
     !!}
     private
     class           (darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_                            => null()
     double precision                                    :: radiusVirialFraction
     logical                                             :: recordMergedSubhaloProperties                            , recordFirstLevelOnly
     integer                                             :: mergedSubhaloIDs             (keplerOrbitCount)          , nodeHierarchyLevelMaximumID
   contains
     !![
     <methods>
       <method description="Compute the radius at which the satellite will be merged." method="radiusMerge" />
     </methods>
     !!]
     final     ::                          satelliteMergingRadiusTriggerDestructor
     procedure :: differentialEvolution => satelliteMergingRadiusTriggerDifferentialEvolution
     procedure :: radiusMerge           => satelliteMergingRadiusTriggerRadiusMerge
  end type nodeOperatorSatelliteMergingRadiusTrigger
  
  interface nodeOperatorSatelliteMergingRadiusTrigger
     !!{
     Constructors for the \refClass{nodeOperatorSatelliteMergingRadiusTrigger} node operator class.
     !!}
     module procedure satelliteMergingRadiusTriggerConstructorParameters
     module procedure satelliteMergingRadiusTriggerConstructorInternal
  end interface nodeOperatorSatelliteMergingRadiusTrigger

  ! Sub-module-scope pointer to self used in callback function.
  class(nodeOperatorSatelliteMergingRadiusTrigger), pointer :: self_
  !$omp threadprivate(self)
  
contains

  function satelliteMergingRadiusTriggerConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorSatelliteMergingRadiusTrigger} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nodeOperatorSatelliteMergingRadiusTrigger)                :: self
    type            (inputParameters                          ), intent(inout) :: parameters
    class           (darkMatterHaloScaleClass                 ), pointer       :: darkMatterHaloScale_
    double precision                                                           :: radiusVirialFraction
    logical                                                                    :: recordMergedSubhaloProperties, recordFirstLevelOnly

    !![
    <inputParameter>
      <name>radiusVirialFraction</name>
      <defaultValue>0.01d0</defaultValue>
      <description>The fraction of the virial radius below which satellites are merged.</description>
      <source>parameters</source>
    </inputParameter>
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
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=nodeOperatorSatelliteMergingRadiusTrigger(radiusVirialFraction,recordMergedSubhaloProperties,recordFirstLevelOnly,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function satelliteMergingRadiusTriggerConstructorParameters

  function satelliteMergingRadiusTriggerConstructorInternal(radiusVirialFraction,recordMergedSubhaloProperties,recordFirstLevelOnly,darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorSatelliteMergingRadiusTrigger} node operator class.
    !!}
    use :: Kepler_Orbits, only : keplerOrbitTimeInitial     , keplerOrbitMassSatellite, keplerOrbitMassHost, keplerOrbitRadius, &
         &                       keplerOrbitRadiusPericenter, keplerOrbitTimeCurrent
    implicit none
    type            (nodeOperatorSatelliteMergingRadiusTrigger)                        :: self
    double precision                                           , intent(in   )         :: radiusVirialFraction
    logical                                                    , intent(in   )         :: recordMergedSubhaloProperties, recordFirstLevelOnly
    class           (darkMatterHaloScaleClass                 ), intent(in   ), target :: darkMatterHaloScale_
    !![
    <constructorAssign variables="radiusVirialFraction, recordMergedSubhaloProperties, recordFirstLevelOnly, *darkMatterHaloScale_"/>
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
  end function satelliteMergingRadiusTriggerConstructorInternal
  
  subroutine satelliteMergingRadiusTriggerDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorSatelliteMergingRadiusTrigger} node operator class.
    !!}
    implicit none
    type(nodeOperatorSatelliteMergingRadiusTrigger), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_"/>
    !!]
    return
  end subroutine satelliteMergingRadiusTriggerDestructor
  
  subroutine satelliteMergingRadiusTriggerDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Trigger merging of a satellite halo based on its orbital radius.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentSatellite
    use :: Vectors         , only : Vector_Magnitude
    implicit none
    class           (nodeOperatorSatelliteMergingRadiusTrigger), intent(inout), target  :: self
    type            (treeNode                                 ), intent(inout), target  :: node
    logical                                                    , intent(inout)          :: interrupt
    procedure       (interruptTask                            ), intent(inout), pointer :: functionInterrupt
    integer                                                    , intent(in   )          :: propertyType
    class           (nodeComponentSatellite                   )               , pointer :: satellite
    double precision                                           , dimension(3)           :: position
    double precision                                                                    :: radius
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
  end subroutine satelliteMergingRadiusTriggerDifferentialEvolution

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
             propertyNew(size(propertyNew))=property
             call basicHost%floatRank1MetaPropertySet(self_%mergedSubhaloIDs(ID),propertyNew)
             deallocate(propertyCurrent)
             deallocate(propertyNew    )
          end do
       end if
    end if
    return
  end subroutine mergerTrigger

  double precision function satelliteMergingRadiusTriggerRadiusMerge(self,node)
    !!{
    Compute the merging radius for a node.
    !!}
    use :: Galacticus_Nodes          , only : treeNode
    use :: Galactic_Structure_Options, only : massTypeGalactic
    use :: Mass_Distributions        , only : massDistributionClass
    implicit none
    class           (nodeOperatorSatelliteMergingRadiusTrigger), intent(inout) :: self
    type            (treeNode                                 ), intent(inout) :: node
    type            (treeNode                                 ), pointer       :: nodeHost
    class           (massDistributionClass                    ), pointer       :: massDistribution_    , massDistributionHost_
    double precision                                                           :: radiusHalfMassCentral, radiusHalfMassSatellite

    ! Find the host node.
    nodeHost => node%mergesWith()
    ! Get mass distributions.
    massDistribution_     => node    %massDistribution(massType=massTypeGalactic)
    massDistributionHost_ => nodeHost%massDistribution(massType=massTypeGalactic)
    ! Get half-mass radii of central and satellite galaxies. We first check that the total mass in the galactic component is
    ! non-zero as we do not want to attempt to find the half-mass radius of the galactic component, if no galactic component
    ! exists. To correctly handle the case that numerical errors lead to a zero-size galactic component (the enclosed mass
    ! within zero radius is non-zero and equals to the total mass of this component), we do a further check that the enclosed
    ! mass within zero radius is smaller than half of the total mass in the galactic component.
    if     (                                                                                    &
         &             massDistributionHost_%massTotal()                                        &
         &   >                                                                                  &
         &   max(                                                                               &
         &       0.0d0,                                                                         &
         &       2.0d0*massDistributionHost_%massEnclosedBySphere(radius=0.0d0)                 &
         &      )                                                                               &
         & ) then
       radiusHalfMassCentral  =massDistributionHost_%radiusEnclosingMass(massFractional=0.5d0)
    else
       radiusHalfMassCentral  =0.0d0
    end if
    if     (                                                                                    &
         &             massDistribution_    %massTotal()                                        &
         &   >                                                                                  &
         &   max(                                                                               &
         &       0.0d0,                                                                         &
         &       2.0d0*massDistribution_    %massEnclosedBySphere(radius=0.0d0)                 &
         &      )                                                                               &
         & ) then
       radiusHalfMassSatellite=massDistribution_    %radiusEnclosingMass(massFractional=0.5d0)
    else
       radiusHalfMassSatellite=0.0d0
    end if
    !![
    <objectDestructor name="massDistribution_"    />
    <objectDestructor name="massDistributionHost_"/>
    !!]
    satelliteMergingRadiusTriggerRadiusMerge=max(                                                              &
         &                                       +                          radiusHalfMassSatellite            &
         &                                       +                          radiusHalfMassCentral            , &
         &                                       +self%                     radiusVirialFraction               &
         &                                       *self%darkMatterHaloScale_%radiusVirial           (nodeHost)  &
         &                                      )
    return
  end function satelliteMergingRadiusTriggerRadiusMerge
