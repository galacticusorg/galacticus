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
  Implements a node operator class that computes the stellar mass-weighted ages of disk, spheroid and nuclear star cluster
  components.
  !!}
  
  use :: Star_Formation_Rates_Disks                , only : starFormationRateDisksClass
  use :: Star_Formation_Rates_Spheroids            , only : starFormationRateSpheroidsClass
  use :: Star_Formation_Rates_Nuclear_Star_Clusters, only : starFormationRateNuclearStarClustersClass
  use :: Satellite_Merging_Mass_Movements          , only : mergerMassMovementsClass

  !![
  <nodeOperator name="nodeOperatorAgesStellarMassWeighted">
    <description>
      A node operator class that computes the stellar mass-weighted ages of disk, spheroid and nuclear star cluster components. Intended to be paired
      with the \refClass{nodePropertyExtractorAgesStellarMassWeighted} class to extract those ages for output.
    </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorAgesStellarMassWeighted
     !!{
     A node operator class that computes the stellar mass-weighted ages of disk and spheroid components.
     !!}
     private
     class  (starFormationRateDisksClass              ), pointer :: starFormationRateDisks_               => null()
     class  (starFormationRateSpheroidsClass          ), pointer :: starFormationRateSpheroids_           => null()
     class  (starFormationRateNuclearStarClustersClass), pointer :: starFormationRateNuclearStarClusters_ => null()
     class  (mergerMassMovementsClass                 ), pointer :: mergerMassMovements_                  => null()
     integer                                                     :: stellarMassFormedDiskID                        , timeStellarMassFormedDiskID    , &
          &                                                         stellarMassFormedSpheroidID                    , timeStellarMassFormedSpheroidID, &
          &                                                         stellarMassFormedNSCID                         , timeStellarMassFormedNSCID
   contains
     final     ::                                   agesStellarMassWeightedDestructor
     procedure :: differentialEvolutionScales    => agesStellarMassWeightedDifferentialEvolutionScales
     procedure :: differentialEvolutionInactives => agesStellarMassWeightedDifferentialEvolutionInactives
     procedure :: differentialEvolution          => agesStellarMassWeightedDifferentialEvolution
     procedure :: galaxiesMerge                  => agesStellarMassWeightedGalaxiesMerge
  end type nodeOperatorAgesStellarMassWeighted
  
  interface nodeOperatorAgesStellarMassWeighted
     !!{
     Constructors for the \refClass{nodeOperatorAgesStellarMassWeighted} node operator class.
     !!}
     module procedure agesStellarMassWeightedConstructorParameters
     module procedure agesStellarMassWeightedConstructorInternal
  end interface nodeOperatorAgesStellarMassWeighted
  
contains

  function agesStellarMassWeightedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorAgesStellarMassWeighted} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorAgesStellarMassWeighted      )                :: self
    type (inputParameters                          ), intent(inout) :: parameters
    class(starFormationRateDisksClass              ), pointer       :: starFormationRateDisks_
    class(starFormationRateSpheroidsClass          ), pointer       :: starFormationRateSpheroids_
    class(starFormationRateNuclearStarClustersClass), pointer       :: starFormationRateNuclearStarClusters_  
    class(mergerMassMovementsClass                 ), pointer       :: mergerMassMovements_
    
    !![
    <objectBuilder class="starFormationRateDisks"               name="starFormationRateDisks_"               source="parameters"/>
    <objectBuilder class="starFormationRateSpheroids"           name="starFormationRateSpheroids_"           source="parameters"/>
    <objectBuilder class="starFormationRateNuclearStarClusters" name="starFormationRateNuclearStarClusters_" source="parameters"/>
    <objectBuilder class="mergerMassMovements"                  name="mergerMassMovements_"                  source="parameters"/>
    !!]
    self=nodeOperatorAgesStellarMassWeighted(starFormationRateDisks_,starFormationRateSpheroids_,starFormationRateNuclearStarClusters_,mergerMassMovements_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="starFormationRateDisks_"              />
    <objectDestructor name="starFormationRateSpheroids_"          />
    <objectDestructor name="starFormationRateNuclearStarClusters_"/>
    <objectDestructor name="mergerMassMovements_"                 />
    !!]
    return
  end function agesStellarMassWeightedConstructorParameters

  function agesStellarMassWeightedConstructorInternal(starFormationRateDisks_,starFormationRateSpheroids_,starFormationRateNuclearStarClusters_,mergerMassMovements_) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorAgesStellarMassWeighted} node operator class.
    !!}
    implicit none
    type (nodeOperatorAgesStellarMassWeighted      )                        :: self
    class(starFormationRateDisksClass              ), intent(in   ), target :: starFormationRateDisks_
    class(starFormationRateSpheroidsClass          ), intent(in   ), target :: starFormationRateSpheroids_
    class(starFormationRateNuclearStarClustersClass), intent(in   ), target :: starFormationRateNuclearStarClusters_
    class(mergerMassMovementsClass                 ), intent(in   ), target :: mergerMassMovements_
    !![
    <constructorAssign variables="*starFormationRateDisks_, *starFormationRateSpheroids_, *starFormationRateNuclearStarClusters_, *mergerMassMovements_"/>
    !!]
    
    !![
    <addMetaProperty component="disk"     name="agesStellarMassFormed"     id="self%stellarMassFormedDiskID"         isEvolvable="yes" isCreator="yes"/>
    <addMetaProperty component="disk"     name="agesTimeStellarMassFormed" id="self%timeStellarMassFormedDiskID"     isEvolvable="yes" isCreator="yes"/>
    <addMetaProperty component="spheroid" name="agesStellarMassFormed"     id="self%stellarMassFormedSpheroidID"     isEvolvable="yes" isCreator="yes"/>
    <addMetaProperty component="spheroid" name="agesTimeStellarMassFormed" id="self%timeStellarMassFormedSpheroidID" isEvolvable="yes" isCreator="yes"/>
    <addMetaProperty component="NSC"      name="agesStellarMassFormed"     id="self%stellarMassFormedNSCID"          isEvolvable="yes" isCreator="yes"/>
    <addMetaProperty component="NSC"      name="agesTimeStellarMassFormed" id="self%timeStellarMassFormedNSCID"      isEvolvable="yes" isCreator="yes"/>
    !!]
    return
  end function agesStellarMassWeightedConstructorInternal

  subroutine agesStellarMassWeightedDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorAgesStellarMassWeighted} node operator class.
    !!}
    implicit none
    type(nodeOperatorAgesStellarMassWeighted), intent(inout) :: self

    !![
    <objectDestructor name="self%starFormationRateDisks_"              />
    <objectDestructor name="self%starFormationRateSpheroids_"          />
    <objectDestructor name="self%starFormationRateNuclearStarClusters_"/>
    <objectDestructor name="self%mergerMassMovements_"                 />
    !!]
    return
  end subroutine agesStellarMassWeightedDestructor

  subroutine agesStellarMassWeightedDifferentialEvolutionInactives(self,node)
    !!{
    Mark disk and spheroid age integrals as inactive for ODE solving.    
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDisk, nodeComponentSpheroid, nodeComponentNSC
    implicit none
    class(nodeOperatorAgesStellarMassWeighted), intent(inout) :: self
    type (treeNode                           ), intent(inout) :: node
    class(nodeComponentDisk                  ), pointer       :: disk
    class(nodeComponentSpheroid              ), pointer       :: spheroid
    class(nodeComponentNSC                   ), pointer       :: nuclearStarCluster
    
    ! Get disk and spheroid components.
    disk               => node%disk    ()
    spheroid           => node%spheroid()
    nuclearStarCluster => node%NSC     ()
    ! Mark as inactive.
    select type (disk    )
    type is (nodeComponentDisk     )
       ! Disk does not yet exist - nothing to do here.
    class default
       call disk              %floatRank0MetaPropertyInactive(self%    stellarMassFormedDiskID    )
       call disk              %floatRank0MetaPropertyInactive(self%timeStellarMassFormedDiskID    )
    end select
    select type (spheroid)
    type is (nodeComponentSpheroid )
       ! Spheroid does not yet exist - nothing to do here.
    class default
       call spheroid          %floatRank0MetaPropertyInactive(self%    stellarMassFormedSpheroidID)
       call spheroid          %floatRank0MetaPropertyInactive(self%timeStellarMassFormedSpheroidID)
    end select
    select type (nuclearStarCluster)
    type is (nodeComponentNSC     )
        ! NSC does not yet exist - nothing to do here.
    class default
       call nuclearStarCluster%floatRank0MetaPropertyInactive(self%    stellarMassFormedNSCID     )
       call nuclearStarCluster%floatRank0MetaPropertyInactive(self%timeStellarMassFormedNSCID     )
    end select
    return
  end subroutine agesStellarMassWeightedDifferentialEvolutionInactives
  
  subroutine agesStellarMassWeightedDifferentialEvolutionScales(self,node)
    !!{
    Set absolute ODE solver scale for the unweighted and time-weighted stellar masses.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDisk, nodeComponentSpheroid, nodeComponentNSC
    implicit none
    class           (nodeOperatorAgesStellarMassWeighted), intent(inout) :: self
    type            (treeNode                           ), intent(inout) :: node
    class           (nodeComponentDisk                  ), pointer       :: disk
    class           (nodeComponentSpheroid              ), pointer       :: spheroid
    class           (nodeComponentNSC                   ), pointer       :: nuclearStarCluster
    double precision                                     , parameter     :: massMinimum       =1.0d+0
    double precision                                     , parameter     :: timeScale         =1.0d-3
    double precision                                                     :: mass
    
    ! Get disk and spheroid components.
    disk               =>  node%disk       ()
    spheroid           =>  node%spheroid   ()
    nuclearStarCluster =>  node%NSC        ()
    ! Set scale for masses.
    mass     =  +disk%massGas    ()+spheroid%massGas    ()+nuclearStarCluster%massGas    () &
         &      +disk%massStellar()+spheroid%massStellar()+nuclearStarCluster%massStellar()        
    ! Set scales.
    select type (disk    )
    type is (nodeComponentDisk    )
       ! Disk does not yet exist - nothing to do here.
    class default
       call disk              %floatRank0MetaPropertyScale(self%    stellarMassFormedDiskID    ,max(mass,massMinimum)          )
       call disk              %floatRank0MetaPropertyScale(self%timeStellarMassFormedDiskID    ,max(mass,massMinimum)*timeScale)
    end select
    select type (spheroid)
    type is (nodeComponentSpheroid)
       ! Spheroid does not yet exist - nothing to do here.
    class default
       call spheroid          %floatRank0MetaPropertyScale(self%    stellarMassFormedSpheroidID,max(mass,massMinimum)          )
       call spheroid          %floatRank0MetaPropertyScale(self%timeStellarMassFormedSpheroidID,max(mass,massMinimum)*timeScale)
    end select
    select type (nuclearStarCluster     )
    type is (nodeComponentNSC     )
       ! NSC does not yet exist - nothing to do here.
    class default
       call nuclearStarCluster%floatRank0MetaPropertyScale(self%    stellarMassFormedNSCID     ,max(mass,massMinimum)          )
       call nuclearStarCluster%floatRank0MetaPropertyScale(self%timeStellarMassFormedNSCID     ,max(mass,massMinimum)*timeScale)
    end select
    return
  end subroutine agesStellarMassWeightedDifferentialEvolutionScales
  
  subroutine agesStellarMassWeightedDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Integrates unweighted and time-weighted star formation rates in disk, spheroid, and  nuclear star cluster components.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDisk, nodeComponentSpheroid, nodeComponentNSC, nodeComponentBasic, &
         &                          propertyActive
    implicit none
    class           (nodeOperatorAgesStellarMassWeighted), intent(inout), target  :: self
    type            (treeNode                           ), intent(inout), target  :: node
    logical                                              , intent(inout)          :: interrupt
    procedure       (interruptTask                      ), intent(inout), pointer :: functionInterrupt
    integer                                              , intent(in   )          :: propertyType
    class           (nodeComponentBasic                 )               , pointer :: basic
    class           (nodeComponentDisk                  )               , pointer :: disk
    class           (nodeComponentSpheroid              )               , pointer :: spheroid
    class           (nodeComponentNSC                   )               , pointer :: nuclearStarCluster
    double precision                                                              :: rateStarFormationDisk              , rateStarFormationSpheroid, &
         &                                                                           rateStarFormationNuclearStarCluster, time
    !$GLC attributes unused :: interrupt, functionInterrupt, propertyType

    ! Return immediately if active variables are requested.
    if (propertyActive(propertyType)) return
    ! Get the star formation rates.
    disk                                => node                                      %disk    (    )
    spheroid                            => node                                      %spheroid(    )
    nuclearStarCluster                  => node                                      %NSC     (    )
    rateStarFormationDisk               =  self%starFormationRateDisks_              %rate    (node)
    rateStarFormationSpheroid           =  self%starFormationRateSpheroids_          %rate    (node)
    rateStarFormationNuclearStarCluster =  self%starFormationRateNuclearStarClusters_%rate    (node)
    ! Find the current cosmic time.
    basic => node %basic()
    time  =  basic%time ()
    ! Accumulate rates.
    select type (disk    )
    type is (nodeComponentDisk    )
       ! Disk does not yet exist  - nothing to do here.
    class default
       call disk              %floatRank0MetaPropertyRate(self%timeStellarMassFormedDiskID    ,rateStarFormationDisk              *time)
       call disk              %floatRank0MetaPropertyRate(self%    stellarMassFormedDiskID    ,rateStarFormationDisk                   )
    end select

    select type (spheroid)
    type is (nodeComponentSpheroid )
       ! Spheroid does not yet exist - nothing to do here.
    class default
       call spheroid          %floatRank0MetaPropertyRate(self%timeStellarMassFormedSpheroidID,rateStarFormationSpheroid          *time)
       call spheroid          %floatRank0MetaPropertyRate(self%    stellarMassFormedSpheroidID,rateStarFormationSpheroid               )
    end select

    select type (nuclearStarCluster)
    type is (nodeComponentNSC     )
       ! NSC does not yet exist - nothing to do here. class default
    class default
       call nuclearStarCluster%floatRank0MetaPropertyRate(self%timeStellarmassFormedNSCID     ,rateStarFormationNuclearStarCluster*time)
       call nuclearStarCluster%floatRank0MetaPropertyRate(self%    stellarMassFormedNSCID     ,rateStarFormationNuclearStarCluster     )
    end select
    return
  end subroutine agesStellarMassWeightedDifferentialEvolution

  subroutine agesStellarMassWeightedGalaxiesMerge(self,node)
    !!{
    Combine integrals of star formation rate when galaxies merge.
    !!}
    use :: Error                           , only : Error_Report
    use :: Galacticus_Nodes                , only : nodeComponentDisk    , nodeComponentSpheroid    , nodeComponentNSC
    use :: Satellite_Merging_Mass_Movements, only : destinationMergerDisk, destinationMergerSpheroid, destinationMergerUnmoved, enumerationDestinationMergerType
    implicit none
    class  (nodeOperatorAgesStellarMassWeighted), intent(inout) :: self
    type   (treeNode                           ), intent(inout) :: node
    type   (treeNode                           ), pointer       :: nodeHost
    class  (nodeComponentDisk                  ), pointer       :: disk                      , diskHost
    class  (nodeComponentSpheroid              ), pointer       :: spheroid                  , spheroidHost
    class  (nodeComponentNSC                   ), pointer       :: nuclearStarCluster        , nuclearStarClusterHost
    type   (enumerationDestinationMergerType   )                :: destinationGasSatellite   , destinationStarsSatellite, &
         &                                                         destinationGasHost        , destinationStarsHost
    logical                                                     :: mergerIsMajor             , haveNuclearStarCluster   , &
         &                                                         haveNuclearStarClusterHost
    double precision                                            :: massNuclearStarCluster    , timeNuclearStarCluster

    ! Find the node to merge with.
    nodeHost               => node    %mergesWith(                 )
    disk                   => node    %disk      (autoCreate=.true.)
    spheroid               => node    %spheroid  (autoCreate=.true.)
    nuclearStarCluster     => node    %NSC       (                 )
    diskHost               => nodeHost%disk      (autoCreate=.true.)
    spheroidHost           => nodeHost%spheroid  (autoCreate=.true.)
    nuclearStarClusterHost => nodeHost%NSC       (                 )
    select type (nuclearStarCluster    )
    type is (nodeComponentNSC)
       haveNuclearStarCluster    =.false.
    class default
       haveNuclearStarCluster    =.true.
    end select
    select type (nuclearStarClusterHost)
    type is (nodeComponentNSC)
       haveNuclearStarClusterHost=.false.
    class default
       haveNuclearStarClusterHost=.true.
    end select
    ! Get mass movement descriptors.
    call self%mergerMassMovements_%get(node,destinationGasSatellite,destinationStarsSatellite,destinationGasHost,destinationStarsHost,mergerIsMajor)
    ! Move star formation rates within the host if necessary.
    if (haveNuclearStarClusterHost) then
       massNuclearStarCluster=nuclearStarClusterHost%floatRank0MetaPropertyGet(self%    stellarMassFormedNSCID)
       timeNuclearStarCluster=nuclearStarClusterHost%floatRank0MetaPropertyGet(self%timeStellarMassFormedNSCID)
    else
       massNuclearStarCluster=0.0d0
       timeNuclearStarCluster=0.0d0
    end if
    select case (destinationStarsHost%ID)
    case (destinationMergerDisk    %ID)
       call diskHost                     %floatRank0MetaPropertySet(                                                  self%timeStellarMassFormedDiskID     , &
            &                                                       +diskHost              %floatRank0MetaPropertyGet(self%timeStellarMassFormedDiskID    )  &
            &                                                       +spheroidHost          %floatRank0MetaPropertyGet(self%timeStellarMassFormedSpheroidID)  &
            &                                                       +                                                      timeNuclearStarCluster            &
            &                                                      )
       call spheroidHost                 %floatRank0MetaPropertySet(                                                  self%timeStellarMassFormedSpheroidID , &
            &                                                       +0.0d0                                                                                   &
            &                                                      )
       if (haveNuclearStarClusterHost)                                                                                                                       &
            & call nuclearStarClusterHost%floatRank0MetaPropertySet(                                                  self%timeStellarMassFormedNSCID      , &
            &                                                       +0.0d0                                                                                   & 
            &                                                      )
       call diskHost                     %floatRank0MetaPropertySet(                                                  self%    stellarMassFormedDiskID     , &
            &                                                       +diskHost              %floatRank0MetaPropertyGet(self%    stellarMassFormedDiskID    )  &
            &                                                       +spheroidHost          %floatRank0MetaPropertyGet(self%    StellarMassFormedSpheroidID)  &
            &                                                       +                                                          massNuclearStarCluster        &
            &                                                      )
       call spheroidHost                 %floatRank0MetaPropertySet(                                                  self%    stellarMassFormedSpheroidID , &
            &                                                       +0.0d0                                                                                   &
            &                                                      )
       if (haveNuclearStarClusterHost)                                                                                                                       &
            & call nuclearStarClusterHost%floatRank0MetaPropertySet(                                                  self%    stellarMassFormedNSCID      , &
            &                                                       +0.0d0                                                                                   &
            &                                                      )
    case (destinationMergerSpheroid%ID)
       call spheroidHost                 %floatRank0MetaPropertySet(                                                  self%timeStellarMassFormedSpheroidID , &
            &                                                       +diskHost              %floatRank0MetaPropertyGet(self%timeStellarMassFormedDiskID    )  &
            &                                                       +spheroidHost          %floatRank0MetaPropertyGet(self%timeStellarMassFormedSpheroidID)  &
            &                                                       +                                                      timeNuclearStarCluster            &
            &                                                      )
       call diskHost                     %floatRank0MetaPropertySet(                                                  self%timeStellarMassFormedDiskID     , &
            &                                                       +0.0d0                                                                                   &
            &                                                      )
       if (haveNuclearStarClusterHost)                                                                                                                       &
            & call nuclearStarClusterHost%floatRank0MetaPropertySet(                                                  self%timeStellarMassFormedNSCID      , &
            &                                                       +0.0d0                                                                                   &
            &                                                      )
       call spheroidHost                 %floatRank0MetaPropertySet(                                                  self%    stellarMassFormedSpheroidID , &
            &                                                       +diskHost              %floatRank0MetaPropertyGet(self%    stellarMassFormedDiskID    )  &
            &                                                       +spheroidHost          %floatRank0MetaPropertyGet(self%    stellarMassFormedSpheroidID)  &
            &                                                       +                                                          massNuclearStarCluster        &
            &                                                      )
       call diskHost                     %floatRank0MetaPropertySet(                                                  self%    stellarMassFormedDiskID     , &
            &                                                       +0.0d0                                                                                   &
            &                                                      )
       if (haveNuclearStarClusterHost)                                                                                                                       &
            & call nuclearStarClusterHost%floatRank0MetaPropertySet(                                                  self%    stellarMassFormedNSCID      , &
            &                                                       +0.0d0                                                                                   &
            &                                                      )
    case (destinationMergerUnmoved%ID)
       ! Do nothing.
    case default
       call Error_Report('unrecognized movesTo descriptor'//{introspection:location})
    end select
    ! Move the star formation rates from secondary to primary.
    if (haveNuclearStarCluster) then
       massNuclearStarCluster=nuclearStarCluster%floatRank0MetaPropertyGet(self%    stellarMassFormedNSCID)
       timeNuclearStarCluster=nuclearStarCluster%floatRank0MetaPropertyGet(self%timeStellarMassFormedNSCID)
    else
       massNuclearStarCluster=0.0d0
       timeNuclearStarCluster=0.0d0
    end if
    select case (destinationStarsSatellite%ID)
    case (destinationMergerDisk    %ID)
       call diskHost    %floatRank0MetaPropertySet(                                              self%timeStellarMassFormedDiskID     , &
            &                                      +diskHost          %floatRank0MetaPropertyGet(self%timeStellarMassFormedDiskID    )  &
            &                                      +disk              %floatRank0MetaPropertyGet(self%timeStellarMassFormedDiskID    )  &
            &                                      +spheroid          %floatRank0MetaPropertyGet(self%timeStellarMassFormedSpheroidID)  &
            &                                      +                                                  timeNuclearStarCluster            &
            &                                     )
       call diskHost    %floatRank0MetaPropertySet(                                              self%    stellarMassFormedDiskID     , &
            &                                      +diskHost          %floatRank0MetaPropertyGet(self%    stellarMassFormedDiskID    )  &
            &                                      +disk              %floatRank0MetaPropertyGet(self%    stellarMassFormedDiskID    )  &
            &                                      +spheroid          %floatRank0MetaPropertyGet(self%    stellarMassFormedSpheroidID)  &
            &                                      +                                                      massNuclearStarCluster        &
            &                                     )
    case (destinationMergerSpheroid%ID)
      call spheroidHost%floatRank0MetaPropertySet(                                               self%timeStellarMassFormedSpheroidID , &
            &                                      +spheroidHost      %floatRank0MetaPropertyGet(self%timeStellarMassFormedSpheroidID)  &
            &                                      +disk              %floatRank0MetaPropertyGet(self%timeStellarMassFormedDiskID    )  &
            &                                      +spheroid          %floatRank0MetaPropertyGet(self%timeStellarMassFormedSpheroidID)  &
            &                                      +                                                  timeNuclearStarCluster            &
            &                                     )
       call spheroidHost%floatRank0MetaPropertySet(                                              self%    stellarMassFormedSpheroidID , &
            &                                      +spheroidHost      %floatRank0MetaPropertyGet(self%    stellarMassFormedSpheroidID)  &
            &                                      +disk              %floatRank0MetaPropertyGet(self%    stellarMassFormedDiskID    )  &
            &                                      +spheroid          %floatRank0MetaPropertyGet(self%    stellarMassFormedSpheroidID)  &
            &                                      +                                                      massNuclearStarCluster        &
            &                                     )
    case default
       call Error_Report('unrecognized movesTo descriptor'//{introspection:location})
    end select
    ! Zero rates in the secondary,
    call    disk                     %floatRank0MetaPropertySet(                                        self%timeStellarMassFormedDiskID     , &
         &                                                      +0.0d0                                                                         &
         &                                                     )
    call    spheroid                 %floatRank0MetaPropertySet(                                        self%timeStellarMassFormedSpheroidID , &
         &                                                      +0.0d0                                                                         &
         &                                                     )
    if (haveNuclearStarCluster)                                                                                                                & 
         & call    nuclearStarCluster%floatRank0MetaPropertySet(                                        self%timeStellarMassFormedNSCID      , &
         &                                                      +0.0d0                                                                         &
         &                                                     )
    call    disk                     %floatRank0MetaPropertySet(                                        self%    stellarMassFormedDiskID     , &
         &                                                      +0.0d0                                                                         &
         &                                                     )
    call    spheroid                 %floatRank0MetaPropertySet(                                        self%    stellarMassFormedSpheroidID , &
         &                                                      +0.0d0                                                                         &
         &                                                     )
    if (haveNuclearStarCluster)                                                                                                                &
         & call    nuclearStarCluster%floatRank0MetaPropertySet(                                        self%    stellarMassFormedNSCID      , &
         &                                                      +0.0d0                                                                         &
         &                                                     )
    return
  end subroutine agesStellarMassWeightedGalaxiesMerge
  
