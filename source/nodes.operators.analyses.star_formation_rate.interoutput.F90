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
  Implements a node operator class that tracks the mean star formation rate between successive outputs.
  !!}

  use :: Output_Times                              , only : outputTimesClass
  use :: Satellite_Merging_Mass_Movements          , only : mergerMassMovementsClass
  use :: Star_Formation_Rates_Disks                , only : starFormationRateDisksClass
  use :: Star_Formation_Rates_Spheroids            , only : starFormationRateSpheroidsClass
  use :: Star_Formation_Rates_Nuclear_Star_Clusters, only : starFormationRateNuclearStarClustersClass

  !![
  <nodeOperator name="nodeOperatorStarFormationRateInterOutput">
    <description>
      A node operator class that tracks the mean star formation rate between successive outputs. Intended to be paired with the
      \refClass{nodePropertyExtractorStarFormationRateInterOutput} class to extract those rates for output.
    </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorStarFormationRateInterOutput
     !!{
     A node operator class that tracks the mean star formation rate between successive outputs.
     !!}
     private
     class  (outputTimesClass                         ), pointer :: outputTimes_                          => null()
     class  (mergerMassMovementsClass                 ), pointer :: mergerMassMovements_                  => null()
     class  (starFormationRateDisksClass              ), pointer :: starFormationRateDisks_               => null()
     class  (starFormationRateSpheroidsClass          ), pointer :: starFormationRateSpheroids_           => null()
     class  (starFormationRateNuclearStarClustersClass), pointer :: starFormationRateNuclearStarClusters_ => null()
     integer                                                     :: starFormationRateDiskInterOutputID             , starFormationRateSpheroidInterOutputID, &
          &                                                         starFormationRateNSCInterOutputID              , starFormationRateInterOutputNextID
   contains
     final     ::                                starFormationRateInterOutputDestructor
     procedure :: galaxiesMerge               => starFormationRateInterOutputGalaxiesMerge
     procedure :: differentialEvolutionPre    => starFormationRateInterOutputDifferentialEvolutionPre
     procedure :: differentialEvolutionScales => starFormationRateInterOutputDifferentialEvolutionScales
     procedure :: differentialEvolution       => starFormationRateInterOutputDifferentialEvolution
  end type nodeOperatorStarFormationRateInterOutput
  
  interface nodeOperatorStarFormationRateInterOutput
     !!{
     Constructors for the {\normalfont \ttfamily starFormationRateInterOutput} node operator class.
     !!}
     module procedure starFormationRateInterOutputConstructorParameters
     module procedure starFormationRateInterOutputConstructorInternal
  end interface nodeOperatorStarFormationRateInterOutput
  
contains

  function starFormationRateInterOutputConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily starFormationRateInterOutput} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorStarFormationRateInterOutput )                :: self
    type (inputParameters                          ), intent(inout) :: parameters
    class(outputTimesClass                         ), pointer       :: outputTimes_
    class(mergerMassMovementsClass                 ), pointer       :: mergerMassMovements_
    class(starFormationRateDisksClass              ), pointer       :: starFormationRateDisks_
    class(starFormationRateSpheroidsClass          ), pointer       :: starFormationRateSpheroids_
    class(starFormationRateNuclearStarClustersClass), pointer       :: starFormationRateNuclearStarClusters_

    !![
    <objectBuilder class="outputTimes"                          name="outputTimes_"                          source="parameters"/>
    <objectBuilder class="mergerMassMovements"                  name="mergerMassMovements_"                  source="parameters"/>
    <objectBuilder class="starFormationRateDisks"               name="starFormationRateDisks_"               source="parameters"/>
    <objectBuilder class="starFormationRateSpheroids"           name="starFormationRateSpheroids_"           source="parameters"/>
    <objectBuilder class="starFormationRateNuclearStarClusters" name="starFormationRateNuclearStarClusters_" source="parameters"/>
    !!]
    self=nodeOperatorStarFormationRateInterOutput(outputTimes_,mergerMassMovements_,starFormationRateDisks_,starFormationRateSpheroids_,starFormationRateNuclearStarClusters_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="outputTimes_"                         />
    <objectDestructor name="mergerMassMovements_"                 />
    <objectDestructor name="starFormationRateDisks_"              />
    <objectDestructor name="starFormationRateSpheroids_"          />
    <objectDestructor name="starFormationRateNuclearStarClusters_"/>
    !!]
    return
  end function starFormationRateInterOutputConstructorParameters

  function starFormationRateInterOutputConstructorInternal(outputTimes_,mergerMassMovements_,starFormationRateDisks_,starFormationRateSpheroids_,starFormationRateNuclearStarClusters_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily starFormationRateInterOutput} node operator class.
    !!}
    implicit none
    type (nodeOperatorStarFormationRateInterOutput )                        :: self
    class(outputTimesClass                         ), intent(in   ), target :: outputTimes_
    class(mergerMassMovementsClass                 ), intent(in   ), target :: mergerMassMovements_
    class(starFormationRateDisksClass              ), intent(in   ), target :: starFormationRateDisks_
    class(starFormationRateSpheroidsClass          ), intent(in   ), target :: starFormationRateSpheroids_
    class(starFormationRateNuclearStarClustersClass), intent(in   ), target :: starFormationRateNuclearStarClusters_

    !![
    <constructorAssign variables="*outputTimes_, *mergerMassMovements_, *starFormationRateDisks_, *starFormationRateSpheroids_, *starFormationRateNuclearStarClusters_"/>
    !!]
    
    !![
    <addMetaProperty component="basic"    name="starFormationRateInterOutputNext"     id="self%starFormationRateInterOutputNextID"     isEvolvable="no"  isCreator="yes"/>
    <addMetaProperty component="disk"     name="starFormationRateDiskInterOutput"     id="self%starFormationRateDiskInterOutputID"     isEvolvable="yes" isCreator="yes"/>
    <addMetaProperty component="spheroid" name="starFormationRateSpheroidInterOutput" id="self%starFormationRateSpheroidInterOutputID" isEvolvable="yes" isCreator="yes"/>
    <addMetaProperty component="NSC"      name="starFormationRateNSCInterOutput"      id="self%starFormationRateNSCInterOutputID"      isEvolvable="yes" isCreator="yes"/>
    !!]
    return
  end function starFormationRateInterOutputConstructorInternal

  subroutine starFormationRateInterOutputDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily starFormationRateInterOutput} node operator class.
    !!}
    implicit none
    type(nodeOperatorStarFormationRateInterOutput), intent(inout) :: self

    !![
    <objectDestructor name="self%outputTimes_"                         />
    <objectDestructor name="self%mergerMassMovements_"                 />
    <objectDestructor name="self%starFormationRateDisks_"              />
    <objectDestructor name="self%starFormationRateSpheroids_"          />
    <objectDestructor name="self%starFormationRateNuclearStarClusters_"/>
    !!]
    return
  end subroutine starFormationRateInterOutputDestructor

  subroutine starFormationRateInterOutputDifferentialEvolutionPre(self,node)
    !!{
    Reset the inter-output mean star formation rate to zero when starting a new output interval.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentDisk, nodeComponentSpheroid, nodeComponentNSC
    implicit none
    class(nodeOperatorStarFormationRateInterOutput), intent(inout) :: self
    type (treeNode                                ), intent(inout) :: node
    class(nodeComponentBasic                      ), pointer       :: basic
    class(nodeComponentDisk                       ), pointer       :: disk
    class(nodeComponentSpheroid                   ), pointer       :: spheroid
    class(nodeComponentNSC                        ), pointer       :: nuclearStarCluster


    basic => node%basic()
    if (self%outputTimes_%timeNext(basic%time()) > basic%floatRank0MetaPropertyGet(self%starFormationRateInterOutputNextID)) then
       call basic                %floatRank0MetaPropertySet(self%starFormationRateInterOutputNextID    ,self%outputTimes_%timeNext(basic%time()))
       disk               =>  node%disk    ()
       spheroid           =>  node%spheroid()
       nuclearStarCluster =>  node%NSC     ()
       select type (disk    )
       type is (nodeComponentDisk     )
          ! Disk does not yet exist - nothing to do here.
       class default
          call disk              %floatRank0MetaPropertySet(self%starFormationRateDiskInterOutputID    ,0.0d0                                   )
       end select
       select type (spheroid)
       type is (nodeComponentSpheroid )
          ! Spheroid does not yet exist - nothing to do here.
       class default
          call spheroid          %floatRank0MetaPropertySet(self%starFormationRateSpheroidInterOutputID,0.0d0                                   )
       end select
       select type (nuclearStarCluster)
       type is (nodeComponentNSC)
          ! Spheroid does not yet exist - nothing to do here.
       class default
          call nuclearStarCluster%floatRank0MetaPropertySet(self%starFormationRateNSCInterOutputID     ,0.0d0                                   )
       end select
    end if
    return
  end subroutine starFormationRateInterOutputDifferentialEvolutionPre

  subroutine starFormationRateInterOutputDifferentialEvolutionScales(self,node)
    !!{
    Set absolute ODE solver scale for the mass cooled out of the \gls{cgm}.    
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDisk, nodeComponentSpheroid, nodeComponentNSC
    implicit none
    class           (nodeOperatorStarFormationRateInterOutput), intent(inout) :: self
    type            (treeNode                                ), intent(inout) :: node
    class           (nodeComponentDisk                       ), pointer       :: disk
    class           (nodeComponentSpheroid                   ), pointer       :: spheroid
    class           (nodeComponentNSC                        ), pointer       :: nuclearStarCluster
    double precision                                          , parameter     :: massMinimum       =1.0d0
    double precision                                          , parameter     :: timeScale         =1.0d0
    double precision                                                          :: mass
    
    disk               =>  node%disk       ()
    spheroid           =>  node%spheroid   ()
    nuclearStarCluster =>  node%NSC        ()
    mass               =  +disk%massGas    ()+spheroid%massGas    ()+nuclearStarCluster%massGas    () &
         &                +disk%massStellar()+spheroid%massStellar()+nuclearStarCluster%massStellar() 

    select type (disk               )
    type is (nodeComponentDisk     )
       ! Disk does not yet exist - nothing to do here.
    class default
       call disk              %floatRank0MetaPropertyScale(self%starFormationRateDiskInterOutputID    ,max(mass,massMinimum)/timeScale)
    end select
    select type (spheroid           )
    type is (nodeComponentSpheroid )
       ! Spheroid does not yet exist - nothing to do here.
    class default
       call spheroid          %floatRank0MetaPropertyScale(self%starFormationRateSpheroidInterOutputID,max(mass,massMinimum)/timeScale)
    end select
    select type (nuclearStarCluster)
    type is (nodeComponentNSC     )
       ! Nuclear star cluster does not yet exist - nothing to do here.
    class default
       call nuclearStarCluster%floatRank0MetaPropertyScale(self%starFormationRateNSCInterOutputID     ,max(mass,massMinimum)/timeScale)
    end select
    return
  end subroutine starFormationRateInterOutputDifferentialEvolutionScales
  
  subroutine starFormationRateInterOutputDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Accumulate the mean rate of star formation between outputs.
    !!}
    use :: Galacticus_Nodes, only : propertyInactive, nodeComponentBasic, nodeComponentDisk, nodeComponentSpheroid, &
       &                            nodeComponentNSC
    implicit none
    class           (nodeOperatorStarFormationRateInterOutput), intent(inout), target  :: self
    type            (treeNode                                ), intent(inout), target  :: node
    logical                                                   , intent(inout)          :: interrupt
    procedure       (interruptTask                           ), intent(inout), pointer :: functionInterrupt
    integer                                                   , intent(in   )          :: propertyType
    class           (nodeComponentBasic                      ), pointer                :: basic
    class           (nodeComponentDisk                       ), pointer                :: disk
    class           (nodeComponentSpheroid                   ), pointer                :: spheroid
    class           (nodeComponentNSC                        ), pointer                :: nuclearStarCluster
    double precision                                                                   :: timeInterval
    
    ! Return immediately if inactive variables are requested.
    if (propertyInactive(propertyType)) return
    ! Get required components.
    basic              => node%basic   ()
    disk               => node%disk    ()
    spheroid           => node%spheroid()
    nuclearStarCluster => node%NSC     ()
    ! Find the time interval for this inter-output.
    timeInterval=+          basic%floatRank0MetaPropertyGet             (self %starFormationRateInterOutputNextID  )  &
         &       -max(0.0d0,self %outputTimes_             %timePrevious(basic%time                              ()))
    if (timeInterval <= 0.0d0) return
    ! Accumulate rates.
    select type (disk              )
    type is (nodeComponentDisk    )
       ! Disk does not yet exist - nothing to do here.
    class default
       call disk              %floatRank0MetaPropertyRate(self%starFormationRateDiskInterOutputID    ,self%starFormationRateDisks_              %rate(node)/timeInterval)
    end select
    select type (spheroid          )
    type is (nodeComponentSpheroid)
       ! Spheroid does not yet exist - nothing to do here.
    class default
       call spheroid          %floatRank0MetaPropertyRate(self%starFormationRateSpheroidInterOutputID,self%starFormationRateSpheroids_          %rate(node)/timeInterval)
    end select
    select type (nuclearStarCluster)
    type is (nodeComponentNSC     )
       ! Nuclear star cluster does not yet exist - nothing to do here.
    class default
       call nuclearStarCluster%floatRank0MetaPropertyRate(self%starFormationRateNSCInterOutputID     ,self%starFormationRateNuclearStarClusters_%rate(node)/timeInterval)
    end select      
    return
  end subroutine starFormationRateInterOutputDifferentialEvolution

  subroutine starFormationRateInterOutputGalaxiesMerge(self,node)
    !!{
    Combine integrals of star formation rate when galaxies merge.
    !!}
    use :: Error                           , only : Error_Report
    use :: Galacticus_Nodes                , only : nodeComponentDisk      , nodeComponentSpheroid    , nodeComponentNSC
    use :: Satellite_Merging_Mass_Movements, only : destinationMergerDisk  , destinationMergerSpheroid, destinationMergerUnmoved, enumerationDestinationMergerType
    implicit none
    class  (nodeOperatorStarFormationRateInterOutput), intent(inout) :: self
    type   (treeNode                                ), intent(inout) :: node
    type   (treeNode                                ), pointer       :: nodeHost
    class  (nodeComponentDisk                       ), pointer       :: disk                   , diskHost
    class  (nodeComponentSpheroid                   ), pointer       :: spheroid               , spheroidHost
    class  (nodeComponentNSC                        ), pointer       :: nuclearStarCluster     , nuclearStarClusterHost

    type   (enumerationDestinationMergerType        )                :: destinationGasSatellite, destinationStarsSatellite, &
         &                                                              destinationGasHost     , destinationStarsHost
    logical                                                          :: mergerIsMajor

    ! Find the node to merge with.
    nodeHost               => node    %mergesWith(                 )
    disk                   => node    %disk      (autoCreate=.true.)
    spheroid               => node    %spheroid  (autoCreate=.true.)
    nuclearStarCluster     => node    %NSC       (autoCreate=.true.)
    diskHost               => nodeHost%disk      (autoCreate=.true.)
    spheroidHost           => nodeHost%spheroid  (autoCreate=.true.)
    nuclearStarClusterHost => nodeHost%NSC       (autoCreate=.true.)
    ! Get mass movement descriptors.
    call self%mergerMassMovements_%get(node,destinationGasSatellite,destinationStarsSatellite,destinationGasHost,destinationStarsHost,mergerIsMajor)
    ! Move the star formation rates from secondary to primary.
    select case (destinationStarsSatellite%ID)
    case (destinationMergerDisk    %ID)
       call     diskHost%floatRank0MetaPropertySet(                                                                                                             &
            &                                       self                  %starFormationRateDiskInterOutputID                                                 , &
            &                                      +    diskHost          %floatRank0MetaPropertyGet             (self%starFormationRateDiskInterOutputID    )  &
            &                                      +    disk              %floatRank0MetaPropertyGet             (self%starFormationRateDiskInterOutputID    )  &
            &                                      +    spheroid          %floatRank0MetaPropertyGet             (self%starFormationRateSpheroidInterOutputID)  &
            &                                      +    nuclearStarCluster%floatRank0MetaPropertyGet             (self%starFormationRateNSCInterOutputID     )  &
            &                                     )

    case (destinationMergerSpheroid%ID)
       call spheroidHost%floatRank0MetaPropertySet(                                                                                                             &
            &                                       self                  %starFormationRateSpheroidInterOutputID                                             , &
            &                                      +spheroidHost          %floatRank0MetaPropertyGet             (self%starFormationRateSpheroidInterOutputID)  &
            &                                      +    disk              %floatRank0MetaPropertyGet             (self%starFormationRateDiskInterOutputID    )  &
            &                                      +    spheroid          %floatRank0MetaPropertyGet             (self%starFormationRateSpheroidInterOutputID)  &
            &                                      +    nuclearStarCluster%floatRank0MetaPropertyGet             (self%starFormationRateNSCInterOutputID     )  &
            &                                     )
    case default
       call Error_Report('unrecognized movesTo descriptor'//{introspection:location})
    end select
    ! Zero rates in the secondary,
    call    disk              %floatRank0MetaPropertySet(                                                                                                       &
         &                                                self        %starFormationRateDiskInterOutputID                                                     , &
         &                                               +0.0d0                                                                                                 &
         &                                              )
    call    spheroid          %floatRank0MetaPropertySet(                                                                                                       &
         &                                                self        %starFormationRateSpheroidInterOutputID                                                 , &
         &                                               +0.0d0                                                                                                 &
         &                                              )
    call    nuclearStarCluster%floatRank0MetaPropertySet(                                                                                                       &
         &                                                self        %starFormationRateNSCInterOutputID                                                      , &
         &                                               +0.0d0                                                                                                 &
         &                                              )  
    ! Move star formation rates within the host if necessary.
    select case (destinationStarsHost%ID)
    case (destinationMergerDisk    %ID)
       call diskHost              %floatRank0MetaPropertySet(                                                                                                                       &
            &                                                 self                  %starFormationRateDiskInterOutputID                                                           , &
            &                                                +              diskHost%floatRank0MetaPropertyGet                       (self%starFormationRateDiskInterOutputID    )  &
            &                                                +          spheroidHost%floatRank0MetaPropertyGet                       (self%starFormationRateSpheroidInterOutputID)  &
            &                                                +nuclearStarClusterHost%floatRank0MetaPropertyGet                       (self%starFormationRateNSCInterOutputID     )  &
            &                                               )
       call spheroidHost          %floatRank0MetaPropertySet(                                                                                                                       &
            &                                                 self                  %starFormationRateSpheroidInterOutputID                                                       , &
            &                                                +0.0d0                                                                                                                 &
            &                                               )
       call nuclearStarClusterHost%floatRank0MetaPropertySet(                                                                                                                       &
            &                                                self                   %starFormationRateNSCInterOutputID                                                            , &
            &                                                +0.0d0                                                                                                                 &
            &                                               )
    case (destinationMergerSpheroid%ID)
       call spheroidHost          %floatRank0MetaPropertySet(                                                                                                                       &
            &                                                 self                  %starFormationRateSpheroidInterOutputID                                                       , &
            &                                                +          spheroidHost%floatRank0MetaPropertyGet                       (self%starFormationRateSpheroidInterOutputID)  &
            &                                                +              diskHost%floatRank0MetaPropertyGet                       (self%starFormationRateDiskInterOutputID    )  &
            &                                                +nuclearStarClusterHost%floatRank0MetaPropertyGet                       (self%starFormationRateNSCInterOutputID     )  &
            &                                               )
       call diskHost              %floatRank0MetaPropertySet(                                                                                                                       &
            &                                                 self                  %starFormationRateDiskInterOutputID                                                           , &
            &                                                +0.0d0                                                                                                                 &
            &                                               )
       call nuclearStarClusterHost%floatRank0MetaPropertySet(                                                                                                                       &
            &                                                self                   %starFormationRateNSCInterOutputID                                                            , & 
            &                                                +0.0d0                                                                                                                 &
            &                                               )
    case (destinationMergerUnmoved%ID)
       ! Do nothing.
    case default
       call Error_Report('unrecognized movesTo descriptor'//{introspection:location})
    end select
    return
  end subroutine starFormationRateInterOutputGalaxiesMerge
