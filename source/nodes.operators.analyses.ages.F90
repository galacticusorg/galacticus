!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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
  Implements a node operator class that computes the stellar mass-weighted ages of disk and spheroid components.
  !!}
  
  use :: Star_Formation_Rates_Disks      , only : starFormationRateDisksClass
  use :: Star_Formation_Rates_Spheroids  , only : starFormationRateSpheroidsClass
  use :: Satellite_Merging_Mass_Movements, only : mergerMassMovementsClass

  !![
  <nodeOperator name="nodeOperatorAgesStellarMassWeighted">
    <description>
      A node operator class that computes the stellar mass-weighted ages of disk and spheroid components. Intended to be paired
      with the \refClass{nodePropertyExtractorAgesStellarMassWeighted} class to extract those ages for output.
    </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorAgesStellarMassWeighted
     !!{
     A node operator class that computes the stellar mass-weighted ages of disk and spheroid components.
     !!}
     private
     class  (starFormationRateDisksClass    ), pointer :: starFormationRateDisks_     => null()
     class  (starFormationRateSpheroidsClass), pointer :: starFormationRateSpheroids_ => null()
     class  (mergerMassMovementsClass       ), pointer :: mergerMassMovements_        => null()
     integer                                           :: stellarMassFormedDiskID              , timeStellarMassFormedDiskID    , &
          &                                               stellarMassFormedSpheroidID          , timeStellarMassFormedSpheroidID
   contains
     final     ::                                   agesStellarMassWeightedDestructor
     procedure :: differentialEvolutionScales    => agesStellarMassWeightedDifferentialEvolutionScales
     procedure :: differentialEvolutionInactives => agesStellarMassWeightedDifferentialEvolutionInactives
     procedure :: differentialEvolution          => agesStellarMassWeightedDifferentialEvolution
     procedure :: galaxiesMerge                  => agesStellarMassWeightedGalaxiesMerge
  end type nodeOperatorAgesStellarMassWeighted
  
  interface nodeOperatorAgesStellarMassWeighted
     !!{
     Constructors for the {\normalfont \ttfamily agesStellarMassWeighted} node operator class.
     !!}
     module procedure agesStellarMassWeightedConstructorParameters
     module procedure agesStellarMassWeightedConstructorInternal
  end interface nodeOperatorAgesStellarMassWeighted
  
contains

  function agesStellarMassWeightedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily agesStellarMassWeighted} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorAgesStellarMassWeighted)                :: self
    type (inputParameters                    ), intent(inout) :: parameters
    class(starFormationRateDisksClass        ), pointer       :: starFormationRateDisks_
    class(starFormationRateSpheroidsClass    ), pointer       :: starFormationRateSpheroids_
    class(mergerMassMovementsClass           ), pointer       :: mergerMassMovements_
    
    !![
    <objectBuilder class="starFormationRateDisks"     name="starFormationRateDisks_"     source="parameters"/>
    <objectBuilder class="starFormationRateSpheroids" name="starFormationRateSpheroids_" source="parameters"/>
    <objectBuilder class="mergerMassMovements"        name="mergerMassMovements_"        source="parameters"/>
    !!]
    self=nodeOperatorAgesStellarMassWeighted(starFormationRateDisks_,starFormationRateSpheroids_,mergerMassMovements_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="starFormationRateDisks_"    />
    <objectDestructor name="starFormationRateSpheroids_"/>
    <objectDestructor name="mergerMassMovements_"       />
    !!]
    return
  end function agesStellarMassWeightedConstructorParameters

  function agesStellarMassWeightedConstructorInternal(starFormationRateDisks_,starFormationRateSpheroids_,mergerMassMovements_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily agesStellarMassWeighted} node operator class.
    !!}
    use :: Galacticus_Nodes, only : defaultDiskComponent, defaultSpheroidComponent
    implicit none
    type (nodeOperatorAgesStellarMassWeighted)                        :: self
    class(starFormationRateDisksClass        ), intent(in   ), target :: starFormationRateDisks_
    class(starFormationRateSpheroidsClass    ), intent(in   ), target :: starFormationRateSpheroids_
    class(mergerMassMovementsClass           ), intent(in   ), target :: mergerMassMovements_
    !![
    <constructorAssign variables="*starFormationRateDisks_, *starFormationRateSpheroids_, *mergerMassMovements_"/>
    !!]
    
    self%stellarMassFormedDiskID        =defaultDiskComponent    %addMetaProperty(var_str('agesStellarMassFormedDisk'        ),'disk:stellarMassFormed'        ,isEvolvable=.true.)  
    self%timeStellarMassFormedDiskID    =defaultDiskComponent    %addMetaProperty(var_str('agesTimeStellarMassFormedDisk'    ),'disk:timeStellarMassFormed'    ,isEvolvable=.true.)  
    self%stellarMassFormedSpheroidID    =defaultSpheroidComponent%addMetaProperty(var_str('agesStellarMassFormedSpheroid'    ),'spheroid:stellarMassFormed'    ,isEvolvable=.true.)  
    self%timeStellarMassFormedSpheroidID=defaultSpheroidComponent%addMetaProperty(var_str('agesTimeStellarMassFormedSpheroid'),'spheroid:timeStellarMassFormed',isEvolvable=.true.)  
    return
  end function agesStellarMassWeightedConstructorInternal

  subroutine agesStellarMassWeightedDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily agesStellarMassWeighted} node operator class.
    !!}
    implicit none
    type(nodeOperatorAgesStellarMassWeighted), intent(inout) :: self

    !![
    <objectDestructor name="self%starFormationRateDisks_"    />
    <objectDestructor name="self%starFormationRateSpheroids_"/>
    <objectDestructor name="self%mergerMassMovements_"       />
    !!]
    return
  end subroutine agesStellarMassWeightedDestructor

  subroutine agesStellarMassWeightedDifferentialEvolutionInactives(self,node)
    !!{
    Mark disk and spheroid age integrals as inactive for ODE solving.    
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDisk, nodeComponentSpheroid
    implicit none
    class(nodeOperatorAgesStellarMassWeighted), intent(inout) :: self
    type (treeNode                           ), intent(inout) :: node
    class(nodeComponentDisk                  ), pointer       :: disk
    class(nodeComponentSpheroid              ), pointer       :: spheroid
    
    ! Get disk and spheroid components.
    disk     => node%disk    ()
    spheroid => node%spheroid()
    ! Mark as inactive.
    select type (disk    )
    type is (nodeComponentDisk    )
       ! Disk does not yet exist - nothing to do here.
    class default
       call disk    %metaPropertyInactive(self%    stellarMassFormedDiskID    )
       call disk    %metaPropertyInactive(self%timeStellarMassFormedDiskID    )
    end select
    select type (spheroid)
    type is (nodeComponentSpheroid)
       ! Spheroid does not yet exist - nothing to do here.
    class default
       call spheroid%metaPropertyInactive(self%    stellarMassFormedSpheroidID)
       call spheroid%metaPropertyInactive(self%timeStellarMassFormedSpheroidID)
    end select
     return
  end subroutine agesStellarMassWeightedDifferentialEvolutionInactives
  
  subroutine agesStellarMassWeightedDifferentialEvolutionScales(self,node)
    !!{
    Set absolute ODE solver scale for the unweighted and time-weighted stellar masses.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDisk, nodeComponentSpheroid
    implicit none
    class           (nodeOperatorAgesStellarMassWeighted), intent(inout) :: self
    type            (treeNode                           ), intent(inout) :: node
    class           (nodeComponentDisk                  ), pointer       :: disk
    class           (nodeComponentSpheroid              ), pointer       :: spheroid
    double precision                                     , parameter     :: massMinimum=1.0d+0
    double precision                                     , parameter     :: timeScale  =1.0d-3
    double precision                                                     :: mass
    
    ! Get disk and spheroid components.
    disk     =>  node%disk       ()
    spheroid =>  node%spheroid   ()
    ! Set scale for masses.
    mass     =  +disk%massGas    ()+spheroid%massGas    () &
         &      +disk%massStellar()+spheroid%massStellar()
    ! Set scales.
    select type (disk    )
    type is (nodeComponentDisk    )
       ! Disk does not yet exist - nothing to do here.
    class default
       call disk    %metaPropertyScale(self%    stellarMassFormedDiskID    ,max(mass,massMinimum)          )
       call disk    %metaPropertyScale(self%timeStellarMassFormedDiskID    ,max(mass,massMinimum)*timeScale)
    end select
    select type (spheroid)
    type is (nodeComponentSpheroid)
       ! Spheroid does not yet exist - nothing to do here.
    class default
       call spheroid%metaPropertyScale(self%    stellarMassFormedSpheroidID,max(mass,massMinimum)          )
       call spheroid%metaPropertyScale(self%timeStellarMassFormedSpheroidID,max(mass,massMinimum)*timeScale)
    end select
    return
  end subroutine agesStellarMassWeightedDifferentialEvolutionScales
  
  subroutine agesStellarMassWeightedDifferentialEvolution(self,node,interrupt,functionInterrupt,propertyType)
    !!{
    Integrates unweighted and time-weighted star formation rates in disk and spheroid components.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDisk, nodeComponentSpheroid, nodeComponentBasic, propertyTypeActive
    implicit none
    class           (nodeOperatorAgesStellarMassWeighted), intent(inout), target  :: self
    type            (treeNode                           ), intent(inout)          :: node
    logical                                              , intent(inout)          :: interrupt
    procedure       (interruptTask                      ), intent(inout), pointer :: functionInterrupt
    integer                                              , intent(in   )          :: propertyType
    class           (nodeComponentBasic                 )               , pointer :: basic
    class           (nodeComponentDisk                  )               , pointer :: disk
    class           (nodeComponentSpheroid              )               , pointer :: spheroid
    double precision                                                              :: rateStarFormationDisk, rateStarFormationSpheroid, &
         &                                                                           time
    !$GLC attributes unused :: interrupt, functionInterrupt, propertyType

    ! Return immediately if active variables are requested.
    if (propertyType == propertyTypeActive) return
    ! Get the star formation rates.
    disk                      => node                            %disk    (    )
    spheroid                  => node                            %spheroid(    )
    rateStarFormationDisk     =  self%starFormationRateDisks_    %rate    (node)
    rateStarFormationSpheroid =  self%starFormationRateSpheroids_%rate    (node)
    ! Find the current cosmic time.
    basic => node %basic()
    time  =  basic%time ()
    ! Accumulate rates.
    select type (disk    )
    type is (nodeComponentDisk    )
       ! Disk does not yet exist - nothing to do here.
    class default
       call disk    %metaPropertyRate(self%timeStellarMassFormedDiskID    ,rateStarFormationDisk    *time)
       call disk    %metaPropertyRate(self%    stellarMassFormedDiskID    ,rateStarFormationDisk         )
    end select
    select type (spheroid)
    type is (nodeComponentSpheroid)
       ! Spheroid does not yet exist - nothing to do here.
    class default
       call spheroid%metaPropertyRate(self%timeStellarMassFormedSpheroidID,rateStarFormationSpheroid*time)
       call spheroid%metaPropertyRate(self%    stellarMassFormedSpheroidID,rateStarFormationSpheroid     )
    end select
    return
  end subroutine agesStellarMassWeightedDifferentialEvolution

  subroutine agesStellarMassWeightedGalaxiesMerge(self,node)
    !!{
    Combine integrals of star formation rate when galaxies merge.
    !!}
    use :: Galacticus_Error                , only : Galacticus_Error_Report
    use :: Galacticus_Nodes                , only : nodeComponentDisk      , nodeComponentSpheroid
    use :: Satellite_Merging_Mass_Movements, only : destinationMergerDisk  , destinationMergerSpheroid, destinationMergerUnmoved
    implicit none
    class  (nodeOperatorAgesStellarMassWeighted), intent(inout) :: self
    type   (treeNode                           ), intent(inout) :: node
    type   (treeNode                           ), pointer       :: nodeHost
    class  (nodeComponentDisk                  ), pointer       :: disk                   , diskHost
    class  (nodeComponentSpheroid              ), pointer       :: spheroid               , spheroidHost
    integer                                                     :: destinationGasSatellite, destinationStarsSatellite, &
         &                                                         destinationGasHost     , destinationStarsHost
    logical                                                     :: mergerIsMajor

    ! Find the node to merge with.
    nodeHost     => node    %mergesWith(                 )
    disk         => node    %disk      (autoCreate=.true.)
    spheroid     => node    %spheroid  (autoCreate=.true.)
    diskHost     => nodeHost%disk      (autoCreate=.true.)
    spheroidHost => nodeHost%spheroid  (autoCreate=.true.)
    ! Get mass movement descriptors.
    call self%mergerMassMovements_%get(node,destinationGasSatellite,destinationStarsSatellite,destinationGasHost,destinationStarsHost,mergerIsMajor)
    ! Move star formation rates within the host if necessary.
    select case (destinationStarsHost)
    case (destinationMergerDisk    )
       call diskHost    %metaPropertySet(                              self%timeStellarMassFormedDiskID     , &
            &                            +diskHost    %metaPropertyGet(self%timeStellarMassFormedDiskID    )  &
            &                            +spheroidHost%metaPropertyGet(self%timeStellarMassFormedSpheroidID)  &
            &                           )
       call spheroidHost%metaPropertySet(                              self%timeStellarMassFormedSpheroidID , &
            &                            +0.0d0                                                               &
            &                           )
       call diskHost    %metaPropertySet(                              self%    stellarMassFormedDiskID     , &
            &                            +diskHost    %metaPropertyGet(self%    stellarMassFormedDiskID    )  &
            &                            +spheroidHost%metaPropertyGet(self%    stellarMassFormedSpheroidID)  &
            &                           )
       call spheroidHost%metaPropertySet(                              self%    stellarMassFormedSpheroidID , &
            &                            +0.0d0                                                               &
            &                           )
    case (destinationMergerSpheroid)
       call spheroidHost%metaPropertySet(                              self%timeStellarMassFormedDiskID     , &
            &                            +diskHost    %metaPropertyGet(self%timeStellarMassFormedDiskID    )  &
            &                            +spheroidHost%metaPropertyGet(self%timeStellarMassFormedSpheroidID)  &
            &                           )
       call diskHost    %metaPropertySet(                              self%timeStellarMassFormedSpheroidID , &
            &                            +0.0d0                                                               &
            &                           )
       call spheroidHost%metaPropertySet(                              self%    stellarMassFormedDiskID     , &
            &                            +diskHost    %metaPropertyGet(self%    stellarMassFormedDiskID    )  &
            &                            +spheroidHost%metaPropertyGet(self%    stellarMassFormedSpheroidID)  &
            &                           )
       call diskHost    %metaPropertySet(                              self%    stellarMassFormedSpheroidID , &
            &                            +0.0d0                                                               &
            &                           )
    case (destinationMergerUnmoved)
       ! Do nothing.
    case default
       call Galacticus_Error_Report('unrecognized movesTo descriptor'//{introspection:location})
    end select
    ! Move the star formation rates from secondary to primary.
    select case (destinationStarsSatellite)
    case (destinationMergerDisk    )
       call diskHost    %metaPropertySet(                              self%timeStellarMassFormedDiskID     , &
            &                            +diskHost    %metaPropertyGet(self%timeStellarMassFormedDiskID    )  &
            &                            +disk        %metaPropertyGet(self%timeStellarMassFormedDiskID    )  &
            &                            +spheroid    %metaPropertyGet(self%timeStellarMassFormedSpheroidID)  &
            &                           )
       call diskHost    %metaPropertySet(                              self%    stellarMassFormedDiskID     , &
            &                            +diskHost    %metaPropertyGet(self%    stellarMassFormedDiskID    )  &
            &                            +disk        %metaPropertyGet(self%    stellarMassFormedDiskID    )  &
            &                            +spheroid    %metaPropertyGet(self%    stellarMassFormedSpheroidID)  &
            &                           )
    case (destinationMergerSpheroid)
      call spheroidHost%metaPropertySet(                              self%timeStellarMassFormedSpheroidID , &
            &                            +spheroidHost%metaPropertyGet(self%timeStellarMassFormedSpheroidID)  &
            &                            +disk        %metaPropertyGet(self%timeStellarMassFormedDiskID    )  &
            &                            +spheroid    %metaPropertyGet(self%timeStellarMassFormedSpheroidID)  &
            &                           )
       call spheroidHost%metaPropertySet(                              self%    stellarMassFormedSpheroidID , &
            &                            +spheroidHost%metaPropertyGet(self%    stellarMassFormedSpheroidID)  &
            &                            +disk        %metaPropertyGet(self%    stellarMassFormedDiskID    )  &
            &                            +spheroid    %metaPropertyGet(self%    stellarMassFormedSpheroidID)  &
            &                           )
    case default
       call Galacticus_Error_Report('unrecognized movesTo descriptor'//{introspection:location})
    end select
    ! Zero rates in the secondary,
    call    disk        %metaPropertySet(                              self%timeStellarMassFormedDiskID     , &
         &                               +0.0d0                                                               &
         &                              )
    call    spheroid    %metaPropertySet(                              self%timeStellarMassFormedSpheroidID , &
         &                               +0.0d0                                                               &
         &                              )
    call    disk        %metaPropertySet(                              self%    stellarMassFormedDiskID     , &
         &                               +0.0d0                                                               &
         &                              )
    call    spheroid    %metaPropertySet(                              self%    stellarMassFormedSpheroidID , &
         &                               +0.0d0                                                               &
         &                              )
    return
  end subroutine agesStellarMassWeightedGalaxiesMerge
  
