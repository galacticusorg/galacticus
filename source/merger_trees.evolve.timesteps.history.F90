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
Implements a merger tree evolution timestepping class which limits the step the next epoch at which to store global history.
!!}

  use :: Cosmology_Functions           , only : cosmologyFunctions             , cosmologyFunctionsClass
  use :: Numerical_Interpolation       , only : interpolator
  use :: Star_Formation_Rates_Disks    , only : starFormationRateDisksClass
  use :: Star_Formation_Rates_Spheroids, only : starFormationRateSpheroidsClass

  !![
  <mergerTreeEvolveTimestep name="mergerTreeEvolveTimestepHistory">
   <description>
    A merger tree evolution timestepping class which records and outputs volume averaged properties of the model universe as a
    function of time. Timesteps are enforced such that:
    \begin{equation}
     \Delta t \le t_{\mathrm{history},i} - t
    \end{equation}
    where $t$ is the current time, $t_{\mathrm{history},i}$ is the $i^\mathrm{th}$ time at which the global history of galaxies
    is to be output and $i$ is chosen to be the smallest $i$ such that $t_{\mathrm{history},i} &gt; t$. If there is no $i$ for
    which $t_{\mathrm{history},i} &gt; t$ this criterion is not applied. If this criterion is the limiting criterion for $\Delta
    t$ then the properties of the galaxy will be accumulated to the global history arrays at the end of the timestep.
  
    Volume-averaged properties are stored to the {\normalfont \ttfamily globalHistory} group of the output file. Currently, the
    properties stored are:
    \begin{description}
     \item[{\normalfont \ttfamily historyTime}] Cosmic time (in Gyr);
     \item[{\normalfont \ttfamily historyExpansion}] Expansion factor;
     \item[{\normalfont \ttfamily historyStarFormationRate}] Volume averaged star formation rate (in $M_\odot/$Gyr/Mpc$^3$).
     \item[{\normalfont \ttfamily historyDiskStarFormationRate}] Volume averaged star formation rate in disks (in $M_\odot/$Gyr/Mpc$^3$).
     \item[{\normalfont \ttfamily historySpheroidStarFormationRate}] Volume averaged star formation rate in spheroids (in $M_\odot/$Gyr/Mpc$^3$).
     \item[{\normalfont \ttfamily historyStellarDensity}] Volume averaged stellar mass density (in $M_\odot/$Mpc$^3$).
     \item[{\normalfont \ttfamily historyDiskStellarDensity}] Volume averaged stellar mass density in disks (in $M_\odot/$Mpc$^3$).
     \item[{\normalfont \ttfamily historySpheroidStellarDensity}] Volume averaged stellar mass density in spheroids (in $M_\odot/$Mpc$^3$).
     \item[{\normalfont \ttfamily historyGasDensity}] Volume averaged cooled gas density (in $M_\odot/$Mpc$^3$).
     \item[{\normalfont \ttfamily historyNodeDensity}] Volume averaged resolved node density (in $M_\odot/$Mpc$^3$).
    \end{description}
    Dimensionful datasets have a {\normalfont \ttfamily unitsInSI} attribute which gives their units\index{units} in the SI
    system.
   </description>
  </mergerTreeEvolveTimestep>
  !!]
  type, extends(mergerTreeEvolveTimestepClass) :: mergerTreeEvolveTimestepHistory
     !!{
     Implementation of a merger tree evolution timestepping class which limits the step the next epoch at which to store global history.
     !!}
     private
     class           (cosmologyFunctionsClass        ), pointer                   :: cosmologyFunctions_         => null()
     class           (starFormationRateDisksClass    ), pointer                   :: starFormationRateDisks_     => null()
     class           (starFormationRateSpheroidsClass), pointer                   :: starFormationRateSpheroids_ => null()
     integer                                                                      :: historyCount
     double precision                                                             :: timeBegin                            , timeEnd
     double precision                                 , allocatable, dimension(:) :: rateStarFormationDisk                , densityStellarDisk    , &
          &                                                                          expansionFactor                      , densityColdGas        , &
          &                                                                          densityHotHaloGas                    , densityNode           , &
          &                                                                          rateStarFormationSpheroid            , densityStellarSpheroid, &
          &                                                                          rateStarFormation                    , densityStellar        , &
          &                                                                          time
     type            (interpolator                   )                            :: interpolator_
   contains
     final     ::                 historyDestructor
     procedure :: timeEvolveTo => historyTimeEvolveTo
     procedure :: autoHook     => historyAutoHook
  end type mergerTreeEvolveTimestepHistory

  interface mergerTreeEvolveTimestepHistory
     !!{
     Constructors for the \refClass{mergerTreeEvolveTimestepHistory} merger tree evolution timestep class.
     !!}
     module procedure historyConstructorParameters
     module procedure historyConstructorInternal
  end interface mergerTreeEvolveTimestepHistory

contains

  function historyConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerTreeEvolveTimestepHistory} merger tree evolution timestep class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (mergerTreeEvolveTimestepHistory)                :: self
    type            (inputParameters                ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass        ), pointer       :: cosmologyFunctions_
    class           (starFormationRateDisksClass    ), pointer       :: starFormationRateDisks_
    class           (starFormationRateSpheroidsClass), pointer       :: starFormationRateSpheroids_
    integer                                                          :: historyCount
    double precision                                                 :: timeBegin                  , timeEnd, &
         &                                                              ageUniverse

    !![
    <objectBuilder class="cosmologyFunctions"         name="cosmologyFunctions_"         source="parameters"/>
    <objectBuilder class="starFormationRateDisks"     name="starFormationRateDisks_"     source="parameters"/>
    <objectBuilder class="starFormationRateSpheroids" name="starFormationRateSpheroids_" source="parameters"/>
    !!]
    ageUniverse=cosmologyFunctions_%cosmicTime(1.0d0)
    !![
    <inputParameter>
      <name>timeBegin</name>
      <defaultValue>0.05d0*ageUniverse</defaultValue>
      <description>The earliest time at which to tabulate the volume averaged history of galaxies (in Gyr).</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>timeEnd</name>
      <defaultValue>ageUniverse</defaultValue>
      <description>The latest time at which to tabulate the volume averaged history of galaxies (in Gyr).</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>historyCount</name>
      <defaultValue>30</defaultValue>
      <description>The number of steps (spaced logarithmically in cosmic time) at which to tabulate the volume averaged history of galaxies.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=mergerTreeEvolveTimestepHistory(historyCount,timeBegin,timeEnd,cosmologyFunctions_,starFormationRateDisks_,starFormationRateSpheroids_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"        />
    <objectDestructor name="starFormationRateDisks_"    />
    <objectDestructor name="starFormationRateSpheroids_"/>
    !!]
    return
  end function historyConstructorParameters

  function historyConstructorInternal(historyCount,timeBegin,timeEnd,cosmologyFunctions_,starFormationRateDisks_,starFormationRateSpheroids_) result(self)
    !!{
    Constructor for the \refClass{mergerTreeEvolveTimestepHistory} merger tree evolution timestep class which takes a parameter set as input.
    !!}
    use, intrinsic :: ISO_C_Binding    , only : c_size_t
    use            :: Numerical_Ranges , only : Make_Range   , rangeTypeLogarithmic
    implicit none
    type            (mergerTreeEvolveTimestepHistory)                        :: self
    integer                                          , intent(in   )         :: historyCount
    double precision                                 , intent(in   )         :: timeBegin                  , timeEnd
    class           (cosmologyFunctionsClass        ), intent(in   ), target :: cosmologyFunctions_
    class           (starFormationRateDisksClass    ), intent(in   ), target :: starFormationRateDisks_
    class           (starFormationRateSpheroidsClass), intent(in   ), target :: starFormationRateSpheroids_
    integer         (c_size_t                       )                        :: timeIndex
    !![
    <constructorAssign variables="historyCount, timeBegin, timeEnd, *cosmologyFunctions_, *starFormationRateDisks_, *starFormationRateSpheroids_"/>
    !!]

    ! Allocate storage arrays.
    allocate(self%time                     (self%historyCount))
    allocate(self%expansionFactor          (self%historyCount))
    allocate(self%rateStarFormation        (self%historyCount))
    allocate(self%rateStarFormationDisk    (self%historyCount))
    allocate(self%rateStarFormationSpheroid(self%historyCount))
    allocate(self%densityStellar           (self%historyCount))
    allocate(self%densityStellarDisk       (self%historyCount))
    allocate(self%densityStellarSpheroid   (self%historyCount))
    allocate(self%densityColdGas           (self%historyCount))
    allocate(self%densityHotHaloGas        (self%historyCount))
    allocate(self%densityNode              (self%historyCount))
    ! Initialize arrays.
    self%time=Make_Range(self%timeBegin,self%timeEnd,self%historyCount,rangeTypeLogarithmic)
    do timeIndex=1,self%historyCount
       self%expansionFactor(timeIndex)=self%cosmologyFunctions_%expansionFactor(self%time(timeIndex))
    end do
    self%rateStarFormation        =0.0d0
    self%rateStarFormationDisk    =0.0d0
    self%rateStarFormationSpheroid=0.0d0
    self%densityStellar           =0.0d0
    self%densityStellarDisk       =0.0d0
    self%densityStellarSpheroid   =0.0d0
    self%densityColdGas           =0.0d0
    self%densityHotHaloGas        =0.0d0
    self%densityNode              =0.0d0
    self%interpolator_            =interpolator(self%time)
    return
  end function historyConstructorInternal

  subroutine historyAutoHook(self)
    !!{
    Create a hook to the HDF5 pre-close event to allow us to finalize and write out our data.
    !!}
    use :: Events_Hooks, only : outputFileCloseEventGlobal
    implicit none
    class(mergerTreeEvolveTimestepHistory), intent(inout) :: self

    call outputFileCloseEventGlobal%attach(self,historyWrite,label='mergerTreeEvolveTimestepHistory')
    return
  end subroutine historyAutoHook

  subroutine historyDestructor(self)
    !!{
    Destructor for the \refClass{mergerTreeEvolveTimestepHistory} merger tree evolution timestep class.
    !!}
    use :: Events_Hooks, only : outputFileCloseEventGlobal
    implicit none
    type(mergerTreeEvolveTimestepHistory), intent(inout) :: self

    if (outputFileCloseEventGlobal%isAttached(self,historyWrite)) call outputFileCloseEventGlobal%detach(self,historyWrite)
    !![
    <objectDestructor name="self%cosmologyFunctions_"        />
    <objectDestructor name="self%starFormationRateDisks_"    />
    <objectDestructor name="self%starFormationRateSpheroids_"/>
    !!]
    return
  end subroutine historyDestructor

  double precision function historyTimeEvolveTo(self,timeEnd,node,task,taskSelf,report,lockNode,lockType)
    !!{
    Determine a suitable timestep for {\normalfont \ttfamily node} using the history method.
    !!}
    use            :: Evolve_To_Time_Reports, only : Evolve_To_Time_Report
    use            :: Galacticus_Nodes      , only : nodeComponentBasic   , treeNode
    use, intrinsic :: ISO_C_Binding         , only : c_size_t
    use            :: ISO_Varying_String    , only : varying_string
    implicit none
    class           (mergerTreeEvolveTimestepHistory), intent(inout), target            :: self
    double precision                                 , intent(in   )                    :: timeEnd
    type            (treeNode                       ), intent(inout), target            :: node
    procedure       (timestepTask                   ), intent(  out), pointer           :: task
    class           (*                              ), intent(  out), pointer           :: taskSelf
    logical                                          , intent(in   )                    :: report
    type            (treeNode                       ), intent(  out), pointer, optional :: lockNode
    type            (varying_string                 ), intent(  out)         , optional :: lockType
    class           (nodeComponentBasic             )               , pointer           :: basic
    integer         (c_size_t                       )                                   :: timeIndex
    double precision                                                                    :: time
    !$GLC attributes unused :: timeEnd

    ! Determine how long until next available timestep.
    basic     => node %basic()
    time      =  basic%time ()
    timeIndex =  self%interpolator_%locate(time)
    if (time < self%time(timeIndex+1)) then
       historyTimeEvolveTo =  self%time(timeIndex+1)
       task                => historyStore
       taskSelf            => self
    else
       historyTimeEvolveTo =  huge(0.0d0)
       task                => null(     )
       taskSelf            => null(     )
    end if
    if (present(lockNode)) lockNode => node
    if (present(lockType)) lockType =  "history"
    if (        report   ) call Evolve_To_Time_Report("history: ",historyTimeEvolveTo)
    return
  end function historyTimeEvolveTo

  subroutine historyStore(self,tree,node,deadlockStatus)
    !!{
    Store various properties in global arrays.
    !!}
    use            :: Galactic_Structure_Options, only : componentTypeDisk    , componentTypeHotHalo, componentTypeSpheroid, massTypeGaseous      , &
          &                                              massTypeStellar
    use            :: Error                     , only : Error_Report
    use            :: Galacticus_Nodes          , only : mergerTree           , nodeComponentBasic  , nodeComponentDisk    , nodeComponentSpheroid, &
         &                                               treeNode
    use            :: Mass_Distributions        , only : massDistributionClass
    use, intrinsic :: ISO_C_Binding             , only : c_size_t
    implicit none
    class           (*                            ), intent(inout)          :: self
    type            (mergerTree                   ), intent(in   )          :: tree
    type            (treeNode                     ), intent(inout), pointer :: node
    type            (enumerationDeadlockStatusType), intent(inout)          :: deadlockStatus
    class           (nodeComponentBasic           )               , pointer :: basic
    class           (nodeComponentDisk            )               , pointer :: disk
    class           (nodeComponentSpheroid        )               , pointer :: spheroid
    class           (massDistributionClass        )               , pointer :: massDistributionHotHalo_    , massDistributionGaseous_        , &
         &                                                                     massDistributionStellarDisk_, massDistributionStellarSpheroid_, &
         &                                                                     massDistributionStellar_
    integer         (c_size_t                     )                         :: timeIndex
    double precision                                                        :: rateStarFormationDisk        , massHotGas, &
         &                                                                     rateStarFormationSpheroid    , time
    !$GLC attributes unused :: deadlockStatus

    select type (self)
    class is (mergerTreeEvolveTimestepHistory)
       ! Get components.
       basic    => node%basic   ()
       disk     => node%disk    ()
       spheroid => node%spheroid()
       ! Get current cosmic time.
       time=basic%time()
       ! Determine how long until next available timestep.
       if (time == self%time(self%historyCount)) then
          timeIndex=self%historyCount
       else
          timeIndex=self%interpolator_%locate(time)
       end if
       ! Extract disk and spheroid star formation rates.
       rateStarFormationDisk    =self%starFormationRateDisks_    %rate(node)
       rateStarFormationSpheroid=self%starFormationRateSpheroids_%rate(node)
       ! Accumulate the properties.
       ! Star formation rate:
       self%rateStarFormation                       (timeIndex) =  +  self                   %rateStarFormation        (timeIndex                                                        ) &
            &                                                      +(+rateStarFormationDisk+rateStarFormationSpheroid)                                                                     &
            &                                                      *  tree                   %volumeWeight
       self%rateStarFormationDisk                   (timeIndex) =  +  self                   %rateStarFormationDisk    (timeIndex                                                        ) &
            &                                                      +  rateStarFormationDisk                                                                                                &
            &                                                      *  tree                   %volumeWeight
       self%rateStarFormationSpheroid               (timeIndex) =  +  self                   %rateStarFormationSpheroid(timeIndex                                                        ) &
            &                                                      +                          rateStarFormationSpheroid                                                                    &
            &                                                      *  tree                   %volumeWeight
       ! Stellar densities.
       massDistributionStellar_                                 =>    node                            %massDistribution      (                                    massType=massTypeStellar)
       massDistributionStellarDisk_                             =>    node                            %massDistribution      (componentType=componentTypeDisk    ,massType=massTypeStellar)
       massDistributionStellarSpheroid_                         =>    node                            %massDistribution      (componentType=componentTypeSpheroid,massType=massTypeStellar)
       self%densityStellar                          (timeIndex) =  +  self                            %densityStellar        (timeIndex                                                   ) &
            &                                                      +  massDistributionStellar_        %massTotal             (                                                            ) &
            &                                                      *  tree                            %volumeWeight
       self%densityStellarDisk                      (timeIndex) =  +  self                            %densityStellarDisk    (timeIndex                                                   ) &
            &                                                      +  massDistributionStellarDisk_    %massTotal             (                                                            ) &
            &                                                      *  tree                            %volumeWeight
       self%densityStellarSpheroid                  (timeIndex) =  +  self                            %densityStellarSpheroid(timeIndex                                                   ) &
            &                                                      +  massDistributionStellarSpheroid_%massTotal             (                                                            ) &
            &                                                      *  tree                            %volumeWeight
       !![
       <objectDestructor name="massDistributionStellar_"        />
       <objectDestructor name="massDistributionStellarDisk_"    />
       <objectDestructor name="massDistributionStellarSpheroid_"/>
       !!]
       ! Hot gas density.
       massDistributionHotHalo_                                 =>    node                            %massDistribution      (componentType=componentTypeHotHalo                          )
       massHotGas                                               =  +  massDistributionHotHalo_        %massTotal             (                                                            )
       self%densityHotHaloGas                       (timeIndex) =  +  self                            %densityHotHaloGas     (timeIndex                                                   ) &
            &                                                      +                                   massHotGas                                                                           &
            &                                                      *  tree                            %volumeWeight
       !![
       <objectDestructor name="massDistributionHotHalo_"/>
       !!]
       ! Galactic gas density.
       massDistributionGaseous_                                 =>    node                            %massDistribution      (                                    massType=massTypeGaseous)
       self%densityColdGas                          (timeIndex) =  +  self                            %densityColdGas        (timeIndex                                                   ) &
            &                                                      +(                                                                                                                       &
            &                                                        +massDistributionGaseous_        %massTotal             (                                                            ) &
            &                                                        -massHotGas                                                                                                            &
            &                                                       )                                                                                                                       &
            &                                                      *  tree                            %volumeWeight
       !![
       <objectDestructor name="massDistributionGaseous_"/>
       !!]
       ! Node density
       if (.not.node%isSatellite()) self%densityNode(timeIndex)=+    self                             %densityNode           (timeIndex                                                   ) &
            &                                                   +    basic                            %mass                  (                                                            ) &
            &                                                   *    tree                             %volumeWeight
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine historyStore

  subroutine historyWrite(self)
    !!{
    Store the global history data to the \glc\ output file.
    !!}
    use :: Error                           , only : Error_Report
    use :: Output_HDF5                     , only : outputFile
    use :: HDF5_Access                     , only : hdf5Access
    use :: IO_HDF5                         , only : hdf5Object
    use :: Numerical_Constants_Astronomical, only : gigaYear    , massSolar, megaParsec
    implicit none
    class           (*         ), intent(inout)               :: self
    double precision            , allocatable  , dimension(:) :: rateStarFormationDisk , rateStarFormationSpheroid, &
         &                                                       densityStellar        , densityStellarDisk       , &
         &                                                       densityStellarSpheroid, densityColdGas           , &
         &                                                       densityHotHaloGas     , densityNode
    type            (hdf5Object)                              :: historyDataset        , historyGroup

    select type (self)
    class is (mergerTreeEvolveTimestepHistory)
      !$ call hdf5Access%set()
       historyGroup=outputFile%openGroup('globalHistory','Global (volume averaged) history for this model.')
       if (.not.historyGroup%hasDataset('time')) then
          ! Write non-cumulative datasets on the first write.
          call historyGroup  %writeDataset  (self%time                       ,"time"             ,"Time [Gyr]"                       ,datasetReturned=historyDataset)
          call historyDataset%writeAttribute(gigaYear                        ,"unitsInSI"                                                                           )
          call historyDataset%close         (                                                                                                                       )
          call historyGroup  %writeDataset  (self%expansionFactor            ,"expansionFactor"  ,"Expansion factor []"                                             )
          call historyGroup  %writeDataset  (self%rateStarFormation          ,"rateStarFormation","Star formation rate [M⊙/Gyr/Mpc³]",datasetReturned=historyDataset)
          call historyDataset%writeAttribute(massSolar/gigaYear/megaParsec**3,"unitsInSI"                                                                           )
          call historyDataset%close         (                                                                                                                       )
       else
          ! Read existing datasets, and accumulate.
          call historyGroup%readDataset('rateStarFormationDisk'    ,rateStarFormationDisk    )
          call historyGroup%readDataset('rateStarFormationSpheroid',rateStarFormationSpheroid)
          call historyGroup%readDataset('densityStellar'           ,densityStellar           )
          call historyGroup%readDataset('densityStellarDisk'       ,densityStellarDisk       )
          call historyGroup%readDataset('densityStellarSpheroid'   ,densityStellarSpheroid   )
          call historyGroup%readDataset('densityColdGas'           ,densityColdGas           )
          call historyGroup%readDataset('densityHotHaloGas'        ,densityHotHaloGas        )
          call historyGroup%readDataset('densityNode'              ,densityNode              )
          self%rateStarFormationDisk    =self%rateStarFormationDisk    +rateStarFormationDisk
          self%rateStarFormationSpheroid=self%rateStarFormationSpheroid+rateStarFormationSpheroid
          self%densityStellar           =self%densityStellar           +densityStellar
          self%densityStellarDisk       =self%densityStellarDisk       +densityStellarDisk
          self%densityStellarSpheroid   =self%densityStellarSpheroid   +densityStellarSpheroid
          self%densityColdGas           =self%densityColdGas           +densityColdGas
          self%densityHotHaloGas        =self%densityHotHaloGas        +densityHotHaloGas
          self%densityNode              =self%densityNode              +densityNode
       end if
       ! Accumulate other datasets. (We are doing this inside an hdf5Access lock, so OpenMP threads will not conflict here.)
       call historyGroup  %writeDataset  (self%rateStarFormationDisk      ,"rateStarFormationDisk"    ,"Star formation rate in disks [M⊙/Gyr/Mpc³]"    ,datasetReturned=historyDataset)
       call historyDataset%writeAttribute(massSolar/gigaYear/megaParsec**3,"unitsInSI"                                                                                                )
       call historyDataset%close         (                                                                                                                                            )
       call historyGroup  %writeDataset  (self%rateStarFormationSpheroid  ,"rateStarFormationSpheroid","Star formation rate in spheroids [M⊙/Gyr/Mpc³]",datasetReturned=historyDataset)
       call historyDataset%writeAttribute(massSolar/gigaYear/megaParsec**3,"unitsInSI"                                                                                                )
       call historyDataset%close         (                                                                                                                                            )
       call historyGroup  %writeDataset  (self%densityStellar             ,"densityStellar"           ,"Stellar mass density [M⊙/Mpc³]"                ,datasetReturned=historyDataset)
       call historyDataset%writeAttribute(massSolar/megaParsec**3         ,"unitsInSI"                                                                                                )
       call historyDataset%close         (                                                                                                                                            )
       call historyGroup  %writeDataset  (self%densityStellarDisk         ,"densityStellarDisk"       ,"Stellar mass density in disks [M⊙/Mpc³]"       ,datasetReturned=historyDataset)
       call historyDataset%writeAttribute(massSolar/megaParsec**3         ,"unitsInSI"                                                                                                )
       call historyDataset%close         (                                                                                                                                            )
       call historyGroup  %writeDataset  (self%densityStellarSpheroid     ,"densityStellarSpheroid"   ,"Stellar mass density in spheroids [M⊙/Mpc³]"   ,datasetReturned=historyDataset)
       call historyDataset%writeAttribute(massSolar/megaParsec**3         ,"unitsInSI"                                                                                                )
       call historyDataset%close         (                                                                                                                                            )
       call historyGroup  %writeDataset  (self%densityColdGas             ,"densityColdGas"           ,"Gas mass density [M⊙/Mpc³]"                    ,datasetReturned=historyDataset)
       call historyDataset%writeAttribute(massSolar/megaParsec**3         ,"unitsInSI"                                                                                                )
       call historyDataset%close         (                                                                                                                                            )
       call historyGroup  %writeDataset  (self%densityHotHaloGas          ,"densityHotHaloGas"        ,"Hot gas mass density [M⊙/Mpc³]"                ,datasetReturned=historyDataset)
       call historyDataset%writeAttribute(massSolar/megaParsec**3         ,"unitsInSI"                                                                                                )
       call historyDataset%close         (                                                                                                                                            )
       call historyGroup  %writeDataset  (self%densityNode                ,"densityNode"              ,"Node mass density [M⊙/Mpc³]"                   ,datasetReturned=historyDataset)
       call historyDataset%writeAttribute(massSolar/megaParsec**3         ,"unitsInSI"                                                                                                )
       call historyDataset%close         (                                                                                                                                            )
       call historyGroup  %close         (                                                                                                                                            )
       !$ call hdf5Access %unset         (                                                                                                                                            )
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine historyWrite
