!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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

!% Contains a module which implements a time-stepping criterion for merger tree evolution which permits global histories to be
!% stored.

module Merger_Tree_Timesteps_History
  !% Implements a simple time-stepping criterion for merger tree evolution.
  use FGSL
  use Galacticus_Nodes
  implicit none
  private
  public :: Merger_Tree_Timestep_History            , Merger_Tree_History_Write                  , &
       &    Merger_Tree_Timestep_History_State_Store, Merger_Tree_Timestep_History_State_Retrieve

  ! Variable inidicating if module is initialized and active.
  logical                                                        :: timestepHistoryInitialized      =.false.
  logical                                                        :: diskActive                              , spheroidActive               , &
       &                                                            outputTimestepHistory

  ! Variables which control the distribution of timesteps.
  integer                                                        :: timestepHistorySteps
  double precision                                               :: timestepHistoryBegin                    , timestepHistoryEnd

  ! Storage arrays.
  double precision                   , allocatable, dimension(:) :: historyDiskStarFormationRate            , historyDiskStellarDensity    , &
       &                                                            historyExpansion                        , historyGasDensity            , &
       &                                                            historyHotGasDensity                    , historyNodeDensity           , &
       &                                                            historySpheroidStarFormationRate        , historySpheroidStellarDensity, &
       &                                                            historyStarFormationRate                , historyStellarDensity        , &
       &                                                            historyTime

  ! Interpolation variables.
  type            (fgsl_interp_accel)                            :: interpolationAccelerator
  !$omp threadprivate(interpolationAccelerator)
  
contains

  !# <timeStepsTask>
  !#  <unitName>Merger_Tree_Timestep_History</unitName>
  !# </timeStepsTask>
  subroutine Merger_Tree_Timestep_History(node,timeStep,End_Of_Timestep_Task,report,nodeLock,lockType)
    !% Determines the timestep to go to the next tabulation point for global history storage.
    use, intrinsic :: ISO_C_Binding
    use Input_Parameters
    use Cosmology_Functions
    use Memory_Management
    use Numerical_Ranges
    use Numerical_Interpolation
    use Merger_Trees_Evolve_Timesteps_Template
    use Evolve_To_Time_Reports
    use ISO_Varying_String
    implicit none
    type            (treeNode                     ), intent(inout)          , pointer :: node
    procedure       (End_Of_Timestep_Task_Template), intent(inout)          , pointer :: End_Of_Timestep_Task
    double precision                               , intent(inout)                    :: timeStep
    logical                                        , intent(in   )                    :: report
    type            (treeNode                     ), intent(inout), optional, pointer :: nodeLock
    type            (varying_string               ), intent(inout), optional          :: lockType
    class           (nodeComponentBasic           )                         , pointer :: basic
    class           (cosmologyFunctionsClass      )                         , pointer :: cosmologyFunctions_
    integer         (c_size_t                     )                                   :: timeIndex
    double precision                                                                  :: ourTimeStep         , time, &
         &                                                                               universeAge

    if (.not.timestepHistoryInitialized) then
       !$omp critical (timestepHistoryInitialize)
       if (.not.timestepHistoryInitialized) then
          ! Get module parameters.
          !# <inputParameter>
          !#   <name>outputTimestepHistory</name>
          !#   <cardinality>1</cardinality>
          !#   <defaultValue>.true.</defaultValue>
          !#   <description>Specifies whether or not accumulated galaxy population history statistics should be output.</description>
          !#   <group>timeStepping</group>
          !#   <source>globalParameters</source>
          !#   <type>real</type>
          !# </inputParameter>
          if (outputTimestepHistory) then
             ! Get the default cosmology functions object.
             cosmologyFunctions_ => cosmologyFunctions()
             ! Determine if we have active components that can provide star formation rates.
             diskActive           =    defaultDiskComponent%starFormationRateIsGettable()
             spheroidActive       =defaultSpheroidComponent%starFormationRateIsGettable()
             ! Get time at present day.
             universeAge=cosmologyFunctions_%cosmicTime(expansionFactor=0.999d0)
             ! Get module parameters.
             !# <inputParameter>
             !#   <name>timestepHistoryBegin</name>
             !#   <cardinality>1</cardinality>
             !#   <defaultValue>0.05d0*universeAge</defaultValue>
             !#   <description>The earliest time at which to tabulate the volume averaged history of galaxies (in Gyr).</description>
             !#   <group>timeStepping</group>
             !#   <source>globalParameters</source>
             !#   <type>real</type>
             !# </inputParameter>
             !# <inputParameter>
             !#   <name>timestepHistoryEnd</name>
             !#   <cardinality>1</cardinality>
             !#   <defaultValue>universeAge</defaultValue>
             !#   <description>The latest time at which to tabulate the volume averaged history of galaxies (in Gyr).</description>
             !#   <group>timeStepping</group>
             !#   <source>globalParameters</source>
             !#   <type>real</type>
             !# </inputParameter>
             !# <inputParameter>
             !#   <name>timestepHistorySteps</name>
             !#   <cardinality>1</cardinality>
             !#   <defaultValue>30</defaultValue>
             !#   <description>The number of steps (spaced logarithmically in cosmic time) at which to tabulate the volume averaged history of galaxies.</description>
             !#   <group>timeStepping</group>
             !#   <source>globalParameters</source>
             !#   <type>integer</type>
             !# </inputParameter>
             ! Allocate storage arrays.
             call allocateArray(historyTime                     ,[timestepHistorySteps])
             call allocateArray(historyExpansion                ,[timestepHistorySteps])
             call allocateArray(historyStarFormationRate        ,[timestepHistorySteps])
             call allocateArray(historyDiskStarFormationRate    ,[timestepHistorySteps])
             call allocateArray(historySpheroidStarFormationRate,[timestepHistorySteps])
             call allocateArray(historyStellarDensity           ,[timestepHistorySteps])
             call allocateArray(historyDiskStellarDensity       ,[timestepHistorySteps])
             call allocateArray(historySpheroidStellarDensity   ,[timestepHistorySteps])
             call allocateArray(historyGasDensity               ,[timestepHistorySteps])
             call allocateArray(historyHotGasDensity            ,[timestepHistorySteps])
             call allocateArray(historyNodeDensity              ,[timestepHistorySteps])
             ! Initialize arrays.
             historyTime=Make_Range(timestepHistoryBegin,timestepHistoryEnd,timestepHistorySteps,rangeTypeLogarithmic)
             do timeIndex=1,timestepHistorySteps
                historyExpansion(timeIndex)=cosmologyFunctions_%expansionFactor(historyTime(timeIndex))
             end do
             historyStarFormationRate        =0.0d0
             historyDiskStarFormationRate    =0.0d0
             historySpheroidStarFormationRate=0.0d0
             historyStellarDensity           =0.0d0
             historyDiskStellarDensity       =0.0d0
             historySpheroidStellarDensity   =0.0d0
             historyGasDensity               =0.0d0
             historyHotGasDensity            =0.0d0
             historyNodeDensity              =0.0d0
          end if
          timestepHistoryInitialized=.true.
       end if
       !$omp end critical (timestepHistoryInitialize)
    end if

    ! Return if history data is not being gathered.
    if (.not.outputTimestepHistory) return
    
    ! Adjust timestep.
    ! Get current cosmic time.
    basic => node%basic()
    time=basic%time()

    ! Determine how long until next available timestep.
    timeIndex=Interpolate_Locate(historyTime,interpolationAccelerator,time)
    if (time < historyTime(timeIndex+1)) then
       ! Find next time for storage.
       ourTimeStep=historyTime(timeIndex+1)-time

       ! Set return value if our timestep is smaller than current one.
       if (ourTimeStep <= timeStep) then
          if (present(nodeLock)) nodeLock => node
          if (present(lockType)) lockType =  "history"
          timeStep=ourTimeStep
          End_Of_Timestep_Task => Merger_Tree_History_Store
       end if
    end if
    if (report) call Evolve_To_Time_Report("history: ",timeStep)
    return
  end subroutine Merger_Tree_Timestep_History

  subroutine Merger_Tree_History_Store(tree,node,deadlockStatus)
    !% Store various properties in global arrays.
    use, intrinsic :: ISO_C_Binding
    use Numerical_Interpolation
    use Galactic_Structure_Options
    use Galactic_Structure_Enclosed_Masses
    implicit none
    type            (mergerTree           ), intent(in   )          :: tree
    type            (treeNode             ), intent(inout), pointer :: node
    integer                                , intent(inout)          :: deadlockStatus
    class           (nodeComponentBasic   )               , pointer :: basic
    class           (nodeComponentDisk    )               , pointer :: disk
    class           (nodeComponentSpheroid)               , pointer :: spheroid
    integer         (c_size_t             )                         :: timeIndex
    double precision                                                :: rateStarFormationDisk    , massHotGas, &
         &                                                             rateStarFormationSpheroid, time
    !GCC$ attributes unused :: deadlockStatus
    
    ! Get components.
    basic    => node%basic   ()
    disk     => node%disk    ()
    spheroid => node%spheroid()
    ! Get current cosmic time.
    time=basic%time()
    ! Determine how long until next available timestep.
    if (time == historyTime(timestepHistorySteps)) then
       timeIndex=timestepHistorySteps
    else
       timeIndex=Interpolate_Locate(historyTime,interpolationAccelerator,time)
    end if
    ! Extract disk and spheroid star formation rates.
    rateStarFormationDisk                        =0.0d0
    rateStarFormationSpheroid                    =0.0d0
    if (diskActive)     rateStarFormationDisk    =disk    %starFormationRate()
    if (spheroidActive) rateStarFormationSpheroid=spheroid%starFormationRate()
    ! Accumulate the properties.
    ! Star formation rate:
    !$omp critical(timestepHistoryAccumulate)
    historyStarFormationRate                       (timeIndex)=+historyStarFormationRate          (timeIndex)                                                         &
         &                                                     +(rateStarFormationDisk+rateStarFormationSpheroid)                                                     &
         &                                                     *tree%volumeWeight
    historyDiskStarFormationRate                   (timeIndex)=+historyDiskStarFormationRate      (timeIndex)                                                         &
         &                                                     + rateStarFormationDisk                                                                                &
         &                                                     *tree%volumeWeight
    historySpheroidStarFormationRate               (timeIndex)= historySpheroidStarFormationRate  (timeIndex)                                                         &
         &                                                     +                       rateStarFormationSpheroid                                                      &
         &                                                     *tree%volumeWeight
    ! Stellar densities.
    historyStellarDensity                          (timeIndex)=+historyStellarDensity             (timeIndex)                                                         &
         &                                                     +  Galactic_Structure_Enclosed_Mass(node,                                    massType=massTypeStellar) &
         &                                                     *tree%volumeWeight
    historyDiskStellarDensity                      (timeIndex)= historyDiskStellarDensity         (timeIndex)                                                         &
         &                                                     +  Galactic_Structure_Enclosed_Mass(node,componentType=componentTypeDisk    ,massType=massTypeStellar) &
         &                                                     *tree%volumeWeight
    historySpheroidStellarDensity                  (timeIndex)= historySpheroidStellarDensity     (timeIndex)                                                         &
         &                                                     +  Galactic_Structure_Enclosed_Mass(node,componentType=componentTypeSpheroid,massType=massTypeStellar) &
         &                                                     *tree%volumeWeight

    ! Hot gas density.
    massHotGas                                                =+  Galactic_Structure_Enclosed_Mass(node,componentType=componentTypeHotHalo                          )
    historyHotGasDensity                           (timeIndex)=+historyHotGasDensity              (timeIndex)                                                         &
         &                                                     +massHotGas                                                                                            &
         &                                                     *tree%volumeWeight
    ! Galactic gas density.
    historyGasDensity                              (timeIndex)=+historyGasDensity                 (timeIndex)                                                         &
         &                                                     +(                                                                                                     &
         &                                                        Galactic_Structure_Enclosed_Mass(node,massType=massTypeGaseous                                    ) &
         &                                                       -massHotGas                                                                                          &
         &                                                      )                                                                                                     &
         &                                                     *tree%volumeWeight
    ! Node density
    if (.not.node%isSatellite()) historyNodeDensity(timeIndex)=+historyNodeDensity                (timeIndex)                                                         &
            &                                                  +basic%mass()                                                                                          &
            &                                                  *tree%volumeWeight
    !$omp end critical(timestepHistoryAccumulate)
    return
  end subroutine Merger_Tree_History_Store

  !# <hdfPreCloseTask>
  !#  <unitName>Merger_Tree_History_Write</unitName>
  !# </hdfPreCloseTask>
  subroutine Merger_Tree_History_Write
    !% Store the global history data to the \glc\ output file.
    use Galacticus_HDF5
    use Numerical_Constants_Astronomical
    implicit none
    type(hdf5Object) :: historyDataset, historyGroup

    ! Output the history data if and only if any has been collated.
    if (outputTimestepHistory) then
       !$omp critical (HDF5_Access)
       !@ <outputType>
       !@   <name>globalHistory</name>
       !@   <description>A set of volume-averaged properites describing the mass content of the Universe.</description>
       !@ </outputType>

       historyGroup=galacticusOutputFile%openGroup('globalHistory','Global (volume averaged) history for this model.')

       !@ <outputProperty>
       !@   <name>historyTime</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>The times at which global history data are stored.</description>
       !@   <label>???</label>
       !@   <outputType>globalHistory</outputType>
       !@ </outputProperty>
       call historyGroup%writeDataset(historyTime                     ,"historyTime"                     ,"Time [Gyr]"                                       ,datasetReturned=historyDataset)
       call historyDataset%writeAttribute(gigaYear                        ,"unitsInSI")
       call historyDataset%close()

       !@ <outputProperty>
       !@   <name>historyExpansion</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>The expansion factors at which global history data are stored.</description>
       !@   <label>???</label>
       !@   <outputType>globalHistory</outputType>
       !@ </outputProperty>
       call historyGroup%writeDataset(historyExpansion                ,"historyExpansion"                ,"Expansion factor []"                                                             )

       !@ <outputProperty>
       !@   <name>historyStarFormationRate</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>The mean rate of star formation in the Universe.</description>
       !@   <label>???</label>
       !@   <outputType>globalHistory</outputType>
       !@ </outputProperty>
       call historyGroup%writeDataset(historyStarFormationRate        ,"historyStarFormationRate"        ,"Star formation rate [M⊙/Gyr/Mpc³]"             ,datasetReturned=historyDataset)
       call historyDataset%writeAttribute(massSolar/gigaYear/megaParsec**3,"unitsInSI")
       call historyDataset%close()

       !@ <outputProperty>
       !@   <name>historyDiskStarFormationRate</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>The mean rate of star formation in disks in the Universe.</description>
       !@   <label>???</label>
       !@   <outputType>globalHistory</outputType>
       !@ </outputProperty>
       call historyGroup%writeDataset(historyDiskStarFormationRate    ,"historyDiskStarFormationRate"    ,"Star formation rate in disks [M⊙/Gyr/Mpc³]"    ,datasetReturned=historyDataset)
       call historyDataset%writeAttribute(massSolar/gigaYear/megaParsec**3,"unitsInSI")
       call historyDataset%close()

       !@ <outputProperty>
       !@   <name>historySpheroidStarFormationRate</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>The mean rate of star formation in spheroids in the Universe.</description>
       !@   <label>???</label>
       !@   <outputType>globalHistory</outputType>
       !@ </outputProperty>
       call historyGroup%writeDataset(historySpheroidStarFormationRate,"historySpheroidStarFormationRate","Star formation rate in spheroids [M⊙/Gyr/Mpc³]",datasetReturned=historyDataset)
       call historyDataset%writeAttribute(massSolar/gigaYear/megaParsec**3,"unitsInSI")
       call historyDataset%close()

       !@ <outputProperty>
       !@   <name>historyStellarDensity</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>The mean density of stars in the Universe.</description>
       !@   <label>???</label>
       !@   <outputType>globalHistory</outputType>
       !@ </outputProperty>
       call historyGroup%writeDataset(historyStellarDensity           ,"historyStellarDensity"           ,"Stellar mass density [M⊙/Mpc³]"                ,datasetReturned=historyDataset)
       call historyDataset%writeAttribute(massSolar/megaParsec**3         ,"unitsInSI")
       call historyDataset%close()

       !@ <outputProperty>
       !@   <name>historyDiskStellarDensity</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>The mean density of stars in disks in the Universe.</description>
       !@   <label>???</label>
       !@   <outputType>globalHistory</outputType>
       !@ </outputProperty>
       call historyGroup%writeDataset(historyDiskStellarDensity       ,"historyDiskStellarDensity"       ,"Stellar mass density in disks [M⊙/Mpc³]"       ,datasetReturned=historyDataset)
       call historyDataset%writeAttribute(massSolar/megaParsec**3         ,"unitsInSI")
       call historyDataset%close()

       !@ <outputProperty>
       !@   <name>historySpheroidStellarDensity</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>The mean density of stars in spheroids in the Universe.</description>
       !@   <label>???</label>
       !@   <outputType>globalHistory</outputType>
       !@ </outputProperty>
       call historyGroup%writeDataset(historySpheroidStellarDensity   ,"historySpheroidStellarDensity"   ,"Stellar mass density in spheroids [M⊙/Mpc³]"   ,datasetReturned=historyDataset)
       call historyDataset%writeAttribute(massSolar/megaParsec**3         ,"unitsInSI")
       call historyDataset%close()

       !@ <outputProperty>
       !@   <name>historyGasDensity</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>The mean density of ISM gas in the Universe.</description>
       !@   <label>???</label>
       !@   <outputType>globalHistory</outputType>
       !@ </outputProperty>
       call historyGroup%writeDataset(historyGasDensity               ,"historyGasDensity"               ,"Gas mass density [M⊙/Mpc³]"                    ,datasetReturned=historyDataset)
       call historyDataset%writeAttribute(massSolar/megaParsec**3         ,"unitsInSI")
       call historyDataset%close()

       !@ <outputProperty>
       !@   <name>historyHotGasDensity</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>The mean density of hot gas in the Universe.</description>
       !@   <label>???</label>
       !@   <outputType>globalHistory</outputType>
       !@ </outputProperty>
       call historyGroup%writeDataset(historyHotGasDensity            ,"historyHotGasDensity"            ,"Hot gas mass density [M⊙/Mpc³]"                ,datasetReturned=historyDataset)
       call historyDataset%writeAttribute(massSolar/megaParsec**3         ,"unitsInSI")
       call historyDataset%close()

       !@ <outputProperty>
       !@   <name>historyNodeDensity</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>The mean density of resolved nodes in the Universe.</description>
       !@   <label>???</label>
       !@   <outputType>globalHistory</outputType>
       !@ </outputProperty>
       call historyGroup%writeDataset(historyNodeDensity              ,"historyNodeDensity"              ,"Node mass density [M⊙/Mpc³]"                   ,datasetReturned=historyDataset)
       call historyDataset%writeAttribute(massSolar/megaParsec**3         ,"unitsInSI")
       call historyDataset%close()

       call historyGroup%close()
       !$omp end critical (HDF5_Access)
    end if
    return
  end subroutine Merger_Tree_History_Write

  !# <galacticusStateStoreTask>
  !#  <unitName>Merger_Tree_Timestep_History_State_Store</unitName>
  !# </galacticusStateStoreTask>
  subroutine Merger_Tree_Timestep_History_State_Store(stateFile,fgslStateFile,stateOperatorID)
    !% Write the tablulation state to file.
    use Galacticus_Display
    use FGSL
    implicit none
    integer           , intent(in   ) :: stateFile    , stateOperatorID
    type   (fgsl_file), intent(in   ) :: fgslStateFile
    !GCC$ attributes unused :: stateOperatorID, fgslStateFile

    call Galacticus_Display_Indent  ('storing state for "merger tree history"',verbosity=verbosityWorking)
    write (stateFile) timestepHistoryInitialized
    if (timestepHistoryInitialized) then
       write (stateFile) outputTimestepHistory
       if (outputTimestepHistory) then
          write (stateFile) diskActive                      , spheroidActive               , &
               &            timestepHistorySteps
          write (stateFile) timestepHistoryBegin            , timestepHistoryEnd           , &
               &            historyTime                     , historyExpansion             , &
               &            historyStarFormationRate        , historyDiskStarFormationRate , &
               &            historySpheroidStarFormationRate, historyStellarDensity        , &
               &            historyDiskStellarDensity       , historySpheroidStellarDensity, &
               &            historyGasDensity               , historyHotGasDensity         , &
               &            historyNodeDensity
       end if
    end if
    call Galacticus_Display_Unindent('done'                                   ,verbosity=verbosityWorking)
    return
  end subroutine Merger_Tree_Timestep_History_State_Store

  !# <galacticusStateRetrieveTask>
  !#  <unitName>Merger_Tree_Timestep_History_State_Retrieve</unitName>
  !# </galacticusStateRetrieveTask>
  subroutine Merger_Tree_Timestep_History_State_Retrieve(stateFile,fgslStateFile,stateOperatorID)
    !% Retrieve the tabulation state from the file.
    use Galacticus_Display
    use FGSL
    implicit none
    integer           , intent(in   ) :: stateFile    , stateOperatorID
    type   (fgsl_file), intent(in   ) :: fgslStateFile
    !GCC$ attributes unused :: stateOperatorID, fgslStateFile

    call Galacticus_Display_Indent  ('restoring state for "merger tree history"',verbosity=verbosityWorking)
    if (allocated(historyTime                     )) deallocate(historyTime                     )
    if (allocated(historyExpansion                )) deallocate(historyExpansion                )
    if (allocated(historyStarFormationRate        )) deallocate(historyStarFormationRate        )
    if (allocated(historyDiskStarFormationRate    )) deallocate(historyDiskStarFormationRate    )
    if (allocated(historySpheroidStarFormationRate)) deallocate(historySpheroidStarFormationRate)
    if (allocated(historyStellarDensity           )) deallocate(historyStellarDensity           )
    if (allocated(historyDiskStellarDensity       )) deallocate(historyDiskStellarDensity       )
    if (allocated(historySpheroidStellarDensity   )) deallocate(historySpheroidStellarDensity   )
    if (allocated(historyGasDensity               )) deallocate(historyGasDensity               )
    if (allocated(historyHotGasDensity            )) deallocate(historyHotGasDensity            )
    if (allocated(historyNodeDensity              )) deallocate(historyNodeDensity              )
    read (stateFile) timestepHistoryInitialized
    if (timestepHistoryInitialized) then
       read (stateFile) outputTimestepHistory
       if (outputTimestepHistory) then
          read (stateFile) diskActive          , spheroidActive, &
               &           timestepHistorySteps
          allocate(historyTime                     (timestepHistorySteps))
          allocate(historyExpansion                (timestepHistorySteps))
          allocate(historyStarFormationRate        (timestepHistorySteps))
          allocate(historyDiskStarFormationRate    (timestepHistorySteps))
          allocate(historySpheroidStarFormationRate(timestepHistorySteps))
          allocate(historyStellarDensity           (timestepHistorySteps))
          allocate(historyDiskStellarDensity       (timestepHistorySteps))
          allocate(historySpheroidStellarDensity   (timestepHistorySteps))
          allocate(historyGasDensity               (timestepHistorySteps))
          allocate(historyHotGasDensity            (timestepHistorySteps))
          allocate(historyNodeDensity              (timestepHistorySteps))
          read (stateFile) timestepHistoryBegin            , timestepHistoryEnd           , &
               &           historyTime                     , historyExpansion             , &
               &           historyStarFormationRate        , historyDiskStarFormationRate , &
               &           historySpheroidStarFormationRate, historyStellarDensity        , &
               &           historyDiskStellarDensity       , historySpheroidStellarDensity, &
               &           historyGasDensity               , historyHotGasDensity         , &
               &           historyNodeDensity
       end if
    end if
    call Galacticus_Display_Unindent('done'                                     ,verbosity=verbosityWorking)
    return
  end subroutine Merger_Tree_Timestep_History_State_Retrieve

end module Merger_Tree_Timesteps_History
