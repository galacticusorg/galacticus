!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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
  public :: Merger_Tree_Timestep_History, Merger_Tree_History_Write

  ! Variable inidicating if module is initialized and active.
  logical                                                        :: timestepHistoryInitialized      =.false.
  logical                                                        :: diskActive                              , spheroidActive

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
  subroutine Merger_Tree_Timestep_History(thisNode,timeStep,End_Of_Timestep_Task,report,lockNode,lockType)
    !% Determines the timestep to go to the next tabulation point for global history storage.
    use Input_Parameters
    use Cosmology_Functions
    use Memory_Management
    use Numerical_Ranges
    use Numerical_Interpolation
    use Merger_Trees_Evolve_Timesteps_Template
    use Evolve_To_Time_Reports
    use ISO_Varying_String
    implicit none
    type            (treeNode                     ), intent(inout)          , pointer :: thisNode
    procedure       (End_Of_Timestep_Task_Template), intent(inout)          , pointer :: End_Of_Timestep_Task
    double precision                               , intent(inout)                    :: timeStep
    logical                                        , intent(in   )                    :: report
    type            (treeNode                     ), intent(inout), optional, pointer :: lockNode
    type            (varying_string               ), intent(inout), optional          :: lockType
    class           (nodeComponentBasic           )                         , pointer :: thisBasicComponent
    integer                                                                           :: timeIndex
    double precision                                                                  :: ourTimeStep         , time

    if (.not.timestepHistoryInitialized) then
       !$omp critical (timestepHistoryInitialize)
       if (.not.timestepHistoryInitialized) then
          ! Determine if we have active components that can provide star formation rates.
          diskActive           =    defaultDiskComponent%starFormationRateIsGettable()
          spheroidActive       =defaultSpheroidComponent%starFormationRateIsGettable()
          ! Get time at present day.
          time=Cosmology_Age(aExpansion=0.999d0)
          ! Get module parameters.
          !@ <inputParameter>
          !@   <name>timestepHistoryBegin</name>
          !@   <defaultValue>5\% of the age of the Universe</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The earliest time at which to tabulate the volume averaged history of galaxies (in Gyr).
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@   <group>timeStepping</group>
          !@ </inputParameter>
          call Get_Input_Parameter('timestepHistoryBegin',timestepHistoryBegin,defaultValue=0.05d0*time)
          !@ <inputParameter>
          !@   <name>timestepHistoryEnd</name>
          !@   <defaultValue>The age of the Universe</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The latest time at which to tabulate the volume averaged history of galaxies (in Gyr).
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@   <group>timeStepping</group>
          !@ </inputParameter>
          call Get_Input_Parameter('timestepHistoryEnd'  ,timestepHistoryEnd  ,defaultValue=       time)
          !@ <inputParameter>
          !@   <name>timestepHistorySteps</name>
          !@   <defaultValue>30</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The number of steps (spaced logarithmically in cosmic time) at which to tabulate the volume averaged history of galaxies.
          !@   </description>
          !@   <type>integer</type>
          !@   <cardinality>1</cardinality>
          !@   <group>timeStepping</group>
          !@ </inputParameter>
          call Get_Input_Parameter('timestepHistorySteps',timestepHistorySteps,defaultValue=30         )
          ! Allocate storage arrays.
          call Alloc_Array(historyTime                     ,[timestepHistorySteps])
          call Alloc_Array(historyExpansion                ,[timestepHistorySteps])
          call Alloc_Array(historyStarFormationRate        ,[timestepHistorySteps])
          call Alloc_Array(historyDiskStarFormationRate    ,[timestepHistorySteps])
          call Alloc_Array(historySpheroidStarFormationRate,[timestepHistorySteps])
          call Alloc_Array(historyStellarDensity           ,[timestepHistorySteps])
          call Alloc_Array(historyDiskStellarDensity       ,[timestepHistorySteps])
          call Alloc_Array(historySpheroidStellarDensity   ,[timestepHistorySteps])
          call Alloc_Array(historyGasDensity               ,[timestepHistorySteps])
          call Alloc_Array(historyHotGasDensity            ,[timestepHistorySteps])
          call Alloc_Array(historyNodeDensity              ,[timestepHistorySteps])
          ! Initialize arrays.
          historyTime=Make_Range(timestepHistoryBegin,timestepHistoryEnd,timestepHistorySteps,rangeTypeLogarithmic)
          do timeIndex=1,timestepHistorySteps
             historyExpansion(timeIndex)=Expansion_Factor(historyTime(timeIndex))
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
          timestepHistoryInitialized      =.true.
       end if
       !$omp end critical (timestepHistoryInitialize)
    end if

    ! Adjust timestep.
    ! Get current cosmic time.
    thisBasicComponent => thisNode%basic()
    time=thisBasicComponent%time()

    ! Determine how long until next available timestep.
    timeIndex=Interpolate_Locate(timestepHistorySteps,historyTime,interpolationAccelerator,time)
    if (time < historyTime(timeIndex+1)) then
       ! Find next time for storage.
       ourTimeStep=historyTime(timeIndex+1)-time

       ! Set return value if our timestep is smaller than current one.
       if (ourTimeStep <= timeStep) then
          if (present(lockNode)) lockNode => thisNode
          if (present(lockType)) lockType =  "history"
          timeStep=ourTimeStep
          End_Of_Timestep_Task => Merger_Tree_History_Store
       end if
    end if
    if (report) call Evolve_To_Time_Report("history: ",timeStep)
    return
  end subroutine Merger_Tree_Timestep_History

  subroutine Merger_Tree_History_Store(thisTree,thisNode,deadlockStatus)
    !% Store various properties in global arrays.
    use Numerical_Interpolation
    use Galactic_Structure_Options
    use Galactic_Structure_Enclosed_Masses
    implicit none
    type            (mergerTree           ), intent(in   )          :: thisTree
    type            (treeNode             ), intent(inout), pointer :: thisNode
    integer                                , intent(inout)          :: deadlockStatus
    class           (nodeComponentBasic   )               , pointer :: thisBasicComponent
    class           (nodeComponentDisk    )               , pointer :: thisDiskComponent
    class           (nodeComponentSpheroid)               , pointer :: thisSpheroidComponent
    integer                                                         :: timeIndex
    double precision                                                :: diskStarFormationRate    , hotGasMass, &
         &                                                             spheroidStarFormationRate, time

    ! Get components.
    thisBasicComponent    => thisNode%basic   ()
    thisDiskComponent     => thisNode%disk    ()
    thisSpheroidComponent => thisNode%spheroid()

    ! Get current cosmic time.
    time=thisBasicComponent%time()

    ! Determine how long until next available timestep.
    if (time == historyTime(timestepHistorySteps)) then
       timeIndex=timestepHistorySteps
    else
       timeIndex=Interpolate_Locate(timestepHistorySteps,historyTime,interpolationAccelerator,time)
    end if

    ! Accumulate the properties.
    ! Star formation rate:
    diskStarFormationRate    =0.0d0
    spheroidStarFormationRate=0.0d0
    if (diskActive)     diskStarFormationRate    =thisDiskComponent    %starFormationRate()
    if (spheroidActive) spheroidStarFormationRate=thisSpheroidComponent%starFormationRate()
    historyStarFormationRate        (timeIndex)= historyStarFormationRate        (timeIndex)                                                               &
         &                                      +(diskStarFormationRate+spheroidStarFormationRate)                                                         &
         &                                      *thisTree%volumeWeight
    historyDiskStarFormationRate    (timeIndex)= historyDiskStarFormationRate    (timeIndex)                                                               &
         &                                      + diskStarFormationRate                                                                                    &
         &                                      *thisTree%volumeWeight
    historySpheroidStarFormationRate(timeIndex)= historySpheroidStarFormationRate(timeIndex)                                                               &
         &                                      +                       spheroidStarFormationRate                                                          &
         &                                      *thisTree%volumeWeight
    ! Stellar densities.
    historyStellarDensity           (timeIndex)= historyStellarDensity           (timeIndex)                                                               &
         &                                      +  Galactic_Structure_Enclosed_Mass(thisNode,massType=massTypeStellar                                    ) &
         &                                      *thisTree%volumeWeight
    historyDiskStellarDensity       (timeIndex)= historyDiskStellarDensity       (timeIndex)                                                               &
         &                                      +  Galactic_Structure_Enclosed_Mass(thisNode,componentType=componentTypeDisk    ,massType=massTypeStellar) &
         &                                      *thisTree%volumeWeight
    historySpheroidStellarDensity   (timeIndex)= historySpheroidStellarDensity   (timeIndex)                                                               &
         &                                        +Galactic_Structure_Enclosed_Mass(thisNode,componentType=componentTypeSpheroid,massType=massTypeStellar) &
         &                                      *thisTree%volumeWeight

    ! Hot gas density.
    hotGasMass                                 =   Galactic_Structure_Enclosed_Mass(thisNode,componentType=componentTypeHotHalo                          )
    historyHotGasDensity            (timeIndex)= historyHotGasDensity            (timeIndex)                                                               &
         &                                      +hotGasMass                                                                                                &
         &                                      *thisTree%volumeWeight
    ! Galactic gas density.
    historyGasDensity               (timeIndex)= historyGasDensity               (timeIndex)                                                               &
         &                                      +(                                                                                                         &
         &                                         Galactic_Structure_Enclosed_Mass(thisNode,massType=massTypeGaseous                                    ) &
         &                                        -hotGasMass                                                                                              &
         &                                       )                                                                                                         &
         &                                      *thisTree%volumeWeight
    ! Node density
    if (.not.thisNode%isSatellite()) then
       historyNodeDensity           (timeIndex)= historyNodeDensity              (timeIndex)                                                               &
            &                                   +thisBasicComponent%mass()                                                                                 &
            &                                   *thisTree%volumeWeight
    end if
    return
  end subroutine Merger_Tree_History_Store

  !# <hdfPreCloseTask>
  !#  <unitName>Merger_Tree_History_Write</unitName>
  !# </hdfPreCloseTask>
  subroutine Merger_Tree_History_Write
    !% Store the global history data to the \glc\ output file.
    use ISO_Varying_String
    use IO_HDF5
    use Galacticus_HDF5
    use Numerical_Constants_Astronomical
    implicit none
    type(hdf5Object) :: historyDataset, historyGroup

    ! Output the history data if and only if any has been collated.
    if (timestepHistoryInitialized) then
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

    end if
    return
  end subroutine Merger_Tree_History_Write

end module Merger_Tree_Timesteps_History
