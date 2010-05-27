!% Contains a module which implements a time-stepping criterion for merger tree evolution which permits global histories to be
!% stored.

module Merger_Tree_Timesteps_History
  !% Implements a simple time-stepping criterion for merger tree evolution.
  use FGSL
  private
  public :: Merger_Tree_Timestep_History, Merger_Tree_History_Write

  ! Variable inidicating if module is initialized and active.
  logical          :: timestepHistoryInitialized=.false.
  logical          :: timestepHistoryActive,diskActive,spheroidActive

  ! Variables which control the distribution of timesteps.
  integer          :: timestepHistorySteps
  double precision :: timestepHistoryBegin,timestepHistoryEnd

  ! Storage arrays.
  double precision, dimension(:), allocatable :: historyTime,historyExpansion,historyStarFormationRate

  ! Interpolation variables.
  type(fgsl_interp_accel)                     :: interpolationAccelerator
  !$omp threadprivate(interpolationAccelerator)

contains

  !# <timeStepsTask>
  !#  <unitName>Merger_Tree_Timestep_History</unitName>
  !# </timeStepsTask>
  subroutine Merger_Tree_Timestep_History(thisNode,timeStep,End_Of_Timestep_Task)
    !% Determines the timestep to go to the next tabulation point for global history storage.
    use Tree_Nodes
    use Tree_Node_Methods
    use Input_Parameters
    use Cosmology_Functions
    use Memory_Management
    use Numerical_Ranges
    use Numerical_Interpolation
    use Merger_Trees_Evolve_Timesteps_Template
    implicit none
    type(treeNode),                           intent(inout), pointer :: thisNode
    procedure(End_Of_Timestep_Task_Template), intent(inout), pointer :: End_Of_Timestep_Task
    double precision,                         intent(inout)          :: timeStep
    integer                                                          :: timeIndex
    double precision                                                 :: time,ourTimeStep
    
    !$omp critical (timestepHistoryInitialize)
    if (.not.timestepHistoryInitialized) then
       ! Determine if we have active components that can provide star formation rates.
       diskActive    =associated(Tree_Node_Disk_SFR)
       spheroidActive=associated(Tree_Node_Spheroid_SFR)
       timestepHistoryActive=diskActive.or.spheroidActive
       if (timestepHistoryActive) then
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
          !@ </inputParameter>
          call Get_Input_Parameter('timestepHistoryBegin',timestepHistoryBegin,defaultValue=0.05d0*time)
          !@ <inputParameter>
          !@   <name>timestepHistoryEnd</name>
          !@   <defaultValue>The age of the Universe</defaultValue>       
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The latest time at which to tabulate the volume averaged history of galaxies (in Gyr).
          !@   </description>
          !@ </inputParameter>
          call Get_Input_Parameter('timestepHistoryEnd'  ,timestepHistoryEnd  ,defaultValue=       time)
          !@ <inputParameter>
          !@   <name>timestepHistorySteps</name>
          !@   <defaultValue>30</defaultValue>       
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The number of steps (spaced logarithmically in cosmic time) at which to tabulate the volume averaged history of galaxies.
          !@   </description>
          !@ </inputParameter>
          call Get_Input_Parameter('timestepHistorySteps',timestepHistorySteps,defaultValue=30         )
          ! Allocate storage arrays.
          call Alloc_Array(historyTime             ,timestepHistorySteps,'historyTime'             )
          call Alloc_Array(historyExpansion        ,timestepHistorySteps,'historyExpansion'        )
          call Alloc_Array(historyStarFormationRate,timestepHistorySteps,'historyStarFormationRate')
          ! Initialize arrays.
          historyTime=Make_Range(timestepHistoryBegin,timestepHistoryEnd,timestepHistorySteps,rangeTypeLogarithmic)
          do timeIndex=1,timestepHistorySteps
             historyExpansion(timeIndex)=Expansion_Factor(historyTime(timeIndex))
          end do
          historyStarFormationRate=0.0d0
       end if
       timestepHistoryInitialized=.true.
    end if
    !$omp end critical (timestepHistoryInitialize)

    ! Adjust timestep if this module is active.
    if (timestepHistoryActive) then
       
       ! Get current cosmic time.
       time=Tree_Node_Time(thisNode)
       
       ! Determine how long until next available timestep.
       timeIndex=Interpolate_Locate(timestepHistorySteps,historyTime,interpolationAccelerator,time)
       if (time < historyTime(timeIndex+1)) then
          ! Find next time for storage.
          ourTimeStep=historyTime(timeIndex+1)-time
          
          ! Set return value if our timestep is smaller than current one.
          if (ourTimeStep <= timeStep) then
             timeStep=ourTimeStep
             End_Of_Timestep_Task => Merger_Tree_History_Store
          end if
       end if
    end if
    return
  end subroutine Merger_Tree_Timestep_History

  subroutine Merger_Tree_History_Store(thisTree,thisNode)
    !% Store various properties in global arrays.
    use Merger_Trees
    use Tree_Nodes
    use Tree_Node_Methods
    use Numerical_Interpolation
    implicit none
    type(mergerTree), intent(in)             :: thisTree
    type(treeNode),   intent(inout), pointer :: thisNode
    integer                                  :: timeIndex
    double precision                         :: time,starFormationRate

    ! Get current cosmic time.
    time=Tree_Node_Time(thisNode)

    ! Determine how long until next available timestep.
    if (time == historyTime(timestepHistorySteps)) then
       timeIndex=timestepHistorySteps
    else
       timeIndex=Interpolate_Locate(timestepHistorySteps,historyTime,interpolationAccelerator,time)
    end if

    ! Accumulate the properties.
    starFormationRate=0.0d0
    if (diskActive) starFormationRate=starFormationRate+Tree_Node_Disk_SFR(thisNode)
    if (spheroidActive) starFormationRate=starFormationRate+Tree_Node_Spheroid_SFR(thisNode)
    historyStarFormationRate(timeIndex)=historyStarFormationRate(timeIndex)+starFormationRate*thisTree%volumeWeight

    return
  end subroutine Merger_Tree_History_Store

  !# <hdfPreCloseTask>
  !#  <unitName>Merger_Tree_History_Write</unitName>
  !# </hdfPreCloseTask>
  subroutine Merger_Tree_History_Write
    !% Store the global history data to the \glc\ output file.
    use ISO_Varying_String
    use Galacticus_HDF5_Groups
    use HDF5
    implicit none
    integer(kind=HID_T)              :: historyGroupID=0,historyDataID=0
    type(varying_string)             :: groupName,groupComment

    if (timestepHistoryActive) then
       groupName='globalHistory'
       groupComment='Global (volume averaged) history for this model.'
       historyGroupID=Galacticus_Output_Make_Group(groupName,groupComment)
       historyDataID=0
       call Galacticus_Output_Dataset(historyGroupID,historyDataID,'historyTime'                 ,'Time [Gyr]'                               ,historyTime                 )
       historyDataID=0
       call Galacticus_Output_Dataset(historyGroupID,historyDataID,'historyExpansion'            ,'Expansion factor []'                      ,historyExpansion            )
       historyDataID=0
       call Galacticus_Output_Dataset(historyGroupID,historyDataID,'historyStarFormationRate','Star formation rate [Msun/Gyr/Mpc^3]',historyStarFormationRate)
    end if
    return
  end subroutine Merger_Tree_History_Write

end module Merger_Tree_Timesteps_History
