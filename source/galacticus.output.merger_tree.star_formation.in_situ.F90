!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
!!    Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which handles computation and output of star formation histories split by in-situ and accreted star
!% formation.

module Star_Formation_Histories_In_Situ
  !% Handles computation and output of star formation histories split by in-situ and accreted star
  implicit none
  private
  public :: Star_Formation_Histories_In_Situ_Initialize, Star_Formation_History_Scales_In_Situ_Satellite_Merging

  ! Record of whether this method is active.
  logical          :: moduleActive=.false.
  
  ! Parameters controlling the tabulation of star formation rate history.
  double precision :: starFormationHistoryFineTime, starFormationHistoryFineTimeStep, &
       &              starFormationHistoryTimeStep

  ! Type used to store timestep range information.
  type timeStepRange
     integer                                  :: count
     double precision                         :: timeBegin, timeEnd
     type            (timeStepRange), pointer :: next
  end type timeStepRange

contains

  !# <starFormationHistoriesMethod>
  !#  <unitName>Star_Formation_Histories_In_Situ_Initialize</unitName>
  !# </starFormationHistoriesMethod>
  subroutine Star_Formation_Histories_In_Situ_Initialize(starFormationHistoriesMethod,Star_Formation_History_Create_Do&
       &,Star_Formation_History_Scales_Do,Star_Formation_History_Record_Do,Star_Formation_History_Output_Do)
    !% Initializes the in-situ split star formation history module.
    use ISO_Varying_String
    use Input_Parameters
    use Numerical_Ranges
    use Memory_Management
    use Galacticus_Error
    implicit none
    type     (varying_string                       ), intent(in   )          :: starFormationHistoriesMethod
    procedure(Star_Formation_History_Create_In_Situ), intent(inout), pointer :: Star_Formation_History_Create_Do
    procedure(Star_Formation_History_Scales_In_Situ), intent(inout), pointer :: Star_Formation_History_Scales_Do
    procedure(Star_Formation_History_Record_In_Situ), intent(inout), pointer :: Star_Formation_History_Record_Do
    procedure(Star_Formation_History_Output_In_Situ), intent(inout), pointer :: Star_Formation_History_Output_Do

    if (starFormationHistoriesMethod == 'inSitu') then
       ! Record that module is active.
       moduleActive=.true.
       ! Associate procedure pointers.
       Star_Formation_History_Create_Do => Star_Formation_History_Create_In_Situ
       Star_Formation_History_Scales_Do => Star_Formation_History_Scales_In_Situ
       Star_Formation_History_Record_Do => Star_Formation_History_Record_In_Situ
       Star_Formation_History_Output_Do => Star_Formation_History_Output_In_Situ
       ! Get controlling parameters.
       !@ <inputParameter>
       !@   <name>starFormationHistoryTimeStep</name>
       !@   <defaultValue>$0.1$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The time step to use in tabulations of star formation histories [Gyr].
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>output</group>
       !@ </inputParameter>
       call Get_Input_Parameter('starFormationHistoryTimeStep'          ,starFormationHistoryTimeStep          ,defaultValue=0.1d0 )
       !@ <inputParameter>
       !@   <name>starFormationHistoryFineTimeStep</name>
       !@   <defaultValue>$0.01$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The fine time step to use in tabulations of star formation histories [Gyr].
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>output</group>
       !@ </inputParameter>
       call Get_Input_Parameter('starFormationHistoryFineTimeStep'      ,starFormationHistoryFineTimeStep      ,defaultValue=0.01d0)
       !@ <inputParameter>
       !@   <name>starFormationHistoryFineTime</name>
       !@   <defaultValue>$0.1$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The period prior to each output for which the fine time step is used in tabulations of star formation histories [Gyr].
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>output</group>
       !@ </inputParameter>
       call Get_Input_Parameter('starFormationHistoryFineTime'          ,starFormationHistoryFineTime          ,defaultValue=0.1d0 )
    end if
    return
  end subroutine Star_Formation_Histories_In_Situ_Initialize

  subroutine Star_Formation_History_Create_In_Situ(thisNode,thisHistory,timeBegin)
    !% Create the history required for storing star formation history.
    use Histories
    use Galacticus_Nodes
    use Galacticus_Output_Times
    implicit none
    type            (treeNode          ), intent(inout), pointer :: thisNode
    type            (history           ), intent(inout)          :: thisHistory
    class           (nodeComponentBasic)               , pointer :: thisBasicComponent
    double precision                    , intent(in   )          :: timeBegin
    double precision                                             :: timeBeginActual   , timeEnd

    ! Find the start and end times for this history.
    thisBasicComponent => thisNode%basic()
    timeBeginActual=min(timeBegin,thisBasicComponent%time())
    timeEnd  =Galacticus_Next_Output_Time(timeBegin)
    call Star_Formation_History_In_Situ_Make_History(thisHistory,timeBeginActual,timeEnd)
    return
  end subroutine Star_Formation_History_Create_In_Situ

  subroutine Star_Formation_History_In_Situ_Make_History(thisHistory,timeBegin,timeEnd,currentTimes)
    !% Create the history required for storing star formation history.
    use Histories
    use Numerical_Ranges
    use Galacticus_Output_Times
    use Galacticus_Error
    implicit none
    type            (history      )              , intent(inout)           :: thisHistory
    double precision                             , intent(in   )           :: timeBegin       , timeEnd
    double precision               , dimension(:), intent(in   ), optional :: currentTimes
    type            (timeStepRange), pointer                               :: firstTimeStep   , nextTimeStep , thisTimeStep
    integer                                                                :: coarseTimeCount , fineTimeCount, timeCount
    logical                                                                :: gotFirstTimeStep
    double precision                                                       :: timeCoarseBegin , timeCoarseEnd, timeFineBegin, timeNext, &
         &                                                                    timeNow

    ! Exit with a null history if it would contain no time.
    if (timeEnd <= timeBegin) then
       call thisHistory%destroy()
       return
    end if
    ! If we have a set of times tabulated already, do some sanity checks.
    if (present(currentTimes)) then
       ! Complain if the beginning time is before the given list of times.
       if (timeBegin < currentTimes(1                 )) call Galacticus_Error_Report('Star_Formation_History_In_Situ_Make_History','requested begin time is before currently tabulated times')
       ! Complain if the end time is less than the maximum tabulated time.
       if (timeEnd   < currentTimes(size(currentTimes))) call Galacticus_Error_Report('Star_Formation_History_In_Situ_Make_History','requested end time is within currently tabulated times')
    end if

    ! Step through time, creating a set of timesteps as needed.
    if (present(currentTimes)) then
       timeNow       =currentTimes(size(currentTimes))
    else
       timeNow       =timeBegin
    end if
    timeCount       =0
    gotFirstTimeStep=.false.
    do while (timeNow < timeEnd)
       ! Get the time of the next output
       timeNext=Galacticus_Next_Output_Time(timeNow)
       ! Unphysical (negative) value indicates no next output.
       if (timeNext < 0.0d0 .or. timeNext > timeEnd) timeNext=timeEnd
       ! Construct coarse and fine timesteps for this output, recording the parameters of each range.
       ! Determine the number of fine timestep bins required and the time at which we begin using fine timesteps.
       if (starFormationHistoryFineTime > 0.0d0) then
          fineTimeCount  =int(min(timeNext-timeNow,starFormationHistoryFineTime)/starFormationHistoryFineTimeStep)+1
          timeFineBegin  =timeNext-starFormationHistoryFineTimeStep*dble(fineTimeCount-1)
          timeCoarseBegin=timeNow      +starFormationHistoryTimeStep
          timeCoarseEnd  =timeFineBegin-starFormationHistoryFineTimeStep
       else
          fineTimeCount  =0
          timeFineBegin  =timeNext
          timeCoarseBegin=timeNow      +starFormationHistoryTimeStep
          timeCoarseEnd  =timeNext
       end if
       ! Determine the number of coarse time bins required for this history.
       if (timeCoarseEnd > timeCoarseBegin) then
          coarseTimeCount=max(int((timeCoarseEnd-timeCoarseBegin)/starFormationHistoryTimeStep)+1,2)
       else if (fineTimeCount == 0) then
          coarseTimeCount=2
          timeCoarseBegin=(timeCoarseEnd-timeNow)/3.0d0+timeNow
       else
          coarseTimeCount=0
       end if
       ! Create the time steps.
       if (gotFirstTimeStep) then
          allocate(thisTimeStep%next)
          thisTimeStep => thisTimeStep%next
       else
          allocate(firstTimeStep)
          thisTimeStep => firstTimeStep
          if (coarseTimeCount > 0) then
             coarseTimeCount=coarseTimeCount+1
             timeCoarseBegin=max(timeCoarseBegin-starFormationHistoryTimeStep    ,0.0d0)
          else
             fineTimeCount  =fineTimeCount  +1
             timeFineBegin  =max(timeFineBegin  -starFormationHistoryFineTimeStep,0.0d0)
          end if
          gotFirstTimeStep=.true.
       end if
       if (coarseTimeCount > 0) then
          thisTimeStep%count    =  coarseTimeCount
          thisTimeStep%timeBegin=  timeCoarseBegin
          thisTimeStep%timeEnd  =  timeCoarseEnd
          allocate(thisTimeStep%next)
          thisTimeStep          => thisTimeStep%next
       end if
       if (fineTimeCount > 0) then
          thisTimeStep%count    =  fineTimeCount
          thisTimeStep%timeBegin=  timeFineBegin
          thisTimeStep%timeEnd  =  timeNext
       end if
       thisTimeStep%next => null()

       ! Increment the total number of steps required.
       timeCount=timeCount+fineTimeCount+coarseTimeCount

       ! Increment the time.
       timeNow=timeNext
    end do
    ! Shift the end point for the final step to the overall end time.
    if (gotFirstTimeStep) thisTimeStep%timeEnd=timeNext

    ! Copy in existing times if necessary.
    if (present(currentTimes)) then
       timeCount=timeCount+size(currentTimes)
       if (gotFirstTimeStep) timeCount=timeCount-1
    end if
    call thisHistory%create(2,timeCount)
    timeCount=0
    if (present(currentTimes)) then
       if (gotFirstTimeStep) then
          thisHistory%time(timeCount+1:timeCount+size(currentTimes)-1)=currentTimes(1:size(currentTimes)-1)
          timeCount=size(currentTimes)-1
       else
          thisHistory%time(timeCount+1:timeCount+size(currentTimes)  )=currentTimes(1:size(currentTimes)  )
          timeCount=size(currentTimes)
       end if
    end if

    ! Create new times if necessary.
    if (gotFirstTimeStep) then
       thisTimeStep => firstTimeStep
       do while (associated(thisTimeStep))
          ! Populate the time array.
          if      (thisTimeStep%count == 1) then
             thisHistory%time(timeCount+1)=thisTimeStep%timeEnd
          else if (thisTimeStep%count >  1) then
             thisHistory%time(timeCount+1:timeCount+thisTimeStep%count)=Make_Range(thisTimeStep%timeBegin,thisTimeStep%timeEnd&
                  &,thisTimeStep%count,rangeTypeLinear)
          end if
          timeCount=timeCount+thisTimeStep%count
          ! Jump to the next time step.
          nextTimeStep => thisTimeStep%next
          deallocate(thisTimeStep)
          thisTimeStep => nextTimeStep
       end do
    end if
    return
  end subroutine Star_Formation_History_In_Situ_Make_History

  subroutine Star_Formation_History_Record_In_Situ(thisNode,thisHistory,fuelAbundances,starFormationRate)
    !% Record the star formation history for {\normalfont \ttfamily thisNode}.
    use, intrinsic :: ISO_C_Binding
    use Histories
    use Galacticus_Nodes
    use Abundances_Structure
    use Arrays_Search
    implicit none
    type            (treeNode          ), intent(inout), pointer :: thisNode
    type            (history           ), intent(inout)          :: thisHistory
    type            (abundances        ), intent(in   )          :: fuelAbundances
    double precision                    , intent(in   )          :: starFormationRate
    class           (nodeComponentBasic)               , pointer :: thisBasicComponent
    integer                                                      :: historyCount
    integer         (c_size_t          )                         :: iHistory
    double precision                                             :: timeNode
    !GCC$ attributes unused :: fuelAbundances

    ! Get the current time for this node.
    thisBasicComponent => thisNode%basic()
    timeNode=thisBasicComponent%time()

    ! Get the number of times at which star formation rate is tabulated for this node.
    historyCount=size(thisHistory%time)

    ! Find the point in the table at which to accumulate the star formation rate.
    iHistory=Search_Array(thisHistory%time,timeNode)+1

    ! Accumulate to the appropriate time.
    thisHistory%data(iHistory,:)=starFormationRate

    return
  end subroutine Star_Formation_History_Record_In_Situ

  subroutine Star_Formation_History_Output_In_Situ(thisNode,nodePassesFilter,thisHistory,iOutput,treeIndex,componentLabel)
    !% Output the star formation history for {\normalfont \ttfamily thisNode}.
    use, intrinsic :: ISO_C_Binding
    use Histories
    use ISO_Varying_String
    use Galacticus_HDF5
    use Galacticus_Nodes
    use String_Handling
    use Kind_Numbers
    use Galacticus_Output_Times
    implicit none
    type            (treeNode          ), intent(inout), pointer :: thisNode
    logical                             , intent(in   )          :: nodePassesFilter
    type            (history           ), intent(inout)          :: thisHistory
    integer         (c_size_t          ), intent(in   )          :: iOutput
    integer         (kind=kind_int8    ), intent(in   )          :: treeIndex
    character       (len=*             ), intent(in   )          :: componentLabel
    class           (nodeComponentBasic)               , pointer :: parentBasicComponent
    type            (treeNode          )               , pointer :: parentNode
    double precision                                             :: timeBegin           , timeEnd
    type            (varying_string    )                         :: groupName
    type            (hdf5Object        )                         :: historyGroup        , outputGroup, treeGroup
    type            (history           )                         :: newHistory

    ! Return if the history does not exist.
    if (.not.thisHistory%exists()) return

    ! Check if the node passes any filtering, and output it if it does.
    if (nodePassesFilter) then
       !$omp critical(HDF5_Access)
       ! Create a group for the profile datasets.
       historyGroup=galacticusOutputFile%openGroup("starFormationHistories","Star formation history data.")
       groupName="Output"
       groupName=groupName//iOutput
       outputGroup=historyGroup%openGroup(char(groupName),"Star formation histories for all trees at each output.")
       groupName="mergerTree"
       groupName=groupName//treeIndex
       treeGroup=outputGroup%openGroup(char(groupName),"Star formation histories for each tree.")
       ! Write dataset to the group.
       groupName=trim(componentLabel)//"Time"
       groupname=groupName//thisNode%index()
       call treeGroup%writeDataset(thisHistory%time,char(groupName),"Star formation history times of the "//trim(componentLabel)//" component.")
       groupName=trim(componentLabel)//"SFH"
       groupname=groupName//thisNode%index()
       call treeGroup%writeDataset(thisHistory%data,char(groupName),"Star formation history stellar masses of the "//trim(componentLabel)//" component.")
       ! Close the star formation history group.
       call treeGroup   %close()
       call outputGroup %close()
       call historyGroup%close()
       !$omp end critical(HDF5_Access)
    end if

    timeBegin=thisHistory%time(1)
    if (iOutput < Galacticus_Output_Time_Count()) then
       timeEnd =Galacticus_Output_Time(iOutput+1)
    else
       parentNode => thisNode
       do while (associated(parentNode%parent))
          parentNode => parentNode%parent
       end do
       parentBasicComponent => parentNode%basic()
       timeEnd=parentBasicComponent%time()
    end if
    call Star_Formation_History_In_Situ_Make_History(newHistory,timeBegin,timeEnd,thisHistory%time)
    newHistory%data(1:size(thisHistory%time),:)=thisHistory%data(:,:)
    call thisHistory%destroy()
    thisHistory=newHistory
    call newHistory %destroy(recordMemory=.false.)
    return
  end subroutine Star_Formation_History_Output_In_Situ

  subroutine Star_Formation_History_Scales_In_Situ(thisHistory,stellarMass,stellarAbundances)
    !% Set the scalings for error control on the absolute values of star formation histories.
    use Histories
    use Abundances_Structure
    use Memory_Management
    implicit none
    double precision            , intent(in   )               :: stellarMass
    type            (abundances), intent(in   )               :: stellarAbundances
    type            (history   ), intent(inout)               :: thisHistory
    double precision            , parameter                   :: stellarMassMinimum=1.0d0
    double precision            , allocatable  , dimension(:) :: timeSteps
    integer                                                   :: i
    !GCC$ attributes unused :: stellarAbundances
    
    ! Return immediately if the history does not exist.
    if (.not.thisHistory%exists()) return

    ! Get timesteps.
    call thisHistory%timeSteps(timeSteps)

    ! Set scaling factors for star formation rate.
    forall(i=1:2)
       thisHistory%data(:,i)=max(stellarMass,stellarMassMinimum)/timeSteps
    end forall
  
    ! Destroy temporary array.
    call Dealloc_Array(timeSteps)

    return
  end subroutine Star_Formation_History_Scales_In_Situ

  !# <satelliteMergerTask>
  !#  <unitName>Star_Formation_History_Scales_In_Situ_Satellite_Merging</unitName>
  !#  <before>re:Node_Component_.+_Satellite_Merging</before>
  !# </satelliteMergerTask>
  subroutine Star_Formation_History_Scales_In_Situ_Satellite_Merging(thisNode)
    !% Zero any in-situ star formation history for galaxy about to merge.
    use Galacticus_Nodes
    use Histories
    implicit none
    type (treeNode             ), intent(inout), pointer :: thisNode
    class(nodeComponentDisk    )               , pointer :: disk
    class(nodeComponentSpheroid)               , pointer :: spheroid
    type (history              )                         :: diskStarFormationHistory, spheroidStarformationHistory
    
    if (moduleActive) then
       disk                                   => thisNode%disk                ()
       spheroid                               => thisNode%spheroid            ()
       diskStarFormationHistory               =  disk    %starFormationHistory()
       spheroidStarformationHistory           =  spheroid%starFormationHistory()
       if (diskStarFormationHistory    %exists()) then
          diskStarFormationHistory    %data(:,1)=0.0d0
          call disk    %starFormationHistorySet(    diskStarFormationHistory)
       end if
       if (spheroidStarFormationHistory%exists()) then
          spheroidStarFormationHistory%data(:,1)=0.0d0
          call spheroid%starFormationHistorySet(spheroidStarFormationHistory)
       end if
    end if
    return
  end subroutine Star_Formation_History_Scales_In_Situ_Satellite_Merging
    
end module Star_Formation_Histories_In_Situ
