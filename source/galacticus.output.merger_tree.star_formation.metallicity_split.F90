!! Copyright 2009, 2010, Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


!% Contains a module which handles computation and output of star formation histories split by metallicity for galaxies.

module Star_Formation_Histories_Metallicity_Split
  !% Handles computation and output of star formation histories split by metallicity for galaxies.
  private
  public :: Star_Formation_Histories_Metallicity_Split_Initialize

  ! Parameters controlling the tabulation of star formation rate history.
  integer                                     :: starFormationHistoryMetallicityCount
  double precision                            :: starFormationHistoryTimeStep,starFormationHistoryFineTimeStep&
       &,starFormationHistoryFineTime,starFormationHistoryMetallicityMinimum,starFormationHistoryMetallicityMaximum

  ! Array of metallicities at which to tabulate.
  double precision, allocatable, dimension(:) :: metallicityTable
  double precision, parameter                 :: metallicityInfinite=1.0d30 ! Effective infinite metallicity.

  ! Flag indicating if metallicity data has been written.
  logical                                     :: metallicityTableWritten=.false.

  ! Type used to store timestep range information.
  type timeStepRange
     integer                      :: count
     double precision             :: timeBegin,timeEnd
     type(timeStepRange), pointer :: next
  end type timeStepRange

contains

  !# <starFormationHistoriesMethod>
  !#  <unitName>Star_Formation_Histories_Metallicity_Split_Initialize</unitName>
  !# </starFormationHistoriesMethod>
  subroutine Star_Formation_Histories_Metallicity_Split_Initialize(starFormationHistoriesMethod,Star_Formation_History_Create_Do&
       &,Star_Formation_History_Scales_Do,Star_Formation_History_Record_Do,Star_Formation_History_Output_Do)
    !% Initializes the metallicity split star formation history module.
    use ISO_Varying_String
    use Input_Parameters
    use Numerical_Ranges
    use Memory_Management
    use Galacticus_Error
    implicit none
    type(varying_string),          intent(in)    :: starFormationHistoriesMethod
    procedure(),          pointer, intent(inout) :: Star_Formation_History_Create_Do,Star_Formation_History_Scales_Do&
         &,Star_Formation_History_Record_Do ,Star_Formation_History_Output_Do
    
    if (starFormationHistoriesMethod == 'metallicity split') then
       ! Associate procedure pointers.
       Star_Formation_History_Create_Do => Star_Formation_History_Create_Metallicity_Split
       Star_Formation_History_Scales_Do => Star_Formation_History_Scales_Metallicity_Split
       Star_Formation_History_Record_Do => Star_Formation_History_Record_Metallicity_Split
       Star_Formation_History_Output_Do => Star_Formation_History_Output_Metallicity_Split

       ! Get controlling parameters.
       !@ <inputParameter>
       !@   <name>starFormationHistoryTimeStep</name>
       !@   <defaultValue>$0.1$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The time step to use in tabulations of star formation histories [Gyr].
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('starFormationHistoryTimeStep'          ,starFormationHistoryTimeStep          ,defaultValue=0.1d0 )
       !@ <inputParameter>
       !@   <name>starFormationHistoryFineTimeStep</name>
       !@   <defaultValue>$0.01$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The fine time step to use in tabulations of star formation histories [Gyr].
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('starFormationHistoryFineTimeStep'      ,starFormationHistoryFineTimeStep      ,defaultValue=0.01d0)
       !@ <inputParameter>
       !@   <name>starFormationHistoryFineTime</name>
       !@   <defaultValue>$0.1$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The period prior to each output for which the fine time step is used in tabulations of star formation histories [Gyr].
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('starFormationHistoryFineTime'          ,starFormationHistoryFineTime          ,defaultValue=0.1d0 )
       !@ <inputParameter>
       !@   <name>starFormationHistoryMetallicityCount</name>
       !@   <defaultValue>10</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The number of bins in metallicity to use when tabulating star formation histories.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('starFormationHistoryMetallicityCount'  ,starFormationHistoryMetallicityCount  ,defaultValue=10    )
       !@ <inputParameter>
       !@   <name>starFormationHistoryMetallicityMinimum</name>
       !@   <defaultValue>$10^{-4}$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The upper limit to the metallicity in the lowest metallicity bin when tabulating star formation histories [Solar units].
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('starFormationHistoryMetallicityMinimum',starFormationHistoryMetallicityMinimum,defaultValue=1.0d-4)
       !@ <inputParameter>
       !@   <name>starFormationHistoryMetallicityMaximum</name>
       !@   <defaultValue>$10^{1}$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The upper limit to the metallicity in the highest metallicity bin when tabulating star formation histories [Solar units].
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('starFormationHistoryMetallicityMaximum',starFormationHistoryMetallicityMaximum,defaultValue=1.0d+1)

       ! Construct a table of metallicities at which to tabulate. Add an extra bin since we want to catch all metallicities,
       ! including those below and above the maximum. A single bin is not allowed, but zero bins implies that no metallicity
       ! resolution is required.
       select case (starFormationHistoryMetallicityCount)
       case (:-1,1)
          call Galacticus_Error_Report('Star_Formation_Histories_Metallicity_Split_Initialize','number of bins must be 0, or greater than 1')
       case default
          call Alloc_Array(metallicityTable,[starFormationHistoryMetallicityCount+1])
          if (starFormationHistoryMetallicityCount > 1) then
             metallicityTable(1:starFormationHistoryMetallicityCount)=Make_Range(starFormationHistoryMetallicityMinimum&
                  &,starFormationHistoryMetallicityMaximum ,starFormationHistoryMetallicityCount,rangeType=rangeTypeLogarithmic)
          end if
          metallicityTable(starFormationHistoryMetallicityCount+1)=metallicityInfinite
       end select
    end if
    return
  end subroutine Star_Formation_Histories_Metallicity_Split_Initialize

  subroutine Star_Formation_History_Create_Metallicity_Split(thisNode,thisHistory)
    !% Create the history required for storing star formation history.
    use Histories
    use Tree_Nodes
    use Galacticus_Output_Times
    implicit none
    type(treeNode),      intent(inout), pointer :: thisNode
    type(history),       intent(inout)          :: thisHistory
    double precision                            :: timeBegin,timeEnd

    ! Find the start and end times for this history.
    timeBegin=Tree_Node_Time(thisNode  )
    timeEnd  =Galacticus_Next_Output_Time(timeBegin)

    call Star_Formation_History_Metallicity_Split_Make_History(thisHistory,timeBegin,timeEnd)
    return
  end subroutine Star_Formation_History_Create_Metallicity_Split
  
  subroutine Star_Formation_History_Metallicity_Split_Make_History(thisHistory,timeBegin,timeEnd)
    !% Create the history required for storing star formation history.
    use Histories
    use Numerical_Ranges
    use Tree_Nodes
    use Galacticus_Output_Times
    implicit none
    type(history),       intent(inout) :: thisHistory
    double precision,    intent(in)    :: timeBegin,timeEnd
    type(timeStepRange), pointer       :: firstTimeStep,thisTimeStep,nextTimeStep
    integer                            :: timeCount,coarseTimeCount,fineTimeCount
    logical                            :: gotFirstTimeStep
    double precision                   :: timeFineBegin,timeCoarseBegin,timeCoarseEnd,timeNow,timeNext

    ! Exit with a null history if it would contain no time.
    if (timeEnd <= timeBegin) then
       call thisHistory%destroy()
       return
    end if

    ! Step through time, creating a set of timesteps as needed.
    timeNow         =timeBegin
    timeCount       =0
    gotFirstTimeStep=.false.
    allocate(firstTimeStep)
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
    thisTimeStep%timeEnd=timeNext

    ! Create the history.
    call thisHistory%create(starFormationHistoryMetallicityCount+1,timeCount)

    thisTimeStep => firstTimeStep
    timeCount=0
    do while (associated(thisTimeStep))
       ! Populate the time array.
       if      (thisTimeStep%count == 1) then
          thisHistory%time(timeCount+1)=thisTimeStep%timeEnd
       else if (thisTimeStep%count >  1) then
          thisHistory%time(timeCount+1:timeCount+thisTimeStep%count)=Make_Range(thisTimeStep%timeBegin,thisTimeStep%timeEnd,thisTimeStep%count,rangeTypeLinear)
       end if
       timeCount=timeCount+thisTimeStep%count
       ! Jump to the next time step.
       nextTimeStep => thisTimeStep%next
       deallocate(thisTimeStep)
       thisTimeStep => nextTimeStep
    end do

    return
  end subroutine Star_Formation_History_Metallicity_Split_Make_History

  subroutine Star_Formation_History_Record_Metallicity_Split(thisNode,thisHistory,fuelAbundances,starFormationRate)
    !% Record the star formation history for {\tt thisNode}.
    use Histories
    use Numerical_Ranges
    use Tree_Nodes
    use Abundances_Structure
    use Arrays_Search
    implicit none
    type(treeNode),            intent(inout), pointer :: thisNode
    type(history),             intent(inout)          :: thisHistory
    type(abundancesStructure), intent(in)             :: fuelAbundances
    double precision,          intent(in)             :: starFormationRate
    integer                                           :: iHistory,iMetallicity,historyCount
    double precision                                  :: timeNode,fuelMetallicity

    ! Get the current time for this node.
    timeNode    =Tree_Node_Time(thisNode)

    ! Get the number of times at which star formation rate is tabulated for this node.
    historyCount=size(thisHistory%time)

    ! Find the point in the table at which to accumulate the star formation rate.
    iHistory=Search_Array(thisHistory%time,timeNode)+1

    ! Find the metallicity bin to accumulate to.
    fuelMetallicity=fuelAbundances%metallicity(metallicityType=linearByMassSolar)
    if (fuelMetallicity < metallicityTable(1) .or. starFormationHistoryMetallicityCount == 0) then
       iMetallicity=1
    else
       iMetallicity=Search_Array(metallicityTable,fuelMetallicity)+1
    end if

    ! Accumulate to all future times.
    thisHistory%rates(iHistory,iMetallicity)=thisHistory%rates(iHistory,iMetallicity)+starFormationRate

    return
  end subroutine Star_Formation_History_Record_Metallicity_Split

  subroutine Star_Formation_History_Output_Metallicity_Split(thisNode,thisHistory,iOutput,treeIndex,componentLabel)
    !% Output the star formation history for {\tt thisNode}.
    use Histories
    use ISO_Varying_String
    use Galacticus_HDF5
    use IO_HDF5
    use Tree_Nodes
    use String_Handling
    use Kind_Numbers
    use Galacticus_Output_Times
    implicit none
    type(treeNode),          intent(inout), pointer :: thisNode
    type(history),           intent(inout)          :: thisHistory
    integer,                 intent(in)             :: iOutput
    integer(kind=kind_int8), intent(in)             :: treeIndex
    character(len=*),        intent(in)             :: componentLabel
    type(treeNode),                         pointer :: parentNode
    double precision                                :: timeBegin,timeEnd
    type(varying_string)                            :: groupName
    type(hdf5Object)                                :: historyGroup,treeGroup,outputGroup
    type(history)                                   :: newHistory
 
    ! Return if the history does not exist.
    if (.not.thisHistory%exists()) return

    ! Write metallicities if not already done.
    if (.not.metallicityTableWritten) then
       ! Open the histories group.
       historyGroup=IO_HDF5_Open_Group(galacticusOutputFile,"starFormationHistories","Star formation history data.")
       
       ! Write the metallicities.
       call historyGroup%writeDataset(metallicityTable,"metallicities","Metallicities at which star formation histories are tabulated.")
       
       ! Close the history group.
       call historyGroup%close()
       
       ! Flag that metallicities have been written.
       metallicityTableWritten=.true.
    end if

    ! Create a group for the profile datasets.
    historyGroup=IO_HDF5_Open_Group(galacticusOutputFile,"starFormationHistories","Star formation history data.")
    groupName="Output"
    groupName=groupName//iOutput
    outputGroup=IO_HDF5_Open_Group(historyGroup,char(groupName),"Star formation histories for all trees at each output.")
    groupName="mergerTree"
    groupName=groupName//treeIndex
    treeGroup=IO_HDF5_Open_Group(outputGroup,char(groupName),"Star formation histories for each tree.")
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

    timeBegin=thisHistory%time(1)
    if (iOutput < Galacticus_Output_Time_Count()) then
       timeEnd  =Galacticus_Output_Time(iOutput+1)
    else
       parentNode => thisNode
       do while (associated(parentNode%parentNode))
          parentNode => parentNode%parentNode
       end do
       timeEnd=Tree_Node_Time(parentNode)
    end if
    call Star_Formation_History_Metallicity_Split_Make_History(newHistory,timeBegin,timeEnd)
    newHistory%data(1:size(thisHistory%time),:)=thisHistory%data(:,:)
    call thisHistory%destroy()
    thisHistory=newHistory
    call newHistory %destroy()
    return
  end subroutine Star_Formation_History_Output_Metallicity_Split

  subroutine Star_Formation_History_Scales_Metallicity_Split(thisHistory,stellarMass,stellarAbundances)
    !% Set the scalings for error control on the absolute values of star formation histories.
    use Histories
    use Stellar_Feedback
    use Abundances_Structure
    use Memory_Management
    implicit none
    double precision,          intent(in)                :: stellarMass
    type(abundancesStructure), intent(in)                :: stellarAbundances
    type(history),             intent(inout)             :: thisHistory
    double precision,          parameter                 :: stellarMassMinimum=1.0d0
    double precision,          dimension(:), allocatable :: timeSteps
    integer                                              :: iMetallicity

    ! Return immediately if the history does not exist.
    if (.not.thisHistory%exists()) return

    ! Get timesteps.
    call thisHistory%timeSteps(timeSteps)

    ! Set scaling factors for recycled mass.
    forall(iMetallicity=1:starFormationHistoryMetallicityCount+1)
       thisHistory%scales(:,iMetallicity)=max(stellarMass,stellarMassMinimum)/timeSteps
    end forall

    ! Destroy temporary array.
    call Dealloc_Array(timeSteps)

    return
  end subroutine Star_Formation_History_Scales_Metallicity_Split
  
end module Star_Formation_Histories_Metallicity_Split
