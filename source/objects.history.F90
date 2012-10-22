!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module defining the history object type.

module Histories
  !% Defines the history object type.
  implicit none
  private
  public :: history, History_Set_Times, Histories_State_Store, Histories_State_Retrieve

  type history
     !% The history object type.
     double precision, allocatable, dimension(:  ) :: time
     double precision, allocatable, dimension(:,:) :: data
     integer                                       :: rangeType
   contains
     ! Operators.
     procedure                 :: add                    => History_Add
     procedure                 :: subtract               => History_Subtract
     procedure                 :: divide                 => History_Divide
     generic                   :: operator(+)            => add
     generic                   :: operator(-)            => subtract
     generic                   :: operator(/)            => divide
     procedure                 :: isZero                 => History_Is_Zero
     !@ <objectMethods>
     !@   <object>history</object>
     !@   <objectMethod>
     !@     <method>create</method>
     !@     <description>Creates a history object with a specified range of times.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>dump</method>
     !@     <description>Dump a history object.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>clone</method>
     !@     <description>Clone a history object.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>destroy</method>
     !@     <description>Destroys a history object.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>trim</method>
     !@     <description>Removes any times in a history which have become outdated.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>increment</method>
     !@     <description>Increments a history.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>addRates</method>
     !@     <description>Adds two histories.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>combine</method>
     !@     <description>Combines two histories.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>extend</method>
     !@     <description>Extends the time range of a history to encompass the specified limits.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>reset</method>
     !@     <description>Resets all entries in a history to zero.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>setToUnity</method>
     !@     <description>Set all entries in a history to unity.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>exists</method>
     !@     <description>Returns true if the given history has been created.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>timeSteps</method>
     !@     <description>Returns an array with the timesteps (i.e. the intervals between successive times) in the given history.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure :: dump       => History_Dump
     procedure :: increment  => History_Increment
     procedure :: create     => History_Create
     procedure :: clone      => History_Clone
     procedure :: destroy    => History_Destroy
     procedure :: trim       => History_Trim
     procedure :: extend     => History_Extend
     procedure :: addRates   => History_Add_Rates
     procedure :: combine    => History_Combine
     procedure :: reset      => History_Reset
     procedure :: setToUnity => History_Set_To_Unity
     procedure :: exists     => History_Exists
     procedure :: timeSteps  => History_Timesteps
     !@   <objectMethod>
     !@     <method>serializeCount</method>
     !@     <description>Return a count of the number of properties in a serialized history object.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>serialize</method>
     !@     <description>Serialize a history object to an array.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>deserialize</method>
     !@     <description>Deserialize a history object from an array.</description>
     !@   </objectMethod>
     procedure                 :: serializeCount         => History_Serialize_Count
     procedure                 :: serialize              => History_Serialize
     procedure                 :: deserialize            => History_Deserialize
  end type history

  ! A null history object.
  type(history),    public            :: nullHistory

  ! Earliest and latest times for history storage.
  double precision, public            :: historyStorageEarliestTime= 0.1d0
  double precision, public            :: historyStorageLatestTime  =15.0d0

  ! Labels for targets when adding to histories.
  integer,          public, parameter :: historyData  =1
  integer,          public, parameter :: historyRates =2
  integer,          public, parameter :: historyScales=3

contains

  subroutine History_Set_Times(timeEarliest,timeLatest)
    !% Extend the range of history times to include the given {\tt timeEarliest} and {\tt timeLatest}.
    implicit none
    double precision, intent(in), optional :: timeEarliest,timeLatest

    if (present(timeEarliest)) historyStorageEarliestTime=min(historyStorageEarliestTime,timeEarliest)
    if (present(timeLatest  )) historyStorageLatestTime  =max(historyStorageLatestTime  ,timeLatest  )
    return
  end subroutine History_Set_Times

  subroutine History_Create(thisHistory,historyCount,timesCount,timeBegin,timeEnd,rangeType)
    !% Create a history object.
    use Memory_Management
    use Numerical_Ranges
    use Galacticus_Error
    implicit none
    class(history),   intent(inout)        :: thisHistory
    integer,          intent(in)           :: historyCount,timesCount
    double precision, intent(in), optional :: timeBegin,timeEnd
    integer,          intent(in), optional :: rangeType
    integer                                :: rangeTypeActual

    if (allocated(thisHistory%time)) then
       call Galacticus_Error_Report('History_Create','this history appears to have been created already')
    else
       allocate(thisHistory%time  (timesCount             ))
       allocate(thisHistory%data  (timesCount,historyCount))
       if (timesCount > 0) then
          call Memory_Usage_Record(sizeof(thisHistory%time)+sizeof(thisHistory%data),memoryType=memoryTypeNodes,blockCount=3)
          if (present(timeBegin)) then
             if (.not.present(timeEnd)) call Galacticus_Error_Report('History_Create','an end time must be given if a begin time is given')
             if (present(rangeType)) then
                rangeTypeActual=rangeType
             else
                rangeTypeActual=rangeTypeLogarithmic
             end if
             thisHistory%time     =Make_Range(timeBegin,timeEnd,timesCount,rangeTypeActual)
             thisHistory%rangeType=rangeTypeActual
          else
             thisHistory%rangeType=rangeTypeUndefined
             thisHistory%time=0.0d0
          end if
          if (historyCount > 0) thisHistory%data=0.0d0
       end if
    end if
    return
  end subroutine History_Create

  subroutine History_Destroy(thisHistory,recordMemory)
    !% Destroy a history.
    use Memory_Management
    implicit none
    class(history), intent(inout)          :: thisHistory
    logical,        intent(in),   optional :: recordMemory
    integer                                :: timesCount,historyCount
    logical                                :: recordMemoryActual

    if (allocated(thisHistory%time)) then
       timesCount  =size(thisHistory%time      )
       historyCount=size(thisHistory%data,dim=2)
       if (present(recordMemory)) then
          recordMemoryActual=recordMemory
       else
          recordMemoryActual=.true.
       end if
       if (recordMemoryActual) call Memory_Usage_Record(                             &
            &                                             sizeof(thisHistory%time  ) &
            &                                            +sizeof(thisHistory%data  ) &
            &                                           ,memoryType=memoryTypeNodes  &
            &                                           ,addRemove =-1               &
            &                                           ,blockCount= 3               &
            &                                          )
       deallocate(thisHistory%time)
       if (allocated(thisHistory%data)) deallocate(thisHistory%data)
    end if
    return
  end subroutine History_Destroy

  subroutine History_Dump(self)
    !% Dumps a history object.
    use Galacticus_Display
    use ISO_Varying_String
    implicit none
    class    (history       ), intent(in   ) :: self
    integer                                  :: i,j
    type     (varying_string)                :: message
    character(len=12        )                :: label

    if (allocated(self%time)) then
       do i=1,size(self%time)
          write (label,'(i3)') i
          message="("//trim(label)//") "
          write (label,'(e12.6)') self%time(i)
          message=message//label//" :"
          do j=1,size(self%data,dim=2)
             write (label,'(e12.6)') self%data(i,j)
             message=message//" "//label
          end do
          call Galacticus_Display_Message(message)
       end do
    end if
    return
  end subroutine History_Dump

  subroutine History_Reset(thisHistory)
    !% Reset a history by zeroing all elements, but leaving the structure (and times) intact.
    use Memory_Management
    implicit none
    class(history), intent(inout) :: thisHistory
    
    if (allocated(thisHistory%time)) thisHistory%data=0.0d0
    return
  end subroutine History_Reset

  subroutine History_Set_To_Unity(thisHistory)
    !% Reset a history by zeroing all elements, but leaving the structure (and times) intact.
    use Memory_Management
    implicit none
    class(history), intent(inout) :: thisHistory
    
    if (allocated(thisHistory%time)) thisHistory%data=1.0d0
    return
  end subroutine History_Set_To_Unity

  logical function History_Exists(thisHistory)
    !% Returns true if the history has been created.
    implicit none
    class(history), intent(in) :: thisHistory
    
    History_Exists=allocated(thisHistory%time)
    return
  end function History_Exists

  subroutine History_Clone(self,historyToClone)
    !% Clone a history object.
    use Memory_Management
    implicit none
    class(history), intent(inout) :: self
    type (history), intent(in   ) :: historyToClone

    if (allocated(self%time)) call Dealloc_Array(self%time,memoryType=memoryTypeNodes)
    if (allocated(self%data)) call Dealloc_Array(self%data,memoryType=memoryTypeNodes)
    if (allocated(historyToClone%time)) then
       call Alloc_Array(self%time,shape(historyToClone%time),memoryType=memoryTypeNodes)
       self%time=historyToClone%time
    end if
    if (allocated(historyToClone%data)) then
       call Alloc_Array(self%data,shape(historyToClone%data),memoryType=memoryTypeNodes)
       self%data=historyToClone%data
    end if
    self%rangeType=historyToClone%rangeType
    return
  end subroutine History_Clone

  logical function History_Is_Zero(self)
    !% Test whether a history object is all zero.
    implicit none
    class(history), intent(in) :: self

    History_Is_Zero=.true.
    if (allocated(self%data)) then
       if (any(self%data /= 0.0d0)) History_Is_Zero=.false.
    end if
    return
  end function History_Is_Zero
  
  function History_Add(history1,history2)
    !% Add two history objects.
    use Galacticus_Error
    implicit none
    type (history)                       :: History_Add
    class(history), intent(in)           :: history1
    class(history), intent(in), optional :: history2

    ! Clone the first history.
    if (.not.(history1%exists().or.history2%exists())) then
       call History_Add%reset()
    else
       select type (history1)
       type is (history)
          History_Add=history1
          if (present(history2)) then
             if (any(shape(History_Add%data) /= shape(history2%data))) call Galacticus_Error_Report('History_Add','mismatch in history object shape')
             History_Add%data=History_Add%data+history2%data
          end if
       end select
    end if
    return
  end function History_Add

  subroutine History_Increment(self,increment)
    !% Increment a history.
    use Galacticus_Error
    implicit none
    class(history), intent(inout) :: self
    class(history), intent(in   ) :: increment

    if (any(shape(self%data) /= shape(increment%data))) call Galacticus_Error_Report('History_Increment','mismatch in history object shape')
    self%data=self%data+increment%data
    return
  end subroutine History_Increment

  function History_Subtract(history1,history2)
    !% Subtract two history objects.
    use Galacticus_Error
    implicit none
    type (history)                       :: History_Subtract
    class(history), intent(in)           :: history1
    class(history), intent(in), optional :: history2

    if (present(history2)) then
       if (any(shape(history1%data) /= shape(history2%data))) call Galacticus_Error_Report('History_Subtract','mismatch in history object shape')
       History_Subtract%data= history1%data-history2%data
    else
       History_Subtract%data=-history1%data
    end if
    return
  end function History_Subtract

  integer function History_Serialize_Count(self)
    !% Return the number of properties required to track a history.
    implicit none
    class(history), intent(in) :: self

    if (allocated(self%data)) then
       History_Serialize_Count=size(self%data)
    else
       History_Serialize_Count=0
    end if
    return
  end function History_Serialize_Count

  subroutine History_Deserialize(self,historyArray)
    !% Pack history from an array into a history structure.
    implicit none
    class(history)  , intent(inout)               :: self
    double precision, intent(in   ), dimension(:) :: historyArray

    ! Extract data from array.
    if (allocated(self%data)) self%data=reshape(historyArray,shape(self%data))
    return
  end subroutine History_Deserialize

  subroutine History_Serialize(self,historyArray)
    !% Pack history from an array into an history structure.
    implicit none
    class(history)  , intent(in   )               :: self
    double precision, intent(  out), dimension(:) :: historyArray(:)

    ! Place data into array.
    if (allocated(self%data)) historyArray(:)=reshape(self%data,shape(historyArray))
    return
  end subroutine History_Serialize

  subroutine History_Trim(thisHistory,currentTime,minimumPointsToRemove)
    !% Removes outdated information from ``future histories'' (i.e. histories that store data for future reference). Removes all
    !% but one entry prior to the given {\tt currentTime} (this allows for interpolation of the history to the current
    !% time). Optionally, the remove is done only if it will remove more than {\tt minimumPointsToRemove} entries (since the
    !% removal can be slow this allows for some optimization).
    use Galacticus_Error
    use Memory_Management
    use, intrinsic :: ISO_C_Binding 
    implicit none
    class(history),   intent(inout)        :: thisHistory
    double precision, intent(in)           :: currentTime
    integer,          intent(in), optional :: minimumPointsToRemove
    type(history)                          :: temporaryHistory
    integer                                :: minimumPointsToRemoveActual,currentPointCount,iTrim,historyCount,newPointCount

    ! Return if no history exists.
    if (.not.allocated(thisHistory%time)) return

    ! Find points to remove.
    currentPointCount=size(thisHistory%time)

    ! Return is nothing to trim.
    if (currentPointCount == 0) return

    ! Decide on the minimum number of points that we will remove.
    if (present(minimumPointsToRemove)) then
       if (minimumPointsToRemove < 1) call Galacticus_Error_Report('History_Trim','minimum number of points to remove must be >= 1')
       minimumPointsToRemoveActual=minimumPointsToRemove
    else
       minimumPointsToRemoveActual=1
    end if

    ! Find how much we can trim. Never trim the final two point as they might be needed to extrapolate beyond the end of the
    ! future history. Having found a point which exceeds the current time, pull back two points, so that we leave one point prior
    ! to the current time.
    iTrim=1
    do while (thisHistory%time(iTrim) < currentTime .and. iTrim <= currentPointCount-2)
       iTrim=iTrim+1
    end do
    iTrim=iTrim-2

    ! Check if there are enough removable points to warrant actually doing the removal.
    if (iTrim >= minimumPointsToRemoveActual) then
       ! Move current history to temporary storage.
       call Move_Alloc(thisHistory%time  ,temporaryHistory%time  )
       call Move_Alloc(thisHistory%data  ,temporaryHistory%data  )
       ! Reallocate the history arrays.
       newPointCount=currentPointCount-iTrim
       historyCount =size(temporaryHistory%data,dim=2)
       allocate(thisHistory%time(newPointCount             ))
       allocate(thisHistory%data(newPointCount,historyCount))
       ! Copy the data back into the new arrays.
       thisHistory%time(:  )=temporaryHistory%time(iTrim+1:currentPointCount  )
       thisHistory%data(:,:)=temporaryHistory%data(iTrim+1:currentPointCount,:)
       ! Deallocate the temporary arrays.
       deallocate(temporaryHistory%time  )
       deallocate(temporaryHistory%data  )
       ! Account for change in memory usage.
       call Memory_Usage_Record(int(iTrim*(1+historyCount),C_SIZE_T)*sizeof(thisHistory%time(1)),memoryType=memoryTypeNodes,addRemove=-1,blockCount=0)
    end if
    return
  end subroutine History_Trim

   subroutine History_Add_Rates(thisHistory,addHistory)
     !% Adds the data in {\tt addHistory} to that in {\tt thisHistory}. This function is designed for histories that track
     !% instantaneous rates. The rates in {\tt addHistory} are interpolated to the times in {\tt thisHistory} and added to the
     !% rates in {\tt thisHistory}.
     use FGSL
     use Numerical_Interpolation
     use Galacticus_Error
     implicit none
     class(history          ), intent(inout) :: thisHistory
     type (history          ), intent(in   ) :: addHistory
     integer                                 :: addHistoryPointCount,iPoint,interpolationPoint,iHistory
     double precision                        :: interpolationFactors(2)
     type (fgsl_interp_accel)                :: interpolationAccelerator
     logical                                 :: interpolationReset

     select type (thisHistory)
     type is (history)
       
        ! Return if addHistory does not exist.
        if (.not.allocated(addHistory%time)) return
       
        ! Get size of addHistory.
        addHistoryPointCount=size(addHistory%time)
       
        ! Return if addHistory has zero size.
        if (addHistoryPointCount == 0) return
       
        ! If thisHistory does not exist, just replace it with addHistory.
        if (.not.allocated(thisHistory%time)) then
           call thisHistory%destroy()
           thisHistory=addHistory
           return
        end if
       
        ! If thisHistory has zero size, just replace it with addHistory.
        if (size(thisHistory%time) == 0) then
           call thisHistory%destroy()
           thisHistory=addHistory
           return
        end if
       
        ! addHistory must have at least two points to permit interpolation.
        if (addHistoryPointCount  < 2) call Galacticus_Error_Report('History_Add','history to add must have at least two points')
       
        ! The two objects must contain the same number of histories.
        if (size(thisHistory%data,dim=2) /= size(addHistory%data,dim=2)) call Galacticus_Error_Report('History_Add','two objects contain differing numbers of histories')
      
        ! Loop over each entry in thisHistory.
        interpolationReset=.true.
        do iPoint=1,size(thisHistory%time)
          
           ! If within range of history spanned by addHistory then....
           if (thisHistory%time(iPoint) >= addHistory%time(1) .and. thisHistory%time(iPoint) <= addHistory%time(addHistoryPointCount)) then
             
             ! Interpolate addHistory to point in thisHistory.
              interpolationPoint  =Interpolate_Locate(addHistoryPointCount,addHistory%time,interpolationAccelerator&
                   &,thisHistory%time(iPoint),interpolationReset)
              interpolationFactors=Interpolate_Linear_Generate_Factors(addHistoryPointCount,addHistory%time,interpolationPoint&
                   &,thisHistory%time(iPoint))
           
              ! Add them.
              forall(iHistory=1:size(thisHistory%data,dim=2))
                 thisHistory%data (iPoint,iHistory)=thisHistory%data (iPoint,iHistory)+addHistory%data(interpolationPoint,iHistory)&
                      &*interpolationFactors(1)+addHistory%data(interpolationPoint+1,iHistory)*interpolationFactors(2)
              end forall
           end if
        
        end do

     end select

     return
   end subroutine History_Add_Rates

   subroutine History_Combine(thisHistory,combineHistory)
     !% Combines the data in {\tt combineHistory} with that in {\tt thisHistory}. This function is designed for histories that
     !% track integrated quantities (such as total mass of stars formed in a time interval for example). {\tt thisHistory} will be
     !% extended if necessary to span the range of {\tt combineHistory}. Then, the data from {\tt combineHistory} will be added to
     !% that in {\tt thisHistory} by finding the fraction of each timestep in {\tt combineHistory} that overlaps with each timestep
     !% in {\tt thisHistory} and assuming that the corresponding fraction of the data value should be added to {\tt thisHistory}.
     use Galacticus_Error
     use Arrays_Search
     use Numerical_Ranges
     implicit none
     class(history)  , intent(inout) :: thisHistory
     type (history)  , intent(in   ) :: combineHistory
     integer                         :: combineHistoryPointCount,iPoint,jPoint ,timeBeginIndex ,timeEndIndex ,combineCount
     double precision                :: timeBegin,timeEnd,fractionContributed

     select type (thisHistory)
     type is (history)
     
        ! Return if combineHistory does not exist.
        if (.not.allocated(combineHistory%time)) return
     
        ! Get size of combineHistory.
        combineHistoryPointCount=size(combineHistory%time)
     
        ! Return if combineHistory has zero size.
        if (combineHistoryPointCount == 0) return
     
        ! If thisHistory does not exist, simply replace it with combineHistory.
        if (.not.allocated(thisHistory%time)) then
           call thisHistory%destroy()
           thisHistory=combineHistory
           return
        end if
     
        ! If thisHistory has zero size, simply replace it with combineHistory.
        if (size(thisHistory%time) == 0) then
           call thisHistory%destroy()
           thisHistory=combineHistory
           return
        end if
     
        ! The two objects must contain the same number of histories.
        if (size(thisHistory%data,dim=2) /= size(combineHistory%data,dim=2)) call Galacticus_Error_Report('History_Combine','two objects contain differing numbers of histories')
     
        ! Determine if we need to extend the time range in thisHistory.
        combineCount=size(combineHistory%time)
        if (thisHistory%rangeType == rangeTypeUndefined) then
           ! The history has no defined range type, so pass the time array of the history being combined to use as a template for new times.
           call thisHistory%extend(times=combineHistory%time)
        else
           ! The history has a defined range type, so simply pass the required extent of the range.
           call thisHistory%extend([combineHistory%time(1),combineHistory%time(combineCount)])
        end if
     
        ! Transfer each entry from combineHistory to thisHistory.
        do iPoint=2,combineCount
           ! Find indices in thisHistory spanned by combineHistory point.
           if (iPoint > 2) then
              ! Reuse the end index from the previous loop iteration if available.
              timeBeginIndex=timeEndIndex
           else
              timeBeginIndex=Search_Array(thisHistory%time,combineHistory%time(iPoint-1))
           end if
           timeEndIndex=min(Search_Array(thisHistory%time,combineHistory%time(iPoint))+1,size(thisHistory%time))    
           ! Loop over all points in thisHistory to which we need to add this contribution.
           do jPoint=timeBeginIndex,timeEndIndex
              if (jPoint == 1) then
                 timeBegin=                               combineHistory%time(iPoint-1)
              else
                 timeBegin=max(thisHistory%time(jPoint-1),combineHistory%time(iPoint-1))
              end if
              timeEnd     =min(thisHistory%time(jPoint  ),combineHistory%time(iPoint  ))
              fractionContributed=(timeEnd-timeBegin)/(combineHistory%time(iPoint)-combineHistory%time(iPoint-1))
              thisHistory%data(jPoint,:)=thisHistory%data(jPoint,:)+combineHistory%data(iPoint,:)*fractionContributed
           end do
        end do
     end select
     return
   end subroutine History_Combine

   function History_Divide(self,divisor)
     !% Divides history data by a double precision {\tt divisor}.
     implicit none
     type (history)               :: History_Divide
     class(history)  , intent(in) :: self
     double precision, intent(in) :: divisor
    
     select type(self)
     type is (history)
        History_Divide=self
        if (allocated(History_Divide%data)) History_Divide%data=History_Divide%data/divisor
     end select
     return
   end function History_Divide

   subroutine History_Extend(thisHistory,timeRange,times)
     !% Extends a history to encompass the given time range.
     use Numerical_Ranges
     use Galacticus_Error
     use ISO_Varying_String
     use String_Handling
     implicit none
     class(history),   intent(inout)                         :: thisHistory
     double precision, intent(in),  dimension(2  ), optional :: timeRange
     double precision, intent(in),  dimension(:  ), optional :: times
     double precision, allocatable, dimension(:  )           :: newTimes
     double precision, allocatable, dimension(:,:)           :: historyDataTemporary
     double precision,              dimension(2  )           :: timeRangeActual
     integer                                                 :: timeCount,timeBeginIndex,timeEndIndex,rangeType,historyCount&
          &,addCount,addCountStart,addCountEnd,newTimesAtStart,newTimesAtEnd
     logical                                                 :: useRange
     double precision                                        :: timeBegin,timeEnd,timeDelta
     type(varying_string)                                    :: message

     ! Determine the range of times that must be covered.
     if (present(timeRange)) then
        timeRangeActual=timeRange
     else
        if (present(times)) then
           timeRangeActual(1)=times(1          )
           timeRangeActual(2)=times(size(times))
        else
           call Galacticus_Error_Report('History_Extend','either timeRange or times must be specified')
        end if
     end if
  
     ! Determine if we need to extend the time range in thisHistory.
     timeCount     =size(thisHistory%time           )
     timeBegin     =     thisHistory%time(1        )
     timeEnd       =     thisHistory%time(timeCount)
     timeBeginIndex=1
     timeEndIndex  =timeCount
     if (.not.present(times)) then
        select case (thisHistory%rangeType)
        case (rangeTypeLinear     )
           timeDelta=(timeEnd-timeBegin)/dble(timeCount-1)
        case (rangeTypeLogarithmic)
           timeDelta=dlog(timeEnd/timeBegin)/dble(timeCount-1)
        case default
           if (thisHistory%rangeType == rangeTypeUndefined) then
              message='undefined range type: '//char(10)
           else
              message='unrecognized range type: '
              message=message//thisHistory%rangeType//char(10)
           end if
           message=message//' -> known types are: '//char(10)//' --> linear      : '//rangeTypeLinear//char(10)//' --> logarithmic : '//rangeTypeLogarithmic
           call Galacticus_Error_Report('History_Extend',message)
        end select
        if (timeRangeActual(1) < timeBegin) then
           select case (thisHistory%rangeType)
           case (rangeTypeLinear)
              addCountStart =int(    (timeBegin-timeRangeActual(1))/timeDelta)+1
              timeBegin=timeBegin      -dble(addCountStart)*timeDelta
           case (rangeTypeLogarithmic)
              addCountStart =int(dlog(timeBegin/timeRangeActual(1))/timeDelta)+1
              timeBegin=timeBegin*dexp(-dble(addCountStart)*timeDelta)
           end select
           timeBeginIndex=timeBeginIndex+addCountStart
           timeEndIndex  =timeEndIndex  +addCountStart
        else
           addCountStart=0
        end if
        if (timeRangeActual(2) > timeEnd  ) then
           select case (thisHistory%rangeType)
           case (rangeTypeLinear)
              addCountEnd =int(    (timeRangeActual(2)-timeEnd)/timeDelta)+1
              timeEnd  =timeEnd        +dble(addCountEnd)*timeDelta
           case (rangeTypeLogarithmic)
              addCountEnd =int(dlog(timeRangeActual(2)/timeEnd)/timeDelta)+1
              timeEnd  =timeEnd  *dexp(+dble(addCountEnd)*timeDelta)
           end select
        else
           addCountEnd=0
        end if
        addCount=addCountStart+addCountEnd
        useRange=.true.
     else
        newTimesAtStart=count(times < thisHistory%time(1                     ))
        newTimesAtEnd  =count(times > thisHistory%time(size(thisHistory%time)))
        timeBeginIndex=timeBeginIndex+newTimesAtStart
        timeEndIndex  =timeEndIndex  +newTimesAtStart
        addCount=newTimesAtStart+newTimesAtEnd
        allocate(newTimes(size(thisHistory%time)+addCount))
        if (newTimesAtStart > 0) newTimes(1             :             newTimesAtStart)=times(1                          :newTimesAtStart)
        if (newTimesAtEnd   > 0) newTimes(timeEndIndex+1:timeEndIndex+newTimesAtEnd  )=times(size(times)-newTimesAtEnd+1:size(times)    )
        newTimes(timeBeginIndex:timeEndIndex)=thisHistory%time
        useRange=.false.
     end if

     ! Create new arrays.
     if (addCount > 0) then
        ! Create copies of current histories.
        call Move_Alloc(thisHistory%data,historyDataTemporary)
        ! Store range type and number of histories.
        rangeType   =thisHistory%rangeType
        historyCount=size(historyDataTemporary,dim=2)
        ! Destroy the history and make a new one.
        call thisHistory%destroy()
        select case (useRange)
        case (.true. )
           call thisHistory%create(historyCount,timeCount+addCount,timeBegin,timeEnd,rangeType)
        case (.false.)
           call thisHistory%create(historyCount,timeCount+addCount)
           thisHistory%time=newTimes
           deallocate(newTimes)
        end select
        ! Copy data back to relevant location.
        thisHistory%data(timeBeginIndex:timeEndIndex,:)=historyDataTemporary
        deallocate(historyDataTemporary)
     end if
     return
   end subroutine History_Extend

  !# <galacticusStateStoreTask>
  !#  <unitName>Histories_State_Store</unitName>
  !# </galacticusStateStoreTask>
  subroutine Histories_State_Store(stateFile,fgslStateFile)
    !% Write the history state to file.
    use FGSL
    implicit none
    integer,         intent(in) :: stateFile
    type(fgsl_file), intent(in) :: fgslStateFile

    write (stateFile) historyStorageEarliestTime,historyStorageLatestTime
    return
  end subroutine Histories_State_Store
  
  !# <galacticusStateRetrieveTask>
  !#  <unitName>Histories_State_Retrieve</unitName>
  !# </galacticusStateRetrieveTask>
  subroutine Histories_State_Retrieve(stateFile,fgslStateFile)
    !% Retrieve the history state from the file.
    use FGSL
    implicit none
    integer,         intent(in) :: stateFile
    type(fgsl_file), intent(in) :: fgslStateFile

    read (stateFile) historyStorageEarliestTime,historyStorageLatestTime
    return
  end subroutine Histories_State_Retrieve

  subroutine History_Timesteps(thisHistory,timeSteps)
    !% Return an array of time intervals in {\tt thisHistory}.
    use Memory_Management
    use Numerical_Ranges
    implicit none
    class(history),   intent(in)                               :: thisHistory
    double precision, intent(inout), allocatable, dimension(:) :: timeSteps
    integer                                                    :: iTime
    double precision                                           :: ratio

    call Alloc_Array(timeSteps,shape(thisHistory%time),memoryType=memoryTypeNodes)
    select case (thisHistory%rangeType)
    case (rangeTypeLogarithmic)
       ratio=thisHistory%time(2)/thisHistory%time(1)
       forall(iTime=1:size(thisHistory%time))
          timeSteps(iTime)=thisHistory%time(1)*(ratio**(dble(iTime)-0.5d0)-ratio**(dble(iTime)-1.5d0))
       end forall
    case (rangeTypeLinear     )
       timeSteps=thisHistory%time(2)-thisHistory%time(1)
    case default
       timeSteps(1)=thisHistory%time(1)
       forall(iTime=2:size(thisHistory%time))
          timeSteps(iTime)=thisHistory%time(iTime)-thisHistory%time(iTime-1)
       end forall
    end select
    return
  end subroutine History_Timesteps
  
end module Histories
