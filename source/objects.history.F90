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


!% Contains a module defining the history object type.

module Histories
  !% Defines the history object type.
  private
  public :: history, History_Set_Times

  type history
     !% The history object type.
     double precision, allocatable, dimension(:  ) :: time
     double precision, allocatable, dimension(:,:) :: data
     double precision, allocatable, dimension(:,:) :: rates

   contains
     
     procedure :: create  => History_Create
     procedure :: destroy => History_Destroy
     procedure :: trim    => History_Trim
     procedure :: add     => History_Add
     procedure :: reset   => History_Reset

  end type history

  ! A null history object.
  type(history),    public :: nullHistory

  ! Earliest and latest times for history storage.
  double precision, public :: historyStorageEarliestTime= 0.1d0
  double precision, public :: historyStorageLatestTime  =15.0d0

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
    type(history),    intent(inout)        :: thisHistory
    integer,          intent(in)           :: historyCount,timesCount
    double precision, intent(in), optional :: timeBegin,timeEnd
    integer,          intent(in), optional :: rangeType
    integer                                :: rangeTypeActual

    if (allocated(thisHistory%time)) then
       call Galacticus_Error_Report('History_Create','this history appears to have been created already')
    else
       allocate(thisHistory%time (timesCount             ))
       allocate(thisHistory%data (timesCount,historyCount))
       allocate(thisHistory%rates(timesCount,historyCount))
       call Memory_Usage_Record(timesCount*(1+2*historyCount),addRemove=+8)
       if (present(timeBegin)) then
          if (.not.present(timeEnd)) call Galacticus_Error_Report('History_Create','an end time must be given if a begin time is given')
          if (present(rangeType)) then
             rangeTypeActual=rangeType
          else
             rangeTypeActual=rangeTypeLogarithmic
          end if
          thisHistory%time=Make_Range(timeBegin,timeEnd,timesCount,rangeTypeActual)
       else
          thisHistory%time=0.0d0
       end if
       thisHistory%data =0.0d0
       thisHistory%rates=0.0d0
    end if

    return
  end subroutine History_Create

  subroutine History_Destroy(thisHistory)
    !% Destroy a history.
    use Memory_Management
    implicit none
    type(history), intent(inout) :: thisHistory
    integer                      :: timesCount,historyCount
    
    if (allocated(thisHistory%time)) then
       timesCount  =size(thisHistory%time      )
       historyCount=size(thisHistory%data,dim=2)
       call Memory_Usage_Record(timesCount*(1+2*historyCount),addRemove=-8)
       deallocate(thisHistory%time )
       deallocate(thisHistory%data )
       deallocate(thisHistory%rates)
    end if
    return
  end subroutine History_Destroy

  subroutine History_Reset(thisHistory)
    !% Reset a history by zeroing all elements, but leaving the structure (and times) intact.
    use Memory_Management
    implicit none
    type(history), intent(inout) :: thisHistory
    
    if (allocated(thisHistory%time)) then
       thisHistory%data =0.0d0
       thisHistory%rates=0.0d0
    end if
    return
  end subroutine History_Reset

  subroutine History_Trim(thisHistory,currentTime,minimumPointsToRemove)
    !% Removes outdated information from ``future histories'' (i.e. histories that store data for future reference). Removes all
    !% but one entry prior to the given {\tt currentTime} (this allows for interpolation of the history to the current
    !% time). Optionally, the remove is done only if it will remove more than {\tt minimumPointsToRemove} entries (since the
    !% removal can be slow this allows for some optimization).
    use Galacticus_Error
    use Memory_Management
    implicit none
    type(history),    intent(inout)        :: thisHistory
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
       call Move_Alloc(thisHistory%time ,temporaryHistory%time )
       call Move_Alloc(thisHistory%data ,temporaryHistory%data )
       call Move_Alloc(thisHistory%rates,temporaryHistory%rates)
       ! Reallocate the history arrays.
       newPointCount=currentPointCount-iTrim
       historyCount =size(temporaryHistory%data,dim=2)
       allocate(thisHistory%time (newPointCount             ))
       allocate(thisHistory%data (newPointCount,historyCount))
       allocate(thisHistory%rates(newPointCount,historyCount))
       ! Copy the data back into the new arrays.
       thisHistory%time (:  )=temporaryHistory%time (iTrim+1:currentPointCount  )
       thisHistory%data (:,:)=temporaryHistory%data (iTrim+1:currentPointCount,:)
       thisHistory%rates(:,:)=temporaryHistory%rates(iTrim+1:currentPointCount,:)
       ! Deallocate the temporary arrays.
       deallocate(temporaryHistory%time )
       deallocate(temporaryHistory%data )
       deallocate(temporaryHistory%rates)
       ! Account for change in memory usage.
       call Memory_Usage_Record(iTrim*(1+2*historyCount),addRemove=-8)
    end if
    return
  end subroutine History_Trim

  subroutine History_Add(thisHistory,addHistory)
    !% Adds the data in {\tt addHistory} to that in {\tt thisHistory}.
    use FGSL
    use Numerical_Interpolation
    use Galacticus_Error
    implicit none
    type(history),           intent(inout) :: thisHistory
    type(history),           intent(in)    :: addHistory
    integer                                :: addHistoryPointCount,iPoint,interpolationPoint,iHistory
    double precision                       :: interpolationFactors(2)
    type(fgsl_interp_accel)                :: interpolationAccelerator
    logical                                :: interpolationReset
    
    ! Return if addHistory does not exist.
    if (.not.allocated(addHistory%time)) return
    
    ! Get size of addHistory.
    addHistoryPointCount=size(addHistory%time)
    
    ! Return if addHistory has zero size.
    if (addHistoryPointCount == 0) return

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
             thisHistory%data(iPoint,iHistory)=thisHistory%data(iPoint,iHistory)+addHistory%data(interpolationPoint,iHistory)&
                  &*interpolationFactors(1)+addHistory%data(interpolationPoint+1,iHistory)*interpolationFactors(2)
          end forall

       end if

    end do

    return
  end subroutine History_Add

end module Histories
