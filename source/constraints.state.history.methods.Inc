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

!% Defines methods for the ``history'' state class.

function stateHistoryConstructor(parameterCount,acceptedStateCount)
  !% Constructor for ``history'' state class.
  implicit none
  type   (stateHistory)               :: stateHistoryConstructor
  integer             , intent(in   ) :: parameterCount         , acceptedStateCount

  allocate(stateHistoryConstructor%current        (    parameterCount))
  allocate(stateHistoryConstructor%accepted       (acceptedStateCount))
  allocate(stateHistoryConstructor%stateSum       (    parameterCount))
  allocate(stateHistoryConstructor%stateSquaredSum(    parameterCount))
  stateHistoryConstructor%parameterCount =parameterCount
  stateHistoryConstructor%stepCount      =0
  stateHistoryConstructor%stateSum       =0.0d0
  stateHistoryConstructor%stateSquaredSum=0.0d0
  stateHistoryConstructor%accepted       =0
  return
end function stateHistoryConstructor

subroutine stateHistoryDestroy(self)
  !% Destructor for ``history'' state class.
  implicit none
  type(stateHistory), intent(inout) :: self

  if (allocated(self%current        )) deallocate(self%current        )
  if (allocated(self%stateSum       )) deallocate(self%stateSum       )
  if (allocated(self%stateSquaredSum)) deallocate(self%stateSquaredSum)
  return
end subroutine stateHistoryDestroy

subroutine stateHistoryUpdate(self,stateNew,logState,isConverged,outlierMask)
  !% Update the current state.
  implicit none
  class           (stateHistory), intent(inout)                         :: self
  double precision              , intent(in   ), dimension(:)           :: stateNew
  logical                       , intent(in   )                         :: logState
  logical                       , intent(in   )                         :: isConverged
  logical                       , intent(in   ), dimension(:), optional :: outlierMask

  call self%stateSimple%update(stateNew,logState,isConverged,outlierMask)
  if (logState) then
     self%stateSum       =self%stateSum       +stateNew
     self%stateSquaredSum=self%stateSquaredSum+stateNew**2
  end if
  return
end subroutine stateHistoryUpdate

function stateHistoryMean(self)
  !% Return the mean over state history.
  use Galacticus_Error
  implicit none
  class           (stateHistory), intent(in   )                  :: self
  double precision              , dimension(self%parameterCount) :: stateHistoryMean
  
  stateHistoryMean=self%stateSum/dble(self%stepCount)
  return
end function stateHistoryMean

function stateHistoryVariance(self)
  !% Return the mean over state history.
  use Galacticus_Error
  implicit none
  class           (stateHistory), intent(in   )                  :: self
  double precision              , dimension(self%parameterCount) :: stateHistoryVariance

  stateHistoryVariance=(self%stateSquaredSum/dble(self%stepCount)-self%mean()**2)*dble(self%stepCount)/dble(self%stepCount-1)
  return
end function stateHistoryVariance

subroutine stateHistoryReset(self)
  !% Reset the state object.
  implicit none
  class(stateHistory), intent(inout) :: self
  
  call self%stateSimple%reset()
  self%stateSum       =0.0d0
  self%stateSquaredSum=0.0d0
  return
end subroutine stateHistoryReset

subroutine stateHistoryRestore(self,stateFileRoot)
  !% Restore the state object from file.
  use MPI_Utilities
  use ISO_Varying_String
  implicit none
  class           (stateHistory  ), intent(inout)               :: self
  character       (len=*         ), intent(in   )               :: stateFileRoot
  double precision                , allocatable  , dimension(:) :: stateVector
  integer                                                       :: mpiRank      , stateCount  , &
       &                                                           logFileUnit  , ioStatus
  double precision                                              :: time         , logPosterior
  type            (varying_string)                              :: logFileName

  ! Reset.
  call self%reset()
  ! Perform restore of parent class.
  call self%stateSimple%restore(stateFileRoot)
  ! Allocate the state vector.
  allocate(stateVector(self%parameterCount))
  ! Read state from the log file.
  logFileName=stateFileRoot//'_'//mpiSelf%rankLabel()//'.log'
  open(newunit=logFileUnit,file=char(logFileName),status='unknown',form='formatted')
  ioStatus=0
  do while (ioStatus == 0)
     read (logFileUnit,*,iostat=ioStatus) stateCount  , &
          &               mpiRank     , &
          &               time        , &
          &               logPosterior, &
          &               stateVector
     ! Accumulate state and state-squared.
     self%stateSum       =self%stateSum       +stateVector
     self%stateSquaredSum=self%stateSquaredSum+stateVector**2
  end do
  close(logFileUnit)
  ! Deallocate the state vector.
  deallocate(stateVector)
  return
end subroutine stateHistoryRestore
