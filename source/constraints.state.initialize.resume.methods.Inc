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

!% Defines methods for a resuming state initializor.

function stateInitializorResumeConstructor(logFileRoot,restoreState)
  !% Constructor for the resume state initializor class.
  implicit none
  type     (stateInitializorResume)                :: stateInitializorResumeConstructor
  character(len=*                 ), intent(in   ) :: logFileRoot
  logical                          , intent(in   ) :: restoreState

  stateInitializorResumeConstructor%logFileRoot =logFileRoot
  stateInitializorResumeConstructor%restoreState=restoreState
  return
end function stateInitializorResumeConstructor

subroutine stateInitializorResumeInitialize(self,simulationState,parameterPriors,parameterMappings,timeEvaluatePrevious)
  !% Initialize simulation state by resuming from a previous simulation.
  use Constraints_Priors
  use Constraints_State
  use Constraints_Constants
  use MPI_Utilities
  use String_Handling
  use Galacticus_Display
  implicit none
  class           (stateInitializorResume), intent(inout)               :: self
  class           (state                 ), intent(inout)               :: simulationState
  type            (prior                 ), intent(inout), dimension(:) :: parameterPriors
  type            (mappingList           ), intent(in   ), dimension(:) :: parameterMappings
  double precision                        , intent(  out)               :: timeEvaluatePrevious
  double precision                        , allocatable  , dimension(:) :: stateVector
  type            (varying_string        )                              :: logFileName    , message
  integer                                                               :: stateCount     , mpiRank     , &
       &                                                                   logFileUnit    , ioStatus    , &
       &                                                                   i
  double precision                                                      :: logPosterior
  logical                                                               :: converged      , first
  character       (len=12                )                              :: label
  
  ! Allocate the state vector.
  allocate(stateVector(simulationState%dimension()))
  ! Read state from the log file.
  logFileName=self%logFileRoot//'_'//mpiSelf%rankLabel()//'.log'
  open(newunit=logFileUnit,file=char(logFileName),status='unknown',form='formatted')
  ioStatus=0
  first   =.true.
  do while (ioStatus == 0)
     read (logFileUnit,*,iostat=ioStatus) stateCount          , &
          &                               mpiRank             , &
          &                               timeEvaluatePrevious, &
          &                               converged           , &
          &                               logPosterior        , &
          &                               stateVector     
     if (ioStatus == 0) then
        ! Map the state.
        do i=1,size(stateVector)
           stateVector(i)=parameterMappings(i)%thisMapping%map(stateVector(i))
        end do
        ! Restore the state object.
        if (self%restoreState) call simulationState%restore(stateVector,first)
        first=.false.
     end if
  end do
  close(logFileUnit)
  ! Set the simulation state.
  call simulationState%update(stateVector,.false.,.false.)
  ! Check for out of range state.
  do i=1,simulationState%dimension()
     if (parameterPriors(i)%logDensity(simulationState) <= logImpossible) then
        write (label,'(e12.6)') parameterMappings(i)%thisMapping%unmap(stateVector(i))
        message='Out of range state for parameter '
        message=message//i//' ['//trim(label)//']'
        call Galacticus_Display_Message(message)
     end if
  end do
  ! Deallocate the state vector.
  deallocate(stateVector)
  return
end subroutine stateInitializorResumeInitialize
