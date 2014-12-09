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

function stateInitializorResumeConstructor(logFileRoot)
  !% Constructor for the resume state initializor class.
  implicit none
  type     (stateInitializorResume)                :: stateInitializorResumeConstructor
  character(len=*                 ), intent(in   ) :: logFileRoot

  stateInitializorResumeConstructor%logFileRoot=logFileRoot
  return
end function stateInitializorResumeConstructor

subroutine stateInitializorResumeInitialize(self,simulationState,parameterPriors,parameterMappings)
  !% Initialize simulation state by resuming from a previous simulation.
  use Constraints_Priors
  use Constraints_State
  use MPI_Utilities
  implicit none
  class           (stateInitializorResume), intent(inout)               :: self
  class           (state                 ), intent(inout)               :: simulationState
  type            (prior                 ), intent(inout), dimension(:) :: parameterPriors
  type            (mappingList           ), intent(in   ), dimension(:) :: parameterMappings
  double precision                        , allocatable  , dimension(:) :: stateVector
  type            (varying_string        )                              :: logFileName
  integer                                                               :: stateCount     , mpiRank     , &
       &                                                                   logFileUnit    , ioStatus    , &
       &                                                                   i
  double precision                                                      :: time           , logPosterior

  ! Restore the state object.
  call simulationState%restore(char(self%logFileRoot))
  ! Allocate the state vector.
  allocate(stateVector(simulationState%dimension()))
  ! Read state from the log file.
  logFileName=self%logFileRoot//'_'//mpiSelf%rankLabel()//'.log'
  open(newunit=logFileUnit,file=char(logFileName),status='unknown',form='formatted')
  ioStatus=0
  do while (ioStatus == 0)
     read (logFileUnit,*,iostat=ioStatus) stateCount  , &
          &               mpiRank     , &
          &               time        , &
          &               logPosterior, &
          &               stateVector     
  end do
  close(logFileUnit)
  ! Map the state.
  do i=1,size(stateVector)
     stateVector(i)=parameterMappings(i)%thisMapping%map(stateVector(i))
  end do
  ! Set the simulation state.
  call simulationState%update(stateVector,.false.,.false.)
  ! Deallocate the state vector.
  deallocate(stateVector)
  return
end subroutine stateInitializorResumeInitialize
