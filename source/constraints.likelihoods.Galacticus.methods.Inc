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

!% Defines methods for a \glc\ likelihood function.

function likelihoodGalacticusConstructor(configFileName)
  !% Constructor for the \glc\ likelihood class.
  implicit none
  type(likelihoodGalacticus)                :: likelihoodGalacticusConstructor
  type(varying_string      ), intent(in   ) :: configFileName

  ! Store the config file name
  likelihoodGalacticusConstructor%configFileName=configFileName
  return
end function likelihoodGalacticusConstructor

double precision function likelihoodGalacticusEvaluate(self,simulationState,simulationConvergence,temperature,logLikelihoodCurrent,logPriorCurrent,logPriorProposed)
  !% Return the log-likelihood for the \glc\ likelihood function. This function runs an external script that drives \glc\ and
  !% writes the resulting likelihood to file. This function then reads that likelihood from file. We make use of {\tt tmpfs} for
  !% this likelihood file so that no disk I/O is required.
  use Constraints_State
  use Constraints_Convergence
  use System_Command
  use String_Handling
  use MPI_Utilities
  implicit none
  class           (likelihoodGalacticus), intent(inout)               :: self
  class           (state               ), intent(in   )               :: simulationState
  class           (convergence         ), intent(inout)               :: simulationConvergence
  double precision                      , intent(in   )               :: temperature    , logLikelihoodCurrent, &
       &                                                                 logPriorCurrent, logPriorProposed
  double precision                      , allocatable  , dimension(:) :: stateVector
  type            (varying_string      )                              :: wrapperCommand , likelihoodFileName
  character       (len=23              )                              :: label
  integer                                                             :: i              , likelihoodFileUnit

  ! Generate the file name for the likelihood file.
  likelihoodFileName='/dev/shm/galacticusLikelihood_'//mpiSelf%rankLabel()//'.dat'
  ! Generate temperature label.
  write (label,'(e12.6)') temperature 
  ! Generate the command to run the Galacticus model.
  wrapperCommand='./constraints/constrainGalacticusWrapper.pl '//self%configFileName//' '//mpiSelf%rankLabel()//' '//likelihoodFileName//' '//trim(label)
  ! Append the current state.
  stateVector=simulationState%get()
  do i=1,size(stateVector)
     write (label,'(1x,e22.15)') stateVector(i)
     wrapperCommand=wrapperCommand//label
  end do
  deallocate(stateVector)
  ! Run the Galacticus model and evaluate likelihood.
  call System_Command_Do(wrapperCommand)
  ! Read the likelihood back from file.
  open(newUnit=likelihoodFileUnit,file=char(likelihoodFileName),status='old',form='formatted')
  read (likelihoodFileUnit,*) likelihoodGalacticusEvaluate
  close(likelihoodFileUnit,status='delete')
  return
end function likelihoodGalacticusEvaluate

subroutine likelihoodGalacticusFunctionChanged(self)
  !% Respond to possible changes in the likelihood function.
  implicit none
  class(likelihoodGalacticus), intent(inout) :: self
  
  return
end subroutine likelihoodGalacticusFunctionChanged
