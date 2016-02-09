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

!% Defines methods for a \glc\ likelihood function.

function likelihoodGalacticusConstructor(configFileName,delayInterval)
  !% Constructor for the \glc\ likelihood class.
  use File_Utilities
  use Galacticus_Error
  implicit none
  type            (likelihoodGalacticus)                :: likelihoodGalacticusConstructor
  type            (varying_string      ), intent(in   ) :: configFileName
  double precision                      , intent(in   ) :: delayInterval

  ! Store the config file name
  likelihoodGalacticusConstructor%configFileName    =configFileName
  likelihoodGalacticusConstructor%delayInterval     =delayInterval
  likelihoodGalacticusConstructor%storeCountPrevious=0
  ! Find the time command.
  likelihoodGalacticusConstructor%timeCommand       =Executable_Find('time')
  if (likelihoodGalacticusConstructor%timeCommand == "") call Galacticus_Error_Report('likelihoodGalacticusConstructor','a working GNU time command is required')
  return
end function likelihoodGalacticusConstructor

double precision function likelihoodGalacticusEvaluate(self,simulationState,parameterMappings,simulationConvergence,temperature,logLikelihoodCurrent,logPriorCurrent,logPriorProposed,timeEvaluate,logLikelihoodVariance)
  !% Return the log-likelihood for the \glc\ likelihood function. This function runs an external script that drives \glc\ and
  !% writes the resulting likelihood to file. This function then reads that likelihood from file. We make use of {\normalfont
  !% \ttfamily tmpfs} for this likelihood file so that no disk I/O is required.
  use Constraints_Constants
  use Constraints_State
  use Constraints_Convergence
  use System_Command
  use String_Handling
  use MPI_Utilities
  use Galacticus_Error
  use File_Utilities
  implicit none
  class           (likelihoodGalacticus), intent(inout)               :: self
  class           (state               ), intent(in   )               :: simulationState
  type            (mappingList         ), intent(in   ), dimension(:) :: parameterMappings
  class           (convergence         ), intent(inout)               :: simulationConvergence
  double precision                      , intent(in   )               :: temperature          , logLikelihoodCurrent, &
       &                                                                 logPriorCurrent      , logPriorProposed
  real                                  , intent(inout)               :: timeEvaluate
  double precision                      , intent(  out), optional     :: logLikelihoodVariance
  double precision                      , allocatable  , dimension(:) :: stateVector
  type            (varying_string      )                              :: wrapperCommand       , likelihoodFileName  , &
       &                                                                 timingFileName       , message
  character       (len=23              )                              :: label
  integer                                                             :: i                    , likelihoodFileUnit  , &
       &                                                                 timingFileUnit       , status
  logical                                                             :: storeModel
  real                                                                :: timeSystem           , timeUser

  ! If prior probability is impossible, then no need to waste time evaluating the likelihood.
  if (logPriorProposed <= logImpossible) then
     likelihoodGalacticusEvaluate=logImpossible
     if (present(logLikelihoodVariance)) logLikelihoodVariance=0.0d0
     return
  end if
  ! Generate the file name for the likelihood and timing files.
  likelihoodFileName='/dev/shm/galacticusLikelihood_'//mpiSelf%rankLabel()//'.dat'
  timingFileName    ='/dev/shm/glcTiming_'           //mpiSelf%rankLabel()//'.txt'
  ! Generate temperature label.
  write (label,'(e12.6)') temperature 
  ! Generate the command to run the Galacticus model.
  wrapperCommand=                                                      &
       & self%timeCommand                                  //    ' '// &
       & '--format "%S %U"'                                //    ' '// &
       & '--output '//timingFileName                       //    ' '// &
       & './constraints/constrainGalacticus.pl '                    // &
       & trim(self%configFileName)                         //    ' '// &
       & mpiSelf%rankLabel()                               //    ' '// &
       & likelihoodFileName                                //    ' '// &
       & trim(label              )
  ! Determine whether to store this model.
  select type (simulationState)
  class is (stateCorrelation)
     storeModel=                                                                                &
          &       simulationConvergence%isConverged              (                            ) &
          & .and.                                                                               &
          &  .not.simulationConvergence%chainIsOutlier           (simulationState%chainIndex()) &
          & .and.                                                                               &
          &       simulationState%count                          (                            ) &
          &      >=                                                                             &
          &       self           %storeCountPrevious                                            &
          &      +                                                                              &
          &       simulationState%correlationLength              (                            )
  class default
     storeModel=.false.
  end select
  if (storeModel) self%storeCountPrevious=simulationState%count()
  if (storeModel) then
     wrapperCommand=wrapperCommand//' '//simulationState%count()
  else
     wrapperCommand=wrapperCommand//' none'
  end if
  ! Append the current state.
  stateVector=simulationState%get()
  do i=1,size(stateVector)
     stateVector(i)=parameterMappings(i)%thisMapping%unmap(stateVector(i))
     write (label,'(1x,e22.15)') stateVector(i)
     wrapperCommand=wrapperCommand//label
  end do
  deallocate(stateVector)
  ! If requested, sleep for a short time to desynchronize launching of the Galacticus wrapper command. On some systems this is
  ! necessary to avoid resource conflicts.
  if (self%delayInterval > 0.0d0) call Sleep(int(self%delayInterval*dble(mpiSelf%rank())))
  ! Run the Galacticus model and evaluate likelihood.
  call System_Command_Do(wrapperCommand,status)
  if (status /= 0) then     
     message="Galacticus wrapper failed with status code "
     message=message//status
     if (File_Exists(likelihoodFileName)) then
        message=message//" (likelihood file '"//char(likelihoodFileName)//"' does exist)"
     else
        message=message//" (likelihood file '"//char(likelihoodFileName)//"' DOES NOT exist)"
     end if
     call Galacticus_Error_Report('likelihoodGalacticusEvaluate',message)
  end if
  ! Read the likelihood back from file.
  open(newUnit=likelihoodFileUnit,file=char(likelihoodFileName),status='old',form='formatted')
  read                                     (likelihoodFileUnit,*) likelihoodGalacticusEvaluate
  if (present(logLikelihoodVariance)) read (likelihoodFileUnit,*) logLikelihoodVariance
  close(likelihoodFileUnit,status='delete')
  ! Read the timing data back from file.
  open(newUnit=timingFileUnit,file=char(timingFileName),status='old',form='formatted')
  read (timingFileUnit,*) timeSystem,timeUser
  close(timingFileUnit,status='delete')
  timeEvaluate=timeSystem+timeUser
  return
end function likelihoodGalacticusEvaluate

subroutine likelihoodGalacticusFunctionChanged(self)
  !% Respond to possible changes in the likelihood function.
  implicit none
  class(likelihoodGalacticus), intent(inout) :: self
  
  return
end subroutine likelihoodGalacticusFunctionChanged

logical function likelihoodGalacticusWillEvaluate(self,simulationState,parameterMappings,simulationConvergence,temperature,logLikelihoodCurrent,logPriorCurrent,logPriorProposed)
  !% Return true if the log-likelihood will be evaluated.
  use Constraints_State
  use Constraints_Convergence
  use Constraints_Constants
  implicit none
  class           (likelihoodGalacticus), intent(inout)               :: self
  class           (state               ), intent(in   )               :: simulationState
  type            (mappingList         ), intent(in   ), dimension(:) :: parameterMappings
  class           (convergence         ), intent(inout)               :: simulationConvergence
  double precision                      , intent(in   )               :: temperature                      , logLikelihoodCurrent, &
       &                                                                 logPriorCurrent                  , logPriorProposed

  ! Likelihood will not be evaluated if the proposed prior is impossible.
  likelihoodGalacticusWillEvaluate=(logPriorProposed > logImpossible)
  return
end function likelihoodGalacticusWillEvaluate
