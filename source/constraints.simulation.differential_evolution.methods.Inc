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

!% Defines methods for the ``differentialEvolution'' simulator class.

function simulatorDifferentialEvolutionConstructor(parameterPriors,randomDistributions,parameterMappings,modelLikelihood,simulationConvergence,simulationStoppingCriterion&
     &,simulationState,simulationStateInitializor,proposalSize,randomJump,stepsMaximum,acceptanceAverageCount,stateSwapCount,logFileRoot,sampleOutliers)
  !% Constructor for ``differentialEvolution'' simulator class.
  implicit none
  type     (simulatorDifferentialEvolution)                                      :: simulatorDifferentialEvolutionConstructor
  type     (prior                         ), intent(in   ), target, dimension(:) :: parameterPriors
  type     (distributionList              ), intent(in   ), target, dimension(:) :: randomDistributions
  type     (mappingList                   ), intent(in   ), target, dimension(:) :: parameterMappings
  class    (likelihood                    ), intent(in   ), target               :: modelLikelihood
  class    (convergence                   ), intent(in   ), target               :: simulationConvergence
  class    (stoppingCriterion             ), intent(in   ), target               :: simulationStoppingCriterion
  class    (state                         ), intent(in   ), target               :: simulationState
  class    (stateInitializor              ), intent(in   ), target               :: simulationStateInitializor
  class    (deProposalSize                ), intent(in   ), target               :: proposalSize
  class    (deRandomJump                  ), intent(in   ), target               :: randomJump
  integer                                  , intent(in   )                       :: stepsMaximum,acceptanceAverageCount,stateSwapCount
  character(len=*                         ), intent(in   )                       :: logFileRoot
  logical                                  , intent(in   )                       :: sampleOutliers

  simulatorDifferentialEvolutionConstructor%parameterCount              =  size(parameterPriors)
  simulatorDifferentialEvolutionConstructor%parameterPriors             => parameterPriors
  simulatorDifferentialEvolutionConstructor%randomDistributions         => randomDistributions
  simulatorDifferentialEvolutionConstructor%parameterMappings           => parameterMappings
  simulatorDifferentialEvolutionConstructor%modelLikelihood             => modelLikelihood
  simulatorDifferentialEvolutionConstructor%simulationConvergence       => simulationConvergence
  simulatorDifferentialEvolutionConstructor%simulationStoppingCriterion => simulationStoppingCriterion
  simulatorDifferentialEvolutionConstructor%simulationState             => simulationState
  simulatorDifferentialEvolutionConstructor%simulationStateInitializor  => simulationStateInitializor
  simulatorDifferentialEvolutionConstructor%proposalSize                => proposalSize
  simulatorDifferentialEvolutionConstructor%randomJump                  => randomJump
  simulatorDifferentialEvolutionConstructor%stepsMaximum                =  stepsMaximum
  simulatorDifferentialEvolutionConstructor%acceptanceAverageCount      =  acceptanceAverageCount
  simulatorDifferentialEvolutionConstructor%stateSwapCount              =  stateSwapCount
  simulatorDifferentialEvolutionConstructor%logFileRoot                 =  logFileRoot
  simulatorDifferentialEvolutionConstructor%sampleOutliers              =  sampleOutliers
  return
end function simulatorDifferentialEvolutionConstructor

subroutine simulatorDifferentialEvolutionSimulate(self)
  !% Perform a differential evolution simulation.
  use MPI_Utilities
  use Pseudo_Random
  use Galacticus_Error
  use Galacticus_Display
  use String_Handling
  use Constraints_Constants
  implicit none
  class           (simulatorDifferentialEvolution), intent(inout)                    :: self
  integer                                         , dimension(                    2) :: chainPair
  double precision                                , dimension(self%parameterCount,2) :: statePair
  double precision                                , dimension(self%parameterCount  ) :: stateVector          , stateVectorProposed
  class           (state                         ), allocatable                      :: stateProposed
  real                                                                               :: timePreEvaluate      , timePostEvaluate     , &
       &                                                                                timeEvaluate
  double precision                                                                   :: logLikelihoodProposed, logPosteriorProposed , &
                                                                                        logLikelihood
  type            (pseudoRandom                  )                                   :: randomNumberGenerator
  type            (varying_string                )                                   :: logFileName          , message
  integer                                                                            :: logFileUnit          , convergedAtStep      , &
       &                                                                                convergenceFileUnit  , i
  
  ! Check that we have sufficient chains for differential evolution.
  if (mpiSelf%count() < 3) call Galacticus_Error_Report('simulatorDifferentialEvolutionSimulate','at least 3 chains are required for differential evolution')
  ! Allocate a simple state object for the proposed state.
  allocate(stateSimple :: stateProposed)
  select type (stateProposed)
  type is (stateSimple)
     stateProposed=stateSimple(self%parameterCount,1)
  end select
  ! Initialize chain to some state vector.
  call self%simulationStateInitializor%initialize(self%simulationState,self%parameterPriors,self%parametermappings)
  ! Evaluate the posterior in the initial state.
  call self%posterior(self%simulationState,logImpossible,logImpossible,self%logPosterior,logLikelihood)
  ! Begin stepping.
  logFileName=self%logFileRoot//'_'//mpiSelf%rankLabel()//'.log'
  open(newunit=logFileUnit,file=char(logFileName),status='unknown',form='formatted')
  self%isConverged=.false. 
  do while (                                                                                                                  &
       &          self%simulationState            %count(                                               ) < self%stepsMaximum &
       &    .and.                                                                                                             &
       &     .not.self%simulationStoppingCriterion%stop (self%simulationState,self%simulationConvergence)                     &
       &   )
     ! Pick two random processes to use for proposal.
     chainPair(1)=simulatorDifferentialChainSelect(self,randomNumberGenerator               )
     chainPair(2)=simulatorDifferentialChainSelect(self,randomNumberGenerator,[chainPair(1)])
     ! Receive states from selected chains.
     stateVector=self%simulationState%get()
     statePair  =mpiSelf%requestData(chainPair,stateVector)
     ! Generate proposal.
     stateVectorProposed= stateVector      &
          &              +self%stepSize()  &
          &              *(                &
          &                +statePair(:,1) &
          &                -statePair(:,2) &
          &               )
     ! Add random perturbations to the proposal.
     stateVectorProposed=stateVectorProposed+self%randomJump%sample(self%randomDistributions,self%simulationState)
     ! Store the proposed state vector.
     call stateProposed%update(stateVectorProposed,.false.,.false.)
     ! Evaluate likelihood.
     call CPU_Time(timePreEvaluate )
     call self%posterior(stateProposed,logLikelihood,self%logPosterior-logLikelihood,logPosteriorProposed,logLikelihoodProposed)
     call CPU_Time(timePostEvaluate)
     timeEvaluate=timePostEvaluate-timePreEvaluate
     ! Decide whether to take step.
     if (self%acceptProposal(self%logPosterior,logPosteriorProposed,randomNumberGenerator)) then
        self%logPosterior =logPosteriorProposed
        logLikelihood     =logLikelihoodProposed
        stateVector       =stateVectorProposed
     end if
     call self%update(stateVector)
     ! Unmap parameters and write to log file.
     do i=1,size(stateVector)
        stateVector(i)=self%parameterMappings(i)%thisMapping%unmap(stateVector(i))
     end do
     if (self%logging()) write (logFileUnit,*) self   %simulationState%count(), &
          &                                    mpiSelf%rank                 (), &
          &                                    timeEvaluate                   , &
          &                                    self%isConverged               , &
          &                                    self%logPosterior              , &
          &                                    stateVector
     ! Repeat.
     call mpiBarrier()
     ! Test for convergence.
     if (.not.self%isConverged) then
        self%isConverged=self%simulationConvergence%isConverged(self%simulationState,self%logPosterior)
        if (self%isConverged) then
           convergedAtStep=self%simulationState%count()
           if (mpiSelf%rank() == 0) then
              message='Converged after '
              message=message//convergedAtStep//' steps'
              call Galacticus_Display_Message(message)
              logFileName=self%logFileRoot//'_'//mpiSelf%rankLabel()//'.convergence.log'
              open(newunit=convergenceFileUnit,file=char(logFileName),status='unknown',form='formatted',access='append')
              write (convergenceFileUnit,'(a,i8)') 'Converged at step: ',convergedAtStep
              call self%simulationConvergence%logReport(convergenceFileUnit)
              close(convergenceFileUnit)
           end if
        end if
     end if
  end do
  close(logFileUnit)
  return
end subroutine simulatorDifferentialEvolutionSimulate

subroutine simulatorDifferentialEvolutionUpdate(self,stateVector)
  !% Update the differential evolution simulator state.
  use MPI_Utilities
  implicit none
  class           (simulatorDifferentialEvolution), intent(inout)                                 :: self
  double precision                                , intent(in   ), dimension(self%parameterCount) :: stateVector
  logical                                         , allocatable  , dimension(:                  ) :: outlierMask
  integer                                                                                         :: i

  allocate(outlierMask(mpiSelf%count()))
  do i=0,mpiSelf%count()-1
     outlierMask(i)=self%simulationConvergence%chainIsOutlier(i)
  end do
  call self%simulationState%update(stateVector,self%logging(),self%simulationConvergence%isConverged(),outlierMask)
  return
end subroutine simulatorDifferentialEvolutionUpdate

integer function simulatorDifferentialChainSelect(self,randomNumberGenerator,blockedChains)
  !% Select a chain at random, optionally excluding blocked chains.
  use Pseudo_Random
  use MPI_Utilities
  implicit none
  class  (simulatorDifferentialEvolution), intent(inout)                         :: self
  type   (pseudoRandom                  ), intent(inout)                         :: randomNumberGenerator
  integer                                , intent(in   ), dimension(:), optional :: blockedChains
  logical                                                                        :: accept
  
  accept=.false.
  do while (.not.accept)
     simulatorDifferentialChainSelect=                                        &
          &      min(                                                         &
          &          int(                                                     &
          &               dble(mpiSelf%count())                               &
          &              *randomNumberGenerator%sample(mpiRankOffset=.true.)  &
          &             )                                                   , &
          &                    mpiSelf%count()                                &
          &          -1                                                       &
          &         )
     accept=.true.
     if (.not.self%sampleOutliers.and.self%isConverged)                                              &
          & accept=                                                                                  &
          &              self%simulationConvergence%chainIsOutlier(                  mpiSelf%rank()) &
          &        .or.                                                                              &
          &         .not.self%simulationConvergence%chainIsOutlier(simulatorDifferentialChainSelect)
     if (present(blockedChains)) accept=accept.and.all(simulatorDifferentialChainSelect /= blockedChains)
  end do
  return
end function simulatorDifferentialChainSelect

logical function simulatorDifferentialEvolutionLogging(self)
  !% Specifies whether or not the current state should be logged to file during differential evolution.
  implicit none
  class(simulatorDifferentialEvolution), intent(inout) :: self

  simulatorDifferentialEvolutionLogging=.true.
  return
end function simulatorDifferentialEvolutionLogging

subroutine simulatorDifferentialEvolutionPosterior(self,simulationState,logLikelihoodCurrent,logPriorCurrent,logPosterior,logLikelihood)
  !% Return the log of the posterior for the current state.
  implicit none
  class           (simulatorDifferentialEvolution), intent(inout) :: self
  class           (state                         ), intent(in   ) :: simulationState
  double precision                                , intent(in   ) :: logLikelihoodCurrent, logPriorCurrent
  double precision                                , intent(  out) :: logLikelihood       , logPosterior
  double precision                                                :: logPrior

  logPrior     =priorsEvaluateLog            (self%parameterPriors,simulationState                                                                 )
  logLikelihood=self%modelLikelihood%evaluate(                     simulationState,self%parameterMappings,self%simulationConvergence,self%temperature(),logLikelihoodCurrent,logPriorCurrent,logPrior)
  logPosterior =logPrior+logLikelihood
  return
end subroutine simulatorDifferentialEvolutionPosterior

double precision function simulatorDifferentialEvolutionStepSize(self)
  !% Return the step size parameter, $\gamma$, for a differential evolution step.
  implicit none
  class(simulatorDifferentialEvolution), intent(inout) :: self
  
  if (mod(self%simulationState%count(),self%stateSwapCount) == 0) then
     ! Every self%stateSwapCount steps, set gamma=1 to allow interchange of chains.
     simulatorDifferentialEvolutionStepSize=1.0d0
  else
     ! Otherwise, use the step-size algorithm.
     simulatorDifferentialEvolutionStepSize=self%                                               &
          &                                      proposalSize%                                  &
          &                                                   gamma(                            &
          &                                                         self%simulationState      , &
          &                                                         self%simulationConvergence  &
          &                                                        )
  end if
  return
end function simulatorDifferentialEvolutionStepSize

logical function simulatorDifferentialEvolutionAcceptProposal(self,logPosterior,logPosteriorProposed,randomNumberGenerator)
  !% Return whether or not to accept a proposal.
  use Pseudo_Random
  use MPI_Utilities
  implicit none
  class           (simulatorDifferentialEvolution), intent(inout) :: self
  double precision                                , intent(in   ) :: logPosterior         , logPosteriorProposed
  type            (pseudoRandom                  ), intent(inout) :: randomNumberGenerator
  double precision                                                :: x

  ! Decide whether to take step.
  x=randomNumberGenerator%sample(mpiRankOffset=.true.)
  simulatorDifferentialEvolutionAcceptProposal=                           &
       &   logPosteriorProposed >      logPosterior                       &
       &  .or.                                                            &
       &   x                    < exp(-logPosterior+logPosteriorProposed)
  return
end function simulatorDifferentialEvolutionAcceptProposal

double precision function simulatorDifferentialEvolutionTemperature(self)
  !% Return the temperature.
  implicit none
  class(simulatorDifferentialEvolution), intent(inout) :: self

  simulatorDifferentialEvolutionTemperature=1.0d0
  return
end function simulatorDifferentialEvolutionTemperature
