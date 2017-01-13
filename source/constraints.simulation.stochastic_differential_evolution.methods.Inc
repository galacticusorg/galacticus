!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
!!    Andrew Benson <abenson@carnegiescience.edu>
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

!% Defines methods for the ``stochasticDifferentialEvolution'' simulator class.

function simulatorStochasticDifferentialEvolutionConstructor(parameterPriors,randomDistributions,parameterMappings,modelLikelihood,simulationConvergence,simulationStoppingCriterion &
          &,simulationState,simulationStateInitializor,proposalSize,randomJump,stepsMaximum&
          &,acceptanceAverageCount,stateSwapCount,logFileRoot,sampleOutliers,logFlushCount,reportCount,temperatureScale,interactionRoot)
  !% Constructor for ``stochasticDifferentialEvolution'' simulator class.
  implicit none
  type            (simulatorStochasticDifferentialEvolution)                                      :: simulatorStochasticDifferentialEvolutionConstructor
  type            (prior                                   ), intent(in   ), target, dimension(:) :: parameterPriors
  type            (distributionList                        ), intent(in   ), target, dimension(:) :: randomDistributions
  type            (mappingList                             ), intent(in   ), target, dimension(:) :: parameterMappings
  class           (likelihood                              ), intent(in   ), target               :: modelLikelihood
  class           (convergence                             ), intent(in   ), target               :: simulationConvergence
  class           (stoppingCriterion                       ), intent(in   ), target               :: simulationStoppingCriterion
  class           (state                                   ), intent(in   ), target               :: simulationState
  class           (stateInitializor                        ), intent(in   ), target               :: simulationStateInitializor
  class           (deProposalSize                          ), intent(in   ), target               :: proposalSize
  class           (deRandomJump                            ), intent(in   ), target               :: randomJump
  integer                                                   , intent(in   )                       :: stepsMaximum                   , stateSwapCount          , &
       &                                                                                             acceptanceAverageCount         , logFlushCount           , &
       &                                                                                             reportCount
  double precision                                          , intent(in   )                       :: temperatureScale
  character       (len=*                                   ), intent(in   )                       :: logFileRoot                    , interactionRoot
  logical                                                   , intent(in   )                       :: sampleOutliers

  simulatorStochasticDifferentialEvolutionConstructor%parameterCount              =  size(parameterPriors)
  simulatorStochasticDifferentialEvolutionConstructor%parameterPriors             => parameterPriors
  simulatorStochasticDifferentialEvolutionConstructor%randomDistributions         => randomDistributions
  simulatorStochasticDifferentialEvolutionConstructor%parameterMappings           => parameterMappings
  simulatorStochasticDifferentialEvolutionConstructor%modelLikelihood             => modelLikelihood
  simulatorStochasticDifferentialEvolutionConstructor%simulationConvergence       => simulationConvergence
  simulatorStochasticDifferentialEvolutionConstructor%simulationStoppingCriterion => simulationStoppingCriterion
  simulatorStochasticDifferentialEvolutionConstructor%simulationState             => simulationState
  simulatorStochasticDifferentialEvolutionConstructor%simulationStateInitializor  => simulationStateInitializor
  simulatorStochasticDifferentialEvolutionConstructor%proposalSize                => proposalSize
  simulatorStochasticDifferentialEvolutionConstructor%randomJump                  => randomJump
  simulatorStochasticDifferentialEvolutionConstructor%stepsMaximum                =  stepsMaximum
  simulatorStochasticDifferentialEvolutionConstructor%acceptanceAverageCount      =  acceptanceAverageCount
  simulatorStochasticDifferentialEvolutionConstructor%stateSwapCount              =  stateSwapCount
  simulatorStochasticDifferentialEvolutionConstructor%logFileRoot                 =  logFileRoot
  simulatorStochasticDifferentialEvolutionConstructor%isInteractive               =  trim(interactionRoot) /= "none"
  simulatorStochasticDifferentialEvolutionConstructor%interactionRoot             =  trim(interactionRoot)
  simulatorStochasticDifferentialEvolutionConstructor%sampleOutliers              =  sampleOutliers
  simulatorStochasticDifferentialEvolutionConstructor%logFlushCount               =  logFlushCount
  simulatorStochasticDifferentialEvolutionConstructor%reportCount                 =  reportCount
  simulatorStochasticDifferentialEvolutionConstructor%temperatureScale            =  temperatureScale
  return
end function simulatorStochasticDifferentialEvolutionConstructor

logical function simulatorStochasticDifferentialEvolutionAcceptProposal(self,logPosterior,logPosteriorProposed,logLikelihoodVariance,logLikelihoodVarianceProposed,randomNumberGenerator)
  !% Return whether or not to accept a proposal.
  use Pseudo_Random
  use MPI_Utilities
  use Galacticus_Error
  use Constraints_Convergence
  implicit none
  class           (simulatorStochasticDifferentialEvolution), intent(inout) :: self
  double precision                                          , intent(in   ) :: logPosterior                , logPosteriorProposed         , &
       &                                                                       logLikelihoodVariance       , logLikelihoodVarianceProposed
  type            (pseudoRandom                            ), intent(inout) :: randomNumberGenerator
  double precision                                                          :: x                           , temperature                  , &
       &                                                                       convergenceMeasure          , convergenceMeasureTarget     , &
       &                                                                       temperatureConvergenceFactor

  ! Find the convergence state of the simulation.
  select type (s => self%simulationConvergence)
  class is (convergenceGelmanRubin)
     convergenceMeasure      =s%convergenceMeasure      ()
     convergenceMeasureTarget=s%convergenceMeasureTarget()
     if (convergenceMeasure <= convergenceMeasureTarget) then
        temperatureConvergenceFactor=0.0d0
     else
        temperatureConvergenceFactor=log(convergenceMeasure/convergenceMeasureTarget)
     end if
  class default
     temperatureConvergenceFactor=0.0d0
     call Galacticus_Error_Report('simulatorStochasticDifferentialEvolutionAcceptProposal','this class requires a Gelman-Rubin convergence criterion')
  end select  
  ! Set the chain temperature to the root-variance of the difference between the current and proposed likelihoods (which is
  ! characteristic amount by which we expect them to differ due to random fluctuations). Multiply by a user-specified factor to
  ! allow control over the temperature.
  temperature=+1.0d0                               &
       &      +temperatureConvergenceFactor        &
       &      *self%temperatureScale               &
       &      *sqrt(                               &
       &            +logLikelihoodVariance         &
       &            +logLikelihoodVarianceProposed &
       &      )
  ! Decide whether to take step.
  x=randomNumberGenerator%sample(mpiRankOffset=.true.)
  simulatorStochasticDifferentialEvolutionAcceptProposal=    &
       &   logPosteriorProposed >      logPosterior          &
       &  .or.                                               &
       &   x                    < exp(                       &
       &                              (                      &
       &                               -logPosterior         &
       &                               +logPosteriorProposed &
       &                              )                      &
       &                              /temperature           &
       &                             )
  return
end function simulatorStochasticDifferentialEvolutionAcceptProposal
