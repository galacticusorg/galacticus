!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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

  !!{
  Implementation of a posterior sampling simulation class which implements a stochastic differential evolution algorithm.
  !!}

  !![
  <posteriorSampleSimulation name="posteriorSampleSimulationStochasticDffrntlEvltn">
   <description>
    This option extends the {\normalfont \ttfamily differentialEvolution} option to run chains at a temperature matched to the
    uncertainty in the log-likelihood. This is designed to work with stochastic likelihood functions where an estimate of the
    uncertainty in the log-likelihood is available, and prevents the chains from becoming trapped in local maxima arising purely from
    random fluctuations.  In addition to the options for the {\normalfont \ttfamily differentialEvolution} algorithm, the details of
    the algorithm are controlled by the following parameters:
    \begin{description}
    \item[{\normalfont \ttfamily temperatureScale}] The temperature scaling factor, $alpha$, described below.
    \end{description}
    In computing the acceptance probability for transitions between states, the chain temperature is set to
    \begin{equation}
      T = 1 + \alpha C \sqrt{\sigma^2_\mathrm{current}+\sigma^2_\mathrm{proposed}},
    \end{equation}
    where $C = \log(\hat{R}/\hat{R}_\mathrm{t})$, $\hat{R}$ is the current maximum (across all parameters) Gelman-Rubin convergence
    statistic, $\hat{R}_\mathrm{t}$ is the target Gelman-Rubin convergence statistic at which convergence will be declared, and
    $\sigma^2_\mathrm{current}$ and $\sigma^2_\mathrm{proposed}$ are the variances in the log-likelihood of the current and proposed
    states respectively. This form ensures that the temperature declines to unity once the chains are converged (such that they will
    sample from the true posterior distribution), while ensuring that $T$ is of order the size of the expected random fluctuations in
    the difference in log-likelihoods between states prior to chain convergence.
   </description>
  </posteriorSampleSimulation>
  !!]
  type, extends(posteriorSampleSimulationDifferentialEvolution) :: posteriorSampleSimulationStochasticDffrntlEvltn
     !!{
     Implementation of a posterior sampling simulation class which implements a stochastic differential evolution algorithm.
     !!}
     private
     double precision :: temperatureScale
   contains
     !![
     <methods>
       <method description="Initialize the object." method="initialize" />
     </methods>
     !!]
     procedure :: acceptProposal => stochasticDifferentialEvolutionAcceptProposal
     procedure :: initialize     => stochasticDifferentialEvolutionInitialize
  end type posteriorSampleSimulationStochasticDffrntlEvltn

  interface posteriorSampleSimulationStochasticDffrntlEvltn
     !!{
     Constructors for the \refClass{posteriorSampleSimulationStochasticDffrntlEvltn} posterior sampling convergence class.
     !!}
     module procedure stochasticDifferentialEvolutionConstructorParameters
     module procedure stochasticDifferentialEvolutionConstructorInternal
  end interface posteriorSampleSimulationStochasticDffrntlEvltn

contains

  function stochasticDifferentialEvolutionConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{posteriorSampleSimulationStochasticDffrntlEvltn} posterior sampling simulation class which builds the object from a
    parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (posteriorSampleSimulationStochasticDffrntlEvltn)                :: self
    type            (inputParameters                                ), intent(inout) :: parameters
    double precision                                                                 :: temperatureScale

    self%posteriorSampleSimulationDifferentialEvolution=posteriorSampleSimulationDifferentialEvolution(parameters)
    !![
    <inputParameter>
      <name>temperatureScale</name>
      <description>The temperature scale.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    call self%initialize(temperatureScale)
    !![
    <inputParametersValidate source="parameters" multiParameters="modelParameter"/>
    !!]
    return
  end function stochasticDifferentialEvolutionConstructorParameters

  function stochasticDifferentialEvolutionConstructorInternal(modelParametersActive_,modelParametersInactive_,posteriorSampleLikelihood_,posteriorSampleConvergence_,posteriorSampleStoppingCriterion_,posteriorSampleState_,posteriorSampleStateInitialize_,posteriorSampleDffrntlEvltnProposalSize_,posteriorSampleDffrntlEvltnRandomJump_,randomNumberGenerator_,stepsMaximum,acceptanceAverageCount,stateSwapCount,slowStepCount,recomputeCount,logFileRoot,sampleOutliers,logFlushCount,reportCount,interactionRoot,appendLogs,loadBalance,ignoreChainNumberAdvice,temperatureScale) result(self)
    !!{
    Internal constructor for the ``stochasticDifferentialEvolution'' simulation class.
    !!}
    implicit none
    type            (posteriorSampleSimulationStochasticDffrntlEvltn)                                      :: self
    type            (modelParameterList                             ), intent(in   ), target, dimension(:) :: modelParametersActive_                  , modelParametersInactive_
    class           (posteriorSampleLikelihoodClass                 ), intent(in   ), target               :: posteriorSampleLikelihood_
    class           (posteriorSampleConvergenceClass                ), intent(in   ), target               :: posteriorSampleConvergence_
    class           (posteriorSampleStoppingCriterionClass          ), intent(in   ), target               :: posteriorSampleStoppingCriterion_
    class           (posteriorSampleStateClass                      ), intent(in   ), target               :: posteriorSampleState_
    class           (posteriorSampleStateInitializeClass            ), intent(in   ), target               :: posteriorSampleStateInitialize_
    class           (posteriorSampleDffrntlEvltnProposalSizeClass   ), intent(in   ), target               :: posteriorSampleDffrntlEvltnProposalSize_
    class           (posteriorSampleDffrntlEvltnRandomJumpClass     ), intent(in   ), target               :: posteriorSampleDffrntlEvltnRandomJump_
    class           (randomNumberGeneratorClass                     ), intent(in   ), target               :: randomNumberGenerator_
    integer                                                          , intent(in   )                       :: stepsMaximum                            , acceptanceAverageCount  , &
         &                                                                                                    stateSwapCount                          , logFlushCount           , &
         &                                                                                                    reportCount                             , recomputeCount          , &
         &                                                                                                    slowStepCount
    character       (len=*                                          ), intent(in   )                       :: logFileRoot                             , interactionRoot
    logical                                                          , intent(in   )                       :: sampleOutliers                          , appendLogs              , &
         &                                                                                                    loadBalance                             , ignoreChainNumberAdvice
    double precision                                                 , intent(in   )                       :: temperatureScale

    self%posteriorSampleSimulationDifferentialEvolution=posteriorSampleSimulationDifferentialEvolution(modelParametersActive_,modelParametersInactive_,posteriorSampleLikelihood_,posteriorSampleConvergence_,posteriorSampleStoppingCriterion_,posteriorSampleState_,posteriorSampleStateInitialize_,posteriorSampleDffrntlEvltnProposalSize_,posteriorSampleDffrntlEvltnRandomJump_,randomNumberGenerator_,stepsMaximum,acceptanceAverageCount,stateSwapCount,slowStepCount,recomputeCount,logFileRoot,sampleOutliers,logFlushCount,reportCount,interactionRoot,appendLogs,loadBalance,ignoreChainNumberAdvice)
    call self%initialize(temperatureScale)
    return
  end function stochasticDifferentialEvolutionConstructorInternal

  subroutine stochasticDifferentialEvolutionInitialize(self,temperatureScale)
    !!{
    Finished initialization of stochastic differential evolution simulation objects during construction.
    !!}
    implicit none
    class           (posteriorSampleSimulationStochasticDffrntlEvltn), intent(inout) :: self
    double precision                                                 , intent(in   ) :: temperatureScale

    self%temperatureScale=temperatureScale
    return
  end subroutine stochasticDifferentialEvolutionInitialize
  
  logical function stochasticDifferentialEvolutionAcceptProposal(self,logPosterior,logPosteriorProposed,logLikelihoodVariance,logLikelihoodVarianceProposed)
    !!{
    Return whether or not to accept a proposal.
    !!}
    use :: Error                         , only : Error_Report
    use :: Posterior_Sampling_Convergence, only : posteriorSampleConvergenceGelmanRubin
    implicit none
    class           (posteriorSampleSimulationStochasticDffrntlEvltn), intent(inout) :: self
    double precision                                                 , intent(in   ) :: logPosterior                , logPosteriorProposed         , &
         &                                                                              logLikelihoodVariance       , logLikelihoodVarianceProposed
    double precision                                                                 :: x                           , temperature                  , &
         &                                                                              convergenceMeasure          , convergenceMeasureTarget     , &
         &                                                                              temperatureConvergenceFactor

    ! Find the convergence state of the simulation.
    select type (posteriorSampleConvergence_ => self%posteriorSampleConvergence_)
       class is (posteriorSampleConvergenceGelmanRubin)
       convergenceMeasure      =posteriorSampleConvergence_%convergenceMeasure      ()
       convergenceMeasureTarget=posteriorSampleConvergence_%convergenceMeasureTarget()
       if (convergenceMeasure <= convergenceMeasureTarget) then
          temperatureConvergenceFactor=0.0d0
       else
          temperatureConvergenceFactor=log(convergenceMeasure/convergenceMeasureTarget)
       end if
       class default
       temperatureConvergenceFactor=0.0d0
       call Error_Report('this class requires a Gelman-Rubin convergence criterion'//{introspection:location})
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
    x=self%randomNumberGenerator_%uniformSample()
    stochasticDifferentialEvolutionAcceptProposal=             &
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
  end function stochasticDifferentialEvolutionAcceptProposal
