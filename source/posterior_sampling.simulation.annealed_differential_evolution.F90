!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
  Implementation of a posterior sampling simulation class which implements an annealed differential evolution algorithm.
  !!}

  !![
  <posteriorSampleSimulation name="posteriorSampleSimulationAnnealedDffrntlEvltn">
   <description>
    This class extends the {\normalfont \ttfamily differentialEvolution} class to include an annealing schedule---the simulation begins at
    high temperature, waits for convergence, lowers the temperature and repeats until convergence at $T=1$ is reached. In addition to
    the options for the {\normalfont \ttfamily differentialEvolution} algorithm, the details of the algorithm are controlled by the
    following sub-parameters:
    \begin{description}
    \item[{\normalfont \ttfamily [temperatureMaximum]}] The maximum temperature to use when tempering.
    \item[{\normalfont \ttfamily [temperatureLevels]}] The number of temperature levels to use.
    \end{description}
    
    The temperature at level $i$ is given by:
    \begin{equation}
    \log T_i = {i-1 \over N-1} \log T_\mathrm{max},
    \end{equation}
    where $T_\mathrm{max}=${\normalfont \ttfamily [temperatureMaximum]} and $N=${\normalfont \ttfamily [temperatureLevels]}.
   </description>
  </posteriorSampleSimulation>
  !!]
  type, extends(posteriorSampleSimulationDifferentialEvolution) :: posteriorSampleSimulationAnnealedDffrntlEvltn
     !!{
     Implementation of a posterior sampling simulation class which implements an annealed differential evolution algorithm.
     !!}
     private
     integer                                     :: temperatureLevelCount, temperatureLevelCurrent
     double precision                            :: temperatureMaximum
     double precision, allocatable, dimension(:) :: temperatures
   contains
     !![
     <methods>
       <method description="Initialize the object." method="initialize" />
     </methods>
     !!]
     procedure :: acceptProposal => annealedDifferentialEvolutionAcceptProposal
     procedure :: update         => annealedDifferentialEvolutionUpdate
     procedure :: temperature    => annealedDifferentialEvolutionTemperature
     procedure :: initialize     => annealedDifferentialEvolutionInitialize
  end type posteriorSampleSimulationAnnealedDffrntlEvltn

  interface posteriorSampleSimulationAnnealedDffrntlEvltn
     !!{
     Constructors for the \refClass{posteriorSampleSimulationAnnealedDffrntlEvltn} posterior sampling convergence class.
     !!}
     module procedure annealedDifferentialEvolutionConstructorParameters
     module procedure annealedDifferentialEvolutionConstructorInternal
  end interface posteriorSampleSimulationAnnealedDffrntlEvltn

contains

  function annealedDifferentialEvolutionConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{posteriorSampleSimulationAnnealedDffrntlEvltn} posterior sampling simulation class which builds the object from a
    parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (posteriorSampleSimulationAnnealedDffrntlEvltn)                :: self
    type            (inputParameters                              ), intent(inout) :: parameters
    integer                                                                        :: temperatureLevelCount
    double precision                                                               :: temperatureMaximum

    self%posteriorSampleSimulationDifferentialEvolution=posteriorSampleSimulationDifferentialEvolution(parameters)
    !![
    <inputParameter>
      <name>temperatureLevelCount</name>
      <defaultValue>10</defaultValue>
      <description>The number temperature levels to use.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>temperatureMaximum</name>
      <description>The maximum temperature to reach.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    call self%initialize(temperatureLevelCount,temperatureMaximum)
    !![
    <inputParametersValidate source="parameters" multiParameters="modelParameter"/>
    !!]
    return
  end function annealedDifferentialEvolutionConstructorParameters

  function annealedDifferentialEvolutionConstructorInternal(modelParametersActive_,modelParametersInactive_,posteriorSampleLikelihood_,posteriorSampleConvergence_,posteriorSampleStoppingCriterion_,posteriorSampleState_,posteriorSampleStateInitialize_,posteriorSampleDffrntlEvltnProposalSize_,posteriorSampleDffrntlEvltnRandomJump_,randomNumberGenerator_,stepsMaximum,acceptanceAverageCount,stateSwapCount,slowStepCount,recomputeCount,logFileRoot,sampleOutliers,logFlushCount,reportCount,interactionRoot,appendLogs,loadBalance,ignoreChainNumberAdvice,temperatureLevelCount,temperatureMaximum) result(self)
    !!{
    Internal constructor for the ``annealedDifferentialEvolution'' simulation class.
    !!}
    implicit none
    type            (posteriorSampleSimulationAnnealedDffrntlEvltn)                                      :: self
    type            (modelParameterList                           ), intent(in   ), target, dimension(:) :: modelParametersActive_                  , modelParametersInactive_
    class           (posteriorSampleLikelihoodClass               ), intent(in   ), target               :: posteriorSampleLikelihood_
    class           (posteriorSampleConvergenceClass              ), intent(in   ), target               :: posteriorSampleConvergence_
    class           (posteriorSampleStoppingCriterionClass        ), intent(in   ), target               :: posteriorSampleStoppingCriterion_
    class           (posteriorSampleStateClass                    ), intent(in   ), target               :: posteriorSampleState_
    class           (posteriorSampleStateInitializeClass          ), intent(in   ), target               :: posteriorSampleStateInitialize_
    class           (posteriorSampleDffrntlEvltnProposalSizeClass ), intent(in   ), target               :: posteriorSampleDffrntlEvltnProposalSize_
    class           (posteriorSampleDffrntlEvltnRandomJumpClass   ), intent(in   ), target               :: posteriorSampleDffrntlEvltnRandomJump_
    class           (randomNumberGeneratorClass                   ), intent(in   ), target               :: randomNumberGenerator_
    integer                                                        , intent(in   )                       :: stepsMaximum                            , acceptanceAverageCount  , &
         &                                                                                                  stateSwapCount                          , logFlushCount           , &
         &                                                                                                  reportCount                             , temperatureLevelCount   , &
         &                                                                                                  recomputeCount                          , slowStepCount
    character       (len=*                                        ), intent(in   )                       :: logFileRoot                             , interactionRoot
    logical                                                        , intent(in   )                       :: sampleOutliers                          , appendLogs              , &
         &                                                                                                  loadBalance                             , ignoreChainNumberAdvice
    double precision                                               , intent(in   )                       :: temperatureMaximum

    self%posteriorSampleSimulationDifferentialEvolution=posteriorSampleSimulationDifferentialEvolution(modelParametersActive_,modelParametersInactive_,posteriorSampleLikelihood_,posteriorSampleConvergence_,posteriorSampleStoppingCriterion_,posteriorSampleState_,posteriorSampleStateInitialize_,posteriorSampleDffrntlEvltnProposalSize_,posteriorSampleDffrntlEvltnRandomJump_,randomNumberGenerator_,stepsMaximum,acceptanceAverageCount,stateSwapCount,slowStepCount,recomputeCount,logFileRoot,sampleOutliers,logFlushCount,reportCount,interactionRoot,appendLogs,loadBalance,ignoreChainNumberAdvice)
    call self%initialize(temperatureLevelCount,temperatureMaximum)
    return
  end function annealedDifferentialEvolutionConstructorInternal

  subroutine annealedDifferentialEvolutionInitialize(self,temperatureLevelCount,temperatureMaximum)
    !!{
    Finished initialization of annealed differential evolution simulation objects during construction.
    !!}
    implicit none
    class           (posteriorSampleSimulationAnnealedDffrntlEvltn), intent(inout) :: self
    integer                                                        , intent(in   ) :: temperatureLevelCount
    double precision                                               , intent(in   ) :: temperatureMaximum
    integer                                                                        :: i

    self%temperatureLevelCurrent=temperatureLevelCount
    self%temperatureLevelCount  =temperatureLevelCount
    self%temperatureMaximum     =temperatureMaximum
    allocate(self%temperatures(temperatureLevelCount))
    if (temperatureLevelCount == 1) then
       self%temperatures(1)=1.0d0
    else
       do i=1,temperatureLevelCount
          self%temperatures(i)                        &
               &  =exp(                               &
               &        log (temperatureMaximum     ) &
               &       *dble(i                    -1) &
               &       /dble(temperatureLevelCount-1) &
               &      )
       end do
    end if
    return
  end subroutine annealedDifferentialEvolutionInitialize

  subroutine annealedDifferentialEvolutionUpdate(self,stateVector)
    !!{
    Update the differential evolution simulator state.
    !!}
    use :: Display           , only : displayMessage
    use :: ISO_Varying_String, only : varying_string
    use :: MPI_Utilities     , only : mpiSelf
    use :: String_Handling   , only : operator(//)
    implicit none
    class           (posteriorSampleSimulationAnnealedDffrntlEvltn), intent(inout)                                 :: self
    double precision                                               , intent(in   ), dimension(self%parameterCount) :: stateVector
    logical                                                        , allocatable  , dimension(:                  ) :: outlierMask
    integer                                                                                                        :: i
    type            (varying_string                               )                                                :: message
    character       (len=7                                        )                                                :: label

    ! Check for convergence.
    if (self%posteriorSampleConvergence_%isConverged(self%posteriorSampleState_,self%logPosterior).and.self%temperatureLevelCurrent > 1) then
       ! Decrease the temperature level.
       self%temperatureLevelCurrent=self%temperatureLevelCurrent-1
       ! Reset state.
       call self%posteriorSampleState_      %reset()
       call self%posteriorSampleConvergence_%reset()
       self%isConverged=.false.
       ! Inform likelihood object that likelihood function may have changed.
       call self%posteriorSampleLikelihood_%functionChanged()
       ! Report.
       if (mpiSelf%isMaster()) then
          write (label,'(f7.2)') self%temperature()
          message="Annealing temperature is now "//label//" (level "
          message=message//self%temperatureLevelCurrent//" of "//self%temperatureLevelCount//")"
          call displayMessage(message)
       end if
    end if
    ! Update the simulation state.
    allocate(outlierMask(0:mpiSelf%count()-1))
    do i=0,mpiSelf%count()-1
       outlierMask(i)=self%posteriorSampleConvergence_%stateIsOutlier(i)
    end do
    call self%posteriorSampleState_%update(stateVector,self%logging(),self%posteriorSampleConvergence_%isConverged(),outlierMask)
    return
  end subroutine annealedDifferentialEvolutionUpdate

  double precision function annealedDifferentialEvolutionTemperature(self)
    !!{
    Return the temperature.
    !!}
    implicit none
    class(posteriorSampleSimulationAnnealedDffrntlEvltn), intent(inout) :: self

    annealedDifferentialEvolutionTemperature=self%temperatures(self%temperatureLevelCurrent)
    return
  end function annealedDifferentialEvolutionTemperature

  logical function annealedDifferentialEvolutionAcceptProposal(self,logPosterior,logPosteriorProposed,logLikelihoodVariance,logLikelihoodVarianceProposed)
    !!{
    Return whether or not to accept a proposal.
    !!}
    implicit none
    class           (posteriorSampleSimulationAnnealedDffrntlEvltn), intent(inout) :: self
    double precision                                               , intent(in   ) :: logPosterior         , logPosteriorProposed         , &
         &                                                                            logLikelihoodVariance, logLikelihoodVarianceProposed
    double precision                                                               :: x
    !$GLC attributes unused :: logLikelihoodVariance, logLikelihoodVarianceProposed

    ! Decide whether to take step.
    x=self%randomNumberGenerator_%uniformSample()
    annealedDifferentialEvolutionAcceptProposal=               &
         &   logPosteriorProposed >      logPosterior          &
         &  .or.                                               &
         &   x                    < exp(                       &
         &                              (                      &
         &                               -logPosterior         &
         &                               +logPosteriorProposed &
         &                              )                      &
         &                              /self%temperature()    &
         &                             )
    return
  end function annealedDifferentialEvolutionAcceptProposal
