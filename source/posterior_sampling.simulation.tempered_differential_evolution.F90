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
  Implementation of a posterior sampling simulation class which implements a tempered differential evolution algorithm.
  !!}

  use :: Posterior_Sampling_Prop_Size_Temp_Exp, only : posteriorSampleDffrntlEvltnPrpslSzTmpExpClass

  !![
  <posteriorSampleSimulation name="posteriorSampleSimulationTemperedDffrntlEvltn">
   <description>
    This class extends the {\normalfont \ttfamily differentialEvolution} option to include tempering during which the likelihood
    function is heated up and cooled down to allow chains to more easily walk through the likelihood landscape. In addition to the
    options for the {\normalfont \ttfamily differentialEvolution} algorithm, the details of the algorithm are controlled by the
    following sub-parameters:
    \begin{description}
    \item[{\normalfont \ttfamily [untemperedStepCount]}] The number of untempered (i.e. $T=1$) steps to take between tempering cycles.
    \item[{\normalfont \ttfamily [temperatureMaximum]}] The maximum temperature to use when tempering.
    \item[{\normalfont \ttfamily [temperedLevels]}] The number of tempered levels to use.
    \item[{\normalfont \ttfamily [stepsPerLevel]}] The number of differential evolution steps to take at each tempering level.
    \item[{\normalfont \ttfamily [logFlushCount]}] The number of steps after which the log file will be flushed to disk.
    \end{description}
    
    In each tempering cycle, the temperature is raised through levels $1$\ldots$N$ (where $N=${\normalfont \ttfamily temperedLevels}),
    and then back down through levels $N-1$\ldots$1$. The temperature at level $i$ is given by:
    \begin{equation}
    \log T_i = {i \over N} \log T_\mathrm{max},
    \end{equation}
    where $T_\mathrm{max}=${\normalfont \ttfamily temperatureMaximum}. During tempered steps, the $\gamma$ parameter of the
    differential evolution algorithm is increased by a factor $T^\alpha$, where $\alpha$ is provided by the {\normalfont \ttfamily
    proposalSizeTemperatureExponent} class. A value of $\alpha=1/2$ is optimal for a Gaussian likelihood.
   </description>
  </posteriorSampleSimulation>
  !!]
  type, extends(posteriorSampleSimulationDifferentialEvolution) :: posteriorSampleSimulationTemperedDffrntlEvltn
     !!{
     Implementation of a posterior sampling simulation class which implements a tempered differential evolution algorithm.
     !!}
     private
     integer                                                                                    :: untemperedStepCount                                , temperingLevelCount    , &
          &                                                                                        stepsPerLevel
     integer                                                                                    :: temperingStep                                      , temperingLevelMonotonic
     double precision                                                                           :: temperatureMaximum
     class           (posteriorSampleDffrntlEvltnPrpslSzTmpExpClass), pointer                   :: posteriorSampleDffrntlEvltnPrpslSzTmpExp_ => null()
     double precision                                               , allocatable, dimension(:) :: temperatures
     class           (posteriorSampleStateClass                    ), allocatable, dimension(:) :: temperedStates
   contains
     !![
     <methods>
       <method method="initialize" description="Return the current tempering level."/>
       <method method="level"      description="Return the current tempering level."/>
     </methods>
     !!]
     final     ::                   temperedDifferentialEvolutionDestructor
     procedure :: logging        => temperedDifferentialEvolutionLogging
     procedure :: acceptProposal => temperedDifferentialEvolutionAcceptProposal
     procedure :: update         => temperedDifferentialEvolutionUpdate
     procedure :: stepSize       => temperedDifferentialEvolutionStepSize
     procedure :: level          => temperedDifferentialEvolutionLevel
     procedure :: temperature    => temperedDifferentialEvolutionTemperature
     procedure :: initialize     => temperedDifferentialEvolutionInitialize
  end type posteriorSampleSimulationTemperedDffrntlEvltn

  interface posteriorSampleSimulationTemperedDffrntlEvltn
     !!{
     Constructors for the \refClass{posteriorSampleSimulationTemperedDffrntlEvltn} posterior sampling convergence class.
     !!}
     module procedure temperedDifferentialEvolutionConstructorParameters
     module procedure temperedDifferentialEvolutionConstructorInternal
  end interface posteriorSampleSimulationTemperedDffrntlEvltn

contains

  function temperedDifferentialEvolutionConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{posteriorSampleSimulationTemperedDffrntlEvltn} posterior sampling simulation class which builds the object from a
    parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (posteriorSampleSimulationTemperedDffrntlEvltn)                :: self
    type            (inputParameters                              ), intent(inout) :: parameters
    class           (posteriorSampleDffrntlEvltnPrpslSzTmpExpClass), pointer       :: posteriorSampleDffrntlEvltnPrpslSzTmpExp_
    integer                                                                        :: temperingLevelCount                      , untemperedStepCount, &
         &                                                                            stepsPerLevel
    double precision                                                               :: temperatureMaximum

    self%posteriorSampleSimulationDifferentialEvolution=posteriorSampleSimulationDifferentialEvolution(parameters)
    !![
    <inputParameter>
      <name>untemperedStepCount</name>
      <defaultValue>10</defaultValue>
      <description>The number of untempered steps to take.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>stepsPerLevel</name>
      <defaultValue>10</defaultValue>
      <description>The number of steps to take at each tempering level.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>temperingLevelCount</name>
      <defaultValue>10</defaultValue>
      <description>The number tempering levels to use.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>temperatureMaximum</name>
      <description>The maximum temperature to reach.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="posteriorSampleDffrntlEvltnPrpslSzTmpExp" name="posteriorSampleDffrntlEvltnPrpslSzTmpExp_" source="parameters"/>
    !!]
    call self%initialize(posteriorSampleDffrntlEvltnPrpslSzTmpExp_,temperingLevelCount,untemperedStepCount,stepsPerLevel,temperatureMaximum)
    !![
    <inputParametersValidate source="parameters" multiParameters="modelParameter"/>
    <objectDestructor name="posteriorSampleDffrntlEvltnPrpslSzTmpExp_"/>
    !!]
    return
  end function temperedDifferentialEvolutionConstructorParameters

  function temperedDifferentialEvolutionConstructorInternal(modelParametersActive_,modelParametersInactive_,posteriorSampleLikelihood_,posteriorSampleConvergence_,posteriorSampleStoppingCriterion_,posteriorSampleState_,posteriorSampleStateInitialize_,posteriorSampleDffrntlEvltnProposalSize_,posteriorSampleDffrntlEvltnRandomJump_,posteriorSampleDffrntlEvltnPrpslSzTmpExp_,randomNumberGenerator_,stepsMaximum,acceptanceAverageCount,stateSwapCount,recomputeCount,logFileRoot,sampleOutliers,logFlushCount,reportCount,interactionRoot,appendLogs,loadBalance,ignoreChainNumberAdvice,temperingLevelCount,untemperedStepCount,stepsPerLevel,temperatureMaximum) result(self)
    !!{
    Internal constructor for the ``temperedDifferentialEvolution'' simulation class.
    !!}
    implicit none
    type            (posteriorSampleSimulationTemperedDffrntlEvltn)                                      :: self
    type            (modelParameterList                           ), intent(in   ), target, dimension(:) :: modelParametersActive_                  , modelParametersInactive_
    class           (posteriorSampleLikelihoodClass               ), intent(in   ), target               :: posteriorSampleLikelihood_
    class           (posteriorSampleConvergenceClass              ), intent(in   ), target               :: posteriorSampleConvergence_
    class           (posteriorSampleStoppingCriterionClass        ), intent(in   ), target               :: posteriorSampleStoppingCriterion_
    class           (posteriorSampleStateClass                    ), intent(in   ), target               :: posteriorSampleState_
    class           (posteriorSampleStateInitializeClass          ), intent(in   ), target               :: posteriorSampleStateInitialize_
    class           (posteriorSampleDffrntlEvltnProposalSizeClass ), intent(in   ), target               :: posteriorSampleDffrntlEvltnProposalSize_
    class           (posteriorSampleDffrntlEvltnRandomJumpClass   ), intent(in   ), target               :: posteriorSampleDffrntlEvltnRandomJump_
    class           (posteriorSampleDffrntlEvltnPrpslSzTmpExpClass), intent(in   ), target               :: posteriorSampleDffrntlEvltnPrpslSzTmpExp_
    class           (randomNumberGeneratorClass                   ), intent(in   ), target               :: randomNumberGenerator_
    integer                                                        , intent(in   )                       :: stepsMaximum                             , acceptanceAverageCount  , &
         &                                                                                                  stateSwapCount                           , logFlushCount           , &
         &                                                                                                  reportCount                              , temperingLevelCount     , &
         &                                                                                                  untemperedStepCount                      , stepsPerLevel           , &
         &                                                                                                  recomputeCount
    character       (len=*                                        ), intent(in   )                       :: logFileRoot                              , interactionRoot
    logical                                                        , intent(in   )                       :: sampleOutliers                           , appendLogs              , &
         &                                                                                                  loadBalance                              , ignoreChainNumberAdvice
    double precision                                               , intent(in   )                       :: temperatureMaximum

    self%posteriorSampleSimulationDifferentialEvolution=posteriorSampleSimulationDifferentialEvolution(modelParametersActive_,modelParametersInactive_,posteriorSampleLikelihood_,posteriorSampleConvergence_,posteriorSampleStoppingCriterion_,posteriorSampleState_,posteriorSampleStateInitialize_,posteriorSampleDffrntlEvltnProposalSize_,posteriorSampleDffrntlEvltnRandomJump_,randomNumberGenerator_,stepsMaximum,acceptanceAverageCount,stateSwapCount,recomputeCount,logFileRoot,sampleOutliers,logFlushCount,reportCount,interactionRoot,appendLogs,loadBalance,ignoreChainNumberAdvice)
    call self%initialize(posteriorSampleDffrntlEvltnPrpslSzTmpExp_,temperingLevelCount,untemperedStepCount,stepsPerLevel,temperatureMaximum)
    return
  end function temperedDifferentialEvolutionConstructorInternal

  subroutine temperedDifferentialEvolutionInitialize(self,posteriorSampleDffrntlEvltnPrpslSzTmpExp_,temperingLevelCount,untemperedStepCount,stepsPerLevel,temperatureMaximum)
    !!{
    Finished initialization of tempered differential evolution simulation objects during construction.
    !!}
    use :: Posterior_Sampling_State, only : posteriorSampleStateSimple
    implicit none
    class           (posteriorSampleSimulationTemperedDffrntlEvltn), intent(inout)         :: self
    class           (posteriorSampleDffrntlEvltnPrpslSzTmpExpClass), intent(in   ), target :: posteriorSampleDffrntlEvltnPrpslSzTmpExp_
    integer                                                        , intent(in   )         :: temperingLevelCount                      , untemperedStepCount, &
         &                                                                                    stepsPerLevel
    double precision                                               , intent(in   )         :: temperatureMaximum
    integer                                                                                :: i

    self%posteriorSampleDffrntlEvltnPrpslSzTmpExp_ => posteriorSampleDffrntlEvltnPrpslSzTmpExp_
    self%temperingLevelCount                       =  temperingLevelCount
    self%untemperedStepCount                       =  untemperedStepCount
    self%stepsPerLevel                             =  stepsPerLevel
    self%temperatureMaximum                        =  temperatureMaximum
    self%temperingLevelMonotonic                   =  0
    self%temperingStep                             =  0
    allocate(self%temperatures(temperingLevelCount))
    allocate(posteriorSampleStateSimple :: self%temperedStates(temperingLevelCount))
    do i=1,temperingLevelCount
       self%temperatures(i) &
            &  =exp(                                                     &
            &        log (temperatureMaximum )                           &
            &       *dble(i                  )                           &
            &       /dble(temperingLevelCount)                           &
            &      )
       select type (posteriorSampleState_ => self%temperedStates(i))
       type is (posteriorSampleStateSimple)
          posteriorSampleState_=posteriorSampleStateSimple(self%acceptanceAverageCount)
          call posteriorSampleState_%parameterCountSet(size(self%modelParametersActive_))
       end select
    end do
    return
  end subroutine temperedDifferentialEvolutionInitialize

  subroutine temperedDifferentialEvolutionDestructor(self)
    !!{
    Destroy a tempered differential evolution simulation object.
    !!}
    implicit none
    type(posteriorSampleSimulationTemperedDffrntlEvltn), intent(inout) :: self

    !![
    <objectDestructor name="self%posteriorSampleDffrntlEvltnPrpslSzTmpExp_"/>
    !!]
    return
  end subroutine temperedDifferentialEvolutionDestructor

  logical function temperedDifferentialEvolutionLogging(self)
    !!{
    Specifies whether or not the current state should be logged to file during differential evolution.
    !!}
    implicit none
    class(posteriorSampleSimulationTemperedDffrntlEvltn), intent(inout) :: self

    temperedDifferentialEvolutionLogging=(self%temperingLevelMonotonic == 0)
    return
  end function temperedDifferentialEvolutionLogging

  subroutine temperedDifferentialEvolutionUpdate(self,stateVector)
    !!{
    Update the differential evolution simulator state.
    !!}
    use :: Display           , only : displayIndent     , displayMessage, displayUnindent, displayVerbosity, &
          &                           verbosityLevelInfo
    use :: ISO_Varying_String, only : varying_string
    use :: MPI_Utilities     , only : mpiSelf
    use :: String_Handling   , only : operator(//)
    implicit none
    class           (posteriorSampleSimulationTemperedDffrntlEvltn), intent(inout)                                 :: self
    double precision                                               , intent(in   ), dimension(self%parameterCount) :: stateVector
    logical                                                        , allocatable  , dimension(:                  ) :: outlierMask
    integer                                                                                                        :: i             , temperingLevelSaved
    logical                                                                                                        :: levelChanged  , forceAcceptance
    double precision                                                                                               :: acceptanceRate, temperature        , &
         &                                                                                                            stepSize
    character       (len=30                                       )                                                :: label
    type            (varying_string                               )                                                :: message

    ! Update the simulation state.
    allocate(outlierMask(0:mpiSelf%count()-1))
    do i=0,mpiSelf%count()-1
       outlierMask(i)=self%posteriorSampleConvergence_%stateIsOutlier(i)
    end do
    call self%posteriorSampleState_%update(stateVector,self%logging(),self%posteriorSampleConvergence_%isConverged(),outlierMask)
    ! Update tempering step count and level as necessary.
    levelChanged=.false.
    self%temperingStep=self%temperingStep+1
    if (self%temperingLevelMonotonic == 0) then
       ! Currently not tempering: check if we've run all required untempered steps.
       if (self%temperingStep > self%untemperedStepCount) then
          ! We have, switch to the first tempered level.
          self%temperingStep          =0
          self%temperingLevelMonotonic=1
          levelChanged                =.true.
       end if
    else
       ! Update the tempered level states.
       call self%temperedStates(self%level())%update(stateVector,.true.,self%posteriorSampleConvergence_%isConverged(),outlierMask)
       ! Currently tempering: check if we've run all required steps at this level.
       if (self%temperingStep > self%stepsPerLevel) then
          ! We have, move to the next tempering level. Note that we run tempering levels from 1 to
          ! 2*temperingLevelCount-1 - levels above temperingLevelCount represent the cooling phase.
          self%temperingStep          =0
          self%temperingLevelMonotonic=self%temperingLevelMonotonic+1
          levelChanged                =.true.
          ! Check if we've finished tempering and switch back to untempered evolution is so.
          if (self%temperingLevelMonotonic > 2*self%temperingLevelCount-1) self%temperingLevelMonotonic=0
       end if
    end if
    ! Check for change in level.
    if (levelChanged) then
       if (mpiSelf%isMaster().and.displayVerbosity() >= verbosityLevelInfo) then
          write (label,'(f8.1)') self%temperature()
          message='Tempering state: level='
          message=message//self%level()//'; temperature='//trim(label)
          call displayMessage(message)
       end if
       if (self%temperingLevelMonotonic == 0) then
          ! We've just returned to the untempered level. Report on acceptance rates in tempered levels.
          if (displayVerbosity() >= verbosityLevelInfo) then
             if (mpiSelf%isMaster()) then
                call displayIndent('Acceptance rates in tempered levels')
                call displayMessage('Level Temperature  Gamma  Rate')
                call displayMessage('------------------------------')
             end if
             ! Store the current tempering level so that we can restore it below.
             temperingLevelSaved=self%temperingLevelMonotonic
             do i=1,self%temperingLevelCount
                self%temperingLevelMonotonic=i
                forceAcceptance=.false.
                acceptanceRate =mpiSelf%average(self%temperedStates(i)%acceptanceRate(               ))
                temperature    =                self                  %temperature   (               )
                stepSize       =                self                  %stepSize      (forceAcceptance)
                if (mpiSelf%isMaster())  then
                   write (label,'(2x,i3,4x,f8.1,1x,f6.3,1x,f5.3)') i,temperature,stepSize,acceptanceRate
                   call displayMessage(label)
                end if
             end do
             if (mpiSelf%isMaster()) call displayUnindent('done')
             ! Restore the tempering level to the original.
             self%temperingLevelMonotonic=temperingLevelSaved
          end if
       else
          ! We're in a tempering level, update the state of the new level to the current state, without logging.
          call self%temperedStates(self%level())%update(stateVector,.false.,self%posteriorSampleConvergence_%isConverged(),outlierMask)
       end if
    end if
    return
  end subroutine temperedDifferentialEvolutionUpdate

  integer function temperedDifferentialEvolutionLevel(self)
    !!{
    Return the actual tempering level.
    !!}
    implicit none
    class(posteriorSampleSimulationTemperedDffrntlEvltn), intent(inout) :: self

    temperedDifferentialEvolutionLevel=self%temperingLevelMonotonic
    if (self%temperingLevelMonotonic > self%temperingLevelCount)            &
         & temperedDifferentialEvolutionLevel=+2                            &
         &                                    *self%temperingLevelCount     &
         &                                    -self%temperingLevelMonotonic
    return
  end function temperedDifferentialEvolutionLevel

  double precision function temperedDifferentialEvolutionStepSize(self,forceAcceptance)
    !!{
    Return the step size parameter, $\gamma$, for a differential evolution step.
    !!}
    implicit none
    class           (posteriorSampleSimulationTemperedDffrntlEvltn), intent(inout) :: self
    logical                                                        , intent(inout) :: forceAcceptance
    double precision                                                               :: gammaBoostFactor

    if (self%recomputeCount > 0 .and. mod(self%posteriorSampleState_%count(),self%recomputeCount) == 0) then
       ! Every self%recomputeCount steps, set γ=0 and force likelihood to be recomputed in the current state.
       temperedDifferentialEvolutionStepSize=0.0d0
       forceAcceptance                      =.true.
    else if (mod(self%posteriorSampleState_%count(),self%stateSwapCount) == 0 .and. self%level() == 0) then
       ! Every self%stateSwapCount steps, set γ=1 to allow interchange of chains.
       temperedDifferentialEvolutionStepSize=1.0d0
    else
       gammaBoostFactor=self%temperature()**self%posteriorSampleDffrntlEvltnPrpslSzTmpExp_%exponent(self%temperedStates,self%temperatures,self%posteriorSampleState_,self%posteriorSampleConvergence_)
       temperedDifferentialEvolutionStepSize=gammaBoostFactor *self%posteriorSampleSimulationDifferentialEvolution%stepSize(forceAcceptance)
    end if
    return
  end function temperedDifferentialEvolutionStepSize

  double precision function temperedDifferentialEvolutionTemperature(self)
    !!{
    Return the temperature.
    !!}
    implicit none
    class(posteriorSampleSimulationTemperedDffrntlEvltn), intent(inout) :: self

    if (self%level() == 0) then
       temperedDifferentialEvolutionTemperature=1.0d0
    else
       temperedDifferentialEvolutionTemperature=self%temperatures(self%level())
    end if
    return
  end function temperedDifferentialEvolutionTemperature

  logical function temperedDifferentialEvolutionAcceptProposal(self,logPosterior,logPosteriorProposed,logLikelihoodVariance,logLikelihoodVarianceProposed)
    !!{
    Return whether or not to accept a proposal.
    !!}
    implicit none
    class           (posteriorSampleSimulationTemperedDffrntlEvltn), intent(inout) :: self
    double precision                                               , intent(in   ) :: logPosterior         , logPosteriorProposed         , &
         &                                                                            logLikelihoodVariance, logLikelihoodVarianceProposed
    double precision                                                               :: x
    !$GLC attributes unused :: logLikelihoodVariance, logLikelihoodVarianceProposed

    ! Decide whether to take step.
    x=self%randomNumberGenerator_%uniformSample()
    temperedDifferentialEvolutionAcceptProposal=               &
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
  end function temperedDifferentialEvolutionAcceptProposal
