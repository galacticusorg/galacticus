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
  Implementation of a posterior sampling simulation class which implements the differential evolution algorithm.
  !!}

  use :: Model_Parameters                           , only : modelParameterList
  use :: Models_Likelihoods                         , only : posteriorSampleLikelihoodClass
  use :: Numerical_Random_Numbers                   , only : randomNumberGeneratorClass
  use :: Posterior_Sample_Differential_Proposal_Size, only : posteriorSampleDffrntlEvltnProposalSizeClass
  use :: Posterior_Sample_Differential_Random_Jump  , only : posteriorSampleDffrntlEvltnRandomJumpClass
  use :: Posterior_Sampling_Convergence             , only : posteriorSampleConvergenceClass
  use :: Posterior_Sampling_State                   , only : posteriorSampleStateClass
  use :: Posterior_Sampling_State_Initialize        , only : posteriorSampleStateInitializeClass
  use :: Posterior_Sampling_Stopping_Criteria       , only : posteriorSampleStoppingCriterionClass

  !![
  <posteriorSampleSimulation name="posteriorSampleSimulationDifferentialEvolution">
   <description>
    This class uses the differential evolution algorithm of \cite{terr_braak_markov_2006}. Multiple, parallel chains are run and
    proposals are constructed by selecting two chains at random, taking a fraction, $\gamma$, of the vector connecting the two chain
    states and adding this to the state of the current chain. The details of the algorithm are controlled by the following parameters:
    \begin{description}
    \item[{\normalfont \ttfamily [stepsMaximum]}] The maximum number of steps to take.
    \item[{\normalfont \ttfamily [acceptanceAverageCount]}] The number of steps over which to average the acceptance rate.
    \item[{\normalfont \ttfamily [stateSwapCount]}] The number of steps after which to set $\gamma=1$ to allow chains to swap states.
    \item[{\normalfont \ttfamily [logFileRoot]}] The full path and root name of a file to log results to. The actual file name will
      have the rank of the \gls{mpi} process appended to it.
    \item[{\normalfont \ttfamily [sampleOutliers]}] If set to {\normalfont \ttfamily false} then proposals for non-outlier chains
      post-convergence are constructed only from other non-outlier chains. Otherwise, proposals for non-outlier chains
      post-convergence are constructed from all other chains.
    \end{description}
   </description>
   <descriptorSpecial>descriptorSpecial</descriptorSpecial>
  </posteriorSampleSimulation>
  !!]
  type, extends(posteriorSampleSimulationClass) :: posteriorSampleSimulationDifferentialEvolution
     !!{
     Implementation of a posterior sampling simulation class which implements the differential evolution algorithm.
     !!}
     private
     integer                                                                                   :: parameterCount                                    , stepsMaximum            , &
          &                                                                                       stateSwapCount                                    , acceptanceAverageCount  , &
          &                                                                                       logFlushCount                                     , reportCount             , &
          &                                                                                       recomputeCount                                    , slowStepCount
     double precision                                                                          :: logPosterior                                      , logPrior
     logical                                                                                   :: isConverged                                       , sampleOutliers          , &
          &                                                                                       isInteractive                                     , appendLogs              , &
          &                                                                                       loadBalance                                       , ignoreChainNumberAdvice
     logical                                                       , allocatable, dimension(:) :: modelParametersActiveIsSlow
     type            (modelParameterList                          ), allocatable, dimension(:) :: modelParametersActive_                            , modelParametersInactive_
     class           (posteriorSampleLikelihoodClass              ), pointer                   :: posteriorSampleLikelihood_               => null()
     class           (posteriorSampleConvergenceClass             ), pointer                   :: posteriorSampleConvergence_              => null()
     class           (posteriorSampleStoppingCriterionClass       ), pointer                   :: posteriorSampleStoppingCriterion_        => null()
     class           (posteriorSampleStateClass                   ), pointer                   :: posteriorSampleState_                    => null()
     class           (posteriorSampleStateInitializeClass         ), pointer                   :: posteriorSampleStateInitialize_          => null()
     class           (posteriorSampleDffrntlEvltnProposalSizeClass), pointer                   :: posteriorSampleDffrntlEvltnProposalSize_ => null()
     class           (posteriorSampleDffrntlEvltnRandomJumpClass  ), pointer                   :: posteriorSampleDffrntlEvltnRandomJump_   => null()
     class           (randomNumberGeneratorClass                  ), pointer                   :: randomNumberGenerator_                   => null()
     type            (varying_string                              )                            :: logFileRoot                                       , interactionRoot
   contains
     !![
     <methods>
       <method method="logging"           description="Return true if the simulator is currently logging state."                                           />
       <method method="posterior"         description="Return the log of posterior probability for the given {\normalfont \ttfamily posteriorSampleState}."/>
       <method method="update"            description="Update the simulator to the new {\normalfont \ttfamily stateVector} after a step."                  />
       <method method="temperature"       description="Return the current temperature."                                                                    />
       <method method="acceptProposal"    description="Return true if the proposed state should be accepted."                                              />
       <method method="stepSize"          description="Return the step size parameter, $\gamma$, for the differential evolution proposal vector."          />
       <method method="chainSelect"       description="Select a chain."                                                                                    />
       <method method="descriptorSpecial" description="Handle adding special parameters to the descriptor."                                                />
     </methods>
     !!]
     final     ::                      differentialEvolutionDestructor
     procedure :: simulate          => differentialEvolutionSimulate
     procedure :: logging           => differentialEvolutionLogging
     procedure :: posterior         => differentialEvolutionPosterior
     procedure :: update            => differentialEvolutionUpdate
     procedure :: stepSize          => differentialEvolutionStepSize
     procedure :: acceptProposal    => differentialEvolutionAcceptProposal
     procedure :: temperature       => differentialEvolutionTemperature
     procedure :: chainSelect       => differentialEvolutionChainSelect
     procedure :: descriptorSpecial => differentialEvolutionDescriptorSpecial
  end type posteriorSampleSimulationDifferentialEvolution

  interface posteriorSampleSimulationDifferentialEvolution
     !!{
     Constructors for the \refClass{posteriorSampleSimulationDifferentialEvolution} posterior sampling convergence class.
     !!}
     module procedure differentialEvolutionConstructorParameters
     module procedure differentialEvolutionConstructorInternal
  end interface posteriorSampleSimulationDifferentialEvolution

contains

  function differentialEvolutionConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{posteriorSampleSimulationDifferentialEvolution} posterior sampling simulation class which builds the object from a
    parameter set.
    !!}
    use :: Display         , only : displayMessage      , displayVerbosity      , verbosityLevelInfo
    use :: Error           , only : Error_Report
    use :: Input_Parameters, only : inputParameter      , inputParameters
    use :: MPI_Utilities   , only : mpiSelf
    use :: Model_Parameters, only : modelParameterActive, modelParameterInactive
    use :: String_Handling , only : operator(//)
    implicit none
    type   (posteriorSampleSimulationDifferentialEvolution)                              :: self
    type   (inputParameters                               ), intent(inout)               :: parameters
    type   (modelParameterList                            ), pointer      , dimension(:) :: modelParametersActive_                  , modelParametersInactive_
    class  (modelParameterClass                           ), pointer                     :: modelParameter_
    class  (posteriorSampleLikelihoodClass                ), pointer                     :: posteriorSampleLikelihood_
    class  (posteriorSampleConvergenceClass               ), pointer                     :: posteriorSampleConvergence_
    class  (posteriorSampleStoppingCriterionClass         ), pointer                     :: posteriorSampleStoppingCriterion_
    class  (posteriorSampleStateClass                     ), pointer                     :: posteriorSampleState_
    class  (posteriorSampleStateInitializeClass           ), pointer                     :: posteriorSampleStateInitialize_
    class  (posteriorSampleDffrntlEvltnProposalSizeClass  ), pointer                     :: posteriorSampleDffrntlEvltnProposalSize_
    class  (posteriorSampleDffrntlEvltnRandomJumpClass    ), pointer                     :: posteriorSampleDffrntlEvltnRandomJump_
    class  (randomNumberGeneratorClass                    ), pointer                     :: randomNumberGenerator_
    integer                                                                              :: stepsMaximum                            , acceptanceAverageCount  , &
         &                                                                                  stateSwapCount                          , logFlushCount           , &
         &                                                                                  reportCount                             , activeParameterCount    , &
         &                                                                                  inactiveParameterCount                  , i                       , &
         &                                                                                  iActive                                 , iInactive               , &
         &                                                                                  recomputeCount                          , slowStepCount
    type   (varying_string                                )                              :: logFileRoot                             , interactionRoot         , &
         &                                                                                  message
    logical                                                                              :: sampleOutliers                          , appendLogs              , &
         &                                                                                  loadBalance                             , ignoreChainNumberAdvice

    !![
    <inputParameter>
      <name>stepsMaximum</name>
      <defaultValue>huge(0)</defaultValue>
      <description>The maximum number of steps to take.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>acceptanceAverageCount</name>
      <defaultValue>10</defaultValue>
      <description>The number of steps over which to average the acceptance rate.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>stateSwapCount</name>
      <defaultValue>10</defaultValue>
      <description>The number of steps between state swap steps.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>slowStepCount</name>
      <defaultValue>1</defaultValue>
      <description>The number of steps between slow parameter update steps.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>recomputeCount</name>
      <defaultValue>-1</defaultValue>
      <description>The number of steps between forced recomputations of the likelihood.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>logFlushCount</name>
      <defaultValue>10</defaultValue>
      <description>The number of steps between flushing the log file.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>reportCount</name>
      <defaultValue>10</defaultValue>
      <description>The number of steps between issuing reports.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>sampleOutliers</name>
      <defaultValue>.true.</defaultValue>
      <description>If true, sample from outlier states.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>interactionRoot</name>
      <defaultValue>var_str('none')</defaultValue>
      <description>Root file name for interaction files, or `{\normalfont \ttfamily none}' if interaction is not required.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>logFileRoot</name>
      <description>Root file name for log files.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>appendLogs</name>
      <description>If true, do not overwrite existing log files, but instead append to them.</description>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
    </inputParameter>
    <inputParameter>
      <name>loadBalance</name>
      <description>If true, attempt to balance the workload across different compute nodes.</description>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
    </inputParameter>
    <inputParameter>
      <name>ignoreChainNumberAdvice</name>
      <description>If true, ignore warnings and errors about not being able to span the full parameter space with the number of chains used.</description>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
    </inputParameter>
    <objectBuilder class="posteriorSampleLikelihood"               name="posteriorSampleLikelihood_"               source="parameters"/>
    <objectBuilder class="posteriorSampleConvergence"              name="posteriorSampleConvergence_"              source="parameters"/>
    <objectBuilder class="posteriorSampleStoppingCriterion"        name="posteriorSampleStoppingCriterion_"        source="parameters"/>
    <objectBuilder class="posteriorSampleState"                    name="posteriorSampleState_"                    source="parameters"/>
    <objectBuilder class="posteriorSampleStateInitialize"          name="posteriorSampleStateInitialize_"          source="parameters"/>
    <objectBuilder class="posteriorSampleDffrntlEvltnProposalSize" name="posteriorSampleDffrntlEvltnProposalSize_" source="parameters"/>
    <objectBuilder class="posteriorSampleDffrntlEvltnRandomJump"   name="posteriorSampleDffrntlEvltnRandomJump_"   source="parameters"/>
    <objectBuilder class="randomNumberGenerator"                   name="randomNumberGenerator_"                   source="parameters"/>
    !!]
    ! Determine the number of parameters.
    activeParameterCount  =0
    inactiveParameterCount=0
    do i=1,parameters%copiesCount("modelParameter")
       !![
       <objectBuilder class="modelParameter" name="modelParameter_" source="parameters" copy="i" />
       !!]
       select type (modelParameter_)
       class is (modelParameterActive  )
          activeParameterCount  =activeParameterCount  +1
       class is (modelParameterInactive)
          inactiveParameterCount=inactiveParameterCount+1
       end select
       !![
       <objectDestructor name="modelParameter_"/>
       !!]
    end do
    if (activeParameterCount < 1) call Error_Report('at least one active parameter must be specified in config file'//{introspection:location})
    if (mpiSelf%isMaster() .and. displayVerbosity() >= verbosityLevelInfo) then
       message='Found '
       message=message//activeParameterCount//' active parameters (and '//inactiveParameterCount//' inactive parameters)'
       call displayMessage(message)
    end if
    ! Initialize priors and random perturbers.
    allocate(modelParametersActive_  (  activeParameterCount))
    allocate(modelParametersInactive_(inactiveParameterCount))
    iActive  =0
    iInactive=0
    do i=1,parameters%copiesCount("modelParameter")
       !![
       <objectBuilder class="modelParameter" name="modelParameter_" source="parameters" copy="i" />
       !!]
       select type (modelParameter_)
       class is (modelParameterInactive)
          iInactive=iInactive+1
          modelParametersInactive_(iInactive)%modelParameter_ => modelParameter_
          !![
          <referenceCountIncrement owner="modelParametersInactive_(iInactive)" object="modelParameter_"/>
          !!]
       class is (modelParameterActive  )
          iActive  =iActive  +1
          modelParametersActive_  (  iActive)%modelParameter_ => modelParameter_
          !![
          <referenceCountIncrement owner="modelParametersActive_  (iActive  )" object="modelParameter_"/>
          !!]
       end select
       !![
       <objectDestructor name="modelParameter_"/>
       !!]
    end do
    self=posteriorSampleSimulationDifferentialEvolution(modelParametersActive_,modelParametersInactive_,posteriorSampleLikelihood_,posteriorSampleConvergence_,posteriorSampleStoppingCriterion_,posteriorSampleState_,posteriorSampleStateInitialize_,posteriorSampleDffrntlEvltnProposalSize_,posteriorSampleDffrntlEvltnRandomJump_,randomNumberGenerator_,stepsMaximum,acceptanceAverageCount,stateSwapCount,slowStepCount,recomputeCount,char(logFileRoot),sampleOutliers,logFlushCount,reportCount,char(interactionRoot),appendLogs,loadBalance,ignoreChainNumberAdvice)
    !![
    <inputParametersValidate source="parameters" multiParameters="modelParameter"/>
    <objectDestructor name="posteriorSampleLikelihood_"              />
    <objectDestructor name="posteriorSampleConvergence_"             />
    <objectDestructor name="posteriorSampleStoppingCriterion_"       />
    <objectDestructor name="posteriorSampleState_"                   />
    <objectDestructor name="posteriorSampleStateInitialize_"         />
    <objectDestructor name="posteriorSampleDffrntlEvltnProposalSize_"/>
    <objectDestructor name="posteriorSampleDffrntlEvltnRandomJump_"  />
    <objectDestructor name="randomNumberGenerator_"                  />
    !!]
    do i=1,  activeParameterCount
       !![
       <objectDestructor name="modelParametersActive_  (i)%modelParameter_"/>
       !!]
    end do
    do i=1,inactiveParameterCount
       !![
       <objectDestructor name="modelParametersInactive_(i)%modelParameter_"/>
       !!]
    end do
    deallocate(modelParametersActive_  )
    deallocate(modelParametersInactive_)
    return
  end function differentialEvolutionConstructorParameters
  
  function differentialEvolutionConstructorInternal(modelParametersActive_,modelParametersInactive_,posteriorSampleLikelihood_,posteriorSampleConvergence_,posteriorSampleStoppingCriterion_,posteriorSampleState_,posteriorSampleStateInitialize_,posteriorSampleDffrntlEvltnProposalSize_,posteriorSampleDffrntlEvltnRandomJump_,randomNumberGenerator_,stepsMaximum,acceptanceAverageCount,stateSwapCount,slowStepCount,recomputeCount,logFileRoot,sampleOutliers,logFlushCount,reportCount,interactionRoot,appendLogs,loadBalance,ignoreChainNumberAdvice) result(self)
    !!{
    Internal constructor for the ``differentialEvolution'' simulation class.
    !!}
    use :: Model_Parameters, only : modelParameterActive
    implicit none
    type     (posteriorSampleSimulationDifferentialEvolution)                              :: self
    type     (modelParameterList                            ), intent(in   ), dimension(:) :: modelParametersActive_                  , modelParametersInactive_
    class    (posteriorSampleLikelihoodClass                ), intent(in   ), target       :: posteriorSampleLikelihood_
    class    (posteriorSampleConvergenceClass               ), intent(in   ), target       :: posteriorSampleConvergence_
    class    (posteriorSampleStoppingCriterionClass         ), intent(in   ), target       :: posteriorSampleStoppingCriterion_
    class    (posteriorSampleStateClass                     ), intent(in   ), target       :: posteriorSampleState_
    class    (posteriorSampleStateInitializeClass           ), intent(in   ), target       :: posteriorSampleStateInitialize_
    class    (posteriorSampleDffrntlEvltnProposalSizeClass  ), intent(in   ), target       :: posteriorSampleDffrntlEvltnProposalSize_
    class    (posteriorSampleDffrntlEvltnRandomJumpClass    ), intent(in   ), target       :: posteriorSampleDffrntlEvltnRandomJump_
    class    (randomNumberGeneratorClass                    ), intent(in   ), target       :: randomNumberGenerator_
    integer                                                  , intent(in   )               :: stepsMaximum                            , acceptanceAverageCount  , &
         &                                                                                    stateSwapCount                          , logFlushCount           , &
         &                                                                                    reportCount                             , recomputeCount          , &
         &                                                                                    slowStepCount
    character(len=*                                         ), intent(in   )               :: logFileRoot                             , interactionRoot
    logical                                                  , intent(in   )               :: sampleOutliers                          , appendLogs              , &
         &                                                                                    loadBalance                             , ignoreChainNumberAdvice
    integer                                                                                :: i
    !![
    <constructorAssign variables="*posteriorSampleLikelihood_, *posteriorSampleConvergence_, *posteriorSampleStoppingCriterion_, *posteriorSampleState_, *posteriorSampleStateInitialize_, *posteriorSampleDffrntlEvltnProposalSize_, *posteriorSampleDffrntlEvltnRandomJump_, *randomNumberGenerator_, stepsMaximum, acceptanceAverageCount, stateSwapCount, slowStepCount, recomputeCount, logFlushCount, reportCount, sampleOutliers, logFileRoot, interactionRoot, appendLogs, loadBalance, ignoreChainNumberAdvice"/>
    !!]

    allocate(self%modelParametersActive_     (size(modelParametersActive_  )))
    allocate(self%modelParametersActiveIsSlow(size(modelParametersActive_  )))
    allocate(self%modelParametersInactive_   (size(modelParametersInactive_)))
    do i=1,size(modelParametersActive_  )
       self%modelParametersActive_  (i)                 =  modelParameterList      ( )
       self%modelParametersActive_  (i)%modelParameter_ => modelParametersActive_  (i)%modelParameter_
       !![
       <referenceCountIncrement owner="self%modelParametersActive_  (i)" object="modelParameter_"/>
       !!]
       select type (modelParameter_ => modelParametersActive_(i)%modelParameter_)
       class is (modelParameterActive)
          self%modelParametersActiveIsSlow(i)=modelParameter_%isSlow()
       end select
    end do
    do i=1,size(modelParametersInactive_)
       self%modelParametersInactive_(i)                 =  modelParameterList      ( )
       self%modelParametersInactive_(i)%modelParameter_ => modelParametersInactive_(i)%modelParameter_
       !![
       <referenceCountIncrement owner="self%modelParametersInactive_(i)" object="modelParameter_"/>
       !!]
    end do
    self%parameterCount=size(modelParametersActive_)
    self%isInteractive =trim(interactionRoot) /= "none"
    call self%posteriorSampleState_%parameterCountSet(self%parameterCount)
    return
  end function differentialEvolutionConstructorInternal

  subroutine differentialEvolutionDestructor(self)
    !!{
    Destroy a differential evolution simulation object.
    !!}
    implicit none
    type   (posteriorSampleSimulationDifferentialEvolution), intent(inout) :: self
    integer                                                                :: i

    !![
    <objectDestructor name="self%posteriorSampleLikelihood_"              />
    <objectDestructor name="self%posteriorSampleConvergence_"             />
    <objectDestructor name="self%posteriorSampleStoppingCriterion_"       />
    <objectDestructor name="self%posteriorSampleState_"                   />
    <objectDestructor name="self%posteriorSampleStateInitialize_"         />
    <objectDestructor name="self%posteriorSampleDffrntlEvltnProposalSize_"/>
    <objectDestructor name="self%posteriorSampleDffrntlEvltnRandomJump_"  />
    <objectDestructor name="self%randomNumberGenerator_"                  />
    !!]
    if (allocated(self%modelParametersActive_  )) then
       do i=1,size(self%modelParametersActive_  )
          !![
	  <objectDestructor name="self%modelParametersActive_  (i)%modelParameter_"/>
          !!]
       end do
       deallocate(self%modelParametersActive_  )
    end if
    if (allocated(self%modelParametersInactive_)) then
       do i=1,size(self%modelParametersInactive_)
          !![
	  <objectDestructor name="self%modelParametersInactive_(i)%modelParameter_"/>
          !!]
       end do
       deallocate(self%modelParametersInactive_)
     end if
     return
  end subroutine differentialEvolutionDestructor

  subroutine differentialEvolutionSimulate(self)
    !!{
    Perform a differential evolution simulation.
    !!}
    use :: Display                     , only : displayIndent             , displayMessage , displayUnindent, displayMagenta, &
         &                                      displayReset
    use :: File_Utilities              , only : File_Exists               , File_Remove
    use :: Error                       , only : Error_Report              , Warn
    use :: MPI_Utilities               , only : mpiBarrier                , mpiSelf
    use :: Models_Likelihoods_Constants, only : logImpossible
    use :: Posterior_Sampling_State    , only : posteriorSampleStateSimple
    use :: String_Handling             , only : operator(//)
    implicit none
    class           (posteriorSampleSimulationDifferentialEvolution), intent(inout)                    :: self
    integer                                                         , dimension(                    2) :: chainPair
    double precision                                                , dimension(self%parameterCount,2) :: statePair
    double precision                                                , dimension(self%parameterCount  ) :: stateVector           , stateVectorProposed          , &
         &                                                                                                stateVectorInteractive, stateVectorPerturbation
    class           (posteriorSampleStateClass                     ), allocatable                      :: stateProposed
    real                                                                                               :: timePreEvaluate       , timePostEvaluate             , &
         &                                                                                                timeEvaluate          , timeEvaluatePrevious
    double precision                                                                                   :: logLikelihoodProposed , logPosteriorProposed         , &
         &                                                                                                logLikelihood         , logLikelihoodVariance        , &
         &                                                                                                timeEvaluateInitial   , logLikelihoodVarianceProposed, &
         &                                                                                                logPosteriorInitial   , logLikelihoodInitial
    type            (varying_string                                )                                   :: logFileName           , message                      , &
         &                                                                                                interactionFileName
    integer                                                                                            :: logFileUnit           , convergedAtStep              , &
         &                                                                                                convergenceFileUnit   , i                            , &
         &                                                                                                ioStatus              , interactionFile
    logical                                                                                               forceAcceptance
    character       (len=32                                        )                                   :: label

    ! Check that we have sufficient chains for differential evolution.
    if (.not.self%ignoreChainNumberAdvice) then
       if      (mpiSelf%count() <    size(self%modelParametersActive_)) then
          call Error_Report('the number of chains should at least equal the number of active parameters, otherwise it may not be possible to sample the full posterior distribution'               //{introspection:location})
       else if (mpiSelf%count() <  2*size(self%modelParametersActive_)) then
          call Warn        ('the number of chains should be at least twice the number of active parameters, otherwise it may not be efficient to sample the full posterior distribution'           //{introspection:location})
       else if (mpiSelf%count() < 10*size(self%modelParametersActive_)) then
          call Warn        ('for non-unimodal/otherwise complicated posteriors, Ter Braak (2006) recommends using a number of chains 10 to 20 times the dimension of the posterior parameter space'//{introspection:location})
       end if
    end if
    ! Check that the random number generator is independent across MPI processes.
    if (.not.self%randomNumberGenerator_%mpiIndependent()) call Error_Report('random number generator produces same sequence on all MPI processes'//{introspection:location})
    ! Write start-up message.
    message="Process "//mpiSelf%rankLabel()//" [PID: "
    message=message//getPID()//"] is running on host '"//mpiSelf%hostAffinity()//"'"
    call displayMessage(message)
    ! Allocate a simple state object for the proposed state.
    allocate(posteriorSampleStateSimple :: stateProposed)
    select type (stateProposed)
    type is (posteriorSampleStateSimple)
       stateProposed=posteriorSampleStateSimple(1)
       call stateProposed%parameterCountSet(self%parameterCount)
    end select
    ! Initialize chain to some state vector.
    call self%posteriorSampleStateInitialize_%initialize(self%posteriorSampleState_,self%modelParametersActive_,self%posteriorSampleLikelihood_,timeEvaluateInitial,logLikelihoodInitial,logPosteriorInitial)
    ! Evaluate the posterior in the initial state if it was not set.
    forceAcceptance     =.false.
    timeEvaluate        =-1.0
    timeEvaluatePrevious=real(timeEvaluateInitial)
    call CPU_Time(timePreEvaluate )
    call self%posterior(self%posteriorSampleState_,logImpossible,logImpossible,self%logPosterior,logLikelihood,logLikelihoodVariance,timeEvaluate,timeEvaluatePrevious,forceAcceptance)    
    call CPU_Time(timePostEvaluate)
    if (timeEvaluate < 0.0) timeEvaluate=timePostEvaluate-timePreEvaluate
    timeEvaluatePrevious=timeEvaluate
    ! Check for impossible state.
    if (self%logPosterior <= logImpossible) call Error_Report('impossible initial state'//{introspection:location})
    ! If the initializer returned a non-impossible likelihood, use it instead.
    if (logLikelihoodInitial > logImpossible) then
       logLikelihood    =logLikelihoodInitial
       self%logPosterior=logPosteriorInitial
    end if
    ! Begin stepping.
    logFileName=self%logFileRoot//'_'//mpiSelf%rankLabel()//'.log'
    if (self%appendLogs) then
       open(newunit=logFileUnit,file=char(logFileName),status='unknown',form='formatted',position='append')
    else
       open(newunit=logFileUnit,file=char(logFileName),status='unknown',form='formatted'                  )
       write (logFileUnit,'(a)') '# Simulation state chain file'
       write (logFileUnit,'(a)') '# Columns:'
       write (logFileUnit,'(a)') '#    1 = Simulation step'
       write (logFileUnit,'(a)') '#    2 = Chain index'
       write (logFileUnit,'(a)') '#    3 = Evaluation time (s)'
       write (logFileUnit,'(a)') '#    4 = Chain is converged? [T/F]'
       write (logFileUnit,'(a)') '#    5 = log posterior'
       write (logFileUnit,'(a)') '#    6 = log likelihood'
       do i=1,size(self%modelParametersActive_)
          write (logFileUnit,'(a,i3,a,a,a)') '#  ',i+6,' = Parameter `',char(self%modelParametersActive_(i)%modelParameter_%name()),'`'
       end do
    end if
    self%isConverged=.false.
    do while (                                                                                                   &
         &          self%posteriorSampleState_            %count(                          ) < self%stepsMaximum &
         &    .and.                                                                                              &
         &     .not.self%posteriorSampleStoppingCriterion_%stop (self%posteriorSampleState_)                     &
         &   )
       ! Initialize likelihood acceptance condition.
       forceAcceptance=.false.
       ! Pick two random processes to use for proposal.
       chainPair(1)=self%chainSelect(              )
       chainPair(2)=self%chainSelect([chainPair(1)])
       ! Receive states from selected chains.
       stateVector=self%posteriorSampleState_%get()
       statePair  =mpiSelf%requestData(chainPair,stateVector)
       ! Generate proposal.
       stateVectorProposed= stateVector                    &
            &              +self%stepSize(forceAcceptance) &
            &              *(                              &
            &                +statePair(:,1)               &
            &                -statePair(:,2)               &
            &               )
       ! Add random perturbations to the proposal.
       if (.not.forceAcceptance) then
          stateVectorPerturbation=self%posteriorSampleDffrntlEvltnRandomJump_%sample(self%modelParametersActive_,self%posteriorSampleState_)
          ! Disable perturbations in slow parameters, except for on slow steps.
          if (mod(self%posteriorSampleState_%count(),self%slowStepCount) /= 0) then
             where (self%modelParametersActiveIsSlow)
                stateVectorPerturbation=0.0d0
             end where
          end if
          stateVectorProposed=+stateVectorProposed     &
               &              +stateVectorPerturbation
       end if
       ! If simulation is interactive, check for any interaction file.
       if (self%isInteractive) then
          ! Check if an interaction file exists.
          interactionFileName=self%interactionRoot//"_"//mpiSelf%rankLabel()
          if (File_Exists(interactionFileName)) then
             ! Read the file and validate.
             open(newUnit=interactionFile,file=char(interactionFileName),status='old',form='formatted',ioStat=ioStatus)
             if (ioStatus == 0) then
                read (interactionFile,*,ioStat=ioStatus) stateVectorInteractive
                if (ioStatus == 0) then
                   ! Copy the state to the proposed state vector.
                   stateVectorProposed=stateVectorInteractive
                   message="Chain "//mpiSelf%rankLabel()//" is being interactively moved to state:"
                   call displayIndent(message)
                   ! Map parameters of interactively proposed state.
                   do i=1,size(stateVector)
                      write (label,*) stateVectorProposed(i)
                      message="State["
                      message=message//i//"] = "//trim(adjustl(label))
                      call displayMessage(message)
                      stateVectorProposed(i)=self%modelParametersActive_(i)%modelParameter_%map(stateVectorProposed(i))
                   end do
                   call displayUnindent('end')
                   ! Force acceptance of this state.
                   forceAcceptance=.true.
                else
                   message=displayMagenta()//"WARNING:"//displayReset()//" state proposed in interaction file '"//interactionFileName//"' cannot be read"
                   call displayMessage(message)
                end if
             else
                message=displayMagenta()//"WARNING:"//displayReset()//" unable to open interaction file '"//interactionFileName//"'"
                call displayMessage(message)
             end if
             close(interactionFile)
             ! Remove the interaction file.
             call File_Remove(interactionFileName)
          end if
       end if
       ! Store the proposed state vector.
       call stateProposed%update(stateVectorProposed,.true.,.false.)
       select type (stateProposed)
       type is (posteriorSampleStateSimple)
          call stateProposed%countSet(self%posteriorSampleState_%count())
       end select
       ! Evaluate likelihood.
       timeEvaluatePrevious=timeEvaluate
       timeEvaluate        =-1.0
       call CPU_Time(timePreEvaluate )
       call self%posterior(stateProposed,logLikelihood,self%logPosterior-logLikelihood,logPosteriorProposed,logLikelihoodProposed,logLikelihoodVarianceProposed,timeEvaluate,timeEvaluatePrevious,forceAcceptance)
       call CPU_Time(timePostEvaluate)
       if (timeEvaluate < 0.0) timeEvaluate=timePostEvaluate-timePreEvaluate
       ! Decide whether to take step.
       if     (                                                                                                                 &
            &   self%acceptProposal(self%logPosterior,logPosteriorProposed,logLikelihoodVariance,logLikelihoodVarianceProposed) &
            &  .or.                                                                                                             &
            &   forceAcceptance                                                                                                 &
            & ) then
          self%logPosterior    =logPosteriorProposed
          logLikelihood        =logLikelihoodProposed
          logLikelihoodVariance=logLikelihoodVarianceProposed
          stateVector          =stateVectorProposed
       else
          ! Step not accepted - retain the old estimate of work time.
          timeEvaluate         =timeEvaluatePrevious
       end if
       call self%update(stateVector)
       ! Unmap parameters and write to log file.
       do i=1,size(stateVector)
          stateVector(i)=self%modelParametersActive_(i)%modelParameter_%unmap(stateVector(i))
       end do
       if (self%logging()) then
          write (logFileUnit,*) self   %posteriorSampleState_%count        (), &
               &                mpiSelf                      %rank         (), &
               &                                              timeEvaluate   , &
               &                self                         %isConverged    , &
               &                self                         %logPosterior   , &
               &                                              logLikelihood  , &
               &                                              stateVector
          if (mod(self%posteriorSampleState_%count(),self%logFlushCount) == 0) call flush(logFileUnit)
       end if
       ! Repeat.
       call mpiBarrier()
       ! Test for convergence.
       if (.not.self%isConverged) then
          self%isConverged=self%posteriorSampleConvergence_%isConverged(self%posteriorSampleState_,self%logPosterior)
          if (self%isConverged) then
             convergedAtStep=self%posteriorSampleState_%count()
             if (mpiSelf%rank() == 0) then
                message='Converged after '
                message=message//convergedAtStep//' steps'
                call displayMessage(message)
                logFileName=self%logFileRoot//'_'//mpiSelf%rankLabel()//'.convergence.log'
                open(newunit=convergenceFileUnit,file=char(logFileName),status='unknown',form='formatted',access='append')
                write (convergenceFileUnit,'(a,i8)') 'Converged at step: ',convergedAtStep
                call self%posteriorSampleConvergence_%logReport(convergenceFileUnit)
                close(convergenceFileUnit)
             end if
          end if
       end if
    end do
    close(logFileUnit)
    return
  end subroutine differentialEvolutionSimulate

  subroutine differentialEvolutionUpdate(self,stateVector)
    !!{
    Update the differential evolution simulator state.
    !!}
    use :: MPI_Utilities, only : mpiSelf
    implicit none
    class           (posteriorSampleSimulationDifferentialEvolution), intent(inout)                                 :: self
    double precision                                                , intent(in   ), dimension(self%parameterCount) :: stateVector
    logical                                                         , allocatable  , dimension(:                  ) :: outlierMask
    integer                                                                                                         :: i

    allocate(outlierMask(0:mpiSelf%count()-1))
    do i=0,mpiSelf%count()-1
       outlierMask(i)=self%posteriorSampleConvergence_%stateIsOutlier(i)
    end do
    call self%posteriorSampleState_%update(stateVector,self%logging(),self%posteriorSampleConvergence_%isConverged(),outlierMask)
    return
  end subroutine differentialEvolutionUpdate

  integer function differentialEvolutionChainSelect(self,blockedChains)
    !!{
    Select a chain at random, optionally excluding blocked chains.
    !!}
    use :: MPI_Utilities, only : mpiSelf
    use :: Error        , only : Error_Report
    implicit none
    class  (posteriorSampleSimulationDifferentialEvolution), intent(inout)                         :: self
    integer                                                , intent(in   ), dimension(:), optional :: blockedChains
    logical                                                                                        :: accept

    if (mpiSelf%count() < 2) call Error_Report('at least two chains are needed for differential evolution'//{introspection:location})
    accept=.false.
    do while (.not.accept)
       differentialEvolutionChainSelect=min(                                                  &
            &                               int(                                              &
            &                                   +dble(mpiSelf%count())                        &
            &                                   *self%randomNumberGenerator_%uniformSample()  &
            &                                  )                                            , &
            &                               +mpiSelf%count()                                  &
            &                               -1                                                &
            &                              )
       accept=.true.
       if (.not.self%sampleOutliers.and.self%isConverged)                                           &
            & accept=                                                                               &
            &              self%posteriorSampleConvergence_%stateIsOutlier(         mpiSelf%rank()) &
            &        .or.                                                                           &
            &         .not.self%posteriorSampleConvergence_%stateIsOutlier(differentialEvolutionChainSelect)
       if (present(blockedChains)) accept=accept.and.all(differentialEvolutionChainSelect /= blockedChains)
    end do
    return
  end function differentialEvolutionChainSelect

  logical function differentialEvolutionLogging(self)
    !!{
    Specifies whether or not the current state should be logged to file during differential evolution.
    !!}
    implicit none
    class(posteriorSampleSimulationDifferentialEvolution), intent(inout) :: self
    !$GLC attributes unused :: self

    differentialEvolutionLogging=.true.
    return
  end function differentialEvolutionLogging

  subroutine differentialEvolutionPosterior(self,posteriorSampleState_,logLikelihoodCurrent,logPriorCurrent,logPosterior,logLikelihood,logLikelihoodVariance,timeEvaluate,timeEvaluatePrevious,forceAcceptance)
    !!{
    Return the log of the posterior for the current state.
    !!}
    use            :: Display         , only : displayIndent             , displayMessage, displayUnindent
    use            :: Error           , only : Error_Report
    use, intrinsic :: ISO_C_Binding   , only : c_size_t
    use            :: Kind_Numbers    , only : kind_int4
    use            :: MPI_Utilities   , only : mpiBarrier                , mpiSelf
    use            :: Model_Parameters, only : modelParameterListLogPrior
    use            :: Sorting         , only : sortIndex
    implicit none
    class           (posteriorSampleSimulationDifferentialEvolution), intent(inout)               :: self
    class           (posteriorSampleStateClass                     ), intent(inout)               :: posteriorSampleState_
    double precision                                                , intent(in   )               :: logLikelihoodCurrent  , logPriorCurrent
    double precision                                                , intent(  out)               :: logLikelihood         , logPosterior            , &
         &                                                                                           logLikelihoodVariance
    real                                                            , intent(inout)               :: timeEvaluate
    real                                                            , intent(in   )               :: timeEvaluatePrevious
    logical                                                         , intent(inout)               :: forceAcceptance
    double precision                                                , dimension(:  ), allocatable :: timesEvaluate         , nodeWork                , &
         &                                                                                           stateVectorSelf       , timesEvaluateActual
    double precision                                                , dimension(:,:), allocatable :: stateVectorWork
    integer                                                         , dimension(:  ), allocatable :: processToProcess      , processFromProcess
    integer         (c_size_t                                      ), dimension(:  ), allocatable :: timesEvaluateOrder    , nodeWorkOrder
    integer                                                         , dimension(1  )              :: chainIndexSelf
    integer                                                         , dimension(1,1)              :: chainIndexWork
    double precision                                                , dimension(1,1)              :: logLikelihoodSelf     , logLikelihoodCurrentWork, &
         &                                                                                           logPriorCurrentWork   , logPriorWork            , &
         &                                                                                           timeEvaluateSelf
    logical                                                         , dimension(1,1)              :: forceAcceptanceSelf   , forceAcceptanceWorkIn
    logical                                                         , dimension(1  )              :: forceAcceptanceWork   , forceAcceptanceSelfIn
    double precision                                                , dimension(1  )              :: logLikelihoodWork     , logLikelihoodCurrentSelf, &
         &                                                                                           logPriorCurrentSelf   , logPriorSelf            , &
         &                                                                                           timeEvaluateWork
    double precision                                                                              :: logPrior              , timeEvaluateEffective
    integer                                                                                       :: i                     , processTrial            , &
         &                                                                                           nodeTrial
    type            (varying_string                                )                              :: message
    character       (len=10                                        )                              :: label

    ! Evaluate the proposed prior.
    logPrior=modelParameterListLogPrior(self%modelParametersActive_,posteriorSampleState_)
    ! Gather timing data from all chains.
    allocate(timesEvaluate      (0:mpiSelf%count()-1))
    allocate(timesEvaluateActual(0:mpiSelf%count()-1))
    allocate(processToProcess   (0:mpiSelf%count()-1))
    allocate(processFromProcess (0:mpiSelf%count()-1))
    allocate(timesEvaluateOrder (0:mpiSelf%count()-1))
    timeEvaluateEffective=timeEvaluatePrevious
    if     (                                                                           &
         &  .not.self%posteriorSampleLikelihood_%willEvaluate(                         &
         &                                         posteriorSampleState_             , &
         &                                         self%modelParametersActive_       , &
         &                                         self%posteriorSampleConvergence_  , &
         &                                         self%temperature                (), &
         &                                         logLikelihoodCurrent              , &
         &                                         logPriorCurrent                   , &
         &                                         logPrior                            &
         &                                        )                                    &
         & )                                                                           &
         & timeEvaluateEffective=0.0d0
    timesEvaluate=mpiSelf%gather(dble(timeEvaluateEffective))
    ! If previous time estimate is negative, don't do load balancing.
    if (mpiSelf%isMaster() .and. mod(self%posteriorSampleState_%count(),self%reportCount) == 0) call displayIndent('Load balancing report')
    if (any(timesEvaluate < 0.0d0) .or. .not.self%loadBalance) then
       forall(i=0:mpiSelf%count()-1)
          processToProcess  (i)=i
          processFromProcess(i)=i
       end forall
       if (self%loadBalance .and. mpiSelf%isMaster() .and. mod(self%posteriorSampleState_%count(),self%reportCount) == 0) call displayMessage('Not performing load balancing - missing work cost data')
    else
       ! Distribute tasks across nodes.
       timesEvaluateOrder=sortIndex(timesEvaluate)-1
       processToProcess=-1
       allocate(nodeWork     (mpiSelf%nodeCount()))
       allocate(nodeWorkOrder(mpiSelf%nodeCount()))
       nodeWork=0.0d0
       do i=mpiSelf%count()-1,0,-1
          nodeWorkOrder=sortIndex(nodeWork)
          do nodeTrial=1,mpiSelf%nodeCount()
             do processTrial=0,mpiSelf%count()-1
                if (mpiSelf%nodeAffinity(processTrial) == nodeWorkOrder(nodeTrial) .and. .not.any(processToProcess == processTrial)) then
                   processToProcess  (timesEvaluateOrder(i))=processTrial
                   processFromProcess(processTrial)=int(timesEvaluateOrder(i),kind_int4)
                   nodeWork(nodeWorkOrder(nodeTrial))=nodeWork(nodeWorkOrder(nodeTrial))+timesEvaluate(timesEvaluateOrder(i))
                   exit
                end if
             end do
             if (processToProcess(timesEvaluateOrder(i)) >= 0) exit
          end do
          if (processToProcess(timesEvaluateOrder(i)) < 0) call Error_Report('failed to assign task to process'//{introspection:location})
       end do
       ! Report.
       if (mpiSelf%isMaster() .and. mod(self%posteriorSampleState_%count(),self%reportCount) == 0) then
          call displayIndent('Chain redistribution:')
          do i=0,mpiSelf%count()-1
             write (label,'(i4.4)') i
             message='Chain '//trim(label)//' -> process/node '
             write (label,'(i4.4)') processToProcess(i)
             message=message//trim(label)//'/'
             write (label,'(i4.4)') mpiSelf%nodeAffinity(processToProcess(i))
             message=message//trim(label)//' (work = '
             write (label,'(f9.2)') timesEvaluate(i)
             message=message//trim(label)//')'
             call displayMessage(message)
          end do
          call displayUnindent('done')
          call displayIndent('Node work loads:')
          do i=1,size(nodeWork)
             write (label,'(i4.4)') i
             message='Node '//trim(label)//': work = '
             write (label,'(f9.2)') nodeWork(i)
             message=message//trim(label)
             call displayMessage(message)
          end do
          call displayUnindent('done')
       end if
    end if
    if (mpiSelf%isMaster() .and. mod(self%posteriorSampleState_%count(),self%reportCount) == 0) call displayUnindent('done')
    ! Get state vector, chain index, current likelihood, current prior and proposed prior.
    allocate(stateVectorSelf(self%parameterCount  ))
    allocate(stateVectorWork(self%parameterCount,1))
    stateVectorSelf         =posteriorSampleState_%get       ()
    chainIndexSelf          =posteriorSampleState_%chainIndex()
    logLikelihoodCurrentSelf=logLikelihoodCurrent
    logPriorCurrentSelf     =logPriorCurrent
    logPriorSelf            =logPrior
    forceAcceptanceSelfIn   =forceAcceptance
    stateVectorWork         =mpiSelf%requestData(processFromProcess(mpiSelf%rank():mpiSelf%rank()),stateVectorSelf         )
    chainIndexWork          =mpiSelf%requestData(processFromProcess(mpiSelf%rank():mpiSelf%rank()),chainIndexSelf          )
    logLikelihoodCurrentWork=mpiSelf%requestData(processFromProcess(mpiSelf%rank():mpiSelf%rank()),logLikelihoodCurrentSelf)
    logPriorCurrentWork     =mpiSelf%requestData(processFromProcess(mpiSelf%rank():mpiSelf%rank()),logPriorCurrentSelf     )
    logPriorWork            =mpiSelf%requestData(processFromProcess(mpiSelf%rank():mpiSelf%rank()),logPriorSelf            )
    forceAcceptanceWorkIn   =mpiSelf%requestData(processFromProcess(mpiSelf%rank():mpiSelf%rank()),forceAcceptanceSelfIn   )
    ! Set state and chain index.
    call posteriorSampleState_%update       (stateVectorWork(:,1),logState=.false.,isConverged=.false.)
    call posteriorSampleState_%chainIndexSet( chainIndexWork(1,1)                                     )
    ! Evaluate the likelihood.
    logLikelihood=self%posteriorSampleLikelihood_%evaluate(posteriorSampleState_,self%modelParametersActive_,self%modelParametersInactive_,self%posteriorSampleConvergence_,self%temperature(),logLikelihoodCurrentWork(1,1),logPriorCurrentWork(1,1),logPriorWork(1,1),timeEvaluate,logLikelihoodVariance=logLikelihoodVariance,forceAcceptance=forceAcceptanceWorkIn(1,1))
    call mpiBarrier()
    ! Distribute likelihoods back to origins.
    logLikelihoodWork    =                                                                         logLikelihood
    logLikelihoodSelf    =mpiSelf%requestData(processToProcess(mpiSelf%rank():mpiSelf%rank()),     logLikelihoodWork        )
    logLikelihood        =                                                                         logLikelihoodSelf  (1,1)
    ! Distribute likelihood variances back to origins.
    logLikelihoodWork    =                                                                         logLikelihoodVariance
    logLikelihoodSelf    =mpiSelf%requestData(processToProcess(mpiSelf%rank():mpiSelf%rank()),     logLikelihoodWork        )
    logLikelihoodVariance=                                                                         logLikelihoodSelf  (1,1)
    ! Distribute evaluation times back to origins.
    timeEvaluateWork     =                                                                    dble(timeEvaluate            )
    timeEvaluateSelf     =mpiSelf%requestData(processToProcess(mpiSelf%rank():mpiSelf%rank()),     timeEvaluateWork         )
    timeEvaluate         =                                                                    real(timeEvaluateSelf   (1,1))
    ! Distribute force acceptances back to origins.
    forceAcceptanceWork  =                                                                         forceAcceptanceWorkIn(1,1)
    forceAcceptanceSelf  =mpiSelf%requestData(processToProcess(mpiSelf%rank():mpiSelf%rank()),     forceAcceptanceWork      )
    forceAcceptance      =                                                                         forceAcceptanceSelf(1,1)
    ! Restore state and chain index.
    call posteriorSampleState_%update       (stateVectorSelf   ,logState=.false.,isConverged=.false.)
    call posteriorSampleState_%chainIndexSet( chainIndexSelf(1)                                     )
    ! Compute the log posterior.
    logPosterior=logPrior+logLikelihood
    ! Gather actual evaluation times and report.
    timesEvaluateActual=mpiSelf%gather(dble(timeEvaluate))
    if (mpiSelf%isMaster() .and. mod(self%posteriorSampleState_%count(),self%reportCount) == 0) then
       call displayIndent('Node work done vs. expected:')
       do i=0,mpiSelf%count()-1
          write (label,'(i4.4)') i
          message='Node '//trim(label)//': work (actual/estimated) = '
          write (label,'(f9.2)') timesEvaluateActual(i)
          message=message//trim(label)//" / "
          write (label,'(f9.2)') timesEvaluate      (i)
          message=message//trim(label)
          call displayMessage(message)
       end do
       call displayUnindent('done')
    end if
    return
  end subroutine differentialEvolutionPosterior

  function differentialEvolutionStepSize(self,forceAcceptance) result(stepSize)
    !!{
    Return the step size parameter, $\gamma$, for a differential evolution step.
    !!}
    implicit none
    class           (posteriorSampleSimulationDifferentialEvolution), intent(inout)                  :: self
    logical                                                         , intent(inout)                  :: forceAcceptance
    double precision                                                , dimension(self%parameterCount) :: stepSize

    if (self%recomputeCount > 0 .and. mod(self%posteriorSampleState_%count(),self%recomputeCount) == 0) then
       ! Every self%recomputeCount steps, set γ=0 and force likelihood to be recomputed in the current state.
       stepSize       =0.0d0
       forceAcceptance=.true.
    else if (mod(self%posteriorSampleState_%count(),self%stateSwapCount) == 0) then
       ! Every self%stateSwapCount steps, set γ=1 to allow interchange of chains.
       stepSize       =1.0d0
    else
       ! Otherwise, use the step-size algorithm.
       stepSize       =self%posteriorSampleDffrntlEvltnProposalSize_%gamma(                                  &
            &                                                              self%posteriorSampleState_      , &
            &                                                              self%posteriorSampleConvergence_  &
            &                                                             )
    end if
    ! Disable steps in slow parameters, except for on slow steps, or state swap steps.
    if     (                                                                  &
         &   mod(self%posteriorSampleState_%count(),self%slowStepCount ) /= 0 &
         &  .and.                                                             &
         &   mod(self%posteriorSampleState_%count(),self%stateSwapCount) /= 0 &
         & ) then
       where (self%modelParametersActiveIsSlow)
          stepSize=0.0d0
       end where
    end if
    return
  end function differentialEvolutionStepSize

  logical function differentialEvolutionAcceptProposal(self,logPosterior,logPosteriorProposed,logLikelihoodVariance,logLikelihoodVarianceProposed)
    !!{
    Return whether or not to accept a proposal.
    !!}
    implicit none
    class           (posteriorSampleSimulationDifferentialEvolution), intent(inout) :: self
    double precision                                                , intent(in   ) :: logPosterior         , logPosteriorProposed         , &
         &                                                                             logLikelihoodVariance, logLikelihoodVarianceProposed
    double precision                                                                :: x
    !$GLC attributes unused :: self, logLikelihoodVariance, logLikelihoodVarianceProposed

    ! Decide whether to take step.
    x=self%randomNumberGenerator_%uniformSample()
    differentialEvolutionAcceptProposal= logPosteriorProposed >      logPosterior                       &
         &                              .or.                                                            &
         &                               x                    < exp(-logPosterior+logPosteriorProposed)
    return
  end function differentialEvolutionAcceptProposal

  double precision function differentialEvolutionTemperature(self)
    !!{
    Return the temperature.
    !!}
    implicit none
    class(posteriorSampleSimulationDifferentialEvolution), intent(inout) :: self
    !$GLC attributes unused :: self

    differentialEvolutionTemperature=1.0d0
    return
  end function differentialEvolutionTemperature

  subroutine differentialEvolutionDescriptorSpecial(self,descriptor)
    !!{
    Add special parameters to the descriptor.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    class  (posteriorSampleSimulationDifferentialEvolution), intent(inout) :: self
    type   (inputParameters                               ), intent(inout) :: descriptor
    integer                                                                :: i
    
    if (allocated(self%modelParametersActive_  )) then
       do i=1,size(self%modelParametersActive_  )
          call self%modelParametersActive_  (i)%modelParameter_%descriptor(descriptor)
       end do
    end if
    if (allocated(self%modelParametersInactive_)) then
       do i=1,size(self%modelParametersInactive_)
          call self%modelParametersInactive_(i)%modelParameter_%descriptor(descriptor)
       end do
    end if
    return
  end subroutine differentialEvolutionDescriptorSpecial
