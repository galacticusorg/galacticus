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

!+    Contributions to this file made by: Claude.

  !!{RST
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
  <posteriorSampleSimulation name="posteriorSampleSimulationDifferentialEvolution" docformat="rst">
   <description>
   This class uses the differential evolution algorithm of :cite:t:`terr_braak_markov_2006`. Multiple, parallel chains are run and proposals are constructed by selecting two chains at random, taking a fraction, :math:`\gamma`, of the vector connecting the two chain states and adding this to the state of the current chain. The details of the algorithm are controlled by the following parameters:

   ``[stepsMaximum]``
      The maximum number of steps to take.

   ``[acceptanceAverageCount]``
      The number of steps over which to average the acceptance rate.

   ``[stateSwapCount]``
      The number of steps after which to set :math:`\gamma=1` to allow chains to swap states.

   ``[logFileRoot]``
      The full path and root name of a file to log results to. The actual file name will have the rank of the :term:`MPI` process appended to it.

   ``[sampleOutliers]``
      If set to ``false`` then proposals for non-outlier chains post-convergence are constructed only from other non-outlier chains. Otherwise, proposals for non-outlier chains post-convergence are constructed from all other chains.

   ``[proposalArchive]``
      If set to ``true`` (the default) then the two states from which each proposal is constructed are drawn from an archive of past states of all chains, rather than from the current states of other chains (the "DE-MCz" algorithm of :cite:t:`ter_braak_differential_2008`). All chains are updated simultaneously, one per :term:`MPI` process, so when proposals are constructed from the current states of other chains those states are themselves being updated, and the sampled distribution is slightly biased (too narrow by a few percent when the number of chains is comparable to the number of parameters). The archive is unchanged during each update, so proposals constructed from it are unbiased. The archive also allows far fewer chains to be used (as few as three), as its states span the parameter space even when the chains do not.

   ``[archiveThinning]``
      The states of all chains are added to the archive every ``[archiveThinning]`` steps.

   ``[archiveInitialCount]``
      The number of states drawn from the prior with which to initialize the archive (in addition to the initial states of the chains). :cite:t:`ter_braak_differential_2008` suggest ten times the number of active parameters.

   ``[archiveSizeMaximum]``
      The maximum number of states held in the archive. Once it is full, each new state replaces the oldest, so that the archive holds the most recent states of the chains (and forgets those from the prior and from burn-in, which no longer help to construct proposals once the chains have converged). Proposals remain unbiased, as the archive is unchanged during each update. It should be much larger than the number of chains, so that proposals are not constructed only from recent, correlated states.
   </description>
   <descriptorSpecial>descriptorSpecial</descriptorSpecial>
  </posteriorSampleSimulation>
  !!]
  type, extends(posteriorSampleSimulationClass) :: posteriorSampleSimulationDifferentialEvolution
     !!{RST
     Implementation of a posterior sampling simulation class which implements the differential evolution algorithm.
     !!}
     private
     integer                                                                                     :: parameterCount                                    , stepsMaximum            , &
          &                                                                                         stateSwapCount                                    , acceptanceAverageCount  , &
          &                                                                                         logFlushCount                                     , reportCount             , &
          &                                                                                         recomputeCount                                    , slowStepCount           , &
          &                                                                                         archiveThinning                                   , archiveInitialCount     , &
          &                                                                                         archiveCount                                      , archiveSizeMaximum      , &
          &                                                                                         archiveNext
     double precision                                                                            :: logPosterior                                      , logPrior
     logical                                                                                     :: isConverged                                       , sampleOutliers          , &
          &                                                                                         isInteractive                                     , appendLogs              , &
          &                                                                                         loadBalance                                       , ignoreChainNumberAdvice , &
          &                                                                                         logProposals                                      , proposalArchive
     logical                                                       , allocatable, dimension(:  ) :: modelParametersActiveIsSlow
     double precision                                              , allocatable, dimension(:,:) :: archive
     type            (modelParameterList                          ), allocatable, dimension(:  ) :: modelParametersActive_                            , modelParametersInactive_
     class           (posteriorSampleLikelihoodClass              ), pointer                     :: posteriorSampleLikelihood_               => null()
     class           (posteriorSampleConvergenceClass             ), pointer                     :: posteriorSampleConvergence_              => null()
     class           (posteriorSampleStoppingCriterionClass       ), pointer                     :: posteriorSampleStoppingCriterion_        => null()
     class           (posteriorSampleStateClass                   ), pointer                     :: posteriorSampleState_                    => null()
     class           (posteriorSampleStateInitializeClass         ), pointer                     :: posteriorSampleStateInitialize_          => null()
     class           (posteriorSampleDffrntlEvltnProposalSizeClass), pointer                     :: posteriorSampleDffrntlEvltnProposalSize_ => null()
     class           (posteriorSampleDffrntlEvltnRandomJumpClass  ), pointer                     :: posteriorSampleDffrntlEvltnRandomJump_   => null()
     class           (randomNumberGeneratorClass                  ), pointer                     :: randomNumberGenerator_                   => null()
     type            (varying_string                              )                              :: logFileRoot                                       , interactionRoot
   contains
     !![
     <methods docformat="rst">
       <method method="logging"           description="Return true if the simulator is currently logging state."                                       />
       <method method="posterior"         description="Return the log of posterior probability for the given ``posteriorSampleState``."                />
       <method method="update"            description="Update the simulator to the new ``stateVector`` after a step."                                  />
       <method method="temperature"       description="Return the current temperature."                                                                />
       <method method="acceptProposal"    description="Return true if the proposed state should be accepted."                                          />
       <method method="stepSize"          description="Return the step size parameter, :math:`\gamma`, for the differential evolution proposal vector."/>
       <method method="chainSelect"       description="Select a chain."                                                                                />
       <method method="archiveAppend"     description="Append states to the archive of past states."                                                   />
       <method method="archiveSelect"     description="Select two distinct states at random from the archive of past states."                          />
       <method method="descriptorSpecial" description="Handle adding special parameters to the descriptor."                                            />
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
     procedure :: archiveAppend     => differentialEvolutionArchiveAppend
     procedure :: archiveSelect     => differentialEvolutionArchiveSelect
     procedure :: descriptorSpecial => differentialEvolutionDescriptorSpecial
  end type posteriorSampleSimulationDifferentialEvolution

  interface posteriorSampleSimulationDifferentialEvolution
     !!{RST
     Constructors for the :galacticus-class:`posteriorSampleSimulationDifferentialEvolution` posterior sampling simulation class.
     !!}
     module procedure differentialEvolutionConstructorParameters
     module procedure differentialEvolutionConstructorInternal
  end interface posteriorSampleSimulationDifferentialEvolution

contains

  function differentialEvolutionConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`posteriorSampleSimulationDifferentialEvolution` posterior sampling simulation class which builds the object from a parameter set.
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
         &                                                                                  recomputeCount                          , slowStepCount           , &
         &                                                                                  archiveThinning                         , archiveInitialCount     , &
         &                                                                                  archiveSizeMaximum
    type   (varying_string                                )                              :: logFileRoot                             , interactionRoot         , &
         &                                                                                  message
    logical                                                                              :: sampleOutliers                          , appendLogs              , &
         &                                                                                  logProposals                            , proposalArchive         , &
         &                                                                                  loadBalance                             , ignoreChainNumberAdvice

    !![
    <inputParameter docformat="rst">
      <name>stepsMaximum</name>
      <defaultValue>huge(0)</defaultValue>
      <description>
      The maximum number of steps to take.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>acceptanceAverageCount</name>
      <defaultValue>10</defaultValue>
      <description>
      The number of steps over which to average the acceptance rate.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>stateSwapCount</name>
      <defaultValue>10</defaultValue>
      <description>
      The number of steps between state swap steps.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>slowStepCount</name>
      <defaultValue>1</defaultValue>
      <description>
      The number of steps between slow parameter update steps.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>recomputeCount</name>
      <defaultValue>-1</defaultValue>
      <description>
      The number of steps between forced recomputations of the likelihood.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>logFlushCount</name>
      <defaultValue>10</defaultValue>
      <description>
      The number of steps between flushing the log file.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>reportCount</name>
      <defaultValue>10</defaultValue>
      <description>
      The number of steps between issuing reports.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>sampleOutliers</name>
      <defaultValue>.true.</defaultValue>
      <description>
      If true, sample from outlier states.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>interactionRoot</name>
      <defaultValue>var_str('none')</defaultValue>
      <description>
      Root file name for interaction files, or `\ ``none``' if interaction is not required.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>logFileRoot</name>
      <description>
      Root file name for log files.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>appendLogs</name>
      <description>
      If true, do not overwrite existing log files, but instead append to them.
      </description>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>logProposals</name>
      <description>
      If true, write a log of every proposed state to
      ``&lt;logFileRoot&gt;Proposals_&lt;rank&gt;.log``, whether or not that proposal was
      accepted. The chain log records only the state retained at each step, so a
      rejected proposal's parameters are otherwise never written anywhere. This
      matters for any analysis which pairs recorded model outputs (for example
      those written by ``posteriorSampleLikelihoodHaloMassFunction`` via
      ``pathSamples``) with the parameters that produced them: without this log
      only accepted steps can be used, discarding the majority of the evaluations
      and, with them, the wider coverage of parameter space that rejected
      proposals provide.
      </description>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>loadBalance</name>
      <description>
      If true, attempt to balance the workload across different compute nodes.
      </description>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>ignoreChainNumberAdvice</name>
      <description>
      If true, ignore warnings and errors about not being able to span the full parameter space with the number of chains used.
      </description>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>proposalArchive</name>
      <description>
      If true, construct proposals from pairs of states drawn from an archive of past states of all chains, rather than from
      the current states of other chains.
      </description>
      <source>parameters</source>
      <defaultValue>.true.</defaultValue>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>archiveThinning</name>
      <description>
      The number of steps between additions of the states of all chains to the archive of past states.
      </description>
      <source>parameters</source>
      <defaultValue>10</defaultValue>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>archiveInitialCount</name>
      <description>
      The number of states drawn from the prior with which to initialize the archive of past states.
      </description>
      <source>parameters</source>
      <defaultValue>0</defaultValue>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>archiveSizeMaximum</name>
      <description>
      The maximum number of states held in the archive of past states. Once it is full, each new state replaces the oldest.
      </description>
      <source>parameters</source>
      <defaultValue>100000</defaultValue>
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
    self=posteriorSampleSimulationDifferentialEvolution(modelParametersActive_,modelParametersInactive_,posteriorSampleLikelihood_,posteriorSampleConvergence_,posteriorSampleStoppingCriterion_,posteriorSampleState_,posteriorSampleStateInitialize_,posteriorSampleDffrntlEvltnProposalSize_,posteriorSampleDffrntlEvltnRandomJump_,randomNumberGenerator_,stepsMaximum,acceptanceAverageCount,stateSwapCount,slowStepCount,recomputeCount,char(logFileRoot),sampleOutliers,logFlushCount,reportCount,char(interactionRoot),appendLogs,loadBalance,ignoreChainNumberAdvice,logProposals,proposalArchive,archiveThinning,archiveInitialCount,archiveSizeMaximum)
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
  
  function differentialEvolutionConstructorInternal(modelParametersActive_,modelParametersInactive_,posteriorSampleLikelihood_,posteriorSampleConvergence_,posteriorSampleStoppingCriterion_,posteriorSampleState_,posteriorSampleStateInitialize_,posteriorSampleDffrntlEvltnProposalSize_,posteriorSampleDffrntlEvltnRandomJump_,randomNumberGenerator_,stepsMaximum,acceptanceAverageCount,stateSwapCount,slowStepCount,recomputeCount,logFileRoot,sampleOutliers,logFlushCount,reportCount,interactionRoot,appendLogs,loadBalance,ignoreChainNumberAdvice,logProposals,proposalArchive,archiveThinning,archiveInitialCount,archiveSizeMaximum) result(self)
    !!{RST
    Internal constructor for the "differentialEvolution" simulation class.
    !!}
    use :: Error           , only : Error_Report
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
         &                                                                                     logProposals                                                    , &
         &                                                                                    loadBalance                             , ignoreChainNumberAdvice
    logical                                                  , intent(in   ), optional     :: proposalArchive
    integer                                                  , intent(in   ), optional     :: archiveThinning                         , archiveInitialCount     , &
         &                                                                                    archiveSizeMaximum
    integer                                                                                :: i
    !![
    <constructorAssign variables="*posteriorSampleLikelihood_, *posteriorSampleConvergence_, *posteriorSampleStoppingCriterion_, *posteriorSampleState_, *posteriorSampleStateInitialize_, *posteriorSampleDffrntlEvltnProposalSize_, *posteriorSampleDffrntlEvltnRandomJump_, *randomNumberGenerator_, stepsMaximum, acceptanceAverageCount, stateSwapCount, slowStepCount, recomputeCount, logFlushCount, reportCount, sampleOutliers, logFileRoot, interactionRoot, appendLogs, loadBalance, ignoreChainNumberAdvice, logProposals"/>
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
    self%proposalArchive    =.true.
    self%archiveThinning    =10
    self%archiveInitialCount= 0
    self%archiveSizeMaximum =100000
    self%archiveCount       = 0
    self%archiveNext        = 1
    if (present(proposalArchive    )) self%proposalArchive    =proposalArchive
    if (present(archiveThinning    )) self%archiveThinning    =archiveThinning
    if (present(archiveInitialCount)) self%archiveInitialCount=archiveInitialCount
    if (present(archiveSizeMaximum )) self%archiveSizeMaximum =archiveSizeMaximum
    if (self%archiveThinning     < 1) call Error_Report('[archiveThinning] must be positive'        //{introspection:location})
    if (self%archiveInitialCount < 0) call Error_Report('[archiveInitialCount] must be non-negative'//{introspection:location})
    if (self%archiveSizeMaximum  < 2) call Error_Report('[archiveSizeMaximum] must be at least 2'   //{introspection:location})
    self%parameterCount=size(modelParametersActive_)
    self%isInteractive =trim(interactionRoot) /= "none"
    call self%posteriorSampleState_%parameterCountSet(self%parameterCount)
    return
  end function differentialEvolutionConstructorInternal

  subroutine differentialEvolutionDestructor(self)
    !!{RST
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
    !!{RST
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
         &                                                                                                interactionFileName   , proposalFileName
    integer                                                                                            :: logFileUnit           , convergedAtStep              , &
         &                                                                                                convergenceFileUnit   , i                            , &
         &                                                                                                ioStatus              , interactionFile              , &
         &                                                                                                proposalFileUnit
    double precision                                               , allocatable, dimension(:)         :: stateVectorProposedLog
    double precision                                               , allocatable, dimension(:,:)       :: archiveInitial
    integer                                                                                            :: j
    logical                                                                                               forceAcceptance       , proposalAccepted
    character       (len=32                                        )                                   :: label

    ! Check that we have sufficient chains for differential evolution. (With an archive of past states, proposals are drawn from
    ! the archive, which spans the parameter space even when the chains do not.)
    if (.not.self%ignoreChainNumberAdvice .and. .not.self%proposalArchive) then
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
    ! Initialize the archive of past states with states drawn from the prior (drawn on one process, so that the archive is the
    ! same on all processes), and the initial states of all chains.
    if (self%proposalArchive) then
       if (self%archiveInitialCount > 0) then
          allocate(archiveInitial(self%parameterCount,self%archiveInitialCount))
          if (mpiSelf%rank() == 0) then
             do j=1,self%archiveInitialCount
                do i=1,self%parameterCount
                   archiveInitial(i,j)=self%modelParametersActive_(i)%modelParameter_%map(self%modelParametersActive_(i)%modelParameter_%priorSample())
                end do
             end do
          end if
          call mpiSelf%broadcastData(0,archiveInitial)
          call self%archiveAppend(archiveInitial)
          deallocate(archiveInitial)
       end if
       call self%archiveAppend(mpiSelf%gather(self%posteriorSampleState_%get()))
    end if
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
    ! Optionally open a log of proposed states. The chain log records only the state
    ! retained at each step, so without this a rejected proposal's parameters are
    ! never written anywhere and any model output recorded during its evaluation
    ! cannot be attributed to a point in parameter space.
    if (self%logProposals) then
       allocate(stateVectorProposedLog(size(self%modelParametersActive_)))
       proposalFileName=self%logFileRoot//'Proposals_'//mpiSelf%rankLabel()//'.log'
       if (self%appendLogs) then
          open(newunit=proposalFileUnit,file=char(proposalFileName),status='unknown',form='formatted',position='append')
       else
          open(newunit=proposalFileUnit,file=char(proposalFileName),status='unknown',form='formatted'                  )
          write (proposalFileUnit,'(a)') '# Proposed state log'
          write (proposalFileUnit,'(a)') '# One row per proposal evaluated, whether or not it was accepted.'
          write (proposalFileUnit,'(a)') '# Columns:'
          write (proposalFileUnit,'(a)') '#    1 = Simulation step'
          write (proposalFileUnit,'(a)') '#    2 = Chain index'
          write (proposalFileUnit,'(a)') '#    3 = Proposal accepted? [T/F]'
          write (proposalFileUnit,'(a)') '#    4 = log posterior of proposal'
          write (proposalFileUnit,'(a)') '#    5 = log likelihood of proposal'
          do i=1,size(self%modelParametersActive_)
             write (proposalFileUnit,'(a,i3,a,a,a)') '#  ',i+5,' = Parameter `',char(self%modelParametersActive_(i)%modelParameter_%name()),'`'
          end do
       end if
    end if
    self%isConverged=.false.
    do while (                                                                                                   &
         &          self%posteriorSampleState_            %count(                          ) < self%stepsMaximum &
         &    .and.                                                                                              &
         &     .not.self%posteriorSampleStoppingCriterion_%stop (self%posteriorSampleState_)                     &
         &   )
       ! Initialize likelihood acceptance condition.
       forceAcceptance=.false.
       stateVector=self%posteriorSampleState_%get()
       if (self%proposalArchive) then
          ! Pick two distinct states from the archive of past states to use for proposal.
          statePair=self%archiveSelect()
       else
          ! Pick two random processes, other than this one, to use for proposal. Allowing this chain to be chosen would make the
          ! proposal asymmetric (an expansion away from, or contraction toward, the other chosen chain), biasing the sampled
          ! distribution.
          chainPair(1)=self%chainSelect([mpiSelf%rank()             ])
          chainPair(2)=self%chainSelect([mpiSelf%rank(),chainPair(1)])
          ! Receive states from selected chains.
          statePair  =mpiSelf%requestData(chainPair,stateVector)
       end if
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
          ! Label the proposal with the step at which it will be recorded in the
          ! chain log. `self%posteriorSampleState_` is not advanced until the
          ! `update()` call following the accept/reject decision below, so its
          ! count here is one behind. Any likelihood which tags its output with
          ! `simulationState%count()` - the sampled model vectors written by
          ! `posteriorSampleLikelihoodHaloMassFunction`, the per-step groups
          ! written by `posteriorSampleLikelihoodGalaxyPopulation` - would
          ! otherwise be offset by one step from the chain log, making the two
          ! impossible to join without knowing the convention.
          call stateProposed%countSet(self%posteriorSampleState_%count()+1)
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
          proposalAccepted     =.true.
       else
          ! Step not accepted - retain the old estimate of work time.
          timeEvaluate         =timeEvaluatePrevious
          proposalAccepted     =.false.
       end if
       call self%update(stateVector)
       ! Add the states of all chains to the archive of past states.
       if (self%proposalArchive .and. mod(self%posteriorSampleState_%count(),self%archiveThinning) == 0) &
            & call self%archiveAppend(mpiSelf%gather(stateVector))
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
          if (self%logProposals) then
             ! Unmap a copy of the proposed state - `stateVectorProposed` itself must
             ! be left mapped for use in subsequent steps.
             do i=1,size(stateVectorProposed)
                stateVectorProposedLog(i)=self%modelParametersActive_(i)%modelParameter_%unmap(stateVectorProposed(i))
             end do
             write (proposalFileUnit,*) self%posteriorSampleState_%count        (), &
                  &                     mpiSelf                   %rank         (), &
                  &                                                proposalAccepted, &
                  &                                                logPosteriorProposed, &
                  &                                                logLikelihoodProposed, &
                  &                                                stateVectorProposedLog
             if (mod(self%posteriorSampleState_%count(),self%logFlushCount) == 0) call flush(proposalFileUnit)
          end if
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
    if (self%logProposals) close(proposalFileUnit)
    return
  end subroutine differentialEvolutionSimulate

  subroutine differentialEvolutionUpdate(self,stateVector)
    !!{RST
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
    !!{RST
    Select a chain at random, optionally excluding blocked chains.
    !!}
    use :: MPI_Utilities, only : mpiSelf
    use :: Error        , only : Error_Report
    implicit none
    class  (posteriorSampleSimulationDifferentialEvolution), intent(inout)                         :: self
    integer                                                , intent(in   ), dimension(:), optional :: blockedChains
    logical                                                                                        :: accept

    if (mpiSelf%count() < 3) call Error_Report('at least three chains are needed for differential evolution'//{introspection:location})
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
       if (.not.self%sampleOutliers.and.self%isConverged)                                                    &
            & accept=                                                                                        &
            &              self%posteriorSampleConvergence_%stateIsOutlier(                  mpiSelf%rank()) &
            &        .or.                                                                                    &
            &         .not.self%posteriorSampleConvergence_%stateIsOutlier(differentialEvolutionChainSelect)
       if (present(blockedChains)) accept=accept.and.all(differentialEvolutionChainSelect /= blockedChains)
    end do
    return
  end function differentialEvolutionChainSelect

  subroutine differentialEvolutionArchiveAppend(self,states)
    !!{RST
    Append states (one per column of ``states``) to the archive of past states. The archive grows as needed up to its maximum
    size, after which each new state replaces the oldest.
    !!}
    implicit none
    class           (posteriorSampleSimulationDifferentialEvolution), intent(inout)                 :: self
    double precision                                                , intent(in   ), dimension(:,:) :: states
    double precision                                                , allocatable  , dimension(:,:) :: archive
    integer                                                                                         :: i      , sizeNew

    do i=1,size(states,dim=2)
       if (self%archiveCount < self%archiveSizeMaximum) then
          ! The archive is not yet full - grow the storage as needed.
          if (.not.allocated(self%archive)) then
             allocate(self%archive(self%parameterCount,min(self%archiveSizeMaximum,max(1024,2*size(states,dim=2)))))
          else if (self%archiveCount == size(self%archive,dim=2)) then
             sizeNew=min(self%archiveSizeMaximum,2*size(self%archive,dim=2))
             call move_alloc(self%archive,archive)
             allocate(self%archive(self%parameterCount,sizeNew))
             self%archive(:,1:self%archiveCount)=archive(:,1:self%archiveCount)
             deallocate(archive)
          end if
          self%archiveCount=self%archiveCount+1
       end if
       ! Store the state, replacing the oldest once the archive is full.
       self%archive(:,self%archiveNext)=states(:,i)
       self%archiveNext                =mod(self%archiveNext,self%archiveSizeMaximum)+1
    end do
    return
  end subroutine differentialEvolutionArchiveAppend

  function differentialEvolutionArchiveSelect(self) result(statePair)
    !!{RST
    Select two distinct states at random from the archive of past states.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (posteriorSampleSimulationDifferentialEvolution), intent(inout)                    :: self
    double precision                                                , dimension(self%parameterCount,2) :: statePair
    integer                                                         , dimension(                    2) :: indices

    if (self%archiveCount < 2) call Error_Report('at least two states are needed in the archive'//{introspection:location})
    indices(1)=min(int(dble(self%archiveCount  )*self%randomNumberGenerator_%uniformSample())+1,self%archiveCount  )
    indices(2)=min(int(dble(self%archiveCount-1)*self%randomNumberGenerator_%uniformSample())+1,self%archiveCount-1)
    ! Draw the second index from the remaining states, so that the two are distinct.
    if (indices(2) >= indices(1)) indices(2)=indices(2)+1
    statePair(:,1)=self%archive(:,indices(1))
    statePair(:,2)=self%archive(:,indices(2))
    return
  end function differentialEvolutionArchiveSelect

  logical function differentialEvolutionLogging(self)
    !!{RST
    Specifies whether or not the current state should be logged to file during differential evolution.
    !!}
    implicit none
    class(posteriorSampleSimulationDifferentialEvolution), intent(inout) :: self
    !$GLC attributes unused :: self

    differentialEvolutionLogging=.true.
    return
  end function differentialEvolutionLogging

  subroutine differentialEvolutionPosterior(self,posteriorSampleState_,logLikelihoodCurrent,logPriorCurrent,logPosterior,logLikelihood,logLikelihoodVariance,timeEvaluate,timeEvaluatePrevious,forceAcceptance)
    !!{RST
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
    !!{RST
    Return the step size parameter, :math:`\gamma`, for a differential evolution step.
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
    !!{RST
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
    !!{RST
    Return the temperature.
    !!}
    implicit none
    class(posteriorSampleSimulationDifferentialEvolution), intent(inout) :: self
    !$GLC attributes unused :: self

    differentialEvolutionTemperature=1.0d0
    return
  end function differentialEvolutionTemperature

  subroutine differentialEvolutionDescriptorSpecial(self,descriptor)
    !!{RST
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
