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
  Implementation of a posterior sampling likelihood class which implements a likelihood for \glc\ models.
  !!}

  use :: Display        , only : enumerationVerbosityLevelType
  use :: Output_Analyses, only : outputAnalysis               , outputAnalysisClass

  !![
  <posteriorSampleLikelihood name="posteriorSampleLikelihoodGalaxyPopulation">
   <description>A posterior sampling likelihood class which implements a likelihood for \glc\ models.</description>
  </posteriorSampleLikelihood>
  !!]
  type, extends(posteriorSampleLikelihoodBaseParameters) :: posteriorSampleLikelihoodGalaxyPopulation
     !!{
     Implementation of a posterior sampling likelihood class which implements a likelihood for \glc\ models.
     !!}
     private
     type   (varying_string               )          :: failedParametersFileName
     logical                                         :: doPing                            , outputAnalyses, &
          &                                             reportEvaluationTimes             , setOutputGroup, &
          &                                             firstComeFirstServed
     integer                                         :: countCollaborativeGroups
     type   (enumerationVerbosityLevelType)          :: evolveForestsVerbosity
     class  (*                            ), pointer :: task_                    => null()
     class  (outputAnalysisClass          ), pointer :: outputAnalysis_          => null()
   contains
     final     ::                    galaxyPopulationDestructor
     procedure :: evaluate        => galaxyPopulationEvaluate
     procedure :: functionChanged => galaxyPopulationFunctionChanged
     procedure :: willEvaluate    => galaxyPopulationWillEvaluate
  end type posteriorSampleLikelihoodGalaxyPopulation

  interface posteriorSampleLikelihoodGalaxyPopulation
     !!{
     Constructors for the \refClass{posteriorSampleLikelihoodGalaxyPopulation} posterior sampling likelihood class.
     !!}
     module procedure galaxyPopulationConstructorParameters
     module procedure galaxyPopulationConstructorInternal
  end interface posteriorSampleLikelihoodGalaxyPopulation

  ! Sub-module-scope pointer to self, used to allow writing of current parameters in case of failures.
  class  (posteriorSampleLikelihoodGalaxyPopulation), pointer :: self_
  integer(c_size_t                                 )          :: iRank_
  !$omp threadprivate(self_,iRank_)
  
contains

  function galaxyPopulationConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{posteriorSampleLikelihoodGalaxyPopulation} posterior sampling likelihood class which builds the object
    from a parameter set.
    !!}
    use :: Display         , only : displayVerbosity, enumerationVerbosityLevelDecode, enumerationVerbosityLevelEncode
    use :: Input_Parameters, only : inputParameter  , inputParameters
    implicit none
    type   (posteriorSampleLikelihoodGalaxyPopulation)                              :: self
    type   (inputParameters                          ), intent(inout)               :: parameters
    type   (varying_string)                                                         :: baseParametersFileName   , failedParametersFileName
    type   (varying_string)                           , allocatable  , dimension(:) :: changeParametersFileNames
    integer                                                                         :: countCollaborativeGroups
    logical                                                                         :: doPing                   , outputAnalyses          , &
         &                                                                             reportEvaluationTimes    , reportFileName          , &
         &                                                                             reportState              , setOutputGroup          , &
         &                                                                             firstComeFirstServed
    type   (varying_string)                                                         :: evolveForestsVerbosity
    type   (inputParameters                          ), pointer                     :: parametersModel

    !![
    <inputParameter>
      <name>baseParametersFileName</name>
      <description>The base set of parameters to use.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>outputAnalyses</name>
      <description>If true, results of the analyses on each step will be stored to the output file.</description>
      <defaultValue>.false.</defaultValue>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>setOutputGroup</name>
      <description>If true, set the primary output group for each step.</description>
      <defaultValue>.false.</defaultValue>
      <source>parameters</source>
    </inputParameter>
     <inputParameter>
     <name>reportEvaluationTimes</name>
      <description>If true, report the time taken to evaluate each model.</description>
      <defaultValue>.false.</defaultValue>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>countCollaborativeGroups</name>
      <description>The number of groups into which MPI processes should be split for the purpose of model evaluation. Each group will be populated with MPI processes. Processes within a group will collaborate on running a single model evaluation. Setting this parameter to a negative value will result in using the maximum possible number of groups (equal to the number of MPI processes).</description>
      <defaultValue>1</defaultValue>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>firstComeFirstServed</name>
      <description>If true, collaborative groups will take the next available evaluation and process it. Otherwise, each group will be assigned a fixed set of evaluations in advance.</description>
      <defaultValue>.false.</defaultValue>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>doPing</name>
      <defaultValue>.false.</defaultValue>
      <description>
        If true, the master MPI process will attach to the {\normalfont \ttfamily calculationReset} event and ping the MPI
        counter. This can help to ensure that the counter updates regularly.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>reportFileName</name>
      <description>If true, report the base parameter file name being evaluated.</description>
      <defaultValue>.false.</defaultValue>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>reportState</name>
      <description>If true, report the state being evaluated.</description>
      <defaultValue>.false.</defaultValue>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>evolveForestsVerbosity</name>
      <description>The verbosity level to use while performing evolve forests tasks.</description>
      <defaultValue>enumerationVerbosityLevelDecode(displayVerbosity(),includePrefix=.false.)</defaultValue>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>failedParametersFileName</name>
      <description>The name of the file to which parameters of failed models should be written.</description>
      <defaultValue>var_str('./failedParameters.xml')</defaultValue>
      <source>parameters</source>
    </inputParameter>
    !!]
    allocate(changeParametersFileNames(parameters%count('changeParametersFileNames',zeroIfNotPresent=.true.)))
    if (size(changeParametersFileNames) > 0) then
       !![
       <inputParameter>
	 <name>changeParametersFileNames</name>
	 <description>The names of files containing parameter changes to be applied.</description>
	 <source>parameters</source>
       </inputParameter>
       !!]
    end if
    allocate(parametersModel)
    parametersModel=inputParameters                          (baseParametersFileName,noOutput=.true.,changeFiles=changeParametersFileNames)
    self           =posteriorSampleLikelihoodGalaxyPopulation(parametersModel,baseParametersFileName,outputAnalyses,setOutputGroup,reportEvaluationTimes,countCollaborativeGroups,firstComeFirstServed,doPing,reportFileName,reportState,enumerationVerbosityLevelEncode(evolveForestsVerbosity,includesPrefix=.false.),failedParametersFileName,changeParametersFileNames)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    nullify(parametersModel)
    return
  end function galaxyPopulationConstructorParameters


  function galaxyPopulationConstructorInternal(parametersModel,baseParametersFileName,outputAnalyses,setOutputGroup,reportEvaluationTimes,countCollaborativeGroups,firstComeFirstServed,doPing,reportFileName,reportState,evolveForestsVerbosity,failedParametersFileName,changeParametersFileNames) result(self)
    !!{
    Constructor for the \refClass{posteriorSampleLikelihoodGalaxyPopulation} posterior sampling likelihood class.
    !!}
    use :: MPI_Utilities     , only : mpiSelf
    use :: Error             , only : Error_Report
    use :: Display           , only : displayGreen, displayReset
    use :: ISO_Varying_String, only : var_str
    use :: String_Handling   , only : operator(//)
    implicit none
    type   (posteriorSampleLikelihoodGalaxyPopulation)                              :: self
    type   (inputParameters                          ), intent(inout), target       :: parametersModel
    logical                                           , intent(in   )               :: doPing                   , outputAnalyses        , &
         &                                                                             reportEvaluationTimes    , reportFileName        , &
         &                                                                             reportState              , setOutputGroup        , &
         &                                                                             firstComeFirstServed
    integer                                           , intent(in   )               :: countCollaborativeGroups
    type   (enumerationVerbosityLevelType            ), intent(in   )               :: evolveForestsVerbosity
    type   (varying_string                           ), intent(in   )               :: failedParametersFileName , baseParametersFileName
    type   (varying_string                           ), intent(in   ), dimension(:) :: changeParametersFileNames
    !![
    <constructorAssign variables="*parametersModel, baseParametersFileName, outputAnalyses, setOutputGroup, reportEvaluationTimes, countCollaborativeGroups, firstComeFirstServed, doPing, reportFileName, reportState, evolveForestsVerbosity, failedParametersFileName, changeParametersFileNames"/>
    !!]

    if (setOutputGroup.and.countCollaborativeGroups < mpiSelf%count()) call Error_Report('[setOutputGroup]=true and [countCollaborativeGroups] less than the number of MPI processes is not recommended'//char(10)//displayGreen()//'  HELP: '//displayReset()//'[setOutputGroup]=true suggests that you want results of each model evaluation written to its own group, but [countCollaborativeGroups]>1 results in each MPI process evolving a subset of trees from a model evaluation, and writing them to its own output file - this will result in a random mix of trees in each output group - it is recommended that you set [countCollaborativeGroups] equal to the number of MPI processes to avoid this problem'//{introspection:location})
    if      (countCollaborativeGroups < 0) then
       self%countCollaborativeGroups=mpiSelf%count()
    else if (countCollaborativeGroups < 1 .or. countCollaborativeGroups > mpiSelf%count()) then
       call Error_Report(var_str('1 ≤ [countCollaborativeGroups] ≤ ')//mpiSelf%count()//' is required '//{introspection:location})
    end if
    return
  end function galaxyPopulationConstructorInternal

  subroutine galaxyPopulationDestructor(self)
    !!{
    Destructor for the \refClass{posteriorSampleLikelihoodGalaxyPopulation} posterior sampling likelihood class.
    !!}
    implicit none
    type(posteriorSampleLikelihoodGalaxyPopulation), intent(inout) :: self

    if (associated(self%parametersModel)) then
       call self%parametersModel%destroy()
       deallocate(self%parametersModel)
    end if
    return
  end subroutine galaxyPopulationDestructor

  double precision function galaxyPopulationEvaluate(self,simulationState,modelParametersActive_,modelParametersInactive_,simulationConvergence,temperature,logLikelihoodCurrent,logPriorCurrent,logPriorProposed,timeEvaluate,logLikelihoodVariance,forceAcceptance)
    !!{
    Return the log-likelihood for the \glc\ likelihood function.
    !!}
    use :: Display                       , only : displayIndent                  , displayMessage               , displayUnindent             , displayVerbosity             , &
          &                                       displayVerbositySet            , verbosityLevelSilent         , verbosityLevelStandard      , enumerationVerbosityLevelType
    use :: Functions_Global              , only : Tasks_Evolve_Forest_Construct_ , Tasks_Evolve_Forest_Destruct_, Tasks_Evolve_Forest_Perform_
    use :: Error                         , only : errorStatusSuccess             , signalHandlerRegister        , signalHandlerDeregister     , signalHandlerInterface
    use :: Events_Hooks                  , only : calculationResetEvent          , openMPThreadBindingAllLevels
    use :: ISO_Varying_String            , only : char                           , operator(//)                 , var_str
    use :: Kind_Numbers                  , only : kind_int8
    use :: MPI_Utilities                 , only : mpiBarrier                     , mpiSelf                      , mpiCounter
    use :: Model_Parameters              , only : modelParameterDerived
    use :: Models_Likelihoods_Constants  , only : logImpossible                  , logImprobable
    use :: Output_HDF5_Open              , only : Output_HDF5_Set_Group
    use :: Numerical_Constants_Prefixes  , only : siFormat
    use :: Posterior_Sampling_Convergence, only : posteriorSampleConvergenceClass
    use :: Posterior_Sampling_State      , only : posteriorSampleStateClass
    use :: String_Handling               , only : String_Count_Words             , String_Join                  , String_Split_Words          , operator(//)
    implicit none
    class           (posteriorSampleLikelihoodGalaxyPopulation), intent(inout), target         :: self
    class           (posteriorSampleStateClass                ), intent(inout)                 :: simulationState
    type            (modelParameterList                       ), intent(inout), dimension(:  ) :: modelParametersActive_, modelParametersInactive_
    class           (posteriorSampleConvergenceClass          ), intent(inout)                 :: simulationConvergence
    double precision                                           , intent(in   )                 :: temperature           , logLikelihoodCurrent    , &
         &                                                                                        logPriorCurrent       , logPriorProposed
    real                                                       , intent(inout)                 :: timeEvaluate
    double precision                                           , intent(  out), optional       :: logLikelihoodVariance
    logical                                                    , intent(inout), optional       :: forceAcceptance
    double precision                                           , allocatable  , dimension(:  ) :: logPriorsProposed     , logLikelihoods
    real                                                       , allocatable  , dimension(:  ) :: timesEvaluate
    integer                                                    , allocatable  , dimension(:  ) :: chainIndex
    double precision                                           , allocatable  , dimension(:,:) :: stateVector
    procedure       (signalHandlerInterface                   ), pointer                       :: handler
    type            (mpiCounter                               )                                :: evaluationCounter
    integer                                                                                    :: iRank                 , status                  , &
         &                                                                                        rankStart             , rankStop                , &
         &                                                                                        countPerGroup         , iGroup                  , &
         &                                                                                        rankOriginal          , countOriginal
    type            (enumerationVerbosityLevelType            )                                :: verbosityLevel
    real                                                                                       :: timeBegin             , timeEnd
    double precision                                                                           :: logLikelihoodProposed
    character       (len=24                                   )                                :: valueText
    type            (varying_string                           )                                :: message               , groupName
    logical                                                                                    :: isActive
    !$GLC attributes unused :: logPriorCurrent, logLikelihoodCurrent, forceAcceptance, temperature, simulationConvergence

    ! Register an error handler.
    self_   => self
    handler => posteriorSampleLikelihoodGalaxyPopulationSignalHandler
    call signalHandlerRegister(handler)
    ! Set the output group if required.
    if (self%setOutputGroup) then
       groupName=var_str("step")//simulationState%count()//":chain"//simulationState%chainIndex()
       call Output_HDF5_Set_Group(groupName)
    end if
    ! Switch verbosity level.
    verbosityLevel=displayVerbosity()
    call displayVerbositySet(self%evolveForestsVerbosity)
    ! Allocate arrays for likelihoods and evaluation times. These will accumulate results, and will then be transferred back to
    ! the relevant process.
    allocate(logLikelihoods(0:mpiSelf%count()-1))
    allocate(timesEvaluate (0:mpiSelf%count()-1))
    logLikelihoods=0.0d0
    timesEvaluate =0.0d0
    ! Get proposed priors for all chains so we can decide which to skip.
    allocate(logPriorsProposed(0:mpiSelf%count()-1))
    logPriorsProposed=mpiSelf%gather(logPriorProposed)
    ! Get states for all chains.
    allocate(stateVector(simulationState%dimension(),0:mpiSelf%count()-1))
    stateVector=mpiSelf%gather(simulationState%get())
    ! Get chain indices for all chains.
    allocate(chainIndex(0:mpiSelf%count()-1))
    chainIndex=mpiSelf%gather(simulationState%chainIndex())
    ! Ensure pointers into the base parameters are initialized.
    call self%initialize(modelParametersActive_,modelParametersInactive_)
    ! Get a counter if needed. Do this before we split communicators so that all processes have access to this counter.
    if (self%firstComeFirstServed) then
       evaluationCounter=mpiCounter()
       ! Attach to the calculation reset event so that we can ping the counter to avoid slow responses.
       if (self%doPing) call calculationResetEvent%attach(self,evaluationCounterPing,openMPThreadBindingAllLevels,label='evaluationCounterPing')
    end if
    ! If more than one collaborative group is to be used, split the MPI communicator here.
    if (mpiSelf%rank() == 0 .and. verbosityLevel >= verbosityLevelStandard) call displayIndent('Begin collaborative group assignment',verbosityLevelSilent)
    rankOriginal =mpiSelf%rank ()
    countOriginal=mpiSelf%count()
    if (self%countCollaborativeGroups == 1) then
       rankStart=              0
       rankStop =countOriginal-1
       iGroup   =              1
       if (mpiSelf%rank() == 0 .and. verbosityLevel >= verbosityLevelStandard) call displayMessage('All processes collaborating in a single group',verbosityLevelSilent)
    else
       countPerGroup=    +     countOriginal             &
            &            /self%countCollaborativeGroups
       iGroup       =min(                                &
            &            +     rankOriginal              &
            &            /     countPerGroup             &
            &            +     1                       , &
            &            +self%countCollaborativeGroups  &
            &           )
       rankStart   =(iGroup-1)*countPerGroup
       if (iGroup < self%countCollaborativeGroups) then
          rankStop = iGroup   *countPerGroup-1
       else
          rankStop =           countOriginal-1
       end if
       if (verbosityLevel >= verbosityLevelStandard) then
          message=var_str('Process ')//mpiSelf%rank()//' belongs to group '//iGroup
          if (.not.self%firstComeFirstServed) then
             message=message//' and will collaborate on '
             if (rankStop > rankStart) then
                message=message//'evaluations '//rankStart//' to '//rankStop
             else
                message=message//'evaluation ' //rankStart
             end if
          end if
          call displayMessage(message,verbosityLevelSilent)
       end if
       call mpiSelf%communicatorPush(color=iGroup)
    end if
    if (rankOriginal == 0 .and. verbosityLevel >= verbosityLevelStandard) call displayUnindent('done',verbosityLevelSilent)
    ! If using a first-come-first-served approach, reset the rank range to the full extent - this is sufficient to ensure that we
    ! make enough requests to the counter such that all evaluations are always performed.
    if (self%firstComeFirstServed) then
       rankStart=              0
       rankStop =countOriginal-1
    end if
    ! Iterate over all chains.
    do iRank=rankStart,rankStop
       ! Determine if this is the active rank for reporting and storing results.
       isActive=mpiSelf%rank() == 0
       ! Find which model to evaluate.
       if (self%firstComeFirstServed) then
          if (isActive) then
             call CPU_Time(timeBegin)
             iRank_=evaluationCounter%increment()
             call CPU_Time(timeEnd  )
          else
             iRank_=0_c_size_t
          end if
          call mpiBarrier()
          iRank_=mpiSelf%sum(iRank_)
          if (iRank_ > rankStop) exit
          if (isActive .and. verbosityLevel >= verbosityLevelStandard) then
             message=var_str('group ')//iGroup//' will collaborate on evaluation '//iRank_//' (wait time for task: '//trim(siFormat(dble(timeEnd-timeBegin),'f10.6,1x'))//'s)'
             call displayMessage(message,verbosityLevelSilent)
          end if
       else
          iRank_=iRank
       end if
       ! Initialize likelihood to impossible.
       if (isActive) then
          logLikelihoods(iRank_)=logImpossible
          if (present(logLikelihoodVariance)) logLikelihoodVariance=0.0d0
       end if
       ! If prior probability is impossible, then no need to waste time evaluating the likelihood.
       if (logPriorsProposed(iRank_) <= logImpossible) cycle
       ! If the likelihood was evaluated for the previous rank, and the current state vector is identical to that of the previous
       ! rank, the proposed likelihood must be unchanged.
       if (.not.self%firstComeFirstServed .and. iRank > rankStart) then
          if (logPriorsProposed(iRank-1) > logImpossible .and. all(stateVector(:,iRank) == stateVector(:,iRank-1))) then
             if (isActive) then
                logLikelihoods(iRank_)=logLikelihoodProposed
                timesEvaluate (iRank_)=0.0d0
                if (verbosityLevel >= verbosityLevelStandard) then
                   write (valueText,'(e12.4)') logLikelihoodProposed
                   message=var_str("Chain ")//chainIndex(iRank_)//" has logℒ="//trim(valueText)
                   call displayMessage(message,verbosityLevelSilent)
                end if
             end if
             cycle
          end if
       end if
       ! Update parameter values.
       call self%update(simulationState,modelParametersActive_,modelParametersInactive_,stateVector(:,iRank_),report=isActive)
       ! Build the task and outputter objects.
       call Tasks_Evolve_Forest_Construct_(self%parametersModel,self%task_)
       !![
       <objectBuilder class="outputAnalysis" name="self%outputAnalysis_" source="self%parametersModel"/>
       !!]
       ! Perform the forest evolution tasks.
       call CPU_Time(timeBegin)
       call Tasks_Evolve_Forest_Perform_(self%task_,status)
       if (mpiSelf%any(status /= errorStatusSuccess)) then
          ! Forest evolution failed - record impossible likelihood.
          if (isActive) then
             ! Dump the failed parameter set to file.
             call self%parametersModel%serializeToXML(self%failedParametersFileName//"."//iRank_//".errCode"//status)
             ! Return impossible likelihood. We use a somewhat-less-than-impossible value to avoid this being rejected as the
             ! initial state.
             logLikelihoodProposed        =logImprobable
             logLikelihoods       (iRank_)=logLikelihoodProposed
          end if
       else
          ! Forest evolution was successful - evaluate the likelihood.
          ! Extract the log-likelihood. This is evaluated by all chains (as they likely need to perform reduction across MPI
          ! processes), but only stored for the chain of this rank.
          logLikelihoodProposed=self%outputAnalysis_%logLikelihood()
          ! For active analysis, return this likelihood.
          if (isActive) then
             logLikelihoods(iRank_)=logLikelihoodProposed
             ! Optionally output results.
             if (self%outputAnalyses) then
                groupName=var_str("step")//simulationState%count()//":chain"//simulationState%chainIndex()
                call self%outputAnalysis_%finalize(groupName)
             end if
          end if
       end if
       if (isActive) then
          ! Record timing information.
          call CPU_Time(timeEnd)
          timesEvaluate(iRank_)=timeEnd-timeBegin
          if (verbosityLevel >= verbosityLevelStandard) then
             write (valueText,'(e12.4)') logLikelihoodProposed
             message=var_str("Chain ")//chainIndex(iRank_)//" has logℒ="//trim(valueText)
             if (self%reportEvaluationTimes) message=message//" (evaluation in "//trim(siFormat(dble(timeEvaluate),'f5.1,1x'))//"s)"
             call displayMessage(message,verbosityLevelSilent)
          end if
       end if
       call mpiBarrier()
       call Tasks_Evolve_Forest_Destruct_(self%task_)
       !![
       <objectDestructor name="self%outputAnalysis_"/>
       !!]
       call self%parametersModel%reset()
    end do
    ! If we use multiple collaborative groups, join the MPI communicator here.
    if (self%countCollaborativeGroups > 1) then
       call mpiBarrier                ()
       call mpiSelf   %communicatorPop()
    end if
    if (self%firstComeFirstServed .and. self%doPing) &
         & call calculationResetEvent%detach(self,evaluationCounterPing)
    ! Retrieve and extract the log-likelihood and evaluation time for this process.
    logLikelihoods          =mpiSelf%sum(logLikelihoods)
    timesEvaluate           =mpiSelf%sum(timesEvaluate )
    galaxyPopulationEvaluate=logLikelihoods(rankOriginal)
    timeEvaluate            =timesEvaluate (rankOriginal)
    ! Restore verbosity level.
    call displayVerbositySet(verbosityLevel)
    ! Deregister our error handler.
    call signalHandlerDeregister(handler)
    return

  contains

    subroutine evaluationCounterPing(self,node,uniqueID)
      !!{
      Return the number of the next forest to process.
      !!}
      use :: Galacticus_Nodes, only : treeNode
      use :: Kind_Numbers    , only : kind_int8
      implicit none
      class  (*        ), intent(inout) :: self
      type   (treeNode ), intent(inout) :: node
      integer(kind_int8), intent(in   ) :: uniqueID
#ifdef USEMPI
      integer(c_size_t )                :: evaluationNumber
#endif
      !$GLC attributes unused :: self, node, uniqueID

#ifdef USEMPI
      !$omp master
      if (mpiSelf%isMaster()) evaluationNumber=evaluationCounter%get()
      !$omp end master
#endif
      return
    end subroutine evaluationCounterPing

  end function galaxyPopulationEvaluate

  subroutine galaxyPopulationFunctionChanged(self)
    !!{
    Respond to possible changes in the likelihood function.
    !!}
    implicit none
    class(posteriorSampleLikelihoodGalaxyPopulation), intent(inout) :: self
    !$GLC attributes unused :: self

    return
  end subroutine galaxyPopulationFunctionChanged

  logical function galaxyPopulationWillEvaluate(self,simulationState,modelParameters_,simulationConvergence,temperature,logLikelihoodCurrent,logPriorCurrent,logPriorProposed)
    !!{
    Return true if the log-likelihood will be evaluated.
    !!}
    use :: Models_Likelihoods_Constants  , only : logImpossible
    use :: Posterior_Sampling_Convergence, only : posteriorSampleConvergenceClass
    use :: Posterior_Sampling_State      , only : posteriorSampleStateClass
    implicit none
    class           (posteriorSampleLikelihoodGalaxyPopulation), intent(inout)               :: self
    class           (posteriorSampleStateClass                ), intent(inout)               :: simulationState
    type            (modelParameterList                       ), intent(in   ), dimension(:) :: modelParameters_
    class           (posteriorSampleConvergenceClass          ), intent(inout)               :: simulationConvergence
    double precision                                           , intent(in   )               :: temperature          , logLikelihoodCurrent, &
         &                                                                                      logPriorCurrent      , logPriorProposed
    !$GLC attributes unused :: self, simulationState, modelParameters_, simulationConvergence, temperature, logLikelihoodCurrent, logPriorCurrent

    ! Likelihood will not be evaluated if the proposed prior is impossible.
    galaxyPopulationWillEvaluate=(logPriorProposed > logImpossible)
    return
  end function galaxyPopulationWillEvaluate

  subroutine posteriorSampleLikelihoodGalaxyPopulationSignalHandler(signal)
    !!{
    Write out current parameters if a signal was caught during model evaluation.
    !!}
    use :: Display           , only : displayMessage         , displayBold   , displayRed, displayReset, &
         &                            verbosityLevelSilent
    use :: ISO_Varying_String, only : operator(//)           , varying_string
    use :: String_Handling   , only : operator(//)
    use :: Error_Utilities   , only : enumerationSignalDecode
    implicit none
    integer                , intent(in   ) :: signal
    type   (varying_string)                :: fileName
    
    ! Dump the failed parameter set to file.
    fileName=self_%failedParametersFileName//"."//iRank_//"."//enumerationSignalDecode(signal,includePrefix=.false.)
    call displayMessage(displayRed()//displayBold()//"Error condition:"//displayReset()//" `posteriorSampleLikelihoodGalaxyPopulation` parameter state will be written to '"//fileName//"'",verbosityLevelSilent)
    call self_%parametersModel%serializeToXML(fileName)
    return
  end subroutine posteriorSampleLikelihoodGalaxyPopulationSignalHandler
