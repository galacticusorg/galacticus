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

  use :: Model_Parameters                    , only : modelParameterList
  use :: Models_Likelihoods                  , only : posteriorSampleLikelihoodClass
  use :: Numerical_Random_Numbers            , only : randomNumberGeneratorClass
  use :: Posterior_Sampling_Convergence      , only : posteriorSampleConvergenceClass
  use :: Posterior_Sampling_State            , only : posteriorSampleStateClass
  use :: Posterior_Sampling_State_Initialize , only : posteriorSampleStateInitializeClass
  use :: Posterior_Sampling_Stopping_Criteria, only : posteriorSampleStoppingCriterionClass

  !![
  <posteriorSampleSimulation name="posteriorSampleSimulationParticleSwarm">
   <description>A posterior sampling simulation class which implements the particle swarm algorithm.</description>
   <descriptorSpecial>descriptorSpecial</descriptorSpecial>
  </posteriorSampleSimulation>
  !!]
  type, extends(posteriorSampleSimulationClass) :: posteriorSampleSimulationParticleSwarm
     !!{
     Implementation of a posterior sampling simulation class which implements the particle swarm algorithm.
     !!}
     private
     type            (modelParameterList                   ), allocatable, dimension(:) :: modelParametersActive_                     , modelParametersInactive_
     class           (posteriorSampleLikelihoodClass       ), pointer                   :: posteriorSampleLikelihood_        => null()
     class           (posteriorSampleConvergenceClass      ), pointer                   :: posteriorSampleConvergence_       => null()
     class           (posteriorSampleStoppingCriterionClass), pointer                   :: posteriorSampleStoppingCriterion_ => null()
     class           (posteriorSampleStateClass            ), pointer                   :: posteriorSampleState_             => null()
     class           (posteriorSampleStateInitializeClass  ), pointer                   :: posteriorSampleStateInitialize_   => null()
     class           (randomNumberGeneratorClass           ), pointer                   :: randomNumberGenerator_            => null()
     integer                                                                            :: parameterCount                             , stepsMaximum                 , &
          &                                                                                reportCount                                , logFlushCount
     double precision                                                                   :: accelerationCoefficientPersonal            , accelerationCoefficientGlobal, &
          &                                                                                inertiaWeight                              , velocityCoefficient          , &
          &                                                                                velocityCoefficientInitial
     logical                                                                            :: isInteractive                              , resume                       , &
          &                                                                                appendLogs
     type            (varying_string                       )                            :: logFileRoot                                , interactionRoot              , &
          &                                                                                logFilePreviousRoot
   contains
     !![
     <methods>
       <method method="posterior"         description="Return the log of posterior probability for the given {\normalfont \ttfamily simulationState}."/>
       <method method="descriptorSpecial" description="Handle adding special parameters to the descriptor."                                           />
     </methods>
     !!]
     final     ::                      particleSwarmDestructor
     procedure :: simulate          => particleSwarmSimulate
     procedure :: posterior         => particleSwarmPosterior
     procedure :: descriptorSpecial => particleSwarmDescriptorSpecial
  end type posteriorSampleSimulationParticleSwarm

  interface posteriorSampleSimulationParticleSwarm
     !!{
     Constructors for the \refClass{posteriorSampleSimulationParticleSwarm} posterior sampling convergence class.
     !!}
     module procedure particleSwarmConstructorParameters
     module procedure particleSwarmConstructorInternal
  end interface posteriorSampleSimulationParticleSwarm

contains

  function particleSwarmConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{posteriorSampleSimulationParticleSwarm} posterior sampling simulation class which builds the object from a
    parameter set.
    !!}
    use :: Display         , only : displayMessage      , displayVerbosity      , verbosityLevelInfo
    use :: Error           , only : Error_Report
    use :: Input_Parameters, only : inputParameter      , inputParameters
    use :: MPI_Utilities   , only : mpiSelf
    use :: Model_Parameters, only : modelParameterActive, modelParameterInactive
    use :: String_Handling , only : operator(//)
    implicit none
    type            (posteriorSampleSimulationParticleSwarm        )                              :: self
    type            (inputParameters                               ), intent(inout)               :: parameters
    type            (modelParameterList                            ), pointer      , dimension(:) :: modelParametersActive_           , modelParametersInactive_
    class           (modelParameterClass                           ), pointer                     :: modelParameter_
    class           (posteriorSampleLikelihoodClass                ), pointer                     :: posteriorSampleLikelihood_
    class           (posteriorSampleConvergenceClass               ), pointer                     :: posteriorSampleConvergence_
    class           (posteriorSampleStoppingCriterionClass         ), pointer                     :: posteriorSampleStoppingCriterion_
    class           (posteriorSampleStateClass                     ), pointer                     :: posteriorSampleState_
    class           (posteriorSampleStateInitializeClass           ), pointer                     :: posteriorSampleStateInitialize_
    class           (randomNumberGeneratorClass                    ), pointer                     :: randomNumberGenerator_
    type            (varying_string                                )                              :: logFileRoot                      , interactionRoot                , &
         &                                                                                           logFilePreviousRoot              , message
    integer                                                                                       :: stepsMaximum                     , reportCount                    , &
         &                                                                                           logFlushCount                    , inactiveParameterCount         , &
         &                                                                                           activeParameterCount             , iActive                        , &
         &                                                                                           iInactive                        , i
    double precision                                                                              :: inertiaWeight                    , accelerationCoefficientPersonal, &
         &                                                                                           accelerationCoefficientGlobal    , velocityCoefficient            , &
         &                                                                                           velocityCoefficientInitial
    logical                                                                                       :: resume                           , appendLogs

    !![
    <inputParameter>
      <name>stepsMaximum</name>
      <defaultValue>huge(0)</defaultValue>
      <description>The maximum number of steps to take.</description>
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
      <name>logFilePreviousRoot</name>
      <defaultValue>var_str('none')</defaultValue>
      <description>Root file name for log files from which to resume.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>resume</name>
      <defaultValue>.false.</defaultValue>
      <description>If true, resume from a previous set of log files.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>appendLogs</name>
      <description>If true, do not overwrite existing log files, but instead append to them.</description>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
    </inputParameter>
    <inputParameter>
      <name>inertiaWeight</name>
      <defaultValue>0.72d0</defaultValue>
      <description>Inertia parameter.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>accelerationCoefficientPersonal</name>
      <defaultValue>1.193d0</defaultValue>
      <description>Personal acceleration parameter.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>accelerationCoefficientGlobal</name>
      <defaultValue>1.193d0</defaultValue>
      <description>Global acceleration parameter.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>velocityCoefficient</name>
      <defaultValue>0.5d0</defaultValue>
      <description>Velocity parameter.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>velocityCoefficientInitial</name>
      <defaultValue>0.0d0</defaultValue>
      <description>Velocity parameter.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="posteriorSampleLikelihood"        name="posteriorSampleLikelihood_"        source="parameters"/>
    <objectBuilder class="posteriorSampleConvergence"       name="posteriorSampleConvergence_"       source="parameters"/>
    <objectBuilder class="posteriorSampleStoppingCriterion" name="posteriorSampleStoppingCriterion_" source="parameters"/>
    <objectBuilder class="posteriorSampleState"             name="posteriorSampleState_"             source="parameters"/>
    <objectBuilder class="posteriorSampleStateInitialize"   name="posteriorSampleStateInitialize_"   source="parameters"/>
    <objectBuilder class="randomNumberGenerator"            name="randomNumberGenerator_"            source="parameters"/>
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
    self=posteriorSampleSimulationParticleSwarm(modelParametersActive_,modelParametersInactive_,posteriorSampleLikelihood_,posteriorSampleConvergence_,posteriorSampleStoppingCriterion_,posteriorSampleState_,posteriorSampleStateInitialize_,randomNumberGenerator_,stepsMaximum,char(logFileRoot),logFlushCount,reportCount,inertiaWeight,accelerationCoefficientPersonal,accelerationCoefficientGlobal,velocityCoefficient,velocityCoefficientInitial,char(interactionRoot),resume,appendLogs,char(logFilePreviousRoot))
    !![
    <inputParametersValidate source="parameters" multiParameters="modelParameter"/>
    <objectDestructor name="posteriorSampleLikelihood_"       />
    <objectDestructor name="posteriorSampleConvergence_"      />
    <objectDestructor name="posteriorSampleStoppingCriterion_"/>
    <objectDestructor name="posteriorSampleState_"            />
    <objectDestructor name="posteriorSampleStateInitialize_"  />
    <objectDestructor name="randomNumberGenerator_"           />
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
  end function particleSwarmConstructorParameters

  function particleSwarmConstructorInternal(modelParametersActive_,modelParametersInactive_,posteriorSampleLikelihood_,posteriorSampleConvergence_,posteriorSampleStoppingCriterion_,posteriorSampleState_,posteriorSampleStateInitialize_,randomNumberGenerator_,stepsMaximum,logFileRoot,logFlushCount,reportCount,inertiaWeight,accelerationCoefficientPersonal,accelerationCoefficientGlobal,velocityCoefficient,velocityCoefficientInitial,interactionRoot,resume,appendLogs,logFilePreviousRoot) result(self)
    !!{
    Internal constructor for the ``particleSwarm'' simulation class.
    !!}
    implicit none
    type            (posteriorSampleSimulationParticleSwarm)                                      :: self
    type            (modelParameterList                    ), intent(in   ), target, dimension(:) :: modelParametersActive_           , modelParametersInactive_
    class           (posteriorSampleLikelihoodClass        ), intent(in   ), target               :: posteriorSampleLikelihood_
    class           (posteriorSampleConvergenceClass       ), intent(in   ), target               :: posteriorSampleConvergence_
    class           (posteriorSampleStoppingCriterionClass ), intent(in   ), target               :: posteriorSampleStoppingCriterion_
    class           (posteriorSampleStateClass             ), intent(in   ), target               :: posteriorSampleState_
    class           (posteriorSampleStateInitializeClass   ), intent(in   ), target               :: posteriorSampleStateInitialize_
    class           (randomNumberGeneratorClass            ), intent(in   ), target               :: randomNumberGenerator_
    character       (len=*                                 ), intent(in   )                       :: logFileRoot                      , interactionRoot                , &
         &                                                                                           logFilePreviousRoot
    integer                                                 , intent(in   )                       :: stepsMaximum                     , reportCount                    , &
         &                                                                                           logFlushCount
    double precision                                        , intent(in   )                       :: inertiaWeight                    , accelerationCoefficientPersonal, &
         &                                                                                           accelerationCoefficientGlobal    , velocityCoefficient            , &
         &                                                                                           velocityCoefficientInitial
    logical                                                 , intent(in   )                       :: resume                           , appendLogs
    integer                                                                                       :: i
    !![
    <constructorAssign variables="*posteriorSampleLikelihood_, *posteriorSampleConvergence_, *posteriorSampleStoppingCriterion_, *posteriorSampleState_, *posteriorSampleStateInitialize_, *randomNumberGenerator_, stepsMaximum, logFileRoot, logFlushCount, reportCount, inertiaWeight, accelerationCoefficientPersonal, accelerationCoefficientGlobal, velocityCoefficient, velocityCoefficientInitial, interactionRoot, resume, appendLogs, logFilePreviousRoot"/>
    !!]

    allocate(self%modelParametersActive_  (size(modelParametersActive_  )))
    allocate(self%modelParametersInactive_(size(modelParametersInactive_)))
    do i=1,size(modelParametersActive_  )
       self%modelParametersActive_  (i)                 =  modelParameterList      ( )
       self%modelParametersActive_  (i)%modelParameter_ => modelParametersActive_  (i)%modelParameter_
       !![
       <referenceCountIncrement owner="self%modelParametersActive_  (i)" object="modelParameter_"/>
       !!]
    end do
    do i=1,size(modelParametersInactive_)
       self%modelParametersInactive_(i)                 =  modelParameterList      ( )
       self%modelParametersInactive_(i)%modelParameter_ => modelParametersInactive_(i)%modelParameter_
       !![
       <referenceCountIncrement owner="self%modelParametersInactive_(i)" object="modelParameter_"/>
       !!]
    end do
    self%parameterCount=size(modelParametersActive_)
    self%isInteractive =trim(interactionRoot ) /= "none"
    call self%posteriorSampleState_%parameterCountSet(self%parameterCount)
    return
  end function particleSwarmConstructorInternal

  subroutine particleSwarmDestructor(self)
    !!{
    Destroy a differential evolution simulation object.
    !!}
    implicit none
    type   (posteriorSampleSimulationParticleSwarm), intent(inout) :: self
    integer                                                        :: i

    !![
    <objectDestructor name="self%posteriorSampleLikelihood_"       />
    <objectDestructor name="self%posteriorSampleConvergence_"      />
    <objectDestructor name="self%posteriorSampleStoppingCriterion_"/>
    <objectDestructor name="self%posteriorSampleState_"            />
    <objectDestructor name="self%posteriorSampleStateInitialize_"  />
    <objectDestructor name="self%randomNumberGenerator_"           />
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
  end subroutine particleSwarmDestructor

  subroutine particleSwarmSimulate(self)
    !!{
    Perform a particle swarm simulation.
    !!}
    use :: Display                     , only : displayIndent  , displayMagenta, displayMessage, displayReset, &
          &                                     displayUnindent
    use :: File_Utilities              , only : File_Exists    , File_Remove
    use :: Error                       , only : Error_Report
    use :: MPI_Utilities               , only : mpiBarrier     , mpiSelf
    use :: Models_Likelihoods_Constants, only : logImpossible
    use :: String_Handling             , only : operator(//)
    implicit none
    class           (posteriorSampleSimulationParticleSwarm), intent(inout)                               :: self
    double precision                                        , dimension(self%parameterCount)              :: stateVector                    , positionMinimum                  , &
         &                                                                                                   positionMaximum                , velocityParticle                 , &
         &                                                                                                   velocityMaximum                , stateBestPersonal                , &
         &                                                                                                   stateBestGlobal                , stateVectorInteractive
    double precision                                        , dimension(:,:)                , allocatable :: stateVectors
    double precision                                        , dimension(:  )                , allocatable :: logPosteriorsAll               , logLikelihoodVariancesAll
    real                                                                                                  :: timePreEvaluate                , timePostEvaluate                 , &
         &                                                                                                   timeEvaluate                   , timeEvaluatePrevious
    double precision                                                                                      :: timeEvaluateInitial            , logPosterior                     , &
         &                                                                                                   logPosteriorBestPersonal       , logPosteriorBestGlobal           , &
         &                                                                                                   logLikelihoodVariance          , logLikelihoodVarianceBestPersonal, &
         &                                                                                                   logLikelihoodVarianceBestGlobal, logLikelihood                    , &
         &                                                                                                   positionStep
    type            (varying_string                        )                                              :: logFileName                    , message                          , &
         &                                                                                                   interactionFileName
    integer                                                                                               :: logFileUnit                    , convergedAtStep                  , &
         &                                                                                                   convergenceFileUnit            , i                                , &
         &                                                                                                   ioStatus                       , interactionFile                  , &
         &                                                                                                   stateCount                     , mpiRank
    logical                                                                                               :: isConverged                    , accept                           , &
         &                                                                                                   forceAcceptance
    character       (len=32                                )                                              :: label

    ! Check that the random number generator is independent across MPI processes.
    if (.not.self%randomNumberGenerator_%mpiIndependent()) call Error_Report('random number generator produces same sequence on all MPI processes'//{introspection:location})
    ! Write start-up message.
    message="Process "//mpiSelf%rankLabel()//" [PID: "
    message=message//getPID()//"] is running on host '"//mpiSelf%hostAffinity()//"'"
    call displayMessage(message)  
    ! Initialize particle to some state vector.
    call self%posteriorSampleStateInitialize_%initialize(self%posteriorSampleState_,self%modelParametersActive_,self%posteriorSampleLikelihood_,timeEvaluateInitial,logLikelihood,logPosterior)
    ! Evaluate the posterior in the initial state.
    forceAcceptance     =.false.
    timeEvaluate        =-1.0
    timeEvaluatePrevious=real(timeEvaluateInitial)
    call CPU_Time(timePreEvaluate )
    call self%posterior(self%posteriorSampleState_,logPosterior,logLikelihood,logLikelihoodVariance,timeEvaluate,timeEvaluatePrevious,forceAcceptance)
    call CPU_Time(timePostEvaluate)
    if (timeEvaluate < 0.0) timeEvaluate=timePostEvaluate-timePreEvaluate
    timeEvaluatePrevious=timeEvaluate    
    ! Set the personal best state to the initial state.
    logPosteriorBestPersonal         =logPosterior
    logLikelihoodVarianceBestPersonal=logLikelihoodVariance
    stateBestPersonal                =self%posteriorSampleState_%get()
    ! Set global best state.
    allocate(stateVectors             (self%parameterCount,mpiSelf%count()))
    allocate(logPosteriorsAll         (                    mpiSelf%count()))
    allocate(logLikelihoodVariancesAll(                    mpiSelf%count()))
    logPosteriorBestGlobal         =logPosterior
    logLikelihoodVarianceBestGlobal=logLikelihoodVariance
    stateBestGlobal                =self%posteriorSampleState_%get()
    logPosteriorsAll               =mpiSelf%gather(logPosteriorBestGlobal         )
    logLikelihoodVariancesAll      =mpiSelf%gather(logLikelihoodVarianceBestGlobal)
    stateVectors                   =mpiSelf%gather(stateBestGlobal                )
    if (mpiSelf%isMaster()) then
       do i=1,mpiSelf%count()
          accept=.false.
          if (i == 1) then
             accept=.true.
          else
             accept=(logPosteriorsAll(i) > logPosteriorBestGlobal)
          end if
          if (accept) then
             stateBestGlobal                =stateVectors             (:,i)
             logPosteriorBestGlobal         =logPosteriorsAll         (  i)
             logLikelihoodVarianceBestGlobal=logLikelihoodVariancesAll(  i)
          end if
       end do
    end if
    call mpiBarrier()
    stateBestGlobal                =reshape(mpiSelf%requestData([0],                 stateBestGlobal ),[self%parameterCount])
    logPosteriorBestGlobal         =sum    (mpiSelf%requestData([0],[         logPosteriorBestGlobal])                      )
    logLikelihoodVarianceBestGlobal=sum    (mpiSelf%requestData([0],[logLikelihoodVarianceBestGlobal])                      )
    ! Compute maximum velocities.
    do i=1,self%parameterCount
       positionMinimum(i)=self%modelParametersActive_(i)%modelParameter_%map(self%modelParametersActive_(i)%modelParameter_%priorMinimum())
       positionMaximum(i)=self%modelParametersActive_(i)%modelParameter_%map(self%modelParametersActive_(i)%modelParameter_%priorMaximum())
       ! Adjust minimum and maximum positions to ensure that, when unmapped, they are within the prior ranges of the
       ! boundaries. (They may not be initially due to finite precision.)
       do while (self%modelParametersActive_(i)%modelParameter_%unmap(positionMinimum(i)) < self%modelParametersActive_(i)%modelParameter_%priorMinimum())
          positionMinimum(i)=nearest(positionMinimum(i),+1.0d0)
       end do
       do while (self%modelParametersActive_(i)%modelParameter_%unmap(positionMaximum(i)) > self%modelParametersActive_(i)%modelParameter_%priorMaximum())
          positionMaximum(i)=nearest(positionMaximum(i),-1.0d0)
       end do
       velocityMaximum(i)=self%velocityCoefficient*(positionMaximum(i)-positionMinimum(i))
    end do
    ! Set initial velocities.
    if (self%resume) then
       ! Simulation is being resumed - retrieve velocities from the previous state files.
       logFileName=self%logFilePreviousRoot//'_'//mpiSelf%rankLabel()//'.log'
       open(newunit=logFileUnit,file=char(logFileName),status='old',form='formatted')
       ioStatus=0
       do while (ioStatus == 0)
          read (logFileUnit,*,iostat=ioStatus) stateCount          , &
               &                               mpiRank             , &
               &                               timeEvaluateInitial , &
               &                               isConverged         , &
               &                               logPosterior        , &
               &                               logLikelihood       , &
               &                               stateVector         , &
               &                               stateVector
          if (ioStatus == 0) velocityParticle=stateVector
       end do
       close(logFileUnit)
    else
       ! Simulation is not being resumed, so set an initial velocity for the particle.
       do i=1,self%parameterCount
          velocityParticle(i)=(2.0d0*self%randomNumberGenerator_%uniformSample()-1.0d0)*velocityMaximum(i)*self%velocityCoefficientInitial
       end do
    end if
    ! Begin the simulation.
    logFileName=self%logFileRoot//'_'//mpiSelf%rankLabel()//'.log'
    if (self%appendLogs) then
       open(newunit=logFileUnit,file=char(logFileName),status='unknown',form='formatted',position='append')
    else
       open(newunit=logFileUnit,file=char(logFileName),status='unknown',form='formatted')
    end if
    isConverged=.false.
    do while (                                                                                                   &
         &          self%posteriorSampleState_            %count(                          ) < self%stepsMaximum &
         &    .and.                                                                                              &
         &     .not.self%posteriorSampleStoppingCriterion_%stop (self%posteriorSampleState_)                     &
         &   )
       ! Get the current particle state.
       stateVector=self%posteriorSampleState_%get()
       ! Update the state vector
       do i=1,self%parameterCount
          positionStep=velocityParticle(i)
          do while (positionStep /= 0.0d0)
             stateVector(i)=stateVector(i)+positionStep
             if (stateVector(i) <= positionMinimum(i)) then
                positionStep       =+positionMinimum (i)-stateVector(i)
                stateVector     (i)=+positionMinimum (i)
                velocityParticle(i)=-velocityParticle(i)
             else if (stateVector(i) >= positionMaximum(i)) then
                positionStep       =+positionMaximum (i)-stateVector(i)
                stateVector     (i)=+positionMaximum (i)
                velocityParticle(i)=-velocityParticle(i)
             else
                positionStep       =+0.0d0
             end if
          end do
       end do
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
                   stateVector=stateVectorInteractive
                   message="Chain "//mpiSelf%rankLabel()//" is being interactively moved to state:"
                   call displayIndent(message)
                   ! Map parameters of interactively proposed state.
                   do i=1,size(stateVector)
                      write (label,*) stateVector(i)
                      message="State["
                      message=message//i//"] = "//trim(adjustl(label))
                      call displayMessage(message)
                      stateVector(i)=self%modelParametersActive_(i)%modelParameter_%map(stateVector(i))
                   end do
                   call displayUnindent('end')
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
       ! Store the state vector.
       call self%posteriorSampleState_%update(stateVector,.true.,self%posteriorSampleConvergence_%isConverged())
       ! Update the velocity.
       velocityParticle=+self%inertiaWeight                                   &
            &           *velocityParticle                                     &
            &           +self%accelerationCoefficientPersonal                 &
            &           *self%randomNumberGenerator_         %uniformSample() &
            &           *(                                                    &
            &             +stateBestPersonal                                  &
            &             -stateVector                                        &
            &            )                                                    &
            &           +self%accelerationCoefficientGlobal                   &
            &           *self%randomNumberGenerator_         %uniformSample() &
            &           *(                                                    &
            &             +stateBestGlobal                                    &
            &             -stateVector                                        &
            &            )
       where (velocityParticle > velocityMaximum)
          velocityParticle=+velocityMaximum
       end where
       where (velocityParticle < -velocityMaximum)
          velocityParticle=-velocityMaximum
       end where
       ! Evaluate posterior.
       timeEvaluatePrevious=timeEvaluate
       timeEvaluate        =-1.0
       forceAcceptance     =.false.
       call CPU_Time(timePreEvaluate )
       call self%posterior(self%posteriorSampleState_,logPosterior,logLikelihood,logLikelihoodVariance,timeEvaluate,timeEvaluatePrevious,forceAcceptance)
       call CPU_Time(timePostEvaluate)
       if (timeEvaluate < 0.0) timeEvaluate=timePostEvaluate-timePreEvaluate
       ! Update personal best state.
       accept= forceAcceptance                         &
            & .or.                                     &
            &  logPosterior > logPosteriorBestPersonal
       if (accept) then
          logPosteriorBestPersonal         =logPosterior
          logLikelihoodVarianceBestPersonal=logLikelihoodVariance
          stateBestPersonal                =stateVector
       end if
       ! Update global best state.
       logPosteriorsAll         =mpiSelf%gather(logPosterior         )
       logLikelihoodVariancesAll=mpiSelf%gather(logLikelihoodVariance)
       stateVectors             =mpiSelf%gather(stateVector          )
       if (mpiSelf%isMaster()) then
          do i=1,mpiSelf%count()
             accept=(logPosteriorsAll(i) > logPosteriorBestGlobal)
             if (accept) then
                stateBestGlobal                =stateVectors             (:,i)
                logPosteriorBestGlobal         =logPosteriorsAll         (  i)
                logLikelihoodVarianceBestGlobal=logLikelihoodVariancesAll(  i)
             end if
          end do
       end if
       call mpiBarrier()
       stateBestGlobal                =reshape(mpiSelf%requestData([0],                 stateBestGlobal ),[self%parameterCount])
       logPosteriorBestGlobal         =sum    (mpiSelf%requestData([0],[         logPosteriorBestGlobal])                      )
       logLikelihoodVarianceBestGlobal=sum    (mpiSelf%requestData([0],[logLikelihoodVarianceBestGlobal])                      )
       ! Unmap parameters and write to log file.
       do i=1,size(stateVector)
          stateVector(i)=self%modelParametersActive_(i)%modelParameter_%unmap(stateVector(i))
       end do
       write (logFileUnit,*) self   %posteriorSampleState_%count(), &
            &                mpiSelf%rank                 (), &
            &                timeEvaluate                   , &
            &                isConverged                    , &
            &                logPosterior                   , &
            &                logLikelihood                  , &
            &                stateVector                    , &
            &                velocityParticle
       if (mod(self%posteriorSampleState_%count(),self%logFlushCount) == 0) call flush(logFileUnit)
       ! Repeat.
       call mpiBarrier()
       ! Test for convergence.
       if (.not.isConverged) then
          isConverged=self%posteriorSampleConvergence_%isConverged(self%posteriorSampleState_,logPosterior)
          if (isConverged) then
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
  end subroutine particleSwarmSimulate

  subroutine particleSwarmPosterior(self,posteriorSampleState_,logPosterior,logLikelihood,logLikelihoodVariance,timeEvaluate,timeEvaluatePrevious,forceAcceptance)
    !!{
    Return the log of the posterior for the current state.
    !!}
    use            :: Display                     , only : displayIndent             , displayMessage, displayUnindent
    use            :: Error                       , only : Error_Report
    use, intrinsic :: ISO_C_Binding               , only : c_size_t
    use            :: Kind_Numbers                , only : kind_int4
    use            :: MPI_Utilities               , only : mpiBarrier                , mpiSelf
    use            :: Model_Parameters            , only : modelParameterListLogPrior
    use            :: Models_Likelihoods_Constants, only : logImpossible
    use            :: Sorting                     , only : sortIndex
    implicit none
    class           (posteriorSampleSimulationParticleSwarm), intent(inout)               :: self
    class           (posteriorSampleStateClass             ), intent(inout)               :: posteriorSampleState_
    double precision                                        , intent(  out)               :: logPosterior              , logLikelihoodVariance   , &
         &                                                                                   logLikelihood
    real                                                    , intent(inout)               :: timeEvaluate
    real                                                    , intent(in   )               :: timeEvaluatePrevious
    logical                                                 , intent(inout)               :: forceAcceptance
    double precision                                        , dimension(:  ), allocatable :: timesEvaluate             , nodeWork                , &
         &                                                                                   stateVectorSelf           , timesEvaluateActual
    double precision                                        , dimension(:,:), allocatable :: stateVectorWork
    integer                                                 , dimension(:  ), allocatable :: processToProcess          , processFromProcess
    integer         (c_size_t                              ), dimension(:  ), allocatable :: timesEvaluateOrder        , nodeWorkOrder
    integer                                                 , dimension(1  )              :: particleIndexSelf
    integer                                                 , dimension(1,1)              :: particleIndexWork
    double precision                                        , dimension(1,1)              :: logLikelihoodSelf         , logPriorWork            , &
         &                                                                                   timeEvaluateSelf
    double precision                                        , dimension(1  )              :: logLikelihoodWork         , logPriorSelf            , &
         &                                                                                   timeEvaluateWork
    double precision                                        , parameter                   :: temperature         =1.0d0
    double precision                                                                      :: logPrior                  , timeEvaluateEffective
    integer                                                                               :: i                         , processTrial            , &
         &                                                                                   nodeTrial
    type            (varying_string                        )                              :: message
    character       (len=10                                )                              :: label

    ! Evaluate the proposed prior.
    logPrior=modelParameterListLogPrior(self%modelParametersActive_,posteriorSampleState_)
    ! Gather timing data from all particles.
    allocate(timesEvaluate      (0:mpiSelf%count()-1))
    allocate(timesEvaluateActual(0:mpiSelf%count()-1))
    allocate(processToProcess   (0:mpiSelf%count()-1))
    allocate(processFromProcess (0:mpiSelf%count()-1))
    allocate(timesEvaluateOrder (0:mpiSelf%count()-1))
    timeEvaluateEffective=timeEvaluatePrevious
    if     (                                                                                    &
         &  .not.self%posteriorSampleLikelihood_%willEvaluate(                                  &
         &                                                    posteriorSampleState_           , &
         &                                                    self%modelParametersActive_     , &
         &                                                    self%posteriorSampleConvergence_, &
         &                                                    temperature                     , &
         &                                                    logImpossible                   , &
         &                                                    logImpossible                   , &
         &                                                    logPrior                          &
         &                                                   )                                  &
         & )                                                                                    &
         & timeEvaluateEffective=0.0d0
    timesEvaluate=mpiSelf%gather(dble(timeEvaluateEffective))
    ! If previous time estimate is negative, don't do load balancing.
    if (mpiSelf%isMaster() .and. mod(self%posteriorSampleState_%count(),self%reportCount) == 0) call displayIndent('Load balancing report')
    if (any(timesEvaluate < 0.0d0)) then
       forall(i=0:mpiSelf%count()-1)
          processToProcess  (i)=i
          processFromProcess(i)=i
       end forall
       if (mpiSelf%isMaster() .and. mod(self%posteriorSampleState_%count(),self%reportCount) == 0) call displayMessage('Not performing load balancing - missing work cost data')
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
             if (processToProcess(timesevaluateorder(i)) >= 0) exit
          end do
          if (processToProcess(timesevaluateorder(i)) < 0) call Error_Report('failed to assign task to process'//{introspection:location})
       end do
       ! Report.
       if (mpiSelf%isMaster() .and. mod(self%posteriorSampleState_%count(),self%reportCount) == 0) then
          call displayIndent('Particle redistribution:')
          do i=0,mpiSelf%count()-1
             write (label,'(i4.4)') i
             message='Particle '//trim(label)//' -> process/node '
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
    ! Get state vector, particle index, and prior.
    allocate(stateVectorSelf(self%parameterCount  ))
    allocate(stateVectorWork(self%parameterCount,1))
    stateVectorSelf         =posteriorSampleState_%get       ()
    particleIndexSelf       =posteriorSampleState_%chainIndex()
    logPriorSelf            =logPrior
    stateVectorWork         =mpiSelf%requestData(processFromProcess(mpiSelf%rank():mpiSelf%rank()),stateVectorSelf  )
    particleIndexWork       =mpiSelf%requestData(processFromProcess(mpiSelf%rank():mpiSelf%rank()),particleIndexSelf)
    logPriorWork            =mpiSelf%requestData(processFromProcess(mpiSelf%rank():mpiSelf%rank()),logPriorSelf     )
    ! Set state and particle index.
    call posteriorSampleState_%update       (  stateVectorWork(:,1),logState=.false.,isConverged=.false.)
    call posteriorSampleState_%chainIndexSet(particleIndexWork(1,1)                                     )
    ! Evaluate the likelihood.
    logLikelihood=self%posteriorSampleLikelihood_%evaluate(                                                             &
         &                                                                       posteriorSampleState_                , &
         &                                                                       self%modelParametersActive_          , &
         &                                                                       self%modelParametersInactive_        , &
         &                                                                       self%posteriorSampleConvergence_     , &
         &                                                                       temperature                          , &
         &                                                                       logImpossible                        , &
         &                                                                       logImpossible                        , &
         &                                                                       logPriorWork                    (1,1), &
         &                                                                       timeEvaluate                         , &
         &                                                 logLikelihoodVariance=logLikelihoodVariance                , &
         &                                                 forceAcceptance      =forceAcceptance                        &
         &                                                )
    call mpiBarrier()
    ! Distribute likelihoods back to origins.
    logLikelihoodWork    =                                                                         logLikelihood
    logLikelihoodSelf    =mpiSelf%requestData(processToProcess(mpiSelf%rank():mpiSelf%rank()),     logLikelihoodWork          )
    logLikelihood        =                                                                         logLikelihoodSelf    (1,1)
    ! Distribute likelihood variances back to origins.
    logLikelihoodWork    =                                                                         logLikelihoodVariance
    logLikelihoodSelf    =mpiSelf%requestData(processToProcess(mpiSelf%rank():mpiSelf%rank()),     logLikelihoodWork          )
    logLikelihoodVariance=                                                                         logLikelihoodSelf    (1,1)
    ! Distribute evaluation times back to origins.
    timeEvaluateWork     =                                                                    dble(timeEvaluate              )
    timeEvaluateSelf     =mpiSelf%requestData(processToProcess(mpiSelf%rank():mpiSelf%rank()),     timeEvaluateWork           )
    timeEvaluate         =                                                                    real(timeEvaluateSelf     (1,1))
    ! Restore state and particle index.
    call posteriorSampleState_%update       (  stateVectorSelf   ,logState=.false.,isConverged=.false.)
    call posteriorSampleState_%chainIndexSet(particleIndexSelf(1)                                     )
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
  end subroutine particleSwarmPosterior

  subroutine particleSwarmDescriptorSpecial(self,descriptor)
    !!{
    Add special parameters to the descriptor.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    class  (posteriorSampleSimulationParticleSwarm), intent(inout) :: self
    type   (inputParameters                       ), intent(inout) :: descriptor
    integer                                                        :: i
    
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
  end subroutine particleSwarmDescriptorSpecial
