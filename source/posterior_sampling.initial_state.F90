!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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
  Implementation of a posterior sampling simulation class which implements evaluation of the model likelihood in the initial state
  only.
  !!}

  use :: Model_Parameters                   , only : modelParameterList
  use :: Models_Likelihoods                 , only : posteriorSampleLikelihoodClass
  use :: Posterior_Sampling_State           , only : posteriorSampleStateClass
  use :: Posterior_Sampling_State_Initialize, only : posteriorSampleStateInitializeClass

  !![
  <posteriorSampleSimulation name="posteriorSampleSimulationInitialState">
   <description>A posterior sampling simulation class which implements evaluation of the model likelihood in the initial state only.</description>
   <descriptorSpecial>descriptorSpecial</descriptorSpecial>
  </posteriorSampleSimulation>
  !!]
  type, extends(posteriorSampleSimulationClass) :: posteriorSampleSimulationInitialState
     !!{
     Implementation of a posterior sampling simulation class which implements evaluation of the model likelihood in the initial state only.
     !!}
     private
     type   (modelParameterList                 ), pointer, dimension(:) :: modelParametersActive_          => null(), modelParametersInactive_ => null()
     class  (posteriorSampleLikelihoodClass     ), pointer               :: posteriorSampleLikelihood_      => null()
     class  (posteriorSampleStateClass          ), pointer               :: posteriorSampleState_           => null()
     class  (posteriorSampleStateInitializeClass), pointer               :: posteriorSampleStateInitialize_ => null()
     integer                                                             :: parameterCount
     logical                                                             :: appendLogs
     type   (varying_string                     )                        :: logFileRoot
   contains
     !![
     <methods>
       <method method="posterior"         description="Return the log of posterior probability for the given {\normalfont \ttfamily simulationState}."/>
       <method method="descriptorSpecial" description="Handle adding special parameters to the descriptor."                                           />
     </methods>
     !!]
     final     ::                      initialStateDestructor
     procedure :: simulate          => initialStateSimulate
     procedure :: posterior         => initialStatePosterior
     procedure :: descriptorSpecial => initialStateDescriptorSpecial
  end type posteriorSampleSimulationInitialState

  interface posteriorSampleSimulationInitialState
     !!{
     Constructors for the {\normalfont \ttfamily initialState} posterior sampling convergence class.
     !!}
     module procedure initialStateConstructorParameters
     module procedure initialStateConstructorInternal
  end interface posteriorSampleSimulationInitialState

contains

  function initialStateConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily initialState} posterior sampling simulation class which builds the object from a
    parameter set.
    !!}
    use :: Display         , only : displayMessage      , displayVerbosity      , verbosityLevelInfo
    use :: Error           , only : Error_Report
    use :: Input_Parameters, only : inputParameter      , inputParameters
    use :: MPI_Utilities   , only : mpiSelf
    use :: Model_Parameters, only : modelParameterActive, modelParameterInactive
    use :: String_Handling , only : operator(//)
    implicit none
    type   (posteriorSampleSimulationInitialState)                              :: self
    type   (inputParameters                     ), intent(inout)               :: parameters
    type   (modelParameterList                  ), pointer      , dimension(:) :: modelParametersActive_         , modelParametersInactive_
    class  (modelParameterClass                 ), pointer                     :: modelParameter_
    class  (posteriorSampleLikelihoodClass      ), pointer                     :: posteriorSampleLikelihood_
    class  (posteriorSampleStateClass           ), pointer                     :: posteriorSampleState_
    class  (posteriorSampleStateInitializeClass ), pointer                     :: posteriorSampleStateInitialize_
    type   (varying_string                      )                              :: logFileRoot                    , message
    integer                                                                    :: inactiveParameterCount         , activeParameterCount    , &
         &                                                                        iInactive                      , iActive                 , &
         &                                                                        i
    logical                                                                    :: appendLogs

    !![
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
    <objectBuilder class="posteriorSampleLikelihood"      name="posteriorSampleLikelihood_"      source="parameters"/>
    <objectBuilder class="posteriorSampleState"           name="posteriorSampleState_"           source="parameters"/>
    <objectBuilder class="posteriorSampleStateInitialize" name="posteriorSampleStateInitialize_" source="parameters"/>
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
    ! Initialize model parameters.
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
    self=posteriorSampleSimulationInitialState(modelParametersActive_,modelParametersInactive_,posteriorSampleLikelihood_,posteriorSampleState_,posteriorSampleStateInitialize_,char(logFileRoot),appendLogs)
    !![
    <inputParametersValidate source="parameters" multiParameters="modelParameter"/>
    <objectDestructor name="posteriorSampleLikelihood_"     />
    <objectDestructor name="posteriorSampleState_"          />
    <objectDestructor name="posteriorSampleStateInitialize_"/>
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
  end function initialStateConstructorParameters

  function initialStateConstructorInternal(modelParametersActive_,modelParametersInactive_,posteriorSampleLikelihood_,posteriorSampleState_,posteriorSampleStateInitialize_,logFileRoot,appendLogs) result(self)
    !!{
    Internal constructor for the ``initialState'' simulation class.
    !!}
    implicit none
    type     (posteriorSampleSimulationInitialState)                                      :: self
    type     (modelParameterList                  ), intent(in   ), target, dimension(:) :: modelParametersActive_    , modelParametersInactive_
    class    (posteriorSampleLikelihoodClass      ), intent(in   ), target               :: posteriorSampleLikelihood_
    class    (posteriorSampleStateClass           ), intent(in   ), target               :: posteriorSampleState_
    class    (posteriorSampleStateInitializeClass ), intent(in   ), target               :: posteriorSampleStateInitialize_
    character(len=*                               ), intent(in   )                       :: logFileRoot
    logical                                        , intent(in   )                       :: appendLogs
    integer                                                                              :: i
    !![
    <constructorAssign variables="*posteriorSampleLikelihood_, *posteriorSampleState_, *posteriorSampleStateInitialize_, logFileRoot, appendLogs"/>
    !!]

    allocate(self%modelParametersActive_  (size(modelParametersActive_  )))
    allocate(self%modelParametersInactive_(size(modelParametersInactive_)))
    self%modelParametersActive_  =modelParametersActive_
    self%modelParametersInactive_=modelParametersInactive_
    do i=1,size(modelParametersActive_  )
       !![
       <referenceCountIncrement owner="self%modelParametersActive_  (i)" object="modelParameter_"/>
       !!]
    end do
    do i=1,size(modelParametersInactive_)
       !![
       <referenceCountIncrement owner="self%modelParametersInactive_(i)" object="modelParameter_"/>
       !!]
    end do
    self%parameterCount=size(modelParametersActive_)
    call self%posteriorSampleState_%parameterCountSet(self%parameterCount)
    return
  end function initialStateConstructorInternal

  subroutine initialStateDestructor(self)
    !!{
    Destroy a differential evolution simulation object.
    !!}
    implicit none
    type   (posteriorSampleSimulationInitialState), intent(inout) :: self
    integer                                                      :: i
 
    !![
    <objectDestructor name="self%posteriorSampleLikelihood_"     />
    <objectDestructor name="self%posteriorSampleState_"          />
    <objectDestructor name="self%posteriorSampleStateInitialize_"/>
    !!]
    if (associated(self%modelParametersActive_  )) then
       do i=1,size(self%modelParametersActive_  )
          !![
	  <objectDestructor name="self%modelParametersActive_  (i)%modelParameter_"/>
          !!]
       end do
    end if
    if (associated(self%modelParametersInactive_)) then
       do i=1,size(self%modelParametersInactive_)
          !![
	  <objectDestructor name="self%modelParametersInactive_(i)%modelParameter_"/>
          !!]
       end do
    end if
    return
  end subroutine initialStateDestructor

  subroutine initialStateSimulate(self)
    !!{
    Perform a initialState simulation.
    !!}
    use :: Display        , only : displayIndent  , displayMagenta, displayMessage, displayReset, &
          &                        displayUnindent
    use :: Error          , only : Error_Report
    use :: MPI_Utilities  , only : mpiBarrier     , mpiSelf
    use :: String_Handling, only : operator(//)
    implicit none
    class           (posteriorSampleSimulationInitialState), intent(inout)                  :: self
    double precision                                      , dimension(self%parameterCount) :: stateVector
    real                                                                                   :: timeEvaluate
    double precision                                                                       :: logPosterior       , logLikelihood, &
         &                                                                                    timeEvaluateInitial
    type            (varying_string                      )                                 :: logFileName        , message 
    integer                                                                                :: logFileUnit
    integer         (c_size_t                            )                                 :: i
    
    ! Write start-up message.
    message="Process "//mpiSelf%rankLabel()//" [PID: "
    message=message//getPID()//"] is running on host '"//mpiSelf%hostAffinity()//"'"
    call displayMessage(message)
    ! Begin the simulation.
    logFileName=self%logFileRoot//'_'//mpiSelf%rankLabel()//'.log'
    if (self%appendLogs) then
       open(newunit=logFileUnit,file=char(logFileName),status='unknown',form='formatted',position='append')
    else
       open(newunit=logFileUnit,file=char(logFileName),status='unknown',form='formatted'                  )
    end if
    ! Initialize chain to some state vector.
    call self%posteriorSampleStateInitialize_%initialize(self%posteriorSampleState_,self%modelParametersActive_,self%posteriorSampleLikelihood_,timeEvaluateInitial,logLikelihood,logPosterior)
    ! Evaluate the posterior.
    call self%posterior(self%posteriorSampleState_,logPosterior,logLikelihood,timeEvaluate)
    ! Unmap parameters and write to log file.
    do i=1,size(stateVector)
       stateVector(i)=self%modelParametersActive_(i)%modelParameter_%unmap(stateVector(i))
    end do
    write (logFileUnit,*) self   %posteriorSampleState_%count(), &
         &                mpiSelf%rank                       (), &
         &                timeEvaluate                         , &
         &                .true.                               , &
         &                logPosterior                         , &
         &                logLikelihood                        , &
         &                stateVector
    call mpiBarrier()
    close(logFileUnit)
    return
  end subroutine initialStateSimulate

  subroutine initialStatePosterior(self,posteriorSampleState_,logPosterior,logLikelihood,timeEvaluate)
    !!{
    Return the log of the posterior for the current state.
    !!}
    use :: Model_Parameters              , only : modelParameterListLogPrior
    use :: Models_Likelihoods_Constants  , only : logImpossible
    use :: Posterior_Sampling_Convergence, only : posteriorSampleConvergenceAlways
    implicit none
    class           (posteriorSampleSimulationInitialState), intent(inout) :: self
    class           (posteriorSampleStateClass           ), intent(inout) :: posteriorSampleState_
    double precision                                      , intent(  out) :: logPosterior                     , logLikelihood
    real                                                  , intent(inout) :: timeEvaluate
    double precision                                      , parameter     :: temperature                =1.0d0
    double precision                                                      :: logPrior
    type            (posteriorSampleConvergenceAlways    )                :: posteriorSampleConvergence_
    
    ! Evaluate the proposed prior.
    logPrior                   =modelParameterListLogPrior      (self%modelParametersActive_,posteriorSampleState_)
    ! Assume always converged.
    posteriorSampleConvergence_=posteriorSampleConvergenceAlways(                                                 )
    ! Evaluate the likelihood.
    logLikelihood              =self%posteriorSampleLikelihood_%evaluate(                               &
         &                                                               posteriorSampleState_        , &
         &                                                               self%modelParametersActive_  , &
         &                                                               self%modelParametersInactive_, &
         &                                                               posteriorSampleConvergence_  , &
         &                                                               temperature                  , &
         &                                                               logImpossible                , &
         &                                                               logImpossible                , &
         &                                                               logPrior                     , &
         &                                                               timeEvaluate                   &
         &                                                              )
    logPosterior               =+logPrior      &
         &                      +logLikelihood
    return
  end subroutine initialStatePosterior

  subroutine initialStateDescriptorSpecial(self,descriptor)
    !!{
    Add special parameters to the descriptor.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    class  (posteriorSampleSimulationInitialState), intent(inout) :: self
    type   (inputParameters                     ), intent(inout) :: descriptor
    integer                                                      :: i
    
    if (associated(self%modelParametersActive_  )) then
       do i=1,size(self%modelParametersActive_  )
          call self%modelParametersActive_  (i)%modelParameter_%descriptor(descriptor)
       end do
    end if
    if (associated(self%modelParametersInactive_)) then
       do i=1,size(self%modelParametersInactive_)
          call self%modelParametersInactive_(i)%modelParameter_%descriptor(descriptor)
       end do
    end if
    return
  end subroutine initialStateDescriptorSpecial
