!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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

  !% Implementation of a posterior sampling likelihood class which implements a likelihood for \glc\ models.

  use :: FoX_DOM         , only : node
  use :: Input_Parameters, only : inputParameter, inputParameters
  use :: Output_Analyses , only : outputAnalysis, outputAnalysisClass

  type :: parameterList
     !% Type used to maintain a list of pointers to parameters to be modified.
     type   (inputParameter), pointer :: parameter_
     integer                          :: indexElement
     type   (varying_string)          :: definition
     logical                          :: resolved
  end type parameterList

  !# <posteriorSampleLikelihood name="posteriorSampleLikelihoodGalaxyPopulation">
  !#  <description>A posterior sampling likelihood class which implements a likelihood for \glc\ models.</description>
  !# </posteriorSampleLikelihood>
  type, extends(posteriorSampleLikelihoodClass) :: posteriorSampleLikelihoodGalaxyPopulation
     !% Implementation of a posterior sampling likelihood class which implements a likelihood for \glc\ models.
     private
     type   (varying_string     )                            :: baseParametersFileName          , failedParametersFileName
     logical                                                 :: randomize
     integer                                                 :: evolveForestsVerbosity
     type   (inputParameters    ), pointer                   :: parametersModel        => null()
     class  (*                  ), pointer                   :: task_
     class  (outputAnalysisClass), pointer                   :: outputAnalysis_        => null()
     type   (parameterList      ), dimension(:), allocatable :: modelParametersActive_          , modelParametersInactive_
   contains
     final     ::                    galaxyPopulationDestructor
     procedure :: evaluate        => galaxyPopulationEvaluate
     procedure :: functionChanged => galaxyPopulationFunctionChanged
     procedure :: willEvaluate    => galaxyPopulationWillEvaluate
  end type posteriorSampleLikelihoodGalaxyPopulation

  interface posteriorSampleLikelihoodGalaxyPopulation
     !% Constructors for the {\normalfont \ttfamily galaxyPopulation} posterior sampling likelihood class.
     module procedure galaxyPopulationConstructorParameters
     module procedure galaxyPopulationConstructorInternal
  end interface posteriorSampleLikelihoodGalaxyPopulation

contains

  function galaxyPopulationConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily galaxyPopulation} posterior sampling likelihood class which builds the object
    !% from a parameter set.
    use :: Display         , only : displayVerbosity
    use :: Input_Parameters, only : inputParameter  , inputParameters
    implicit none
    type   (posteriorSampleLikelihoodGalaxyPopulation)                :: self
    type   (inputParameters                          ), intent(inout) :: parameters
    type   (varying_string)                                           :: baseParametersFileName, failedParametersFileName
    logical                                                           :: randomize
    integer                                                           :: evolveForestsVerbosity
    type   (inputParameters                          ), pointer       :: parametersModel

    !# <inputParameter>
    !#   <name>baseParametersFileName</name>
    !#   <description>The base set of parameters to use.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>randomize</name>
    !#   <description>If true, randomize models (i.e. change the random seed).</description>
    !#   <defaultValue>.false.</defaultValue>
    !#   <source>parameters</source>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>evolveForestsVerbosity</name>
    !#   <description>The verbosity level to use while performing evolve forests tasks.</description>
    !#   <defaultValue>displayVerbosity()</defaultValue>
    !#   <source>parameters</source>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>failedParametersFileName</name>
    !#   <description>The filename to which parameters of failed models should be written.</description>
    !#   <defaultValue>var_str('./failedParameters.xml')</defaultValue>
    !#   <source>parameters</source>
    !# </inputParameter>
    allocate(parametersModel)
    parametersModel=inputParameters                          (baseParametersFileName,noOutput=.true.)
    self           =posteriorSampleLikelihoodGalaxyPopulation(parametersModel,randomize,evolveForestsVerbosity,failedParametersFileName)
    !# <inputParametersValidate source="parameters"/>
    nullify(parametersModel)
    return
  end function galaxyPopulationConstructorParameters

  function galaxyPopulationConstructorInternal(parametersModel,randomize,evolveForestsVerbosity,failedParametersFileName) result(self)
    !% Constructor for ``galaxyPopulation'' posterior sampling likelihood class.
    implicit none
    type   (posteriorSampleLikelihoodGalaxyPopulation)                        :: self
    type   (inputParameters                          ), intent(inout), target :: parametersModel
    logical                                           , intent(in   )         :: randomize
    integer                                           , intent(in   )         :: evolveForestsVerbosity
    type   (varying_string                           ), intent(in   )         :: failedParametersFileName
    !# <constructorAssign variables="*parametersModel, randomize, evolveForestsVerbosity, failedParametersFileName"/>

    return
  end function galaxyPopulationConstructorInternal

  subroutine galaxyPopulationDestructor(self)
    !% Destructor for the {\normalfont \ttfamily galaxyPopulation} posterior sampling likelihood class.
    implicit none
    type(posteriorSampleLikelihoodGalaxyPopulation), intent(inout) :: self

    call self%parametersModel%destroy()
    deallocate(self%parametersModel)
    return
  end subroutine galaxyPopulationDestructor

  double precision function galaxyPopulationEvaluate(self,simulationState,modelParametersActive_,modelParametersInactive_,simulationConvergence,temperature,logLikelihoodCurrent,logPriorCurrent,logPriorProposed,timeEvaluate,logLikelihoodVariance,forceAcceptance)
    !% Return the log-likelihood for the \glc\ likelihood function.
    use :: Display                       , only : displayIndent                  , displayMessage               , displayUnindent             , displayVerbosity, &
          &                                       displayVerbositySet            , verbosityLevelSilent         , verbosityLevelStandard
    use :: Functions_Global              , only : Tasks_Evolve_Forest_Construct_ , Tasks_Evolve_Forest_Destruct_, Tasks_Evolve_Forest_Perform_
    use :: Galacticus_Error              , only : Galacticus_Error_Report        , errorStatusSuccess
    use :: ISO_Varying_String            , only : char                           , operator(//)                 , var_str
    use :: Kind_Numbers                  , only : kind_int8
    use :: MPI_Utilities                 , only : mpiBarrier                     , mpiSelf
    use :: Model_Parameters              , only : modelParameterDerived
    use :: Models_Likelihoods_Constants  , only : logImpossible                  , logImprobable
    use :: Posterior_Sampling_Convergence, only : posteriorSampleConvergenceClass
    use :: Posterior_Sampling_State      , only : posteriorSampleStateClass
    use :: String_Handling               , only : String_Count_Words             , String_Join                  , String_Split_Words          , operator(//)
    implicit none
    class           (posteriorSampleLikelihoodGalaxyPopulation), intent(inout)                 :: self
    class           (posteriorSampleStateClass                ), intent(inout)                 :: simulationState
    type            (modelParameterList                       ), intent(in   ), dimension(:  ) :: modelParametersActive_, modelParametersInactive_
    class           (posteriorSampleConvergenceClass          ), intent(inout)                 :: simulationConvergence
    double precision                                           , intent(in   )                 :: temperature           , logLikelihoodCurrent    , &
         &                                                                                        logPriorCurrent       , logPriorProposed
    real                                                       , intent(inout)                 :: timeEvaluate
    double precision                                           , intent(  out), optional       :: logLikelihoodVariance
    logical                                                    , intent(inout), optional       :: forceAcceptance
    double precision                                           , allocatable  , dimension(:  ) :: logPriorsProposed
    double precision                                           , allocatable  , dimension(:,:) :: stateVector
    type            (varying_string                           ), allocatable  , dimension(:  ) :: parameterNames
    integer                                                                                    :: iRank                 , parameterCount          , &
         &                                                                                        i                     , j                       , &
         &                                                                                        instance              , indexElement            , &
         &                                                                                        verbosityLevel        , status
#ifdef MATHEVALAVAIL
    integer         (kind_int8                                )                                :: evaluator
#endif
    real                                                                                       :: timeBegin             , timeEnd
    double precision                                                                           :: logLikelihoodProposed , valueDerived
    type            (inputParameters                          )                                :: parameters_
    character       (len=10                                   )                                :: labelIndex
    character       (len=24                                   )                                :: valueText
    logical                                                                                    :: firstIteration        , dependenciesResolved    , &
         &                                                                                        dependenciesUpdated
    type            (varying_string                           )                                :: parameterText         , message
    ! Declarations of GNU libmatheval procedures used.
    integer         (kind_int8                                ), external                      :: Evaluator_Create_
    double precision                                           , external                      :: Evaluator_Evaluate_
    external                                                                                   :: Evaluator_Destroy_
    !$GLC attributes unused :: logPriorCurrent, logLikelihoodCurrent, forceAcceptance, temperature, simulationConvergence

    ! Switch verbosity level.
    verbosityLevel=displayVerbosity()
    call displayVerbositySet(self%evolveForestsVerbosity)
    ! Initialize likelihood to impossible.
    galaxyPopulationEvaluate=logImpossible
    if (present(logLikelihoodVariance)) logLikelihoodVariance=0.0d0
    ! Get proposed priors for all chains so we can decide which to skip.
    allocate(logPriorsProposed(0:mpiSelf%count()-1))
    logPriorsProposed=mpiSelf%gather(logPriorProposed)
    ! Get states for all chains.
    allocate(stateVector(simulationState%dimension(),0:mpiSelf%count()-1))
    stateVector=mpiSelf%gather(simulationState%get())
    ! On first call we must build pointers to all parameter nodes which will be modified as a function of chain state.
    if (.not.allocated(self%modelParametersActive_)) then
       allocate(self%modelParametersActive_(size(modelParametersActive_)))
       do i=1,size(modelParametersActive_)
          parameterCount=String_Count_Words(char(modelParametersActive_(i)%modelParameter_%name()),"::")
          allocate(parameterNames(parameterCount))
          call String_Split_Words(parameterNames,char(modelParametersActive_(i)%modelParameter_%name()),"::")
          parameters_=self%parametersModel
          do j=1,parameterCount
             instance    =1
             indexElement=0
             if (index(parameterNames(j),"[") /= 0) then
                labelIndex       =extract(parameterNames(j),index(parameterNames(j),"[")+1,index(parameterNames(j),"]")-1)
                parameterNames(j)=extract(parameterNames(j),                             1,index(parameterNames(j),"[")-1)
                read (labelIndex,*) instance
                instance=instance+1
             else if (index(parameterNames(j),"{") /= 0) then
                labelIndex       =extract(parameterNames(j),index(parameterNames(j),"{")+1,index(parameterNames(j),"}")-1)
                parameterNames(j)=extract(parameterNames(j),                             1,index(parameterNames(j),"{")-1)
                read (labelIndex,*) indexElement
                indexElement=indexElement+1
             end if
             if (j == parameterCount) then
                ! This is the final parameter - so get and store a pointer to its node.
                self%modelParametersActive_(i)%parameter_   => parameters_%node         (char(parameterNames(j)),requireValue=.true. ,copyInstance=instance)
                self%modelParametersActive_(i)%indexElement =  indexElement
             else
                ! This is an intermediate parameter, get the appropriate sub-parameters.
                parameters_                                 =  parameters_%subParameters(char(parameterNames(j)),requireValue=.false.,copyInstance=instance)
             end if
          end do
          deallocate(parameterNames)
       end do
    end if
    if (.not.allocated(self%modelParametersInactive_)) then
       allocate(self%modelParametersInactive_(size(modelParametersInactive_)))
       do i=1,size(modelParametersInactive_)
          parameterCount=String_Count_Words(char(modelParametersInactive_(i)%modelParameter_%name()),"::")
          allocate(parameterNames(parameterCount))
          call String_Split_Words(parameterNames,char(modelParametersInactive_(i)%modelParameter_%name()),"::")
          parameters_=self%parametersModel
          do j=1,parameterCount
             instance    =1
             indexElement=0
              if (index(parameterNames(j),"[") /= 0) then
                labelIndex       =extract(parameterNames(j),index(parameterNames(j),"[")+1,index(parameterNames(j),"]")-1)
                parameterNames(j)=extract(parameterNames(j),                             1,index(parameterNames(j),"[")-1)
                read (labelIndex,*) instance
                instance=instance+1
             else if (index(parameterNames(j),"{") /= 0) then
                labelIndex       =extract(parameterNames(j),index(parameterNames(j),"{")+1,index(parameterNames(j),"}")-1)
                parameterNames(j)=extract(parameterNames(j),                             1,index(parameterNames(j),"{")-1)
                read (labelIndex,*) indexElement
                indexElement=indexElement+1
             end if
             if (j == parameterCount) then
                ! This is the final parameter - so get and store a pointer to its node.
                self%modelParametersInactive_(i)%parameter_   => parameters_%node         (char(parameterNames(j)),requireValue=.true. ,copyInstance=instance)
                self%modelParametersInactive_(i)%indexElement =  indexElement
             else
                ! This is an intermediate parameter, get the appropriate sub-parameters.
                parameters_                                   =  parameters_%subParameters(char(parameterNames(j)),requireValue=.false.,copyInstance=instance)
             end if
          end do
          deallocate(parameterNames)
       end do
    end if
    ! Iterate over all chains.
    do iRank=0,mpiSelf%count()-1
       ! If prior probability is impossible, then no need to waste time evaluating the likelihood.
       if (logPriorsProposed(iRank) <= logImpossible) cycle
       ! Update parameter values.
       do i=1,size(modelParametersActive_)
          if (self%modelParametersActive_(i)%indexElement == 0) then
             ! Simply overwrite the parameter.
             call self%modelParametersActive_(i)%parameter_%set(modelParametersActive_(i)%modelParameter_%unmap(stateVector(i,iRank)))
          else
             ! Overwrite only the indexed parameter in the list.
             parameterText =self%modelParametersActive_(i)%parameter_%get()
             parameterCount=String_Count_Words(char(parameterText))
             if (self%modelParametersActive_(i)%indexElement > parameterCount)                                      &
                  & call Galacticus_Error_Report(                                                                   &
                  &                              var_str('attempt to access non-existant element {')             // &
                  &                                     (self%modelParametersActive_(i)%indexElement          -1)// &
                  &                                      '} of parameter "'                                      // &
                  &                              char   (     modelParametersActive_(i)%modelParameter_%name()  )// &
                  &                                      '"'                                                     // &
                  &                              {introspection:location}                                        &
                  &                             )
             allocate(parameterNames(parameterCount))
             call String_Split_Words(parameterNames,char(parameterText))
             write (valueText,'(e24.16)') modelParametersActive_(i)%modelParameter_%unmap(stateVector(i,iRank))
             parameterNames(self%modelParametersActive_(i)%indexElement)=trim(valueText)
             call self%modelParametersActive_(i)%parameter_%set(String_Join(parameterNames," "))
             deallocate(parameterNames)
          endif
       end do
       ! Resolve dependencies in derived parameters.
       if (size(modelParametersInactive_) > 0) then
          do i=1,size(modelParametersInactive_)
             select type (modelParameter_ => modelParametersInactive_(i)%modelParameter_)
                class is (modelParameterDerived)
                self%modelParametersInactive_(i)%definition=modelParameter_%definition()
                self%modelParametersInactive_(i)%resolved  =.false.
             end select
          end do
          firstIteration      =.true.
          dependenciesResolved=.false.
          do while (.not.dependenciesResolved)
             dependenciesResolved=.true.
             dependenciesUpdated =.false.
             do i=1,size(modelParametersInactive_)
                select type (modelParameter_ => modelParametersInactive_(i)%modelParameter_)
                   class is (modelParameterDerived)
                   if (index(self%modelParametersInactive_(i)%definition,"%[") /= 0) then
                      ! The expression contains dependencies on other variables. Substitute the actual values where possible.
                      !! For active parameters we only need to consider substitution on the first iteration (since they are fully defined immediately).
                      !! Also handle special parameters here:
                      !!  * %[posteriorSimulationStep] - this is the current step number in the simulation.
                      if (firstIteration) then
                         if (index(self%modelParametersInactive_(i)%definition,"%[posteriorSimulationStep]") /= 0) then
                            write (valueText,'(i10)') simulationState%count()
                            self%modelParametersInactive_(i)%definition=replace(                                                                     &
                                 &                                                    self% modelParametersInactive_(i)%definition                 , &
                                 &                                                    "%[posteriorSimulationStep]"                                 , &
                                 &                                                    valueText                                                    , &
                                 &                                              every=.true.                                                         &
                                 &                                             )
                         end if
                         do j=1,size(modelParametersActive_)
                            if (index(self%modelParametersInactive_(i)%definition,"%["//modelParametersActive_  (j)%modelParameter_%name()//"]") /= 0) dependenciesUpdated=.true.
                            self%modelParametersInactive_(i)%definition=replace(                                                                     &
                                 &                                                    self% modelParametersInactive_(i)%definition                 , &
                                 &                                                    "%["//modelParametersActive_  (j)%modelParameter_%name()//"]", &
                                 &                                                    self% modelParametersActive_  (j)%parameter_     %get ()     , &
                                 &                                              every=.true.                                                         &
                                 &                                             )
                         end do
                      end if
                      !! For inactive parameters we must consider them each iteration as they become resolved.
                      do j=1,size(modelParametersInactive_)
                         if (i /= j .and. self%modelParametersInactive_(j)%resolved) then
                            if (index(self%modelParametersInactive_(i)%definition,"%["//modelParametersInactive_(j)%modelParameter_%name()//"]") /= 0) dependenciesUpdated=.true.
                            self%modelParametersInactive_(i)%definition=replace(                                                                     &
                                 &                                                    self% modelParametersInactive_(i)%definition                 , &
                                 &                                                    "%["//modelParametersInactive_(j)%modelParameter_%name()//"]", &
                                 &                                                    self% modelParametersInactive_(j)%parameter_     %get ()     , &
                                 &                                              every=.true.                                                         &
                                 &                                             )
                         end if
                      end do
                   end if
                   if (index(self%modelParametersInactive_(i)%definition,"%[") == 0) then
                      ! No dependencies remain, the expression can be evaluated.
                      self%modelParametersInactive_(i)%resolved=.true.
                      dependenciesUpdated                      =.true.
#ifdef MATHEVALAVAIL
                      evaluator   =Evaluator_Create_(char(self%modelParametersInactive_(i)%definition))
                      valueDerived=Evaluator_Evaluate_(evaluator,0,"",0.0d0)
                      call Evaluator_Destroy_(evaluator)
#else
                      call Galacticus_Error_Report('derived parameters require libmatheval, but it is not installed'//{introspection:location})
#endif
                      if (self%modelParametersInactive_(i)%indexElement == 0) then
                         ! Simply overwrite the parameter.
                         call self%modelParametersInactive_(i)%parameter_%set(valueDerived)
                      else
                         ! Overwrite only the indexed parameter in the list.
                         parameterText =self%modelParametersInactive_(i)%parameter_%get()
                         parameterCount=String_Count_Words(char(parameterText))
                         allocate(parameterNames(parameterCount))
                         call String_Split_Words(parameterNames,char(parameterText))
                         write (valueText,'(e24.16)') valueDerived
                         parameterNames(self%modelParametersInactive_(i)%indexElement)=trim(valueText)
                         call self%modelParametersInactive_(i)%parameter_%set(String_Join(parameterNames," "))
                         deallocate(parameterNames)
                      end if
                   else
                      dependenciesResolved=.false.
                   end if
                   class default
                   call Galacticus_Error_Report('support for this parameter type is not implemented'//{introspection:location})
                end select
             end do
             if (.not.dependenciesUpdated) then
                call displayVerbositySet(verbosityLevelStandard)
                call displayIndent('unresolved parameters')
                do i=1,size(modelParametersInactive_)
                   select type (modelParameter_ => modelParametersInactive_(i)%modelParameter_)
                      class is (modelParameterDerived)
                      if (index(self%modelParametersInactive_(i)%definition,"%[") /= 0) call displayMessage(modelParametersInactive_(i)%modelParameter_%name()//" : "//self%modelParametersInactive_(i)%definition)
                   end select
                end do
                call displayUnindent('unresolved parameters')
                call Galacticus_Error_Report('can not resolve parameter dependencies'//{introspection:location})
             end if
             firstIteration=.false.
          end do
       end if
       ! Build the task and outputter objects.
       call Tasks_Evolve_Forest_Construct_(self%parametersModel,self%task_)
       !# <objectBuilder class="outputAnalysis" name="self%outputAnalysis_" source="self%parametersModel"/>
       ! Perform the forest evolution tasks.
       call CPU_Time(timeBegin)
       call Tasks_Evolve_Forest_Perform_(self%task_,status)
       if (mpiSelf%any(status /= errorStatusSuccess)) then
          ! Forest evolution failed - record impossible likelihood.
          if (iRank == mpiSelf%rank()) then
             ! Dump the failed parameter set to file.
             call self%parametersModel%serializeToXML(self%failedParametersFileName//"."//iRank//".errCode"//status)
             ! Return impossible likelihood. We use a somewhat-less-than-impossible value to avoid this being rejected as the
             ! initial state.
             logLikelihoodProposed   =logImprobable
             galaxyPopulationEvaluate=logLikelihoodProposed
          end if
       else
          ! Forest evolution was successful - evaluate the likelihood.
          ! Extract the log-likelihood. This is evaluated by all chains (as they likely need to perform reduction across MPI
          ! processes), but only stored for the chain of this rank.
          logLikelihoodProposed=self%outputAnalysis_%logLikelihood()
          if (iRank == mpiSelf%rank()) galaxyPopulationEvaluate=logLikelihoodProposed
       end if
       if (iRank == mpiSelf%rank()) then
          ! Record timing information.
          call CPU_Time(timeEnd)
          timeEvaluate=timeEnd-timeBegin
          if (verbosityLevel >= verbosityLevelStandard) then
             write (valueText,'(e12.4)') logLikelihoodProposed
             message=var_str("Chain ")//simulationState%chainIndex()//" has logℒ="//trim(valueText)
             call displayMessage(message,verbosityLevelSilent)
          end if
       end if
       call mpiBarrier()
       call Tasks_Evolve_Forest_Destruct_(self%task_)
       !# <objectDestructor name="self%outputAnalysis_"/>
       call self%parametersModel%reset()
    end do
    ! Restore verbosity level.
    call displayVerbositySet(verbosityLevel)
    return
  end function galaxyPopulationEvaluate

  subroutine galaxyPopulationFunctionChanged(self)
    !% Respond to possible changes in the likelihood function.
    implicit none
    class(posteriorSampleLikelihoodGalaxyPopulation), intent(inout) :: self
    !$GLC attributes unused :: self

    return
  end subroutine galaxyPopulationFunctionChanged

  logical function galaxyPopulationWillEvaluate(self,simulationState,modelParameters_,simulationConvergence,temperature,logLikelihoodCurrent,logPriorCurrent,logPriorProposed)
    !% Return true if the log-likelihood will be evaluated.
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
