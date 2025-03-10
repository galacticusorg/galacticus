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
  Implementation of a model likelihood class which combines other likelihoods assumed to be independent.
  !!}

  use :: Posterior_Sampling_State, only : posteriorSampleStateSimple

  type, public :: posteriorSampleLikelihoodList
     class  (posteriorSampleLikelihoodClass), pointer                   :: modelLikelihood_        => null()
     integer                                , dimension(:), allocatable :: parameterMap                     , parameterMapInactive
     type   (modelParameterList            ), dimension(:), allocatable :: modelParametersActive_           , modelParametersInactive_
     type   (varying_string                ), dimension(:), allocatable :: parameterMapNames                , parameterMapNamesInactive
     type   (posteriorSampleLikelihoodList ), pointer                   :: next                    => null()
     type   (posteriorSampleStateSimple    )                            :: simulationState
     logical                                                            :: parameterMapInitialized          , report
  end type posteriorSampleLikelihoodList

  !![
  <enumeration>
   <name>orderRotation</name>
   <description>Specifies how to rotate the order of likelihood evaluation by process number.</description>
   <validator>yes</validator>
   <encodeFunction>yes</encodeFunction>
   <entry label="none"        />
   <entry label="byRank"      />
   <entry label="byRankOnNode"/>
  </enumeration>
  !!]

  !![
  <posteriorSampleLikelihood name="posteriorSampleLikelihoodIndependentLikelihoods">
    <description>
      A posterior sampling likelihood class which combines likelihoods from one or more other \refClass{posteriorSampleLikelihoodClass}
      classes that are assumed to be independent (i.e. the $\log \mathcal{L}$ of the models are simply summed to find the final
      likelihood).
      
      Since each \refClass{posteriorSampleLikelihoodClass} class may require a different set of parameters a {\normalfont \ttfamily
      [parameterMap]} parameter may be specified. If present, the number of {\normalfont \ttfamily [parameterMap]} parameters must
      equal the number of {\normalfont \ttfamily [posteriorSampleLikelihood]} parameters. Each such parameter should give a
      (space-separated) list of the names of parameters (as defined in the \refClass{modelParameterActive} class) which should be
      passed to the corresponding {\normalfont \ttfamily [posteriorSampleLikelihood]}. If no {\normalfont \ttfamily
      [parameterMap]} parameters are given then all parameters are passed to each \refClass{posteriorSampleLikelihoodClass} class.
      
      Similarly, a set of {\normalfont \ttfamily parameterInactiveMap} parameters may be given, to specify which (if any, an empty
      {\normalfont \ttfamily value} is permissible) of the inactive parameters specified by \refClass{modelParameterInactive}
      should be passed to the corresponding {\normalfont \ttfamily [posteriorSampleLikelihood]}. If no {\normalfont \ttfamily
      [parameterInactiveMap]} then no inactive parameters are passed to any of the {\normalfont \ttfamily
      [posteriorSampleLikelihood]} classes.

      Optionally, a parameter {\normalfont \ttfamily [logLikelihoodAccept]} may be specified. Once the likelihood of a chain
      reaches this value, no further evaluations of the likelihood will be made - the chain is assumed to be sufficiently likely
      that it is ``acceptable''.
    </description>
   <linkedList type="posteriorSampleLikelihoodList" variable="modelLikelihoods" next="next" object="modelLikelihood_" objectType="posteriorSampleLikelihoodClass"/>
  </posteriorSampleLikelihood>
  !!]
  type, extends(posteriorSampleLikelihoodClass) :: posteriorSampleLikelihoodIndependentLikelihoods
     !!{
     Implementation of a posterior sampling likelihood class which combines other likelihoods assumed to be independent.
     !!}
     private
     type            (enumerationOrderRotationType )          :: orderRotation
     type            (posteriorSampleLikelihoodList), pointer :: modelLikelihoods    => null()
     double precision                                         :: logLikelihoodAccept
     logical                                                  :: report                       , parameterMapIdentity
   contains
     final     ::                    independentLikelihoodsDestructor
     procedure :: evaluate        => independentLikelihoodsEvaluate
     procedure :: functionChanged => independentLikelihoodsFunctionChanged
  end type posteriorSampleLikelihoodIndependentLikelihoods

  interface posteriorSampleLikelihoodIndependentLikelihoods
     !!{
     Constructors for the {\normalfont \ttfamily independentLikelihoods} posterior sampling convergence class.
     !!}
     module procedure independentLikelihoodsConstructorParameters
     module procedure independentLikelihoodsConstructorInternal
  end interface posteriorSampleLikelihoodIndependentLikelihoods

contains

  function independentLikelihoodsConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily independentLikelihoods} posterior sampling convergence class which builds the object from a
    parameter set.
    !!}
    use :: Error           , only : Error_Report
    use :: Input_Parameters, only : inputParameter                         , inputParameterErrorStatusEmptyValue, inputParameterErrorStatusSuccess, inputParameters, &
         &                          enumerationInputParameterErrorStatusType
    use :: String_Handling , only : String_Count_Words                      , String_Split_Words                 , char
    implicit none
    type   (posteriorSampleLikelihoodIndependentLikelihoods)                :: self
    type   (inputParameters                                ), intent(inout) :: parameters
    type   (posteriorSampleLikelihoodList                  ), pointer       :: modelLikelihood_  , modelLikelihoodLast
    integer                                                                 :: i                 , parameterMapCount  , &
         &                                                                     countRotation     , countLikelihoods
    type   (enumerationInputParameterErrorStatusType       )                :: errorStatus
    type   (varying_string                                 )                :: parameterMapJoined, orderRotation

    !![
    <inputParameter>
      <name>orderRotation</name>
      <source>parameters</source>
      <defaultValue>var_str('none')</defaultValue>
      <description>The order in which evaluation of likelihoods should be rotated as a function of process number.</description>
    </inputParameter>
    <inputParameter>
      <name>logLikelihoodAccept</name>
      <variable>self%logLikelihoodAccept</variable>
      <defaultValue>huge(0.0d0)</defaultValue>
      <description>The log-likelihood which should be ``accepted''---once the log-likelihood reaches this value (or larger) no further updates to the chain will be made.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>report</name>
      <variable>self%report</variable>
      <defaultValue>.false.</defaultValue>
      <description>If true, report on the log-likelihood obtained.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self%parameterMapIdentity=.false.
    if     (                                                                             &
         &   parameters%copiesCount('posteriorSampleLikelihood',zeroIfNotPresent=.true.) &
         &  /=                                                                           &
         &   parameters%copiesCount('parameterMap'             ,zeroIfNotPresent=.true.) &
         & ) then
       if (parameters%copiesCount('parameterMap',zeroIfNotPresent=.true.) == 0) then
          self%parameterMapIdentity=.true.
       else
          call Error_Report('number of parameter maps must match number of likelihoods'//{introspection:location})
       end if
    end if
    self            %modelLikelihoods => null()
    modelLikelihood_                  => null()
    countLikelihoods                  =  parameters%copiesCount('posteriorSampleLikelihood',zeroIfNotPresent=.true.)
    do i=1,countLikelihoods
       if (associated(modelLikelihood_)) then
          allocate(modelLikelihood_%next)
          modelLikelihood_ => modelLikelihood_%next
       else
          allocate(self%modelLikelihoods)
          modelLikelihood_ => self%modelLikelihoods
       end if
       !![
       <objectBuilder class="posteriorSampleLikelihood" name="modelLikelihood_%modelLikelihood_" source="parameters" copy="i"/>
       !!]
       modelLikelihood_%simulationState        =posteriorSampleStateSimple(1)
       modelLikelihood_%parameterMapInitialized=.false.
       if (.not.self%parameterMapIdentity) then
          call parameters%value('parameterMap',parameterMapJoined,copyInstance=i)
          parameterMapCount=String_Count_Words(char(parameterMapJoined)," ")
          allocate(modelLikelihood_%parameterMap          (parameterMapCount))
          allocate(modelLikelihood_%parameterMapNames     (parameterMapCount))
          allocate(modelLikelihood_%modelParametersActive_(parameterMapCount))
          call String_Split_Words(modelLikelihood_%parameterMapNames,char(parameterMapJoined)," ")
          call modelLikelihood_%simulationState%parameterCountSet(parameterMapCount)
       end if
       if (parameters%copiesCount('parameterInactiveMap',zeroIfNotPresent=.true.) == 0) then
          ! No inactive parameter map is acceptable.
          allocate   (modelLikelihood_%modelParametersInactive_ (                0))
       else
          call parameters%value('parameterInactiveMap',parameterMapJoined,copyInstance=i,errorStatus=errorStatus)
          if      (errorStatus == inputParameterErrorStatusSuccess   ) then
             parameterMapCount=String_Count_Words(char(parameterMapJoined)," ")
             allocate(modelLikelihood_%parameterMapInactive     (parameterMapCount))
             allocate(modelLikelihood_%parameterMapNamesInactive(parameterMapCount))
             allocate(modelLikelihood_%modelParametersInactive_ (parameterMapCount))
             call String_Split_Words(modelLikelihood_%parameterMapNamesInactive,char(parameterMapJoined)," ")
          else if (errorStatus == inputParameterErrorStatusEmptyValue) then
             ! Empty value is acceptable.
             allocate(modelLikelihood_%modelParametersInactive_ (                0))
          else
             call Error_Report('invalid parameter'//{introspection:location})
          end if
       end if
    end do
    !! Perform rotation of likelihoods.
    self%orderRotation=enumerationOrderRotationEncode(char(orderRotation),includesPrefix=.false.)
    countRotation     =0
    select case (self%orderRotation%ID)
    case (orderRotationNone        %ID)
       countRotation=0
    case (orderRotationByRank      %ID)
       countRotation=mpiSelf%rank      ()
    case (orderRotationByRankOnNode%ID)
       countRotation=mpiSelf%rankOnNode()
    end select
    if (countRotation > 0 .and. countLikelihoods > 1) then
       do i=1,countRotation
          modelLikelihood_    => self%modelLikelihoods
          modelLikelihoodLast => self%modelLikelihoods
          do while (associated(modelLikelihoodLast%next))
             modelLikelihoodLast => modelLikelihoodLast%next
          end do
          self               %modelLikelihoods => modelLikelihood_%next
          modelLikelihoodLast%next             => modelLikelihood_
          modelLikelihood_   %next             => null()
       end do
    end if
    !![
    <inputParametersValidate source="parameters" multiParameters="posteriorSampleLikelihood, parameterMap, parameterInactiveMap" extraAllowedNames="parameterMap parameterInactiveMap"/>
    !!]
    return
  end function independentLikelihoodsConstructorParameters

  function independentLikelihoodsConstructorInternal(modelLikelihoods,logLikelihoodAccept,report,orderRotation) result(self)
    !!{
    Constructor for ``independentLikelihoods'' posterior sampling likelihood class.
    !!}
    implicit none
    type            (posteriorSampleLikelihoodIndependentLikelihoods)                        :: self
    type            (posteriorSampleLikelihoodList                  ), target, intent(in   ) :: modelLikelihoods
    double precision                                                         , intent(in   ) :: logLikelihoodAccept
    logical                                                                  , intent(in   ) :: report
    type            (enumerationOrderRotationType                   )        , intent(in   ) :: orderRotation
    !![
    <constructorAssign variables="*modelLikelihoods, logLikelihoodAccept, report, orderRotation"/>
    !!]

    return
  end function independentLikelihoodsConstructorInternal

  subroutine independentLikelihoodsDestructor(self)
    !!{
    Destructor for ``independentLikelihoods'' posterior sampling likelihood class.
    !!}
    implicit none
    type   (posteriorSampleLikelihoodIndependentLikelihoods), intent(inout) :: self
    type   (posteriorSampleLikelihoodList                  ), pointer       :: modelLikelihood_, modelLikelihoodNext
    integer                                                                 :: i
    
    if (associated(self%modelLikelihoods)) then
       modelLikelihood_ => self%modelLikelihoods
       do while (associated(modelLikelihood_))
          modelLikelihoodNext => modelLikelihood_%next
          !![
          <objectDestructor name="modelLikelihood_%modelLikelihood_"/>
          !!]
          do i=1,size(modelLikelihood_%modelParametersActive_  )
             !![
             <objectDestructor name="modelLikelihood_%modelParametersActive_  (i)%modelParameter_"/>
             !!]
          end do
          do i=1,size(modelLikelihood_%modelParametersInactive_)
             !![
             <objectDestructor name="modelLikelihood_%modelParametersInactive_(i)%modelParameter_"/>
             !!]
          end do
          deallocate(modelLikelihood_)
          modelLikelihood_ => modelLikelihoodNext
       end do
    end if
    return
  end subroutine independentLikelihoodsDestructor

  double precision function independentLikelihoodsEvaluate(self,simulationState,modelParametersActive_,modelParametersInactive_,simulationConvergence,temperature,logLikelihoodCurrent,logPriorCurrent,logPriorProposed,timeEvaluate,logLikelihoodVariance,forceAcceptance)
    !!{
    Return the log-likelihood for the halo mass function likelihood function.
    !!}
    use :: Display                     , only : displayMessage
    use :: Error                       , only : Error_Report
    use :: Models_Likelihoods_Constants, only : logImpossible
    implicit none
    class           (posteriorSampleLikelihoodIndependentLikelihoods), intent(inout), target       :: self
    class           (posteriorSampleStateClass                      ), intent(inout)               :: simulationState
    type            (modelParameterList                             ), intent(inout), dimension(:) :: modelParametersActive_, modelParametersInactive_
    class           (posteriorSampleConvergenceClass                ), intent(inout)               :: simulationConvergence
    double precision                                                 , intent(in   )               :: temperature           , logLikelihoodCurrent   , &
         &                                                                                            logPriorCurrent       , logPriorProposed
    real                                                             , intent(inout)               :: timeEvaluate
    double precision                                                 , intent(  out), optional     :: logLikelihoodVariance
    logical                                                          , intent(inout), optional     :: forceAcceptance
    type            (posteriorSampleLikelihoodList                  ), pointer                     :: modelLikelihood_
    double precision                                                 , allocatable  , dimension(:) :: stateVector           , stateVectorMapped
    real                                                                                           :: timeEvaluate_
    double precision                                                                               :: logLikelihoodVariance_, logPriorProposed_
    integer                                                                                        :: i                     , j
    character       (len=16                                         )                              :: label
    !$GLC attributes unused :: forceAcceptance

    allocate(stateVector      (simulationState%dimension()))
    allocate(stateVectorMapped(simulationState%dimension()))
    stateVector                                                =  simulationState%get()
    independentLikelihoodsEvaluate                             =  0.0d0
    modelLikelihood_                                           => self%modelLikelihoods
    timeEvaluate                                               =  0.0
    if (present(logLikelihoodVariance)) logLikelihoodVariance  =  0.0d0
    do while (associated(modelLikelihood_))
       if (.not.modelLikelihood_%parameterMapInitialized) then
          if (self%parameterMapIdentity) then
             allocate(modelLikelihood_%parameterMap          (size(modelParametersActive_)))
             allocate(modelLikelihood_%modelParametersActive_(size(modelParametersActive_)))
             call modelLikelihood_%simulationState%parameterCountSet(size(modelParametersActive_))
          end if
          do i=1,size(modelLikelihood_%parameterMap)
             ! Determine the mapping of the simulation state vector to this likelihood.
             if (self%parameterMapIdentity) then
                modelLikelihood_%parameterMap(i)=i
             else
                modelLikelihood_%parameterMap(i)=-1
                do j=1,size(modelParametersActive_)
                   if (modelParametersActive_(j)%modelParameter_%name() == modelLikelihood_%parameterMapNames(i)) then
                      modelLikelihood_%parameterMap(i)=j
                      exit
                   end if
                end do
                if (modelLikelihood_%parameterMap(i) == -1) call Error_Report('failed to find matching parameter ['//char(modelLikelihood_%parameterMapNames(i))//']'//{introspection:location})
             end if
             ! Copy the model parameter definition.
             allocate(modelLikelihood_%modelParametersActive_(i)%modelParameter_,mold=modelParametersActive_(modelLikelihood_%parameterMap(i))%modelParameter_)
             !![
             <deepCopyReset variables="modelParametersActive_(modelLikelihood_%parameterMap(i))%modelParameter_"/>
             <deepCopy source="modelParametersActive_(modelLikelihood_%parameterMap(i))%modelParameter_" destination="modelLikelihood_%modelParametersActive_(i)%modelParameter_"/>
             <deepCopyFinalize variables="modelLikelihood_%modelParametersActive_(i)%modelParameter_"/>
             !!]
          end do
          if (allocated(modelLikelihood_%parameterMapInactive)) then
             do i=1,size(modelLikelihood_%parameterMapInactive)
                ! Determine the mapping of the inactive parameters to this likelihood.
                modelLikelihood_%parameterMapInactive(i)=-1
                do j=1,size(modelParametersInActive_)
                   if (modelParametersInactive_(j)%modelParameter_%name() == modelLikelihood_%parameterMapNamesInactive(i)) then
                      modelLikelihood_%parameterMapInactive(i)=j
                      exit
                   end if
                end do
                if (modelLikelihood_%parameterMapInactive(i) == -1) call Error_Report('failed to find matching parameter ['//char(modelLikelihood_%parameterMapNamesInactive(i))//']'//{introspection:location})
                ! Copy the model parameter definition.
                allocate(modelLikelihood_%modelParametersInactive_(i)%modelParameter_,mold=modelParametersInactive_(modelLikelihood_%parameterMapInactive(i))%modelParameter_)
                !![
                <deepCopyReset variables="modelParametersInactive_(modelLikelihood_%parameterMapInactive(i))%modelParameter_"/>
                <deepCopy source="modelParametersInactive_(modelLikelihood_%parameterMapInactive(i))%modelParameter_" destination="modelLikelihood_%modelParametersInactive_(i)%modelParameter_"/>
                <deepCopyFinalize variables="modelLikelihood_%modelParametersInactive_(i)%modelParameter_"/>
                !!]
             end do
          end if
          ! Mark the likelihood as initialized.
          modelLikelihood_%parameterMapInitialized=.true.
       end if
       ! Map the overall simulation state to the state for this likelihood.
       forall(i=1:size(modelLikelihood_%parameterMap))
          stateVectorMapped(i)=stateVector(modelLikelihood_%parameterMap(i))
       end forall
       call modelLikelihood_%simulationState%update       (stateVectorMapped(1:size(modelLikelihood_%parameterMap)),logState=.false.,isConverged=.false.)
       call modelLikelihood_%simulationState%chainIndexSet(simulationState%chainIndex())
       call modelLikelihood_%simulationState%countSet     (simulationState%count     ())
       ! Determine if the chain is already accepted - if it is we set the proposed prior to be impossible so that the model will not actually be evaluated.
       if (logLikelihoodCurrent > self%logLikelihoodAccept) then
          logPriorProposed_=logImpossible
       else
          logPriorProposed_=logPriorProposed
       end if
       ! Evaluate this likelihood
       timeEvaluate_                                              =  -1.0
       independentLikelihoodsEvaluate                             =  +independentLikelihoodsEvaluate                                                        &
            &                                                        +modelLikelihood_%modelLikelihood_%evaluate(                                           &
            &                                                                                                    modelLikelihood_%simulationState         , &
            &                                                                                                    modelLikelihood_%modelParametersActive_  , &
            &                                                                                                    modelLikelihood_%modelParametersInactive_, &
            &                                                                                                                     simulationConvergence   , &
            &                                                                                                                     temperature             , &
            &                                                                                                                     logLikelihoodCurrent    , &
            &                                                                                                                     logPriorCurrent         , &
            &                                                                                                                     logPriorProposed_       , &
            &                                                                                                                     timeEvaluate_           , &
            &                                                                                                                     logLikelihoodVariance_    &
            &                                                                                                                    )
       if (present(logLikelihoodVariance)) logLikelihoodVariance  =  +logLikelihoodVariance_ &
            &                                                        +logLikelihoodVariance
       if (timeEvaluate_ >= 0.0d0        ) timeEvaluate           =  +timeEvaluate_          &
            &                                                        +timeEvaluate
       modelLikelihood_                                           =>  modelLikelihood_%next
    end do
    if (self%report) then
       write (label,'(e16.10)') independentLikelihoodsEvaluate
       call displayMessage("logℒ (total) = "//trim(label))
    end if
    return
  end function independentLikelihoodsEvaluate

  subroutine independentLikelihoodsFunctionChanged(self)
    !!{
    Respond to possible changes in the likelihood function.
    !!}
    implicit none
    class(posteriorSampleLikelihoodIndependentLikelihoods), intent(inout) :: self
    type (posteriorSampleLikelihoodList                  ), pointer       :: modelLikelihood_

    modelLikelihood_ => self%modelLikelihoods
    do while (associated(modelLikelihood_))
       call modelLikelihood_%modelLikelihood_%functionChanged()
       modelLikelihood_ => modelLikelihood_%next
    end do
    return
  end subroutine independentLikelihoodsFunctionChanged
