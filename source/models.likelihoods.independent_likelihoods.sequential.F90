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

  !![
  <posteriorSampleLikelihood name="posteriorSampleLikelihoodIndpndntLklhdsSqntl">
   <description>
    A posterior sampling likelihood class which sequentially combines other likelihoods assumed to be independent. This class
    begins by evaluating the first likelihood. If the likelihood is negative, then it is immediately returned, without
    evaluation of any further likelihoods. If it is positive, then the next likelihood is evaluated and the same conditions
    applied. This process repeats until either a negative likelihood is found, or all likelihoods are evaluated. Once a given
    likelihood has been evaluated it will be evaluated on all subsequent calls. Additionally, when a new likelihood is
    evaluated for the first time, acceptance of the proposed state will be forced. This class therefore allows a sequence of
    likelihoods to be specified which must be sequentially made sufficiently ``good'' before evaluating the next. The approach
    is intended to allow crude, but rapid constraints to be placed on parameters before progressing to more detailed, but slow
    to evaluate constraints.
   </description>
  </posteriorSampleLikelihood>
  !!]
  type, extends(posteriorSampleLikelihoodIndependentLikelihoods) :: posteriorSampleLikelihoodIndpndntLklhdsSqntl
     !!{
     Implementation of a posterior sampling likelihood class which sequentially combines other likelihoods assumed to be independent.
     !!}
     private
     integer                                     :: evaluateCount                , evaluateCountGlobal, &
          &                                         forceCount
     logical                                     :: finalLikelihoodFullEvaluation, restored           , &
          &                                         restoreLevels
     double precision, dimension(:), allocatable :: likelihoodMultiplier         , likelihoodAccept
   contains
     procedure :: evaluate => independentLikelihoodsSequentialEvaluate
     procedure :: restore  => independentLikelihoodsSequentialRestore
  end type posteriorSampleLikelihoodIndpndntLklhdsSqntl

  interface posteriorSampleLikelihoodIndpndntLklhdsSqntl
     !!{
     Constructors for the {\normalfont \ttfamily indpndntLklhdsSqntl} posterior sampling convergence class.
     !!}
     module procedure independentLikelihoodsSequentialConstructorParameters
     module procedure independentLikelihoodsSequentialConstructorInternal
  end interface posteriorSampleLikelihoodIndpndntLklhdsSqntl

  double precision, parameter :: logLikelihoodIncrement=1.0d2

contains

  function independentLikelihoodsSequentialConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily indpndntLklhdsSqntl} posterior sampling convergence class which builds the object from a
    parameter set.
    !!}
    use :: Error           , only : Error_Report
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (posteriorSampleLikelihoodIndpndntLklhdsSqntl)                :: self
    type   (inputParameters                             ), intent(inout) :: parameters
    integer                                                              :: i

    !![
    <inputParameter>
      <name>finalLikelihoodFullEvaluation</name>
      <variable>self%finalLikelihoodFullEvaluation</variable>
      <defaultValue>.true.</defaultValue>
      <description>If true the final likelihood is evaluated fully, and not treated as a ``lock in'' likelihood.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>restoreLevels</name>
      <variable>self%restoreLevels</variable>
      <defaultValue>.true.</defaultValue>
      <description>If true the level reached by each chain is restored on restarts. Otherwise, the level is initialized to zero (which may be useful for stochastic likelihoods).</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    ! Initialize the parent class.
    self%posteriorSampleLikelihoodIndependentLikelihoods=posteriorSampleLikelihoodIndependentLikelihoods(parameters)
    ! Get likelihood multipliers and acceptances.
    if     (                                                                             &
         &   parameters%copiesCount('posteriorSampleLikelihood',zeroIfNotPresent=.true.) &
         &  /=                                                                           &
         &   parameters%copiesCount('likelihoodMultiplier'     ,zeroIfNotPresent=.true.) &
         & ) call Error_Report('number of likelihood multipliers must match number of likelihoods'//{introspection:location})
    if     (                                                                             &
         &   parameters%copiesCount('posteriorSampleLikelihood',zeroIfNotPresent=.true.) &
         &  /=                                                                           &
         &   parameters%copiesCount('likelihoodAccept'         ,zeroIfNotPresent=.true.) &
         & ) call Error_Report('number of likelihood accepts must match number of likelihoods'    //{introspection:location})
    allocate(self%likelihoodMultiplier(parameters%copiesCount('posteriorSampleLikelihood',zeroIfNotPresent=.true.)))
    allocate(self%likelihoodAccept    (parameters%copiesCount('posteriorSampleLikelihood',zeroIfNotPresent=.true.)))
    do i=1,parameters%copiesCount('posteriorSampleLikelihood',zeroIfNotPresent=.true.)
       call parameters%value('likelihoodMultiplier',self%likelihoodMultiplier(i),copyInstance=i)
       call parameters%value('likelihoodAccept'    ,self%likelihoodAccept    (i),copyInstance=i)
    end do
    ! The evaluateCount value is the likelihood number at which we have achieved success so far.
    self%evaluateCount                                  =0
    self%evaluateCountGlobal                            =0
    self%forceCount                                     =0
    !![
    <inputParametersValidate source="parameters" multiParameters="likelihoodMultiplier, likelihoodAccept"/>
    !!]
    return
  end function independentLikelihoodsSequentialConstructorParameters

  function independentLikelihoodsSequentialConstructorInternal(modelLikelihoods,finalLikelihoodFullEvaluation,restoreLevels,likelihoodMultiplier,likelihoodAccept) result(self)
    !!{
    Constructor for {\normalfont \ttfamily indpndntLklhdsSqntl} posterior sampling likelihood class.
    !!}
    implicit none
    type            (posteriorSampleLikelihoodIndpndntLklhdsSqntl)                              :: self
    type            (posteriorSampleLikelihoodList               ), intent(in   ), target       :: modelLikelihoods
    logical                                                       , intent(in   )               :: finalLikelihoodFullEvaluation, restoreLevels
    double precision                                              , intent(in   ), dimension(:) :: likelihoodMultiplier         , likelihoodAccept
    !![
    <constructorAssign variables="*modelLikelihoods, finalLikelihoodFullEvaluation, restoreLevels, likelihoodMultiplier, likelihoodAccept"/>
    !!]

    ! The evaluateCount value is the likelihood number at which we have achieved success so far.
    self%evaluateCount      =0
    self%evaluateCountGlobal=0
    self%forceCount         =0
    self%restored           =.false.
    return
  end function independentLikelihoodsSequentialConstructorInternal

  double precision function independentLikelihoodsSequentialEvaluate(self,simulationState,modelParametersActive_,modelParametersInactive_,simulationConvergence,temperature,logLikelihoodCurrent,logPriorCurrent,logPriorProposed,timeEvaluate,logLikelihoodVariance,forceAcceptance)
    !!{
    Return the log-likelihood for the halo mass function likelihood function.
    !!}
    use :: Display                     , only : displayMessage
    use :: Error                       , only : Error_Report
    use :: ISO_Varying_String          , only : varying_string
    use :: MPI_Utilities               , only : mpiSelf
    use :: Models_Likelihoods_Constants, only : logImpossible , logImprobable
    use :: String_Handling             , only : operator(//)
    implicit none
    class           (posteriorSampleLikelihoodIndpndntLklhdsSqntl), intent(inout), target       :: self
    class           (posteriorSampleStateClass                   ), intent(inout)               :: simulationState
    type            (modelParameterList                          ), intent(inout), dimension(:) :: modelParametersActive_, modelParametersInactive_
    class           (posteriorSampleConvergenceClass             ), intent(inout)               :: simulationConvergence
    double precision                                              , intent(in   )               :: temperature           , logLikelihoodCurrent   , &
         &                                                                                         logPriorCurrent       , logPriorProposed
    real                                                          , intent(inout)               :: timeEvaluate
    double precision                                              , intent(  out), optional     :: logLikelihoodVariance
    logical                                                       , intent(inout), optional     :: forceAcceptance
    type            (posteriorSampleLikelihoodList               ), pointer                     :: modelLikelihood_
    integer                                                       , allocatable  , dimension(:) :: evaluateCounts        , chainIndices           , &
         &                                                                                         forceCounts
    double precision                                              , allocatable  , dimension(:) :: stateVector           , stateVectorMapped
    real                                                                                        :: timeEvaluate_
    double precision                                                                            :: logLikelihoodVariance_, logPriorProposed_      , &
         &                                                                                         logLikelihood
    integer                                                                                     :: i                     , j                      , &
         &                                                                                         likelihoodCount       , evaluateCount          , &
         &                                                                                         forceCount
    type            (varying_string                              )                              :: message
    logical                                                                                     :: finalLikelihood

    allocate(stateVector      (  simulationState%dimension()  ))
    allocate(stateVectorMapped(  simulationState%dimension()  ))
    allocate(evaluateCounts   (0:mpiSelf        %count    ()-1))
    allocate(forceCounts      (0:mpiSelf        %count    ()-1))
    allocate(chainIndices     (0:mpiSelf        %count    ()-1))
    stateVector                                               =  simulationState%get()
    independentLikelihoodsSequentialEvaluate                  =  0.0d0
    likelihoodCount                                           =  0
    evaluateCounts                                            =  mpiSelf%gather(self%evaluateCount)
    forceCounts                                               =  mpiSelf%gather(self%forceCount   )
    evaluateCount                                             =  minval(evaluateCounts)
    modelLikelihood_                                          => self%modelLikelihoods
    timeEvaluate                                              =  0.0d0
    if (present(logLikelihoodVariance)) logLikelihoodVariance =  0.0d0
    ! If a new global likelihood has been reached, report on it.
    if (mpiSelf%isMaster().and.evaluateCount > self%evaluateCountGlobal) then
       self%evaluateCountGlobal=evaluateCount
       message="sequential likelihood number "
       message=message//evaluateCount//" has been reached globally"
       call displayMessage(message)
    end if
    ! Initialize a local copy of the proposed log prior. In order to ensure MPI synchronization between processes we must always
    ! call evaluate on each independent likelihood. For example, the "Galacticus" likelihood class runs each chain under MPI
    ! across all processes, one at a time, and the MPI processes must coordinate to ensure only one such Galacticus is spawned at
    ! a time. To avoid unnecessary calculation we will force this local proposed log prior to an impossible value if we don't need
    ! to evaluate it. The called model likelihood functions are expected to instantly return in such cases without actually
    ! computing their likelihood.
    logPriorProposed_=logPriorProposed
    ! If this chain has already reached beyond the current evaluate level then we do not want to evaluate it again. Set the
    ! proposed prior to an impossible value to prevent the model likelihood being evaluated. This will result in the chain remaining
    ! locked into its current state.
    if (evaluateCounts(simulationState%chainIndex()) > evaluateCount .and. .not.self%restored) then
       logPriorProposed_                       =logImpossible
       independentLikelihoodsSequentialEvaluate=logImpossible
    end if
    ! Iterate through likelihoods, until none remain, or until we reach the maximum to currently evaluate (we do not evaluate
    ! further than the least likelihood reached over all chains).
    do while (associated(modelLikelihood_).and.likelihoodCount <= evaluateCount)
       likelihoodCount=likelihoodCount+1
       finalLikelihood=.not.associated(modelLikelihood_%next)
       if (.not.modelLikelihood_%parameterMapInitialized) then
          do i=1,size(modelLikelihood_%parameterMap)
             ! Determine the mapping of the simulation state vector to this likelihood.
             modelLikelihood_%parameterMap(i)=-1
             do j=1,size(modelParametersActive_)
                if (modelParametersActive_(j)%modelParameter_%name() == modelLikelihood_%parameterMapNames(i)) then
                   modelLikelihood_%parameterMap(i)=j
                   exit
                end if
             end do
             if (modelLikelihood_%parameterMap(i) == -1) call Error_Report('failed to find matching parameter ['//char(modelLikelihood_%parameterMapNames(i))//']'//{introspection:location})
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
       ! Evaluate this likelihood
       timeEvaluate_=-1.0
       logLikelihood                                             =+modelLikelihood_%modelLikelihood_%evaluate(                                           &
            &                                                                                                 modelLikelihood_%simulationState         , &
            &                                                                                                 modelLikelihood_%modelParametersActive_  , &
            &                                                                                                 modelLikelihood_%modelParametersInactive_, &
            &                                                                                                                  simulationConvergence   , &
            &                                                                                                                  temperature             , &
            &                                                                                                                  logLikelihoodCurrent    , &
            &                                                                                                                  logPriorCurrent         , &
            &                                                                                                                  logPriorProposed_       , &
            &                                                                                                                  timeEvaluate_           , &
            &                                                                                                                  logLikelihoodVariance_    &
            &                                                                                                                 )
       ! Modify the likelihood, unless it is improbable, in which case we simply leave it alone.
       if (logLikelihood > logImprobable) then
          !! First shift it so that acceptable likelihoods are always positive.
          logLikelihood=logLikelihood-self%likelihoodAccept    (likelihoodCount)
          !! Scale by a multiplying factor to steepen the likelihood function.
          logLikelihood=logLikelihood*self%likelihoodMultiplier(likelihoodCount)
          !! If likelihood is positive, fix it to our likelihood increment value.
          if (logLikelihood >= 0.0d0) logLikelihood=logLikelihoodIncrement
       end if
       ! Accumulate the likelihood unless our local proposed prior is impossible (in which case we did not actually need to
       ! compute this likelihood, but were calling the evaluate method just to ensure MPI synchronization).
       if (logPriorProposed_ >  logImpossible) independentLikelihoodsSequentialEvaluate=+independentLikelihoodsSequentialEvaluate &
            &                                                                           +logLikelihood
       if (present(logLikelihoodVariance)    ) logLikelihoodVariance                   =+logLikelihoodVariance                    &
            &                                                                           +logLikelihoodVariance_
       if (timeEvaluate_     >= 0.0d0        ) timeEvaluate                            =+timeEvaluate                             &
            &                                                                           +timeEvaluate_
       if (.not.(self%finalLikelihoodFullEvaluation.and.finalLikelihood)) then
          if (likelihoodCount > evaluateCounts(simulationState%chainIndex()) .or. self%restored) then
             ! We have matched or exceeded the previous number of likelihoods evaluated. (Or this is the first evaluation after being restored.)
             ! Force acceptance if this is the first time we have reached this far (unless the proposed prior is impossible - we
             ! do not want to accept states outside of the prior bounds).
             if (forceCounts(simulationState%chainIndex()) < evaluateCounts(simulationState%chainIndex()) .and. logPriorProposed > logImpossible) then
                if (.not.present(forceAcceptance)) call Error_Report('"forceAcceptance" argument must be present'//{introspection:location})
                forceCounts    (simulationState%chainIndex())=evaluateCounts(simulationState%chainIndex())
                forceAcceptance                              =.true.
             end if
             ! Check for acceptable likelihood.
             if (logLikelihood >= 0.0d0) then
                ! We have achieved a match at a new likelihood.
                evaluateCounts(simulationState%chainIndex())=likelihoodCount
             else
                ! We have evaluated at least as many likelihoods as previously, but now have a negative likelihood. Do not evaluate
                ! further.
                logPriorProposed_=logImpossible
             end if
          else if (logLikelihood < 0.0d0) then
             ! We have a negative likelihood and have not yet exceeded the maximum number of likelihoods previously
             ! achieved. Return an impossible likelihood to prevent this state from being accepted.
             independentLikelihoodsSequentialEvaluate=logImpossible
             logPriorProposed_                       =logImpossible
          end if
       end if
       modelLikelihood_ => modelLikelihood_%next
    end do
    ! Retrieve the evaluation count for this process back from whichever process ran our chain.
    evaluateCount =evaluateCounts(simulationState%chainIndex())
    forceCount    =forceCounts   (simulationState%chainIndex())
    evaluateCounts=mpiSelf%gather(evaluateCount)
    forceCounts   =mpiSelf%gather(forceCount   )
    chainIndices  =mpiSelf%gather(simulationState%chainIndex())
    do i=0,mpiSelf%count()-1
       if (chainIndices(i) == mpiSelf%rank()) then
          self%evaluateCount=evaluateCounts(i)
          self%forceCount   =forceCounts   (i)
          exit
       end if
    end do
    ! Unset restored status - the likelihood has now been evaluated.
    self%restored=.false.
    return
  end function independentLikelihoodsSequentialEvaluate

  subroutine independentLikelihoodsSequentialRestore(self,simulationState,logLikelihood)
    !!{
    Process a previous state to restore progress state.
    !!}
    implicit none
    class           (posteriorSampleLikelihoodIndpndntLklhdsSqntl), intent(inout)               :: self
    double precision                                              , intent(in   ), dimension(:) :: simulationState
    double precision                                              , intent(in   )               :: logLikelihood
    !$GLC attributes unused :: simulationState

    ! Detect the sequential state jumping to the next level.
    if (logLikelihood-dble(self%evaluateCount)*logLikelihoodIncrement > 0.0d0 .and. self%restoreLevels) then
       ! Increment the record of the level achieved.
       self%evaluateCount=self%evaluateCount+1
       self%   forceCount=self%   forceCount+1
    end if
    ! Set restored status - likelihood must always be evaluated after restoration.
    self%restored=.true.
    return
  end subroutine independentLikelihoodsSequentialRestore

