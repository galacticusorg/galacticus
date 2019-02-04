!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

  !% Implementation of a model likelihood class which combines other likelihoods assumed to be independent.

  type, public :: posteriorSampleLikelihoodList
     class  (posteriorSampleLikelihoodClass), pointer                   :: modelLikelihood_
     integer                                , dimension(:), allocatable :: parameterMap                     , parameterMapInactive
     type   (modelParameterList            ), dimension(:), allocatable :: modelParametersActive_           , modelParametersInactive_
     type   (varying_string                ), dimension(:), allocatable :: parameterMapNames                , parameterMapNamesInactive
     type   (posteriorSampleLikelihoodList ), pointer                   :: next                    => null()
     type   (posteriorSampleStateSimple    )                            :: simulationState
     logical                                                            :: parameterMapInitialized
  end type posteriorSampleLikelihoodList

  !# <posteriorSampleLikelihood name="posteriorSampleLikelihoodIndependentLikelihoods">
  !#  <description>A posterior sampling likelihood class which combines other likelihoods assumed to be independent.</description>
  !# </posteriorSampleLikelihood>
  type, extends(posteriorSampleLikelihoodClass) :: posteriorSampleLikelihoodIndependentLikelihoods
     !% Implementation of a posterior sampling likelihood class which combines other likelihoods assumed to be independent.
     private
     type(posteriorSampleLikelihoodList), pointer :: modelLikelihoods
   contains
     final     ::                    independentLikelihoodsDestructor
     procedure :: evaluate        => independentLikelihoodsEvaluate
     procedure :: functionChanged => independentLikelihoodsFunctionChanged
  end type posteriorSampleLikelihoodIndependentLikelihoods

  interface posteriorSampleLikelihoodIndependentLikelihoods
     !% Constructors for the {\normalfont \ttfamily independentLikelihoods} posterior sampling convergence class.
     module procedure independentLikelihoodsConstructorParameters
     module procedure independentLikelihoodsConstructorInternal
  end interface posteriorSampleLikelihoodIndependentLikelihoods
  
contains

  function independentLikelihoodsConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily independentLikelihoods} posterior sampling convergence class which builds the object from a
    !% parameter set.
    use Input_Parameters
    use String_Handling
    use Galacticus_Error
    implicit none
    type   (posteriorSampleLikelihoodIndependentLikelihoods)                :: self
    type   (inputParameters                                ), intent(inout) :: parameters
    type   (posteriorSampleLikelihoodList                  ), pointer       :: modelLikelihood_
    integer                                                                 :: i                 , parameterMapCount, &
         &                                                                     errorStatus
    type   (varying_string                                 )                :: parameterMapJoined
    
    if     (                                                                                   &
         &   parameters%copiesCount('posteriorSampleLikelihoodMethod',zeroIfNotPresent=.true.) &
         &  /=                                                                                 &
         &   parameters%copiesCount('parameterMap'                   ,zeroIfNotPresent=.true.) &
         & ) call Galacticus_Error_Report('number of parameter maps must match number of likelihoods'//{introspection:location})
    self            %modelLikelihoods => null()
    modelLikelihood_                  => null()
    do i=1,parameters%copiesCount('posteriorSampleLikelihoodMethod',zeroIfNotPresent=.true.)
       if (associated(modelLikelihood_)) then
          allocate(modelLikelihood_%next)
          modelLikelihood_ => modelLikelihood_%next
       else
          allocate(self%modelLikelihoods)
          modelLikelihood_ => self%modelLikelihoods
       end if
       modelLikelihood_%modelLikelihood_        => posteriorSampleLikelihood (parameters,i)
       modelLikelihood_%simulationState         =  posteriorSampleStateSimple(           1)
       modelLikelihood_%parameterMapInitialized =  .false.
       call parameters%value('parameterMap',parameterMapJoined,copyInstance=i)
       parameterMapCount=String_Count_Words(char(parameterMapJoined)," ")
       allocate(modelLikelihood_%parameterMap          (parameterMapCount))
       allocate(modelLikelihood_%parameterMapNames     (parameterMapCount))
       allocate(modelLikelihood_%modelParametersActive_(parameterMapCount))
       call String_Split_Words(modelLikelihood_%parameterMapNames,char(parameterMapJoined)," ")
       call modelLikelihood_%simulationState%parameterCountSet(parameterMapCount)
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
          call Galacticus_Error_Report('invalid parameter'//{introspection:location})
       end if
    end do
    return
  end function independentLikelihoodsConstructorParameters

  function independentLikelihoodsConstructorInternal(modelLikelihoods) result(self)
    !% Constructor for ``independentLikelihoods'' posterior sampling likelihood class.
    implicit none
    type(posteriorSampleLikelihoodIndependentLikelihoods)                        :: self
    type(posteriorSampleLikelihoodList                  ), target, intent(in   ) :: modelLikelihoods
    !# <constructorAssign variables="*modelLikelihoods"/>

    return
  end function independentLikelihoodsConstructorInternal
  
  elemental subroutine independentLikelihoodsDestructor(self)
    !% Destructor for ``independentLikelihoods'' posterior sampling likelihood class.
    implicit none
    type(posteriorSampleLikelihoodIndependentLikelihoods), intent(inout) :: self
    type(posteriorSampleLikelihoodList                  ), pointer       :: modelLikelihood_, modelLikelihoodNext

    if (associated(self%modelLikelihoods)) then
       modelLikelihood_ => self%modelLikelihoods
       do while (associated(modelLikelihood_))
          modelLikelihoodNext => modelLikelihood_%next
          deallocate(modelLikelihood_%modelLikelihood_)
          deallocate(modelLikelihood_          )
          modelLikelihood_ => modelLikelihoodNext
       end do
    end if
    return
  end subroutine independentLikelihoodsDestructor
  
  double precision function independentLikelihoodsEvaluate(self,simulationState,modelParametersActive_,modelParametersInactive_,simulationConvergence,temperature,logLikelihoodCurrent,logPriorCurrent,logPriorProposed,timeEvaluate,logLikelihoodVariance,forceAcceptance)
    !% Return the log-likelihood for the halo mass function likelihood function.
    use Galacticus_Error
    implicit none
    class           (posteriorSampleLikelihoodIndependentLikelihoods), intent(inout)               :: self
    class           (posteriorSampleStateClass                      ), intent(inout)               :: simulationState
    type            (modelParameterList                             ), intent(in   ), dimension(:) :: modelParametersActive_, modelParametersInactive_
    class           (posteriorSampleConvergenceClass                ), intent(inout)               :: simulationConvergence
    double precision                                                 , intent(in   )               :: temperature           , logLikelihoodCurrent   , &
         &                                                                                            logPriorCurrent       , logPriorProposed
    real                                                             , intent(inout)               :: timeEvaluate
    double precision                                                 , intent(  out), optional     :: logLikelihoodVariance
    logical                                                          , intent(inout), optional     :: forceAcceptance
    type            (posteriorSampleLikelihoodList                  ), pointer                     :: modelLikelihood_
    double precision                                                 , allocatable  , dimension(:) :: stateVector           , stateVectorMapped
    double precision                                                                               :: logLikelihoodVariance_, timeEvaluate_
    integer                                                                                        :: i                     , j
    !GCC$ attributes unused :: forceAcceptance

    allocate(stateVector      (simulationState%dimension()))
    allocate(stateVectorMapped(simulationState%dimension()))
    stateVector                                                =  simulationState%get()
    independentLikelihoodsEvaluate                             =  0.0d0
    modelLikelihood_                                           => self%modelLikelihoods
    timeEvaluate_                                              =  0.0d0
    if (present(logLikelihoodVariance)) logLikelihoodVariance_ =  0.0d0
    do while (associated(modelLikelihood_))
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
             if (modelLikelihood_%parameterMap(i) == -1) call Galacticus_Error_Report('failed to find matching parameter ['//char(modelLikelihood_%parameterMapNames(i))//']'//{introspection:location})             
             ! Copy the model parameter definition.
             allocate(modelLikelihood_%modelParametersActive_(i)%modelParameter_,mold=modelParametersActive_(modelLikelihood_%parameterMap(i))%modelParameter_)
             call modelParametersActive_(modelLikelihood_%parameterMap(i))%modelParameter_%deepCopy(modelLikelihood_%modelParametersActive_(i)%modelParameter_)
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
                if (modelLikelihood_%parameterMapInactive(i) == -1) call Galacticus_Error_Report('failed to find matching parameter ['//char(modelLikelihood_%parameterMapNamesInactive(i))//']'//{introspection:location})
                ! Copy the model parameter definition.
                allocate(modelLikelihood_%modelParametersInactive_(i)%modelParameter_,mold=modelParametersInactive_(modelLikelihood_%parameterMapInactive(i))%modelParameter_)
                call modelParametersInactive_(modelLikelihood_%parameterMapInactive(i))%modelParameter_%deepCopy(modelLikelihood_%modelParametersInactive_(i)%modelParameter_)
             end do
          end if
          ! Mark the likelihood as initialized.
          modelLikelihood_%parameterMapInitialized=.true.
       end if
       ! Map the overall simulation state to the state for this likelihood.
       forall(i=1:size(modelLikelihood_%parameterMap))
          stateVectorMapped(i)=stateVector(modelLikelihood_%parameterMap(i))
       end forall
       call modelLikelihood_%simulationState%update(stateVectorMapped(1:size(modelLikelihood_%parameterMap)),logState=.false.,isConverged=.false.)
       ! Evaluate this likelihood
       independentLikelihoodsEvaluate                             =  +independentLikelihoodsEvaluate                                                        &
            &                                                        +modelLikelihood_%modelLikelihood_%evaluate(                                           &
            &                                                                                                    modelLikelihood_%simulationState         , &
            &                                                                                                    modelLikelihood_%modelParametersActive_  , &
            &                                                                                                    modelLikelihood_%modelParametersInactive_, &
            &                                                                                                                     simulationConvergence   , &
            &                                                                                                                     temperature             , &
            &                                                                                                                     logLikelihoodCurrent    , &
            &                                                                                                                     logPriorCurrent         , &
            &                                                                                                                     logPriorProposed        , &
            &                                                                                                                     timeEvaluate            , &
            &                                                                                                                     logLikelihoodVariance     &
            &                                                                                                                    )
       if (present(logLikelihoodVariance)) logLikelihoodVariance_ =  +logLikelihoodVariance_ &
            &                                                        +logLikelihoodVariance
       timeEvaluate_                                              =  +timeEvaluate_          &
            &                                                        +timeEvaluate
       modelLikelihood_                                           =>  modelLikelihood_%next
    end do
    return    
  end function independentLikelihoodsEvaluate
  
  subroutine independentLikelihoodsFunctionChanged(self)
    !% Respond to possible changes in the likelihood function.
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
