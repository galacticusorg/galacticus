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
  Implementation of a posterior sampling convergence class which switches between two other options.
  !!}

  !![
  <posteriorSampleStateInitialize name="posteriorSampleStateInitializeSwitched">
   <description>A posterior sampling state initialization class which sets initial state by switching between two other options.</description>
  </posteriorSampleStateInitialize>
  !!]
  type, extends(posteriorSampleStateInitializeClass) :: posteriorSampleStateInitializeSwitched
     !!{
     Implementation of a posterior sampling state initialization class which sets initial state by switching between two other options.
     !!}
     private
     class(posteriorSampleStateInitializeClass), pointer                   :: stateInitializeMethod1 => null(), stateInitializeMethod2 => null()
     type(varying_string                      ), allocatable, dimension(:) :: modelParameterName1             , modelParameterName2
   contains
     final     ::                switchedDestructor
     procedure :: initialize  => switchedInitialize
  end type posteriorSampleStateInitializeSwitched

  interface posteriorSampleStateInitializeSwitched
     !!{
     Constructors for the \refClass{posteriorSampleStateInitializeSwitched} posterior sampling state initialization class.
     !!}
     module procedure switchedConstructorParameters
     module procedure switchedConstructorInternal
  end interface posteriorSampleStateInitializeSwitched

contains

  function switchedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{posteriorSampleStateInitializeSwitched} posterior sampling state initialization class.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (posteriorSampleStateInitializeSwitched)                            :: self
    type (inputParameters                       ), intent(inout)             :: parameters
    class(posteriorSampleStateInitializeClass   ), pointer                   :: stateInitializeMethod1, stateInitializeMethod2
    type (varying_string                        ), allocatable, dimension(:) :: modelParameterName1   , modelParameterName2

    allocate(modelParameterName1(parameters%count('modelParameterName1')))
    allocate(modelParameterName2(parameters%count('modelParameterName2')))
    !![
    <inputParameter>
      <name>modelParameterName1</name>
      <source>parameters</source>
      <description>Names of parameters to be initialized by initializer number 1.</description>
    </inputParameter>
    <inputParameter>
      <name>modelParameterName2</name>
      <source>parameters</source>
      <description>Names of parameters to be initialized by initializer number 2.</description>
    </inputParameter>
    <objectBuilder class="posteriorSampleStateInitialize" name="stateInitializeMethod1" parameterName="posteriorSampleStateInitializeMethod1" source="parameters"/>
    <objectBuilder class="posteriorSampleStateInitialize" name="stateInitializeMethod2" parameterName="posteriorSampleStateInitializeMethod2" source="parameters"/>
    !!]
    self=posteriorSampleStateInitializeSwitched(modelParameterName1,modelParameterName2,stateInitializeMethod1,stateInitializeMethod2)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="stateInitializeMethod1"/>
    <objectDestructor name="stateInitializeMethod2"/>
    !!]
   return
  end function switchedConstructorParameters

  function switchedConstructorInternal(modelParameterName1,modelParameterName2,stateInitializeMethod1,stateInitializeMethod2) result(self)
    !!{
    Constructor for the \refClass{posteriorSampleStateInitializeSwitched} posterior sampling state initialization class.
    !!}
    implicit none
    type (posteriorSampleStateInitializeSwitched)                              :: self
    type (varying_string                        ), dimension(:), intent(in   ) :: modelParameterName1   , modelParameterName2
    class(posteriorSampleStateInitializeClass   ), target      , intent(in   ) :: stateInitializeMethod1, stateInitializeMethod2
    !![
    <constructorAssign variables="modelParameterName1, modelParameterName2, *stateInitializeMethod1, *stateInitializeMethod2"/>
    !!]

    return
  end function switchedConstructorInternal

  subroutine switchedDestructor(self)
    !!{
    Destructor for the \refClass{posteriorSampleStateInitializeSwitched} posterior sampling state initialization class.
    !!}
    implicit none
    type (posteriorSampleStateInitializeSwitched), intent(inout) :: self

    !![
    <objectDestructor name="self%stateInitializeMethod1"/>
    <objectDestructor name="self%stateInitializeMethod2"/>
    !!]
    return
  end subroutine switchedDestructor

  subroutine switchedInitialize(self,simulationState,modelParameters_,modelLikelihood,timeEvaluatePrevious,logLikelihood,logPosterior)
    !!{
    Initialize simulation state by drawing at random from the parameter priors.
    !!}
    use :: Error                       , only : Error_Report
    use :: Models_Likelihoods_Constants, only : logImpossible
    use :: Posterior_Sampling_State    , only : posteriorSampleStateSimple
    implicit none
    class           (posteriorSampleStateInitializeSwitched), intent(inout)               :: self
    class           (posteriorSampleStateClass             ), intent(inout)               :: simulationState
    class           (posteriorSampleLikelihoodClass        ), intent(inout)               :: modelLikelihood
    type            (modelParameterList                    ), intent(inout), dimension(:) :: modelParameters_
    double precision                                        , intent(  out)               :: timeEvaluatePrevious, logLikelihood, &
         &                                                                                   logPosterior
    class           (posteriorSampleStateInitializeClass   ), pointer                     :: stateInitializor_
    type            (varying_string                        ), allocatable  , dimension(:) :: modelParameterNames
    type            (posteriorSampleStateSimple            ), allocatable                 :: simulationState__
    type            (modelParameterList                    ), allocatable  , dimension(:) :: modelParameters__
    double precision                                        , allocatable  , dimension(:) :: stateVector         , stateVector__
    integer                                                 , allocatable  , dimension(:) :: mapping__
    integer                                                                               :: i                   , j            , &
         &                                                                                   initializer
    logical                                                                               :: matched

    ! Validate that all parameters are in one of our lists.
    do i=1,size(modelParameters_)
       matched=.false.
       if (.not.matched) then
          do j=1,size(self%modelParameterName1)
             if (self%modelParameterName1(j) == modelParameters_(i)%modelParameter_%name()) then
                matched=.true.
                exit
             end if
          end do
       end if
       if (.not.matched) then
          do j=1,size(self%modelParameterName2)
             if (self%modelParameterName2(j) == modelParameters_(i)%modelParameter_%name()) then
                matched=.true.
                exit
             end if
          end do
       end if
       if (.not.matched) call Error_Report('parameter "'//modelParameters_(i)%modelParameter_%name()//'" is not listed so would not be initialized'//{introspection:location})
    end do
    ! Iterate over both initializers.
    allocate(stateVector(size(modelParameters_)))
    do initializer=1,2
       select case (initializer)
       case (1)
          allocate(modelParameterNames(size(self%modelParameterName1)))
          modelParameterNames =  self%modelParameterName1
          stateInitializor_   => self%stateInitializeMethod1
       case (2)
          allocate(modelParameterNames(size(self%modelParameterName2)))
          modelParameterNames =  self%modelParameterName2
          stateInitializor_   => self%stateInitializeMethod2
       end select
       ! Construct a state object and set of model parameters for this initializer.
       allocate(simulationState__                           )
       allocate(stateVector__    (size(modelParameterNames)))
       allocate(mapping__        (size(modelParameterNames)))
       allocate(modelParameters__(size(modelParameterNames)))
       simulationState__=posteriorSampleStateSimple(1)
       call simulationState__%parameterCountSet(size(modelParameterNames))
       do i=1,size(modelParameterNames)
          matched=.false.
          do j=1,size(modelParameters_)
             if (modelParameters_(j)%modelParameter_%name() == modelParameterNames(i)) then
                mapping__(i)=j
                matched=.true.
                exit
             end if
          end do
          if (.not.matched) call Error_Report('named parameter "'//modelParameterNames(i)//'" does not appear in active parameters'//{introspection:location})
          allocate(modelParameters__(i)%modelParameter_,mold=modelParameters_(j)%modelParameter_)
          !![
          <deepCopyReset variables="modelParameters_(j)%modelParameter_"/>
          <deepCopy source="modelParameters_(j)%modelParameter_" destination="modelParameters__(i)%modelParameter_"/>
          <deepCopyFinalize variables="modelParameters__(i)%modelParameter_"/>
          !!]
       end do
       ! Apply the initializer
       call stateInitializor_%initialize(simulationState__,modelParameters__,modelLikelihood,timeEvaluatePrevious,logLikelihood,logPosterior)
       ! Combine states into the final state.
       stateVector__=simulationState__%get()
       do i=1,size(modelParameterNames)
          stateVector(mapping__(i))=stateVector__(i)
       end do
       ! Clean up.
       do i=1,size(modelParameterNames)
          !![
          <objectDestructor name="modelParameters__(i)%modelParameter_"/>
          !!]
       end do
       deallocate(modelParameterNames)
       deallocate(modelParameters__  )
       deallocate(simulationState__  )
       deallocate(stateVector__      )
       deallocate(mapping__          )
    end do
    ! Set the initial state.
    call simulationState%update(stateVector,.false.,.false.)
    deallocate(stateVector)
    ! We have no information about evaluation time.
    timeEvaluatePrevious=-1.0d0
    ! We have no information about the likelihood of this state.
    logLikelihood=logImpossible
    logPosterior =logImpossible
    return
  end subroutine switchedInitialize
