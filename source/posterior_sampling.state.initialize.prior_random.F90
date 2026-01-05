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
  Implementation of a posterior sampling convergence class which priorRandom converges.
  !!}

  !![
  <posteriorSampleStateInitialize name="posteriorSampleStateInitializePriorRandom">
   <description>
    A posterior sampling state initialization class which samples the initial state at random from the prior distribution(s).
   </description>
  </posteriorSampleStateInitialize>
  !!]
  type, extends(posteriorSampleStateInitializeClass) :: posteriorSampleStateInitializePriorRandom
     !!{
     Implementation of a posterior sampling state initialization class which samples the initial state at random from the priors.
     !!}
     private
   contains
     procedure :: initialize  => priorRandomInitialize
  end type posteriorSampleStateInitializePriorRandom

  interface posteriorSampleStateInitializePriorRandom
     !!{
     Constructors for the \refClass{posteriorSampleStateInitializePriorRandom} posterior sampling state initialization class.
     !!}
     module procedure priorRandomConstructorParameters
  end interface posteriorSampleStateInitializePriorRandom

contains

  function priorRandomConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{posteriorSampleStateInitializePriorRandom} posterior sampling state initialization class.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(posteriorSampleStateInitializePriorRandom)                :: self
    type(inputParameters                          ), intent(inout) :: parameters

    self=posteriorSampleStateInitializePriorRandom()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function priorRandomConstructorParameters

  subroutine priorRandomInitialize(self,simulationState,modelParameters_,modelLikelihood,timeEvaluatePrevious,logLikelihood,logPosterior)
    !!{
    Initialize simulation state by drawing at random from the parameter priors.
    !!}
    use :: Models_Likelihoods_Constants, only : logImpossible
    implicit none
    class           (posteriorSampleStateInitializePriorRandom), intent(inout)                                    :: self
    class           (posteriorSampleStateClass                ), intent(inout)                                    :: simulationState
    class           (posteriorSampleLikelihoodClass           ), intent(inout)                                    :: modelLikelihood
    type            (modelParameterList                       ), intent(inout), dimension(:                     ) :: modelParameters_
    double precision                                           , intent(  out)                                    :: timeEvaluatePrevious, logLikelihood, &
         &                                                                                                           logPosterior
    double precision                                                          , dimension(size(modelParameters_)) :: state
    integer                                                                                                       :: j
    !$GLC attributes unused ::  self, modelLikelihood

    ! No knowledge of evaluation time.
    timeEvaluatePrevious=-1.0d0
    ! We have no information about the likelihood of this state.
    logLikelihood=logImpossible
    logPosterior =logImpossible
    ! Initialize chain to some state vector.
    do j=1,simulationState%dimension()
       state(j)=modelParameters_(j)%modelParameter_%map        (  &
            &   modelParameters_(j)%modelParameter_%priorSample ( &
            &                                                   ) &
            &                                                  )
    end do
    call simulationState%update(state,.false.,.false.)
  return
  end subroutine priorRandomInitialize
