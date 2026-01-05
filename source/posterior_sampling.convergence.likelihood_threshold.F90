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
  Implementation of a posterior sampling convergence class which declares convergence once all likelihoods are above a threshold.
  !!}

  !![
  <posteriorSampleConvergence name="posteriorSampleConvergenceLikelihoodThreshold">
   <description>A posterior sampling convergence class which declares convergence once all likelihoods are above a threshold.</description>
  </posteriorSampleConvergence>
  !!]
  type, extends(posteriorSampleConvergenceClass) :: posteriorSampleConvergenceLikelihoodThreshold
     !!{
     Implementation of a posterior sampling convergence class which declares convergence once all likelihoods are above a threshold.
     !!}
     private
     integer                                     :: convergedAtStepCount
     logical                                     :: converged
     double precision                            :: likelihoodThreshold
     logical         , allocatable, dimension(:) :: isOutlier
   contains
     procedure :: isConverged     => likelihoodThresholdIsConverged
     procedure :: convergedAtStep => likelihoodThresholdConvergedAtStep
     procedure :: reset           => likelihoodThresholdReset
     procedure :: logReport       => likelihoodThresholdLogReport
     procedure :: stateIsOutlier  => likelihoodThresholdStateIsOutlier
  end type posteriorSampleConvergenceLikelihoodThreshold

  interface posteriorSampleConvergenceLikelihoodThreshold
     !!{
     Constructors for the \refClass{posteriorSampleConvergenceLikelihoodThreshold} posterior sampling convergence class.
     !!}
     module procedure likelihoodThresholdConstructorParameters
     module procedure likelihoodThresholdConstructorInternal
  end interface posteriorSampleConvergenceLikelihoodThreshold

contains

  function likelihoodThresholdConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{posteriorSampleConvergenceLikelihoodThreshold} posterior sampling convergence class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (posteriorSampleConvergenceLikelihoodThreshold)                :: self
    type            (inputParameters                              ), intent(inout) :: parameters
    double precision                                                               :: likelihoodThreshold

    !![
    <inputParameter>
      <name>likelihoodThreshold</name>
      <defaultValue>0.0d0</defaultValue>
      <description>The threshold log-likelihood above which convergence is declared.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=posteriorSampleConvergenceLikelihoodThreshold(likelihoodThreshold)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function likelihoodThresholdConstructorParameters

  function likelihoodThresholdConstructorInternal(likelihoodThreshold) result(self)
    !!{
    Internal constructor for the \refClass{posteriorSampleConvergenceLikelihoodThreshold} posterior sampling convergence class.
    !!}
    use :: MPI_Utilities    , only : mpiSelf
    implicit none
    type            (posteriorSampleConvergenceLikelihoodThreshold)                :: self
    double precision                                               , intent(in   ) :: likelihoodThreshold
    !![
    <constructorAssign variables="likelihoodThreshold"/>
    !!]

    self%converged           =.false.
    self%convergedAtStepCount=huge(0)
    allocate(self%isOutlier(mpiSelf%count()))
    return
  end function likelihoodThresholdConstructorInternal

  logical function likelihoodThresholdIsConverged(self,simulationState,logLikelihood)
    !!{
    Returns true if the posterior sampling is converged (which it likelihoodThreshold is).
    !!}
    use :: MPI_Utilities, only : mpiSelf
    implicit none
    class           (posteriorSampleConvergenceLikelihoodThreshold), intent(inout)               :: self
    class           (posteriorSampleStateClass                    ), intent(inout), optional     :: simulationState
    double precision                                               , intent(in   ), optional     :: logLikelihood
    double precision                                               , allocatable  , dimension(:) :: logLikelihoods
    
    ! If no arguments were provided, return current convergence status without updating.
    if (.not.(present(simulationState).and.present(logLikelihood))) then
       likelihoodThresholdIsConverged=self%converged
       return
    end if
    ! Convergence requires all chains to have a likelihood above the threshold.
    logLikelihoods=mpiSelf%gather(logLikelihood)
    if (.not.self%converged.and.all(logLikelihoods >= self%likelihoodThreshold)) then
       self%converged           =.true.
       self%convergedAtStepCount=simulationState%count()
    else
       ! Mark states that exceed the likelihood threshold as outliers so that they can be ignored in some calculations.
       self%isOutlier=logLikelihoods >= self%likelihoodThreshold
    end if
    likelihoodThresholdIsConverged=self%converged
    return
  end function likelihoodThresholdIsConverged

  integer function likelihoodThresholdConvergedAtStep(self)
    !!{
    Return the step at which the simulation converged.
    !!}
    implicit none
    class(posteriorSampleConvergenceLikelihoodThreshold), intent(inout) :: self

    likelihoodThresholdConvergedAtStep=self%convergedAtStepCount
    return
  end function likelihoodThresholdConvergedAtStep

  subroutine likelihoodThresholdReset(self)
    !!{
    Reset the convergence object.
    !!}
    implicit none
    class(posteriorSampleConvergenceLikelihoodThreshold), intent(inout) :: self

    self%converged           =.false.
    self%convergedAtStepCount=-1
    return
  end subroutine likelihoodThresholdReset

  subroutine likelihoodThresholdLogReport(self,fileUnit)
    !!{
    Write a convergence report to the given {\normalfont \ttfamily fileUnit}.
    !!}
    implicit none
    class  (posteriorSampleConvergenceLikelihoodThreshold), intent(inout) :: self
    integer                                               , intent(in   ) :: fileUnit

    if (self%converged) then
       write (fileUnit,*) 'Convergence: converged'
    else
       write (fileUnit,*) 'Convergence: unconverged'
    end if
    return
  end subroutine likelihoodThresholdLogReport

  logical function likelihoodThresholdStateIsOutlier(self,stateIndex)
    !!{
    Return true if the specified chain is deemed to be an outlier.
    !!}
    implicit none
    class  (posteriorSampleConvergenceLikelihoodThreshold), intent(inout) :: self
    integer                                               , intent(in   ) :: stateIndex

    likelihoodThresholdStateIsOutlier=self%isOutlier(stateIndex+1)
    return
  end function likelihoodThresholdStateIsOutlier
