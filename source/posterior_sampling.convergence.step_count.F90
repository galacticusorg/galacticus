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
  Implementation of a posterior sampling convergence class which declares convergence after a fixed number of steps.
  !!}

  !![
  <posteriorSampleConvergence name="posteriorSampleConvergenceStepCount">
   <description>A posterior sampling convergence class which declares convergence after a fixed number of steps.</description>
  </posteriorSampleConvergence>
  !!]
  type, extends(posteriorSampleConvergenceClass) :: posteriorSampleConvergenceStepCount
     !!{
     Implementation of a posterior sampling convergence class which declares convergence after a fixed number of steps.
     !!}
     private
     integer :: countSteps  , stepStart
     logical :: isConverged_
   contains
     procedure :: isConverged     => stepCountIsConverged
     procedure :: convergedAtStep => stepCountConvergedAtStep
     procedure :: reset           => stepCountReset
     procedure :: logReport       => stepCountLogReport
     procedure :: stateIsOutlier  => stepCountStateIsOutlier
  end type posteriorSampleConvergenceStepCount

  interface posteriorSampleConvergenceStepCount
     !!{
     Constructors for the {\normalfont \ttfamily stepCount} posterior sampling convergence class.
     !!}
     module procedure stepCountConstructorParameters
     module procedure stepCountConstructorInternal
  end interface posteriorSampleConvergenceStepCount

contains

  function stepCountConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily stepCount} convergence class.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (posteriorSampleConvergenceStepCount)                :: self
    type   (inputParameters                    ), intent(inout) :: parameters
    integer                                                     :: countSteps
    
    !![
    <inputParameter>
      <name>countSteps</name>
      <description>The number of steps after which to declare the simulation as converged.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=posteriorSampleConvergenceStepCount(countSteps)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function stepCountConstructorParameters

  function stepCountConstructorInternal(countSteps) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily stepCount} convergence class.
    !!}
    implicit none
    type   (posteriorSampleConvergenceStepCount)                :: self
    integer                                     , intent(in   ) :: countSteps
    !![
    <constructorAssign variables="countSteps"/>
    !!]

    self%stepStart   =0
    self%isConverged_=.false.
    return
  end function stepCountConstructorInternal

  logical function stepCountIsConverged(self,simulationState,logLikelihood)
    !!{
    Returns true if the posterior sampling is converged (which it stepCount is).
    !!}
    implicit none
    class           (posteriorSampleConvergenceStepCount), intent(inout)           :: self
    class           (posteriorSampleStateClass          ), intent(inout), optional :: simulationState
    double precision                                     , intent(in   ), optional :: logLikelihood
    !$GLC attributes unused :: logLikelihood

    if (present(simulationState)) self%isConverged_=simulationState%count() >= self%stepStart+self%countSteps
    stepCountIsConverged=self%isConverged_
    return
  end function stepCountIsConverged

  integer function stepCountConvergedAtStep(self)
    !!{
    Return the step at which the simulation converged.
    !!}
    implicit none
    class(posteriorSampleConvergenceStepCount), intent(inout) :: self

    stepCountConvergedAtStep=self%stepStart+self%countSteps
    return
  end function stepCountConvergedAtStep

  subroutine stepCountReset(self)
    !!{
    Reset the convergence object.
    !!}
    implicit none
    class(posteriorSampleConvergenceStepCount), intent(inout) :: self

    self%isConverged_=.false.
    self%stepStart   =self%stepStart+self%countSteps
    return
  end subroutine stepCountReset

  subroutine stepCountLogReport(self,fileUnit)
    !!{
    Write a convergence report to the given {\normalfont \ttfamily fileUnit}.
    !!}
    implicit none
    class  (posteriorSampleConvergenceStepCount), intent(inout) :: self
    integer                                     , intent(in   ) :: fileUnit

    if (self%isConverged_) then
       write (fileUnit,*) 'Convergence: converged'
    else
       write (fileUnit,*) 'Convergence: not converged'
    end if
    return
  end subroutine stepCountLogReport

  logical function stepCountStateIsOutlier(self,stateIndex)
    !!{
    Return true if the specified chain is deemed to be an outlier. In this case, chains are never outliers.
    !!}
    implicit none
    class  (posteriorSampleConvergenceStepCount), intent(inout) :: self
    integer                                     , intent(in   ) :: stateIndex
    !$GLC attributes unused :: self, stateIndex

    stepCountStateIsOutlier=.false.
    return
  end function stepCountStateIsOutlier
