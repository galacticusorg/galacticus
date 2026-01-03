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
  Implementation of a posterior sampling stopping class which stepCount stops.
  !!}

  use :: Posterior_Sampling_Convergence, only : posteriorSampleConvergenceClass

  !![
  <posteriorSampleStoppingCriterion name="posteriorSampleStoppingCriterionStepCount">
   <description>
    This type will cause the simulation to stop when at least a number of steps (as specified by {\normalfont \ttfamily
    [stopAfterCount]}) have accrued post-convergence.
   </description>
  </posteriorSampleStoppingCriterion>
  !!]
  type, extends(posteriorSampleStoppingCriterionClass) :: posteriorSampleStoppingCriterionStepCount
     !!{
     Implementation of a posterior sampling convergence class which stepCount converges.
     !!}
     private
     integer                                           :: stopAfterCount
     class  (posteriorSampleConvergenceClass), pointer :: posteriorSampleConvergence_ => null()
   contains
     final     ::         stepCountDestructor
     procedure :: stop => stepCountStop
  end type posteriorSampleStoppingCriterionStepCount

  interface posteriorSampleStoppingCriterionStepCount
     !!{
     Constructors for the \refClass{posteriorSampleStoppingCriterionStepCount} posterior sampling convergence class.
     !!}
     module procedure stepCountConstructorParameters
     module procedure stepCountConstructorInternal
  end interface posteriorSampleStoppingCriterionStepCount

contains

  function stepCountConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{posteriorSampleStoppingCriterionStepCount} posterior sampling stopping class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (posteriorSampleStoppingCriterionStepCount)                :: self
    type   (inputParameters                          ), intent(inout) :: parameters
    class  (posteriorSampleConvergenceClass          ), pointer       :: posteriorSampleConvergence_
    integer                                                           :: stopAfterCount

    !![
    <inputParameter>
      <name>stopAfterCount</name>
      <description>The number of steps to continue after convergence before stopping.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="posteriorSampleConvergence" name="posteriorSampleConvergence_" source="parameters"/>
    !!]
    self=posteriorSampleStoppingCriterionStepCount(stopAfterCount,posteriorSampleConvergence_)
    !![
    <inputParametersValidate source="parameters"/>
   <objectDestructor name="posteriorSampleConvergence_"/>
   !!]
   return
  end function stepCountConstructorParameters

  function stepCountConstructorInternal(stopAfterCount,posteriorSampleConvergence_) result(self)
    !!{
    Internal constructor for the \refClass{posteriorSampleStoppingCriterionStepCount} posterior sampling stopping class.
    !!}
    implicit none
    type   (posteriorSampleStoppingCriterionStepCount)                        :: self
    integer                                           , intent(in   )         :: stopAfterCount
    class  (posteriorSampleConvergenceClass          ), intent(in   ), target :: posteriorSampleConvergence_
    !![
    <constructorAssign variables="stopAfterCount, *posteriorSampleConvergence_"/>
    !!]

    return
 end function stepCountConstructorInternal

  subroutine stepCountDestructor(self)
    !!{
    Destructor for the \refClass{posteriorSampleStoppingCriterionStepCount} posterior sampling stopping class.
    !!}
    implicit none
    type(posteriorSampleStoppingCriterionStepCount), intent(inout) :: self

    !![
    <objectDestructor name="self%posteriorSampleConvergence_"/>
    !!]
    return
  end subroutine stepCountDestructor

  logical function stepCountStop(self,simulationState)
    !!{
    Returns true if the posterior sampling should stop.
    !!}
    implicit none
    class(posteriorSampleStoppingCriterionStepCount), intent(inout) :: self
    class(posteriorSampleStateClass                ), intent(inout) :: simulationState

    stepCountStop=   self%posteriorSampleConvergence_%isConverged    () &
         &        .and.                                                 &
         &                simulationState            %count          () &
         &         >                                                    &
         &          +self%posteriorSampleConvergence_%convergedAtStep() &
         &          +self%stopAfterCount
    return
  end function stepCountStop
