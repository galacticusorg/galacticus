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
  Implementation of a posterior sampling convergence class which is always converged.
  !!}

  !![
  <posteriorSampleConvergence name="posteriorSampleConvergenceAlways">
   <description>A posterior sampling convergence class which is always converged.</description>
  </posteriorSampleConvergence>
  !!]
  type, extends(posteriorSampleConvergenceClass) :: posteriorSampleConvergenceAlways
     !!{
     Implementation of a posterior sampling convergence class which is always converged.
     !!}
     private
   contains
     procedure :: isConverged     => alwaysIsConverged
     procedure :: convergedAtStep => alwaysConvergedAtStep
     procedure :: reset           => alwaysReset
     procedure :: logReport       => alwaysLogReport
     procedure :: stateIsOutlier  => alwaysStateIsOutlier
  end type posteriorSampleConvergenceAlways

  interface posteriorSampleConvergenceAlways
     !!{
     Constructors for the \refClass{posteriorSampleConvergenceAlways} posterior sampling convergence class.
     !!}
     module procedure alwaysConstructorParameters
  end interface posteriorSampleConvergenceAlways

contains

  function alwaysConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{posteriorSampleConvergenceAlways} merger tree halo mass function sampling class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(posteriorSampleConvergenceAlways)                :: self
    type(inputParameters                 ), intent(inout) :: parameters

    self=posteriorSampleConvergenceAlways()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function alwaysConstructorParameters

  logical function alwaysIsConverged(self,simulationState,logLikelihood)
    !!{
    Returns true if the posterior sampling is converged (which it always is).
    !!}
    implicit none
    class           (posteriorSampleConvergenceAlways), intent(inout)           :: self
    class           (posteriorSampleStateClass       ), intent(inout), optional :: simulationState
    double precision                                  , intent(in   ), optional :: logLikelihood
    !$GLC attributes unused :: self, simulationState, logLikelihood

    alwaysIsConverged=.true.
    return
  end function alwaysIsConverged

  integer function alwaysConvergedAtStep(self)
    !!{
    Return the step at which the simulation converged.
    !!}
    implicit none
    class(posteriorSampleConvergenceAlways), intent(inout) :: self
    !$GLC attributes unused :: self

    alwaysConvergedAtStep=1
    return
  end function alwaysConvergedAtStep

  subroutine alwaysReset(self)
    !!{
    Reset the convergence object.
    !!}
    implicit none
    class(posteriorSampleConvergenceAlways), intent(inout) :: self
    !$GLC attributes unused :: self

    return
  end subroutine alwaysReset

  subroutine alwaysLogReport(self,fileUnit)
    !!{
    Write a convergence report to the given {\normalfont \ttfamily fileUnit}.
    !!}
    implicit none
    class  (posteriorSampleConvergenceAlways), intent(inout) :: self
    integer                                  , intent(in   ) :: fileUnit
    !$GLC attributes unused :: self

    write (fileUnit,*) 'Convergence: converged'
    return
  end subroutine alwaysLogReport

  logical function alwaysStateIsOutlier(self,stateIndex)
    !!{
    Return true if the specified chain is deemed to be an outlier. In this case, chains are never outliers.
    !!}
    implicit none
    class  (posteriorSampleConvergenceAlways), intent(inout) :: self
    integer                                  , intent(in   ) :: stateIndex
    !$GLC attributes unused :: self, stateIndex

    alwaysStateIsOutlier=.false.
    return
  end function alwaysStateIsOutlier
