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
  Implementation of a posterior sampling convergence class which never converges.
  !!}

  !![
  <posteriorSampleConvergence name="posteriorSampleConvergenceNever">
   <description>
    This option assumes that the simulation never converges, and so the calculation will run indefinitely. It is intended primarily
    for testing purposes.
   </description>
  </posteriorSampleConvergence>
  !!]
  type, extends(posteriorSampleConvergenceClass) :: posteriorSampleConvergenceNever
     !!{
     Implementation of a posterior sampling convergence class which never converges.
     !!}
     private
   contains
     procedure :: isConverged     => neverIsConverged
     procedure :: convergedAtStep => neverConvergedAtStep
     procedure :: reset           => neverReset
     procedure :: logReport       => neverLogReport
     procedure :: stateIsOutlier  => neverStateIsOutlier
  end type posteriorSampleConvergenceNever

  interface posteriorSampleConvergenceNever
     !!{
     Constructors for the \refClass{posteriorSampleConvergenceNever} posterior sampling convergence class.
     !!}
     module procedure neverConstructorParameters
  end interface posteriorSampleConvergenceNever

contains

  function neverConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{posteriorSampleConvergenceNever} merger tree halo mass function sampling class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(posteriorSampleConvergenceNever)                :: self
    type(inputParameters                ), intent(inout) :: parameters

    self=posteriorSampleConvergenceNever()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function neverConstructorParameters

  logical function neverIsConverged(self,simulationState,logLikelihood)
    !!{
    Returns true if the posterior sampling is converged (which it never is).
    !!}
    implicit none
    class           (posteriorSampleConvergenceNever), intent(inout)           :: self
    class           (posteriorSampleStateClass      ), intent(inout), optional :: simulationState
    double precision                                 , intent(in   ), optional :: logLikelihood
    !$GLC attributes unused :: self, simulationState, logLikelihood

    neverIsConverged=.false.
    return
  end function neverIsConverged

  integer function neverConvergedAtStep(self)
    !!{
    Return the step at which the simulation converged.
    !!}
    implicit none
    class(posteriorSampleConvergenceNever), intent(inout) :: self
    !$GLC attributes unused :: self

    neverConvergedAtStep=-1
    return
  end function neverConvergedAtStep

  subroutine neverReset(self)
    !!{
    Reset the convergence object.
    !!}
    implicit none
    class(posteriorSampleConvergenceNever), intent(inout) :: self
    !$GLC attributes unused :: self

    return
  end subroutine neverReset

  subroutine neverLogReport(self,fileUnit)
    !!{
    Write a convergence report to the given {\normalfont \ttfamily fileUnit}.
    !!}
    implicit none
    class  (posteriorSampleConvergenceNever), intent(inout) :: self
    integer                                 , intent(in   ) :: fileUnit
    !$GLC attributes unused :: self

    write (fileUnit,*) 'Convergence: unconverged'
    return
  end subroutine neverLogReport

  logical function neverStateIsOutlier(self,stateIndex)
    !!{
    Return true if the specified chain is deemed to be an outlier. In this case, chains are never outliers.
    !!}
    implicit none
    class  (posteriorSampleConvergenceNever), intent(inout) :: self
    integer                                 , intent(in   ) :: stateIndex
    !$GLC attributes unused :: self, stateIndex

    neverStateIsOutlier=.false.
    return
  end function neverStateIsOutlier
