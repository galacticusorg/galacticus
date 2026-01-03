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
  Implementation of a posterior sampling stopping class which never stops.
  !!}

  !![
  <posteriorSampleStoppingCriterion name="posteriorSampleStoppingCriterionNever">
   <description>A posterior sampling stopping class which never stops.</description>
  </posteriorSampleStoppingCriterion>
  !!]
  type, extends(posteriorSampleStoppingCriterionClass) :: posteriorSampleStoppingCriterionNever
     !!{
     Implementation of a posterior sampling convergence class which never converges.
     !!}
     private
   contains
     procedure :: stop => neverStop
  end type posteriorSampleStoppingCriterionNever

  interface posteriorSampleStoppingCriterionNever
     !!{
     Constructors for the \refClass{posteriorSampleStoppingCriterionNever} posterior sampling convergence class.
     !!}
     module procedure neverConstructorParameters
  end interface posteriorSampleStoppingCriterionNever

contains

  function neverConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{posteriorSampleStoppingCriterionNever} posterior sampling stopping class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(posteriorSampleStoppingCriterionNever)                :: self
    type(inputParameters                      ), intent(inout) :: parameters

    self=posteriorSampleStoppingCriterionNever()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function neverConstructorParameters

  logical function neverStop(self,simulationState)
    !!{
    Returns true if the posterior sampling should stop (which it never should).
    !!}
    implicit none
    class(posteriorSampleStoppingCriterionNever), intent(inout) :: self
    class(posteriorSampleStateClass            ), intent(inout) :: simulationState
    !$GLC attributes unused :: self, simulationState

    neverStop=.false.
    return
  end function neverStop
