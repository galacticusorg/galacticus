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
  Implementation of a zero rate stellar feedback model.
  !!}

  !![
  <stellarFeedbackOutflows name="stellarFeedbackOutflowsZero">
   <description>A zero rate stellar feedback model.</description>
  </stellarFeedbackOutflows>
  !!]
  type, extends(stellarFeedbackOutflowsClass) :: stellarFeedbackOutflowsZero
     !!{
     Implementation of a zero rate stellar feedback model.
     !!}
     private
   contains
     procedure :: outflowRate => zeroOutflowRate
  end type stellarFeedbackOutflowsZero

  interface stellarFeedbackOutflowsZero
     !!{
     Constructors for the zero rate stellar feedback model.
     !!}
     module procedure zeroConstructorParameters
  end interface stellarFeedbackOutflowsZero

contains

  function zeroConstructorParameters(parameters) result(self)
    !!{
    Constructor for the stellar feedback class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(stellarFeedbackOutflowsZero)                :: self
    type(inputParameters            ), intent(inout) :: parameters

    self=stellarFeedbackOutflowsZero()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function zeroConstructorParameters

  subroutine zeroOutflowRate(self,component,rateStarFormation,rateEnergyInput,rateOutflowEjective,rateOutflowExpulsive)
    !!{
    Returns a zero outflow rate from disks
    !!}
    implicit none
    class           (stellarFeedbackOutflowsZero), intent(inout) :: self
    class           (nodeComponent              ), intent(inout) :: component
    double precision                             , intent(in   ) :: rateEnergyInput    , rateStarFormation
    double precision                             , intent(  out) :: rateOutflowEjective, rateOutflowExpulsive
    !$GLC attributes unused :: self, component, rateEnergyInput, rateStarFormation

    rateOutflowEjective=0.0d0
    rateOutflowExpulsive=0.0d0
    return
  end subroutine zeroOutflowRate
