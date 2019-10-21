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

  !% Implementation of a zero expulsive outflow rate due to star formation feedback in galactic spheroids.

  !# <starFormationExpulsiveFeedbackSpheroids name="starFormationExpulsiveFeedbackSpheroidsZero">
  !#  <description>A zero expulsive outflow rate due to star formation feedback in galactic spheroids.</description>
  !# </starFormationExpulsiveFeedbackSpheroids>
  type, extends(starFormationExpulsiveFeedbackSpheroidsClass) :: starFormationExpulsiveFeedbackSpheroidsZero
     !% Implementation of a zero expulsive outflow rate due to star formation feedback in galactic spheroids.
     private
   contains
     procedure :: outflowRate => zeroOutflowRate
  end type starFormationExpulsiveFeedbackSpheroidsZero

  interface starFormationExpulsiveFeedbackSpheroidsZero
     !% Constructors for the superwind expulsive star formation feedback in spheroids class.
     module procedure zeroConstructorParameters
  end interface starFormationExpulsiveFeedbackSpheroidsZero

contains

  function zeroConstructorParameters(parameters) result(self)
    !% Constructor for the superwind expulsive star formation feedback in spheroids class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(starFormationExpulsiveFeedbackSpheroidsZero)                :: self
    type(inputParameters                            ), intent(inout) :: parameters
    !GCC$ attributes unused :: parameters

    self=starFormationExpulsiveFeedbackSpheroidsZero()
    return
  end function zeroConstructorParameters

  double precision function zeroOutflowRate(self,node,rateEnergyInput,rateStarFormation)
    !% Returns a zero expulsive outflow rate from spheroids
    implicit none
    class           (starFormationExpulsiveFeedbackSpheroidsZero), intent(inout) :: self
    type            (treeNode                                   ), intent(inout) :: node
    double precision                                             , intent(in   ) :: rateEnergyInput, rateStarFormation
    !GCC$ attributes unused :: self, node, rateEnergyInput, rateStarFormation

    zeroOutflowRate=0.0d0
    return
  end function zeroOutflowRate
