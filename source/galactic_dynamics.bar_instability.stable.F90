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
  Implementation of a perfectly stable model for galactic disk bar instability.
  !!}

  !![
  <galacticDynamicsBarInstability name="galacticDynamicsBarInstabilityStable">
   <description>
    A galactic dynamics bar instability class which assumes perfect stability for galactic disks and so returns an infinite
    timescale, and no external driving torque.
   </description>
  </galacticDynamicsBarInstability>
  !!]
  type, extends(galacticDynamicsBarInstabilityClass) :: galacticDynamicsBarInstabilityStable
     !!{
     Implementation of a perfectly stable model for galactic disk bar instability.
     !!}
     private
   contains
     procedure :: timescale => stableTimescale
  end type galacticDynamicsBarInstabilityStable

  interface galacticDynamicsBarInstabilityStable
     !!{
     Constructors for the \refClass{galacticDynamicsBarInstabilityStable} model for galactic disk bar instability class.
     !!}
     module procedure stableConstructorParameters
  end interface galacticDynamicsBarInstabilityStable

contains

  function stableConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{galacticDynamicsBarInstabilityStable} model for galactic disk bar instability class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(galacticDynamicsBarInstabilityStable)                :: self
    type(inputParameters                     ), intent(inout) :: parameters

    self=galacticDynamicsBarInstabilityStable()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function stableConstructorParameters

  subroutine stableTimescale(self,node,timescale,externalDrivingSpecificTorque,fractionAngularMomentumRetainedDisk,fractionAngularMomentumRetainedSpheroid)
    !!{
    Computes a timescale for depletion of a disk to a pseudo-bulge via bar instability based on the criterion of
    \cite{efstathiou_stability_1982}.
    !!}
    implicit none
    class           (galacticDynamicsBarInstabilityStable), intent(inout) :: self
    type            (treeNode                            ), intent(inout) :: node
    double precision                                      , intent(  out) :: externalDrivingSpecificTorque      , timescale                              , &
         &                                                                   fractionAngularMomentumRetainedDisk, fractionAngularMomentumRetainedSpheroid
    !$GLC attributes unused :: self, node

    ! Assume infinite timescale (i.e. no instability).
    timescale                      =-1.0d0
    ! Also assume no torque.
    externalDrivingSpecificTorque  =+0.0d0
    ! Fractions of angular momentum retained are arbitrary.
    fractionAngularMomentumRetainedDisk    =+1.0d0
    fractionAngularMomentumRetainedSpheroid=+1.0d0
    return
  end subroutine stableTimescale
