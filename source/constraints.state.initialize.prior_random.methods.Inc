!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Defines methods for a prior random state initializor.

function stateInitializorPriorRandomConstructor()
  !% Constructor for the prior random state initializor class.
  implicit none
  type(stateInitializorPriorRandom) :: stateInitializorPriorRandomConstructor

  return
end function stateInitializorPriorRandomConstructor

subroutine stateInitializorPriorRandomInitialize(self,simulationState,parameterPriors)
  !% Initialize simulation state by drawing at random from the parameter priors.
  implicit none
  class(stateInitializorPriorRandom), intent(inout)               :: self
  class(state                      ), intent(inout)               :: simulationState
  type (prior                      ), intent(inout), dimension(:) :: parameterPriors

  ! Initialize chain to some state vector.
  call priorsSample(parameterPriors,simulationState)
  return
end subroutine stateInitializorPriorRandomInitialize
