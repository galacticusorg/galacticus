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

!% Defines methods for the ``fixed'' differential evolution proposal size temperature exponent class.

function dePropSizeTempExpFixedConstructor(exponent)
  !% Constructor for the ``fixed'' differential evolution proposal size temperature exponent class.
  implicit none
  type            (dePropSizeTempExpFixed)                :: dePropSizeTempExpFixedConstructor
  double precision                        , intent(in   ) :: exponent

  dePropSizeTempExpFixedConstructor%exponentValue=exponent
  return
end function dePropSizeTempExpFixedConstructor

double precision function dePropSizeTempExpFixedExponent(self,temperedStates,temperatures,simulationState,simulationConvergence)
  !% Return the fixed differential evolution proposal size temperature exponent.
  implicit none
  class           (dePropSizeTempExpFixed), intent(inout)               :: self
  class           (state                 ), intent(in   )               :: simulationState
  class           (state                 ), intent(in   ), dimension(:) :: temperedStates
  double precision                        , intent(in   ), dimension(:) :: temperatures
  class           (convergence           ), intent(inout)               :: simulationConvergence

  dePropSizeTempExpFixedExponent=self%exponentValue
  return
end function dePropSizeTempExpFixedExponent
