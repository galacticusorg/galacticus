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

!% Defines methods for the ``correlationLengthCount'' stopping criterion class.

function stoppingCriterionCorrelationLengthCountConstructor(stopAfterCount)
  !% Constructor for ``correlationLengthCount'' stopping criterion class.
  implicit none
  type   (stoppingCriterionCorrelationLengthCount)                :: stoppingCriterionCorrelationLengthCountConstructor
  integer                                         , intent(in   ) :: stopAfterCount
  
  stoppingCriterionCorrelationLengthCountConstructor%stopAfterCount=stopAfterCount
  return
end function stoppingCriterionCorrelationLengthCountConstructor

logical function stoppingCriterionCorrelationLengthCountStop(self,simulationState,simulationConvergence)
  !% Return a Boolean indicating if the simulation should stop.
  use Galacticus_Error
  implicit none
  class(stoppingCriterionCorrelationLengthCount), intent(in   ) :: self
  class(state                                  ), intent(inout) :: simulationState
  class(convergence                            ), intent(inout) :: simulationConvergence
  
  select type (simulationState)
  class is (stateCorrelation)
     stoppingCriterionCorrelationLengthCountStop=                     &
          &  simulationConvergence%isConverged                     () &
          & .and.                                                     &
          &   simulationState      %postConvergenceCorrelationCount() &
          &  >                                                        &
          &   self%stopAfterCount
  class default
     call Galacticus_Error_Report('stoppingCriterionCorrelationLengthCountStop','state object does not support correlation length count')
  end select
  return
end function stoppingCriterionCorrelationLengthCountStop
