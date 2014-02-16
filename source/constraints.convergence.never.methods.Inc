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

!% Defines methods for the ``never'' convergence class.

function convergenceNeverConstructor()
  !% Constructor for ``never'' convergence class.
  implicit none
  type(convergenceNever) :: convergenceNeverConstructor

  return
end function convergenceNeverConstructor

logical function convergenceNeverIsConverged(self,simulationState,logLikelihood)
  !% Return whether the simulation is converged.
  use Constraints_State
  implicit none
  class           (convergenceNever), intent(inout)           :: self
  class           (state           ), intent(in   ), optional :: simulationState
  double precision                  , intent(in   ), optional :: logLikelihood

  convergenceNeverIsConverged=.false.
  return
end function convergenceNeverIsConverged

subroutine convergenceNeverReset(self)
  !% Reset the convergence object.
  implicit none
  class(convergenceNever), intent(inout) :: self

  return
end subroutine convergenceNeverReset

subroutine convergenceNeverLogReport(self,fileUnit)
  !% Write a convergence report to the given {\tt fileUnit}.
  implicit none
  class  (convergenceNever), intent(inout) :: self
  integer                  , intent(in   ) :: fileUnit

  write (fileUnit,*) 'Convergence: unconverged'
  return
end subroutine convergenceNeverLogReport
