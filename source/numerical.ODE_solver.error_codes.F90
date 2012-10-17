!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which defines internal error codes for the \glc\ ODE solver.

module ODE_Solver_Error_Codes
  !% Defines internal error codes for the \glc\ ODE solver.
  use FGSL
  implicit none
  public

  ! An interrupt has been triggered.
  integer,         parameter :: odeSolverInterrupt=FGSL_EBADFUNC

  ! Point during the integration at which an interrupt occurred.
  double precision           :: interruptedAtX
  !$omp threadprivate(interruptedAtX)

end module ODE_Solver_Error_Codes
