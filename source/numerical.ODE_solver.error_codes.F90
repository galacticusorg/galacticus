!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
Contains a module which defines internal error codes for the \glc\ ODE solver.
!!}

module ODE_Solver_Error_Codes
  !!{
  Defines internal error codes for the \glc\ ODE solver.
  !!}
  use :: Interface_GSL, only : GSL_EBadFunc
  implicit none
  public

  ! An interrupt has been triggered.
  integer         , parameter :: odeSolverInterrupt=GSL_EBadFunc

  ! Point during the integration at which an interrupt occurred.
  double precision            :: interruptedAtX
  !$omp threadprivate(interruptedAtX)
end module ODE_Solver_Error_Codes
