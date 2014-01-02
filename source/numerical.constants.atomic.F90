!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module of useful atomic constants.

module Numerical_Constants_Atomic
  !% Contains various useful atomic constants.
  use FGSL
  implicit none
  public

  ! Atomic mass unit (in kg).
  double precision, parameter :: atomicMassUnit              =FGSL_CONST_MKSA_UNIFIED_ATOMIC_MASS

  ! Atomic masses.
  double precision, parameter :: atomicMassHydrogen          =1.007825d0
  double precision, parameter :: atomicMassHelium            =4.002602d0

  ! Mass of hydrogen atom (in kg).
  double precision, parameter :: massHydrogenAtom            =atomicMassHydrogen*atomicMassUnit

  ! Ionization energies/wavelengths (in eV/Angstroms).
  double precision, parameter :: ionizationWavelengthHydrogen=911.862d0                           !  (Allen's Astrophysical Quantities, p. 36).

end module Numerical_Constants_Atomic
