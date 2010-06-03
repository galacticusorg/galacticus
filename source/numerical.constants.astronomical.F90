!! Copyright 2009, 2010, Andrew Benson <abenson@caltech.edu>
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






!% Contains a module of useful astronomical constants.

module Numerical_Constants_Astronomical
  !% Contains various useful astronomical constants.
  use FGSL
  use Numerical_Constants_Prefixes
  use Numerical_Constants_Atomic
  public
  
  ! Solar mass (in kg).
  double precision, parameter :: massSolar=FGSL_CONST_MKSA_SOLAR_MASS

  ! Solar luminosity (in W; Allen's Astrophysical Quantities, page 340).
  double precision, parameter :: luminositySolar=3.845d26

  ! Solar composition (Allen's Atrophysical Quantities, page 28).
  double precision, parameter :: hydrogenByMassSolar=0.707d0
  double precision, parameter :: heliumByMassSolar  =0.274d0
  double precision, parameter :: metallicitySolar   =0.0188d0

  ! Primordial composition.
  double precision, parameter :: hydrogenByMassPrimordial=0.778d0
  double precision, parameter :: heliumByMassPrimordial  =0.222d0
  double precision, parameter :: metallicityPrimordial   =5.36d-10
  double precision, parameter :: meanAtomicMassPrimordial=1.0d0/(2.0d0*hydrogenByMassPrimordial/atomicMassHydrogen+3.0d0&
       &*heliumByMassPrimordial/atomicMassHelium)

  ! Megaparsec (in m).
  double precision, parameter :: parsec    =FGSL_CONST_MKSA_PARSEC
  double precision, parameter :: megaParsec=mega*parsec

  ! Years and related quantities (in s).
  double precision, parameter :: year=31558149.8d0 ! Sidereal year.
  double precision, parameter :: gigaYear=giga*year

end module Numerical_Constants_Astronomical
