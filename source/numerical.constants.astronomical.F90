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

!% Contains a module of useful astronomical constants.

module Numerical_Constants_Astronomical
  !% Contains various useful astronomical constants.
  use Numerical_Constants_Atomic
  use Numerical_Constants_Units
  implicit none
  public

  ! Solar mass (in kg).
  double precision, parameter :: massSolar                 =FGSL_CONST_MKSA_SOLAR_MASS

  ! Solar radius (in m; Allen's Astrophysical Quantities, page 340).
  double precision, parameter :: radiusSolar               =6.95508d8

  ! Solar luminosity (in W; Allen's Astrophysical Quantities, page 340).
  double precision, parameter :: luminositySolar           =3.845d26

  ! Solar composition (Allen's Atrophysical Quantities, page 28).
  double precision, parameter :: hydrogenByMassSolar       =0.707d0
  double precision, parameter :: heliumByMassSolar         =0.274d0
  double precision, parameter :: metallicitySolar          =0.0188d0

  ! Primordial composition.
  double precision, parameter :: hydrogenByMassPrimordial  =0.778d0
  double precision, parameter :: heliumByMassPrimordial    =0.222d0
  double precision, parameter :: metallicityPrimordial     =5.36d-10
  double precision, parameter :: meanAtomicMassPrimordial  =1.0d0/(2.0d0*hydrogenByMassPrimordial/atomicMassHydrogen+3.0d0*heliumByMassPrimordial/atomicMassHelium)

  ! Parsec and related quantities (in m).
  double precision, parameter :: parsec                    =FGSL_CONST_MKSA_PARSEC
  double precision, parameter :: kiloParsec                =kilo*parsec
  double precision, parameter :: megaParsec                =mega*parsec

  ! Years and related quantities (in s).
  double precision, parameter :: year                      =31558149.8d0                                                                                            !   Sidereal year.
  double precision, parameter :: gigaYear                  =giga*year

  ! Conversion from Mpc/(km/s) to Gyr.
  double precision, parameter :: Mpc_per_km_per_s_To_Gyr   =megaParsec/kilo/gigaYear

  ! Conversion from magnitudes to optical depth.
  double precision, parameter :: magnitudesPerOpticalDepth = 2.5d0/log(10.0d0)

  ! AB magnitude system:
  ! The AB magnitude system is defined such that:
  !   m = -2.5log10(F_nu/[ergs/s/cm^2/Hz])-48.57
  ! Computing the flux at 10pc gives us the zero point for the absolute magnitude scale (units of W/Hz).
  double precision, parameter :: offsetAB                =48.57d0
  double precision, parameter :: luminosityZeroPointAB   =(10.0d0**(-offsetAB/2.5d0))*4.0d0*Pi*((10.0d0*parsec*hecto)**2)*ergs

  ! Anglular conversions.
  double precision, parameter :: hoursToRadians          =2.0d0*Pi/ 24.0d0
  double precision, parameter :: degreesToRadians        =2.0d0*Pi/360.0d0

end module Numerical_Constants_Astronomical
