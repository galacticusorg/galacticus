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






!% Contains a module of useful unit conversions.

module Numerical_Constants_Units
  !% Contains various useful unit conversions.
  use Numerical_Constants_Astronomical
  use Numerical_Constants_Prefixes
  public
  
  ! Ergs in Joules.
  double precision, parameter :: ergs=1.0d-7

  ! Conversion from Mpc/(km/s) to Gyr.
  double precision, parameter :: Mpc_per_km_per_s_To_Gyr=megaParsec/kilo/gigaYear

end module Numerical_Constants_Units
