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






!% Contains a module which provides various internal option codes for the galactic structure functions.

module Galactic_Structure_Options
  !% Provides various internal option codes for the galactic structure functions.
  public

  ! Values used to represent different mass types.
  integer,          parameter :: massTypeAll     =0
  integer,          parameter :: massTypeDark    =1
  integer,          parameter :: massTypeBaryonic=2
  integer,          parameter :: massTypeGalactic=3
  integer,          parameter :: massTypeGaseous =4
  integer,          parameter :: massTypeStellar =5

  ! Values used to represent different component types.
  integer,          parameter :: componentTypeAll     =0
  integer,          parameter :: componentTypeDisk    =1
  integer,          parameter :: componentTypeSpheroid=2
  integer,          parameter :: componentTypeHotHalo =3

  ! Coordinate system options.
  integer,          parameter  :: coordinateSystemSpherical  =1
  integer,          parameter  :: coordinateSystemCylindrical=2
  integer,          parameter  :: coordinateSystemCartesian  =3

  ! Suitably large value to represent infinite radius.
  double precision, parameter :: radiusLarge=1.0d10

end module Galactic_Structure_Options
