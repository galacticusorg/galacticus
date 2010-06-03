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






!% Contains a module which defines a structure to hold tables of simple stellar population luminosities

module Stellar_Population_Luminosities_Table_Type
  !% Defines a structure to hold tables of simple stellar population luminosities
  use FGSL
  private
  public :: luminosityTable
  
  type luminosityTable
     !% Structure for holding tables of simple stellar population luminosities.
     integer                                         :: agesCount,metallicitiesCount
     logical,          allocatable, dimension(:)     :: isTabulated
     double precision, allocatable, dimension(:)     :: age,metallicity
     double precision, allocatable, dimension(:,:,:) :: luminosity
     ! Interpolation structures.
     logical                                         :: resetAge=.true., resetMetallicity=.true.
     type(fgsl_interp_accel)                         :: interpolationAcceleratorAge,interpolationAcceleratorMetallicity
    end type luminosityTable
  
end module Stellar_Population_Luminosities_Table_Type
