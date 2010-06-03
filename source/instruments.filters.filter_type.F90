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






!% Contains a module which defines a structure to hold filter response curves.

module Instruments_Filters_Type
  !% Defines a structure to hold filter response curves.
  use ISO_Varying_String
  use FGSL
  private
  public :: filterType
  
  type filterType
     !% A structure which holds filter response curves.
     integer                                         :: nPoints
     double precision,     allocatable, dimension(:) :: wavelength,response
     type(varying_string)                            :: name
     ! Interpolation structures.
     logical                                         :: reset=.true.
     type(fgsl_interp_accel)                         :: interpolationAccelerator
     type(fgsl_interp)                               :: interpolationObject
  end type filterType
  
end module Instruments_Filters_Type
