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






!% Contains a module which defines the data structure for simple stellar population spectral tables.

module Stellar_Population_Spectra_Table
  !% Defines the data structure for simple stellar population spectral tables.
  use FGSL
  private
  public :: spectralTable

  type spectralTable
     !% Structure to hold spectral data.
     ! The spectra tables.
     integer                                         :: stellarPopulationSpectraAgesNumberPoints &
          &,stellarPopulationSpectraMetallicityNumberPoints,stellarPopulationSpectraWavelengthsNumberPoints
     double precision, allocatable, dimension(:)     :: stellarPopulationSpectraMetallicities,stellarPopulationSpectraAges &
          &,stellarPopulationSpectraWavelengths
     double precision, allocatable, dimension(:,:,:) :: stellarPopulationSpectraTable
     
     ! Interpolation structures.
     logical                                         :: resetAge=.true., resetMetallicity=.true., resetWavelength=.true.
     type(fgsl_interp_accel)                         :: interpolationAcceleratorAge,interpolationAcceleratorMetallicity &
          &,interpolationAcceleratorWavelength
  end type spectralTable

end module Stellar_Population_Spectra_Table
