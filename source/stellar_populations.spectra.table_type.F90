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
