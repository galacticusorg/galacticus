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
