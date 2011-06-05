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
