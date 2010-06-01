!% Contains a module which defines a data type for storing atomic data.

module Atomic_Data_Type
  !% Defines a data type for storing atomic data.
  private
  public :: atomicData
  
  type atomicData
     !% Data type for storing atomic data.
     integer                                     :: atomicNumber
     double precision                            :: atomicMass
     double precision, allocatable, dimension(:) :: abundanceByMass
     character(len=3)                            :: shortLabel
     character(len=20)                           :: name
  end type atomicData

end module Atomic_Data_Type
