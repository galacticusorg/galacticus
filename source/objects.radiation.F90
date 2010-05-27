!% Contains a module which defines the radiation structure data type, used to describe radiation fields. (Currently a dummy
!% implementation.)

module Radiation_Structure
  !% Defines the radiation structure data type, used to describe radiation fields. (Currently a dummy implementation.)
  private
  public :: radiationStructure

  type radiationStructure
     !% The radiation structure data type, used to describe radiation fields. (Currently a dummy implementation.)
     logical, private :: dummy
     contains
       procedure :: set => Radiation_Set
  end type radiationStructure

  ! Option labels.
  integer, public, parameter :: noRadiation=0

contains

  subroutine Radiation_Set(radiation,setOption)
    !% Set the {\tt radiation} field as specified.
    implicit none
    type(radiationStructure), intent(inout)          :: radiation
    integer,                  intent(in),   optional :: setOption

    ! AJB:: Currently does nothing as we don't support radiation structures yet.

    return
  end subroutine Radiation_Set

end module Radiation_Structure
