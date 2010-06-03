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
