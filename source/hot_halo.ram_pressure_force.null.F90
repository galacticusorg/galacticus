!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements a null ram pressure force for hot halos.

module Hot_Halo_Ram_Pressure_Force_Null
  !% Implements a null ram pressure force for hot halos.
  implicit none
  private
  public :: Hot_Halo_Ram_Pressure_Force_Null_Initialize

contains

  !# <hotHaloRamPressureForceMethod>
  !#  <unitName>Hot_Halo_Ram_Pressure_Force_Null_Initialize</unitName>
  !# </hotHaloRamPressureForceMethod>
  subroutine Hot_Halo_Ram_Pressure_Force_Null_Initialize(hotHaloRamPressureForceMethod,Hot_Halo_Ram_Pressure_Force_Get)
    !% Initializes the ``Null'' hot halo ram pressure stripping module.
    use ISO_Varying_String
    implicit none
    type     (varying_string                      ), intent(in   )          :: hotHaloRamPressureForceMethod
    procedure(Hot_Halo_Ram_Pressure_Force_Null_Get), intent(inout), pointer :: Hot_Halo_Ram_Pressure_Force_Get

    if (hotHaloRamPressureForceMethod == 'null') Hot_Halo_Ram_Pressure_Force_Get => Hot_Halo_Ram_Pressure_Force_Null_Get
    return
  end subroutine Hot_Halo_Ram_Pressure_Force_Null_Initialize

  double precision function Hot_Halo_Ram_Pressure_Force_Null_Get(thisNode)
    !% Computes the ram pressure force from the hot halo in the {\tt null} implementation. Always returns zero.
    use Galacticus_Nodes
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    Hot_Halo_Ram_Pressure_Force_Null_Get=0.0d0
    return
  end function Hot_Halo_Ram_Pressure_Force_Null_Get

end module Hot_Halo_Ram_Pressure_Force_Null
