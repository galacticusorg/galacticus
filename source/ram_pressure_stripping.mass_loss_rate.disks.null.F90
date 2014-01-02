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

!% Contains a module which implements a null mass loss rates from disks due to ram pressure stripping.

module Ram_Pressure_Stripping_Mass_Loss_Rate_Disks_Null
  !% Implements a null mass loss rates from disks due to ram pressure stripping.
  implicit none
  private
  public :: Ram_Pressure_Stripping_Mass_Loss_Rate_Disks_Null_Init

contains

  !# <ramPressureStrippingMassLossRateDisksMethod>
  !#  <unitName>Ram_Pressure_Stripping_Mass_Loss_Rate_Disks_Null_Init</unitName>
  !# </ramPressureStrippingMassLossRateDisksMethod>
  subroutine Ram_Pressure_Stripping_Mass_Loss_Rate_Disks_Null_Init(ramPressureStrippingMassLossRateDisksMethod,Ram_Pressure_Stripping_Mass_Loss_Rate_Disk_Get)
    !% Initializes the ``null'' ram pressure stripping mass loss rate from disks module.
    use ISO_Varying_String
    implicit none
    type     (varying_string                                 ), intent(in   )          :: ramPressureStrippingMassLossRateDisksMethod
    procedure(Ram_Pressure_Stripping_Mass_Loss_Rate_Disk_Null), intent(inout), pointer :: Ram_Pressure_Stripping_Mass_Loss_Rate_Disk_Get

    if (ramPressureStrippingMassLossRateDisksMethod == 'null') Ram_Pressure_Stripping_Mass_Loss_Rate_Disk_Get => Ram_Pressure_Stripping_Mass_Loss_Rate_Disk_Null
    return
  end subroutine Ram_Pressure_Stripping_Mass_Loss_Rate_Disks_Null_Init

  double precision function Ram_Pressure_Stripping_Mass_Loss_Rate_Disk_Null(thisNode)
    !% Computes the mass loss rate from disks due to ram pressure stripping. Always returns zero.
    use Galacticus_Nodes
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    Ram_Pressure_Stripping_Mass_Loss_Rate_Disk_Null=0.0d0
    return
  end function Ram_Pressure_Stripping_Mass_Loss_Rate_Disk_Null

end module Ram_Pressure_Stripping_Mass_Loss_Rate_Disks_Null
