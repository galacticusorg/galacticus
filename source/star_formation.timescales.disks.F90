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






!% Contains a module which implements calculations of star formation timescales for galactic disks.

module Star_Formation_Timescales_Disks
  !% Implements calculations of star formation timescales for galactic disks.
  use ISO_Varying_String
  use Tree_Nodes
  private
  public :: Star_Formation_Timescale_Disk
  
  ! Flag to indicate if this module has been initialized.  
  logical              :: starFormationTimescaleDisksInitialized=.false.

  ! Name of cooling rate available method used.
  type(varying_string) :: starFormationTimescaleDisksMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Star_Formation_Timescale_Disk_Template), pointer :: Star_Formation_Timescale_Disk_Get => null()
  abstract interface
     double precision function Star_Formation_Timescale_Disk_Template(thisNode)
       import treeNode
       type(treeNode), intent(inout), pointer :: thisNode
     end function Star_Formation_Timescale_Disk_Template
  end interface

contains

  subroutine Star_Formation_Timescale_Disks_Initialize
    !% Initialize the disk star formation timecale module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="starFormationTimescaleDisksMethod" type="moduleUse">
    include 'star_formation.timescales.disks.modules.inc'
    !# </include>
    implicit none

    !$omp critical(Star_Formation_Timescale_Disks_Initialization) 
    ! Initialize if necessary.
    if (.not.starFormationTimescaleDisksInitialized) then
       ! Get the disk star formation timescale method parameter.
       !@ <inputParameter>
       !@   <name>starFormationTimescaleDisksMethod</name>
       !@   <defaultValue>dynamical time</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of the method to be used for computing star formation timescales in disks.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('starFormationTimescaleDisksMethod',starFormationTimescaleDisksMethod,defaultValue='dynamical time')
       ! Include file that makes calls to all available method initialization routines.
       !# <include directive="starFormationTimescaleDisksMethod" type="code" action="subroutine">
       !#  <subroutineArgs>starFormationTimescaleDisksMethod,Star_Formation_Timescale_Disk_Get</subroutineArgs>
       include 'star_formation.timescales.disks.inc'
       !# </include>
       if (.not.associated(Star_Formation_Timescale_Disk_Get)) call Galacticus_Error_Report('Star_Formation_Timescale_Disks'&
            &,'method ' //char(starFormationTimescaleDisksMethod)//' is unrecognized')
       starFormationTimescaleDisksInitialized=.true.
    end if
    !$omp end critical(Star_Formation_Timescale_Disks_Initialization) 

    return
  end subroutine Star_Formation_Timescale_Disks_Initialize

  double precision function Star_Formation_Timescale_Disk(thisNode)
    !% Returns the timescale (in Gyr) for star formation in the disk component of {\tt thisNode}.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode

    ! Initialize the module.
    call Star_Formation_Timescale_Disks_Initialize

    ! Get the energy using the selected method.
    Star_Formation_Timescale_Disk=Star_Formation_Timescale_Disk_Get(thisNode)

    return
  end function Star_Formation_Timescale_Disk
  
end module Star_Formation_Timescales_Disks
