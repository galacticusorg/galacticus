!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements a fixed star formation timescale for galactic disks.

module Star_Formation_Timescale_Disks_Fixed
  !% Implements a fixed star formation timescale for galactic disks.
  implicit none
  private
  public :: Star_Formation_Timescale_Disks_Fixed_Initialize

  ! Parameters of the timescale model.
  double precision :: starFormationTimescaleDisksFixedTimescale  
                                                              
contains

  !# <starFormationTimescaleDisksMethod>
  !#  <unitName>Star_Formation_Timescale_Disks_Fixed_Initialize</unitName>
  !# </starFormationTimescaleDisksMethod>
  subroutine Star_Formation_Timescale_Disks_Fixed_Initialize(starFormationTimescaleDisksMethod,Star_Formation_Timescale_Disk_Get)
    !% Initializes the ``fixed'' disk star formation timescale module.
    use ISO_Varying_String
    use Input_Parameters
    use Galacticus_Error
    use Galacticus_Nodes
    implicit none
    type     (varying_string                     ), intent(in   )          :: starFormationTimescaleDisksMethod  
    procedure(Star_Formation_Timescale_Disk_Fixed), intent(inout), pointer :: Star_Formation_Timescale_Disk_Get  
                                                                                                              
    if (starFormationTimescaleDisksMethod == 'fixed') then
       Star_Formation_Timescale_Disk_Get => Star_Formation_Timescale_Disk_Fixed
       ! Get parameters of for the timescale calculation.
       !@ <inputParameter>
       !@   <name>starFormationTimescaleDisksFixedTimescale</name>
       !@   <defaultValue>1 Gyr</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The timescale for star formation in the fixed timescale model for disks.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@   <group>starFormation</group>
       !@ </inputParameter>
       call Get_Input_Parameter('starFormationTimescaleDisksFixedTimescale',starFormationTimescaleDisksFixedTimescale,defaultValue=1.0d0)
    end if
    return
  end subroutine Star_Formation_Timescale_Disks_Fixed_Initialize

  double precision function Star_Formation_Timescale_Disk_Fixed(thisNode)
    !% Returns the timescale (in Gyr) for star formation in the galactic disk of {\tt thisNode}, assuming a fixed timecale.
    use Galacticus_Nodes
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode  
    
    ! Return the timescale.                                                 
    Star_Formation_Timescale_Disk_Fixed=starFormationTimescaleDisksFixedTimescale
    return
  end function Star_Formation_Timescale_Disk_Fixed
  
end module Star_Formation_Timescale_Disks_Fixed
