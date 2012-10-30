!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements selection of stellar {\gls{imf}}s with one \gls{imf} for disks and another for spheroids.

module Star_Formation_IMF_Select_Disk_Spheroid
  !% Implements selection of stellar {\gls{imf}}s with one \gls{imf} for disks and another for spheroids.
  implicit none
  private
  public :: IMF_Select_Disk_Spheroid_Initialize

  ! Store the indices of the selected IMFs.
  integer :: imfSelectedDiskIndex,imfSelectedSpheroidIndex

contains

  !# <imfSelectionMethod>
  !#  <unitName>IMF_Select_Disk_Spheroid_Initialize</unitName>
  !# </imfSelectionMethod>
  subroutine IMF_Select_Disk_Spheroid_Initialize(imfSelectionMethod,IMF_Select_Do,imfNames)
    !% Initializes the ``diskSpheroid'' IMF selection module.
    use ISO_Varying_String
    use Input_Parameters
    use Star_Formation_IMF_Utilities
    implicit none
    type(varying_string),          intent(in)    :: imfSelectionMethod,imfNames(:)
    procedure(integer),   pointer, intent(inout) :: IMF_Select_Do
    type(varying_string)                         :: imfSelectionDisk,imfSelectionSpheroid
    
    if (imfSelectionMethod == 'diskSpheroid') then
       IMF_Select_Do => IMF_Select_Disk_Spheroid
       ! Get IMF choices.
       !@ <inputParameter>
       !@   <name>imfSelectionDisk</name>
       !@   <defaultValue>Salpeter</defaultValue>       
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of the initial mass function to use in the ``diskSpheroid initial mass function'' module for star formation in disks.
       !@   </description>
       !@   <type>string</type>
       !@   <cardinality>1</cardinality>
       !@   <group>initialMassFunction</group>
       !@ </inputParameter>
       call Get_Input_Parameter('imfSelectionDisk',imfSelectionDisk,defaultValue='Salpeter')
       imfSelectedDiskIndex=IMF_Index_Lookup(imfSelectionDisk,imfNames)
       !@ <inputParameter>
       !@   <name>imfSelectionSpheroid</name>
       !@   <defaultValue>Salpeter</defaultValue>       
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of the initial mass function to use in the ``diskSpheroid initial mass function'' module for star formation in spheroids.
       !@   </description>
       !@   <type>string</type>
       !@   <cardinality>1</cardinality>
       !@   <group>initialMassFunction</group>
       !@ </inputParameter>
       call Get_Input_Parameter('imfSelectionSpheroid',imfSelectionSpheroid,defaultValue='Salpeter')
       imfSelectedSpheroidIndex=IMF_Index_Lookup(imfSelectionSpheroid,imfNames)
    end if
    return
  end subroutine IMF_Select_Disk_Spheroid_Initialize

  integer function IMF_Select_Disk_Spheroid(starFormationRate,fuelAbundances,component)
    !% Return our selection of stellar initial mass function.
    use Abundances_Structure
    use Galacticus_Error
    use Galactic_Structure_Options
    implicit none
    double precision,          intent(in) :: starFormationRate
    type(abundances), intent(in) :: fuelAbundances
    integer,                   intent(in) :: component

    ! Select between disk and spheroid IMFs.
    select case (component)
    case (componentTypeDisk)
       IMF_Select_Disk_Spheroid=imfSelectedDiskIndex
    case (componentTypeSpheroid)
       IMF_Select_Disk_Spheroid=imfSelectedSpheroidIndex
    case default
       call Galacticus_Error_Report('IMF_Select_Disk_Spheroid','only disk and spheroid components are allowed')
    end select
    return
  end function IMF_Select_Disk_Spheroid
  
end module Star_Formation_IMF_Select_Disk_Spheroid
