!! Copyright 2009, 2010, 2011 Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


!% Contains a module which implements selection of stellar {\IMF}s with one \IMF\ for disks and another for spheroids.

module Star_Formation_IMF_Select_Disk_Spheroid
  !% Implements selection of stellar {\IMF}s with one \IMF\ for disks and another for spheroids.
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
    type(abundancesStructure), intent(in) :: fuelAbundances
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
