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


!% Contains a module which implements calculations of the core radius in cored isothermal hot halo profiles.

module Hot_Halo_Density_Cored_Isothermal_Core_Radii
  !% Implements calculations of the core radius in cored isothermal hot halo profiles.
  use ISO_Varying_String
  use Tree_Nodes
  implicit none
  private
  public :: Hot_Halo_Density_Cored_Isothermal_Core_Radius
  
  ! Flag to indicate if this module has been initialized.  
  logical              :: hotHaloCoredIsothermalCoreRadiiInitialized=.false.

  ! Name of cooling rate available method used.
  type(varying_string) :: hotHaloCoredIsothermalCoreRadiiMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Hot_Halo_Density_Cored_Isothermal_Core_Radius_Template), pointer :: Hot_Halo_Density_Cored_Isothermal_Core_Radius_Get => null()
  abstract interface
     double precision function Hot_Halo_Density_Cored_Isothermal_Core_Radius_Template(thisNode)
       import treeNode
       type(treeNode), intent(inout), pointer :: thisNode
     end function Hot_Halo_Density_Cored_Isothermal_Core_Radius_Template
  end interface

contains

  subroutine Hot_Halo_Density_Cored_Isothermal_Core_Radii_Initialize
    !% Initialize the cored isothermal hot halo profile core radius module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="hotHaloCoredIsothermalCoreRadiiMethod" type="moduleUse">
    include 'hot_halo.density_profile.cored_isothermal.core_radius.modules.inc'
    !# </include>
    implicit none

    !$omp critical(Hot_Halo_Density_Cored_Isothermal_Core_Radius_Initialization) 
    ! Initialize if necessary.
    if (.not.hotHaloCoredIsothermalCoreRadiiInitialized) then
       ! Get the spheroid star formation timescale method parameter.
       !@ <inputParameter>
       !@   <name>hotHaloCoredIsothermalCoreRadiiMethod</name>
       !@   <defaultValue>virialRadiusFraction</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of the method to be used for computing the core radii of cored isothermal hot halo profiles.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('hotHaloCoredIsothermalCoreRadiiMethod',hotHaloCoredIsothermalCoreRadiiMethod,defaultValue='virialRadiusFraction')
       ! Include file that makes calls to all available method initialization routines.
       !# <include directive="hotHaloCoredIsothermalCoreRadiiMethod" type="code" action="subroutine">
       !#  <subroutineArgs>hotHaloCoredIsothermalCoreRadiiMethod,Hot_Halo_Density_Cored_Isothermal_Core_Radius_Get</subroutineArgs>
       include 'hot_halo.density_profile.cored_isothermal.core_radius.inc'
       !# </include>
       if (.not.associated(Hot_Halo_Density_Cored_Isothermal_Core_Radius_Get)) call&
            & Galacticus_Error_Report('Hot_Halo_Density_Cored_Isothermal_Core_Radii' ,'method ' &
            &//char(hotHaloCoredIsothermalCoreRadiiMethod)//' is unrecognized')
       hotHaloCoredIsothermalCoreRadiiInitialized=.true.
    end if
    !$omp end critical(Hot_Halo_Density_Cored_Isothermal_Core_Radius_Initialization) 

    return
  end subroutine Hot_Halo_Density_Cored_Isothermal_Core_Radii_Initialize

  double precision function Hot_Halo_Density_Cored_Isothermal_Core_Radius(thisNode)
    !% Returns the radius (in Mpc) of the core in a cored isothermal hot halo density profile for {\tt thisNode}.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode

    ! Initialize the module.
    call Hot_Halo_Density_Cored_Isothermal_Core_Radii_Initialize

    ! Get the energy using the selected method.
    Hot_Halo_Density_Cored_Isothermal_Core_Radius=Hot_Halo_Density_Cored_Isothermal_Core_Radius_Get(thisNode)

    return
  end function Hot_Halo_Density_Cored_Isothermal_Core_Radius
  
end module Hot_Halo_Density_Cored_Isothermal_Core_Radii
