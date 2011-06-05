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


!% Contains a module that implements calculations of the hot halo gas density profile.

module Hot_Halo_Density_Profile
  use ISO_Varying_String
  use Tree_Nodes
  !# <include directive="hotHaloDensityMethod" type="moduleUse">
  include 'hot_halo.density_profile.modules.inc'
  !# </include>
  private
  public :: Hot_Halo_Density, Hot_Halo_Density_Log_Slope

  ! Flag to indicate if this module has been initialized.  
  logical              :: hotHaloDensityInitialized=.false.

  ! Name of cooling time available method used.
  type(varying_string) :: hotHaloDensityMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Hot_Halo_Density_Get_Template), pointer :: Hot_Halo_Density_Get => null()
  procedure(Hot_Halo_Density_Get_Template), pointer :: Hot_Halo_Density_Log_Slope_Get => null()
  abstract interface
     double precision function Hot_Halo_Density_Get_Template(thisNode,radius)
       import treeNode
       type(treeNode),   intent(inout), pointer :: thisNode
       double precision, intent(in)             :: radius
     end function Hot_Halo_Density_Get_Template
  end interface
  
contains

  subroutine Hot_Halo_Density_Initialize
    !% Initialize the hot halo density profile module.
    use Input_Parameters
    use Galacticus_Error
    implicit none

    !$omp critical(Hot_Halo_Density_Initialization) 
    ! Initialize if necessary.
    if (.not.hotHaloDensityInitialized) then
       ! Get the cooling time available method parameter.
       !@ <inputParameter>
       !@   <name>hotHaloDensityMethod</name>
       !@   <defaultValue>cored isothermal</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of the method to be used for calculations of the hot halo density profile.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('hotHaloDensityMethod',hotHaloDensityMethod,defaultValue='cored isothermal')
       ! Include file that makes calls to all available method initialization routines.
       !# <include directive="hotHaloDensityMethod" type="code" action="subroutine">
       !#  <subroutineArgs>hotHaloDensityMethod,Hot_Halo_Density_Get,Hot_Halo_Density_Log_Slope_Get</subroutineArgs>
       include 'hot_halo.density_profile.inc'
       !# </include>
       if (.not.(associated(Hot_Halo_Density_Get).and.associated(Hot_Halo_Density_Log_Slope_Get))) call&
            & Galacticus_Error_Report('Hot_Halo_Density_Initialize','method ' //char(hotHaloDensityMethod)//' is unrecognized')
       hotHaloDensityInitialized=.true.
    end if
    !$omp end critical(Hot_Halo_Density_Initialization) 
    return
  end subroutine Hot_Halo_Density_Initialize

  double precision function Hot_Halo_Density(thisNode,radius)
    !% Return the density of the hot halo in {\tt thisNode} at radius {\tt radius}.
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, intent(in)             :: radius

    ! Initialize the module if necessary.
    call Hot_Halo_Density_Initialize

    ! Get the cooling time using the selected method.
    Hot_Halo_Density=Hot_Halo_Density_Get(thisNode,radius)

    return
  end function Hot_Halo_Density

  double precision function Hot_Halo_Density_Log_Slope(thisNode,radius)
    !% Return the density of the hot halo in {\tt thisNode} at radius {\tt radius}.
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, intent(in)             :: radius

    ! Initialize the module if necessary.
    call Hot_Halo_Density_Initialize

    ! Get the cooling time using the selected method.
    Hot_Halo_Density_Log_Slope=Hot_Halo_Density_Log_Slope_Get(thisNode,radius)

    return
  end function Hot_Halo_Density_Log_Slope

end module Hot_Halo_Density_Profile
