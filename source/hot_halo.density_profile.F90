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
