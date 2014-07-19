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

!% Contains a module which implements calculations of the core radius in cored isothermal cold mode hot halo profiles.

module Hot_Halo_Cold_Mode_Density_Cored_Isothermal_Core_Radii
  !% Implements calculations of the core radius in cored isothermal cold mode hot halo profiles.
  use ISO_Varying_String
  use Galacticus_Nodes
  implicit none
  private
  public :: Hot_Halo_Cold_Mode_Density_Cored_Isothermal_Core_Radius

  ! Flag to indicate if this module has been initialized.
  logical                                                   :: hotHaloColdModeCoredIsothermalCoreRadiiInitialized       =.false.

  ! Name of cooling rate available method used.
  type     (varying_string                       )          :: hotHaloColdModeCoredIsothermalCoreRadiiMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Hot_Halo_Cold_Mode_Density_Cored_Isothermal_Core_Radius), pointer :: Hot_Halo_Cold_Mode_Density_CIso_CoreR_Get => null()

contains

  subroutine Hot_Halo_Cold_Mode_Density_CIso_CoreR_Init
    !% Initialize the cored isothermal cold mode hot halo profile core radius module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="hotHaloColdModeCoredIsothermalCoreRadiiMethod" type="moduleUse">
    include 'hot_halo.cold_mode.density_profile.cored_isothermal.core_radius.modules.inc'
    !# </include>
    implicit none

    ! Initialize if necessary.
    if (.not.hotHaloColdModeCoredIsothermalCoreRadiiInitialized) then
       !$omp critical(Hot_Halo_Cold_Mode_Density_CIso_CoreR_Initialization)
       if (.not.hotHaloColdModeCoredIsothermalCoreRadiiInitialized) then
          ! Get the cold mode core radius method parameter.
          !@ <inputParameter>
          !@   <name>hotHaloColdModeCoredIsothermalCoreRadiiMethod</name>
          !@   <defaultValue>virialRadiusFraction</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be used for computing the core radii of cored isothermal cold mode hot halo profiles.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('hotHaloColdModeCoredIsothermalCoreRadiiMethod',hotHaloColdModeCoredIsothermalCoreRadiiMethod,defaultValue='virialRadiusFraction')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="hotHaloColdModeCoredIsothermalCoreRadiiMethod" type="functionCall" functionType="void">
          !#  <functionArgs>hotHaloColdModeCoredIsothermalCoreRadiiMethod,Hot_Halo_Cold_Mode_Density_CIso_CoreR_Get</functionArgs>
          include 'hot_halo.cold_mode.density_profile.cored_isothermal.core_radius.inc'
          !# </include>
          if (.not.associated(Hot_Halo_Cold_Mode_Density_CIso_CoreR_Get)) call&
               & Galacticus_Error_Report('Hot_Halo_Cold_Mode_Density_Cored_Isothermal_Core_Radii' ,'method ' &
               &//char(hotHaloColdModeCoredIsothermalCoreRadiiMethod)//' is unrecognized')
          hotHaloColdModeCoredIsothermalCoreRadiiInitialized=.true.
       end if
       !$omp end critical(Hot_Halo_Cold_Mode_Density_CIso_CoreR_Initialization)
    end if
    return
  end subroutine Hot_Halo_Cold_Mode_Density_CIso_CoreR_Init

  double precision function Hot_Halo_Cold_Mode_Density_Cored_Isothermal_Core_Radius(thisNode)
    !% Returns the radius (in Mpc) of the core in a cored isothermal hot halo density profile for {\tt thisNode}.
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    ! Initialize the module.
    call Hot_Halo_Cold_Mode_Density_CIso_CoreR_Init()
    ! Get the energy using the selected method.
    Hot_Halo_Cold_Mode_Density_Cored_Isothermal_Core_Radius=Hot_Halo_Cold_Mode_Density_CIso_CoreR_Get(thisNode)
    return
  end function Hot_Halo_Cold_Mode_Density_Cored_Isothermal_Core_Radius

end module Hot_Halo_Cold_Mode_Density_Cored_Isothermal_Core_Radii
