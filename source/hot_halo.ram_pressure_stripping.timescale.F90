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

!% Contains a module that implements calculations of ram pressure stripping timescales for hot halos.

module Hot_Halo_Ram_Pressure_Stripping_Timescales
  !% Implements calculations of ram pressure stripping timescales for hot halos.
  use ISO_Varying_String
  implicit none
  private
  public :: Hot_Halo_Ram_Pressure_Stripping_Timescale

  ! Flag to indicate if this module has been initialized.
  logical                                                    :: hotHaloRamPressureStrippingTimescaleInitialized    =.false.

  ! Name of cooling rate available method used.
  type     (varying_string                        )          :: hotHaloRamPressureStrippingTimescaleMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Hot_Halo_Ram_Pressure_Stripping_Timescale), pointer :: Hot_Halo_Ram_Pressure_Stripping_Timescale_Get=>null()

contains

  double precision function Hot_Halo_Ram_Pressure_Stripping_Timescale(thisNode)
    !% Return the ram pressure stripping radius for the hot halo of {\tt thisNode} (in units of Mpc).
    use Galacticus_Nodes
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="hotHaloRamPressureStrippingTimescaleMethod" type="moduleUse">
    include 'hot_halo.ram_pressure_stripping.timescale.modules.inc'
    !# </include>
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    ! Initialize if necessary.
    if (.not.hotHaloRamPressureStrippingTimescaleInitialized) then
       !$omp critical(Hot_Halo_Ram_Pressure_Stripping_Initialization)
       if (.not.hotHaloRamPressureStrippingTimescaleInitialized) then
          ! Get the cooling rate method parameter.
          !@ <inputParameter>
          !@   <name>hotHaloRamPressureStrippingTimescaleMethod</name>
          !@   <defaultValue>ramPressureAcceleration</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be used when computing ram pressure stripping timescales for hot halos.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('hotHaloRamPressureStrippingTimescaleMethod',hotHaloRamPressureStrippingTimescaleMethod,defaultValue='ramPressureAcceleration')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="hotHaloRamPressureStrippingTimescaleMethod" type="functionCall" functionType="void">
          !#  <functionArgs>hotHaloRamPressureStrippingTimescaleMethod,Hot_Halo_Ram_Pressure_Stripping_Timescale_Get</functionArgs>
          include 'hot_halo.ram_pressure_stripping.timescales.inc'
          !# </include>
          if (.not.associated(Hot_Halo_Ram_Pressure_Stripping_Timescale_Get)) call Galacticus_Error_Report('Hot_Halo_Ram_Pressure_Stripping_Timescale','method ' &
               &//char(hotHaloRamPressureStrippingTimescaleMethod)//' is unrecognized')
          hotHaloRamPressureStrippingTimescaleInitialized=.true.
       end if
       !$omp end critical(Hot_Halo_Ram_Pressure_Stripping_Initialization)
    end if

    ! Get the cooling rate using the selected method.
    Hot_Halo_Ram_Pressure_Stripping_Timescale=Hot_Halo_Ram_Pressure_Stripping_Timescale_Get(thisNode)

    return
  end function Hot_Halo_Ram_Pressure_Stripping_Timescale
  
end module Hot_Halo_Ram_Pressure_Stripping_Timescales
