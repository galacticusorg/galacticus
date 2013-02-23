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

!% Contains a module that implements calculations of ram pressure stripping of hot halos.

module Hot_Halo_Ram_Pressure_Forces
  !% Implements calculations of ram pressure force
  use ISO_Varying_String
  implicit none
  private
  public :: Hot_Halo_Ram_Pressure_Force

  ! Flag to indicate if this module has been initialized.  
  logical              :: hotHaloRamPressureForceInitialized=.false.

  ! Name of ram pressure force method used.
  type(varying_string) :: hotHaloRamPressureForceMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Hot_Halo_Ram_Pressure_Force), pointer :: Hot_Halo_Ram_Pressure_Force_Get => null()

contains

  double precision function Hot_Halo_Ram_Pressure_Force(thisNode)
    !% Return the ram pressure force for the hot halo of {\tt thisNode}.
    use Galacticus_Nodes
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="hotHaloRamPressureForceMethod" type="moduleUse">
    include 'hot_halo.ram_pressure_force.modules.inc'
    !# </include>
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    ! Initialize if necessary.
    if (.not.hotHaloRamPressureForceInitialized) then
       !$omp critical(Hot_Halo_Ram_Pressure_Force_Initialization) 
       if (.not.hotHaloRamPressureForceInitialized) then
          ! Get the cooling rate method parameter.
          !@ <inputParameter>
          !@   <name>hotHaloRamPressureForceMethod</name>
          !@   <defaultValue>Font2008</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be used when computing ram pressure force on hot halos.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('hotHaloRamPressureForceMethod',hotHaloRamPressureForceMethod,defaultValue='Font2008')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="hotHaloRamPressureForceMethod" type="functionCall" functionType="void">
          !#  <functionArgs>hotHaloRamPressureForceMethod,Hot_Halo_Ram_Pressure_Force_Get</functionArgs>
          include 'hot_halo.ram_pressure_force.inc'
          !# </include>
          if (.not.associated(Hot_Halo_Ram_Pressure_Force_Get)) call Galacticus_Error_Report('Hot_Halo_Ram_Pressure_Force','method ' &
               &//char(hotHaloRamPressureForceMethod)//' is unrecognized')
          hotHaloRamPressureForceInitialized=.true.
       end if
       !$omp end critical(Hot_Halo_Ram_Pressure_Force_Initialization) 
    end if

    ! Get the cooling rate using the selected method.
    Hot_Halo_Ram_Pressure_Force=Hot_Halo_Ram_Pressure_Force_Get(thisNode)

    return
  end function Hot_Halo_Ram_Pressure_Force

end module Hot_Halo_Ram_Pressure_Forces
