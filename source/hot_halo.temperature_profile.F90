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

!% Contains a module that implements calculations of the hot halo gas temperature profile.

module Hot_Halo_Temperature_Profile
  use ISO_Varying_String
  use Galacticus_Nodes
  !# <include directive="hotHaloTemperatureMethod" type="moduleUse">
  include 'hot_halo.temperature_profile.modules.inc'
  !# </include>
  implicit none
  private
  public :: Hot_Halo_Temperature, Hot_Halo_Temperature_Logarithmic_Slope

  ! Flag to indicate if this module has been initialized.  
  logical                                               :: hotHaloTemperatureInitialized             =.false.  
  
  ! Name of cooling time available method used.                                                                                                          
  type     (varying_string                   )          :: hotHaloTemperatureMethod                            
  
  ! Pointer to the function that actually does the calculation.                                                                                                          
  procedure(Hot_Halo_Temperature_Get_Template), pointer :: Hot_Halo_Temperature_Get                  =>null()  
  procedure(Hot_Halo_Temperature_Get_Template), pointer :: Hot_Halo_Temperature_Logarithmic_Slope_Get=>null()  
  interface Hot_Halo_Temperature_Get_Template
     double precision function Hot_Halo_Temperature_Get_Template(thisNode,radius)
       import treeNode
       type            (treeNode), intent(inout), pointer :: thisNode  
       double precision          , intent(in   )          :: radius    
     end function Hot_Halo_Temperature_Get_Template
  end interface
  
contains

  subroutine Hot_Halo_Temperature_Initialize
    !% Initialize the hot halo temperature module.
    use Galacticus_Error
    use Input_Parameters
    implicit none

    ! Initialize if necessary.
    if (.not.hotHaloTemperatureInitialized) then
       !$omp critical(Hot_Halo_Temperature_Initialization) 
       if (.not.hotHaloTemperatureInitialized) then
          ! Get the cooling time available method parameter.
          !@ <inputParameter>
          !@   <name>hotHaloTemperatureMethod</name>
          !@   <defaultValue>virial</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be used for computing hot halo temperature profiles.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('hotHaloTemperatureMethod',hotHaloTemperatureMethod,defaultValue='virial')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="hotHaloTemperatureMethod" type="functionCall" functionType="void">
          !#  <functionArgs>hotHaloTemperatureMethod,Hot_Halo_Temperature_Get,Hot_Halo_Temperature_Logarithmic_Slope_Get</functionArgs>
          include 'hot_halo.temperature_profile.inc'
          !# </include>
          if (.not.associated(Hot_Halo_Temperature_Get)) call Galacticus_Error_Report('Hot_Halo_Temperature','method '&
               &//char(hotHaloTemperatureMethod)//' is unrecognized')
          hotHaloTemperatureInitialized=.true.
       end if
       !$omp end critical(Hot_Halo_Temperature_Initialization) 
    end if
    return
  end subroutine Hot_Halo_Temperature_Initialize

  double precision function Hot_Halo_Temperature(thisNode,radius)
    !% Return the temperature of the hot halo in {\tt thisNode} at radius {\tt radius}.
    implicit none
    type            (treeNode), intent(inout), pointer :: thisNode  
    double precision          , intent(in   )          :: radius    
    
    ! Initialize the module.                                                             
    call Hot_Halo_Temperature_Initialize

    ! Get the cooling time using the selected method.
    Hot_Halo_Temperature=Hot_Halo_Temperature_Get(thisNode,radius)

    return
  end function Hot_Halo_Temperature

  double precision function Hot_Halo_Temperature_Logarithmic_Slope(thisNode,radius)
    !% Return the temperature of the hot halo in {\tt thisNode} at radius {\tt radius}.
    implicit none
    type            (treeNode), intent(inout), pointer :: thisNode  
    double precision          , intent(in   )          :: radius    
    
    ! Initialize the module.                                                             
    call Hot_Halo_Temperature_Initialize

    ! Get the cooling time using the selected method.
    Hot_Halo_Temperature_Logarithmic_Slope=Hot_Halo_Temperature_Logarithmic_Slope_Get(thisNode,radius)

    return
  end function Hot_Halo_Temperature_Logarithmic_Slope

end module Hot_Halo_Temperature_Profile
