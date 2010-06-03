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






!% Contains a module that implements calculations of the cooling function.

module Cooling_Functions
  !% Implements calculations of the cooling function.
  use Abundances_Structure
  use Radiation_Structure
  use ISO_Varying_String 
  !# <include directive="coolingFunctionMethod" type="moduleUse">
  include 'cooling.cooling_function.modules.inc'
  !# </include>
  private
  public :: Cooling_Function, Cooling_Function_Density_Log_Slope, Cooling_Function_Temperature_Log_Slope

  ! Flag to indicate if this module has been initialized.  
  logical              :: coolingFunctionInitialized=.false.

  ! Name of cooling time available method used.
  type(varying_string) :: coolingFunctionMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Cooling_Function_Get_Template), pointer :: Cooling_Function_Get => null()
  procedure(Cooling_Function_Get_Template), pointer :: Cooling_Function_Density_Log_Slope_Get => null()
  procedure(Cooling_Function_Get_Template), pointer :: Cooling_Function_Temperature_Log_Slope_Get => null()
  interface Cooling_Function_Get_Template
     double precision function Cooling_Function_Get_Template(temperature,numberDensityHydrogen,abundances,radiation)
       import abundancesStructure, radiationStructure
       double precision,          intent(in) :: temperature,numberDensityHydrogen
       type(abundancesStructure), intent(in) :: abundances
       type(radiationStructure),  intent(in) :: radiation
     end function Cooling_Function_Get_Template
  end interface
  
contains

  subroutine Cooling_Function_Initialize
    !% Initialize the cooling function module
    use Galacticus_Error
    use Input_Parameters
    implicit none
    
    !$omp critical(Cooling_Function_Initialization) 
    ! Initialize if necessary.
    if (.not.coolingFunctionInitialized) then
       ! Get the cooling function method parameter.
       !@ <inputParameter>
       !@   <name>coolingFunctionMethod</name>
       !@   <defaultValue>atomic CIE Cloudy</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of the method to be used for computing the cooling function.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('coolingFunctionMethod',coolingFunctionMethod,defaultValue='atomic CIE Cloudy')
       ! Include file that makes calls to all available method initialization routines.
       !# <include directive="coolingFunctionMethod" type="code" action="subroutine">
       !#  <subroutineArgs>coolingFunctionMethod,Cooling_Function_Get,Cooling_Function_Density_Log_Slope_Get,Cooling_Function_Temperature_Log_Slope_Get</subroutineArgs>
       include 'cooling.cooling_function.inc'
       !# </include>
       if (.not.associated(Cooling_Function_Get)) call Galacticus_Error_Report('Cooling_Function','method ' &
            &//char(coolingFunctionMethod)//' is unrecognized')
       coolingFunctionInitialized=.true.
    end if
    !$omp end critical(Cooling_Function_Initialization) 
    return
  end subroutine Cooling_Function_Initialize

  double precision function Cooling_Function(temperature,numberDensityHydrogen,abundances,radiation)
    !% Return the cooling function at the given temperature and hydrogen density for the specified set of abundances and radiation
    !% field. Units of the returned cooling function are the traditional ergs cm$^3$ s$^{-1}$.
    implicit none
    double precision,          intent(in) :: temperature,numberDensityHydrogen
    type(abundancesStructure), intent(in) :: abundances
    type(radiationStructure),  intent(in) :: radiation

    ! Initialize the module.
    call Cooling_Function_Initialize
    
    ! Get the cooling time using the selected method.
    Cooling_Function=Cooling_Function_Get(temperature,numberDensityHydrogen,abundances,radiation)
    
    return
  end function Cooling_Function
  
  double precision function Cooling_Function_Density_Log_Slope(temperature,numberDensityHydrogen,abundances,radiation)
    !% Return $\d\ln\Lambda/\d\ln\rho$ for a cooling function at the given temperature and hydrogen density for the specified set
    !% of abundances and radiation field.
    implicit none
    double precision,          intent(in) :: temperature,numberDensityHydrogen
    type(abundancesStructure), intent(in) :: abundances
    type(radiationStructure),  intent(in) :: radiation

    ! Initialize the module.
    call Cooling_Function_Initialize
    
    ! Get the cooling time using the selected method.
    Cooling_Function_Density_Log_Slope=Cooling_Function_Density_Log_Slope_Get(temperature,numberDensityHydrogen,abundances,radiation)
    
    return
  end function Cooling_Function_Density_Log_Slope

  double precision function Cooling_Function_Temperature_Log_Slope(temperature,numberDensityHydrogen,abundances,radiation)
    !% Return $\d\ln\Lambda/\d\ln T$ for a cooling function at the given temperature and hydrogen density for the specified set
    !% of abundances and radiation field.
    implicit none
    double precision,          intent(in) :: temperature,numberDensityHydrogen
    type(abundancesStructure), intent(in) :: abundances
    type(radiationStructure),  intent(in) :: radiation

    ! Initialize the module.
    call Cooling_Function_Initialize
    
    ! Get the cooling time using the selected method.
    Cooling_Function_Temperature_Log_Slope=Cooling_Function_Temperature_Log_Slope_Get(temperature,numberDensityHydrogen,abundances,radiation)
    
    return
  end function Cooling_Function_Temperature_Log_Slope

end module Cooling_Functions
