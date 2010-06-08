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
  private
  public :: Cooling_Function, Cooling_Function_Density_Log_Slope, Cooling_Function_Temperature_Log_Slope

  ! Flag to indicate if this module has been initialized.  
  logical                                         :: coolingFunctionInitialized=.false.

  ! Name of cooling time available method used.
  type(varying_string), allocatable, dimension(:) :: coolingFunctionMethods

contains

  subroutine Cooling_Function_Initialize
    !% Initialize the cooling function module.
    use Galacticus_Error
    use Input_Parameters
    use Memory_Management
    !# <include directive="coolingFunctionMethods" type="moduleUse">
    include 'cooling.cooling_function.modules.inc'
    !# </include>
    implicit none
    integer :: coolingFunctionsCount
    
    !$omp critical(Cooling_Function_Initialization) 
    ! Initialize if necessary.
    if (.not.coolingFunctionInitialized) then
       ! Get the cooling function method parameter.
       !@ <inputParameter>
       !@   <name>coolingFunctionMethods</name>
       !@   <defaultValue>atomic CIE Cloudy, CMB Compton</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The names of the methods to be used for computing the cooling function.
       !@   </description>
       !@ </inputParameter>
       coolingFunctionsCount=max(1,Get_Input_Parameter_Array_Size('coolingFunctionMethods'))
       call Alloc_Array(coolingFunctionMethods,coolingFunctionsCount,'coolingFunctionMethods')
       call Get_Input_Parameter('coolingFunctionMethods',coolingFunctionMethods,defaultValue=['atomic_CIE_Cloudy'])

       ! Include file that makes calls to all available method initialization routines.
       !# <include directive="coolingFunctionMethods" type="code" action="subroutine">
       !#  <subroutineArgs>coolingFunctionMethods</subroutineArgs>
       include 'cooling.cooling_function.inc'
       !# </include>
       coolingFunctionInitialized=.true.
    end if
    !$omp end critical(Cooling_Function_Initialization) 
    return
  end subroutine Cooling_Function_Initialize

  double precision function Cooling_Function(temperature,numberDensityHydrogen,abundances,radiation)
    !% Return the cooling function at the given temperature and hydrogen density for the specified set of abundances and radiation
    !% field. Units of the returned cooling function are the traditional ergs cm$^-3$ s$^{-1}$.
    !# <include directive="coolingFunctionCompute" type="moduleUse">
    include 'cooling.cooling_function.compute.modules.inc'
    !# </include>
    implicit none
    double precision,          intent(in) :: temperature,numberDensityHydrogen
    type(abundancesStructure), intent(in) :: abundances
    type(radiationStructure),  intent(in) :: radiation
    double precision                      :: thisCoolingFunction

    ! Initialize the module.
    call Cooling_Function_Initialize
  
    Cooling_Function=0.0d0
    !# <include directive="coolingFunctionCompute" type="code" action="subroutine">
    !#  <subroutineArgs>thisCoolingFunction,temperature,numberDensityHydrogen,abundances,radiation</subroutineArgs>
    !#  <subroutineAction>Cooling_Function=Cooling_Function+thisCoolingFunction</subroutineAction>
    include 'cooling.cooling_function.compute.inc'
    !# </include>
    
    return
  end function Cooling_Function
  
  double precision function Cooling_Function_Density_Log_Slope(temperature,numberDensityHydrogen,abundances,radiation)
    !% Return $\d\ln\Lambda/\d\ln\rho$ for a cooling function at the given temperature and hydrogen density for the specified set
    !% of abundances and radiation field.
    !# <include directive="coolingFunctionDensitySlopeCompute" type="moduleUse">
    include 'cooling.cooling_function.computeDensitySlope.moduleUse.inc'
    !# </include>
    implicit none
    double precision,          intent(in) :: temperature,numberDensityHydrogen
    type(abundancesStructure), intent(in) :: abundances
    type(radiationStructure),  intent(in) :: radiation
    double precision                      :: thisCoolingFunctionDensitySlope,coolingFunction

    ! Initialize the module.
    call Cooling_Function_Initialize
  
    Cooling_Function_Density_Log_Slope=0.0d0
    !# <include directive="coolingFunctionDensitySlopeCompute" type="code" action="subroutine">
    !#  <subroutineArgs>thisCoolingFunctionDensitySlope,temperature,numberDensityHydrogen,abundances,radiation</subroutineArgs>
    !#  <subroutineAction>Cooling_Function_Density_Log_Slope=Cooling_Function_Density_Log_Slope+thisCoolingFunctionDensitySlope</subroutineAction>
    include 'cooling.cooling_function.computeDensitySlope.inc'
    !# </include>
    
    ! Get cooling function.
    coolingFunction=Cooling_Function(temperature,numberDensityHydrogen,abundances,radiation)

    ! Convert to logarithmic slope.
    if (coolingFunction > 0.0d0) then
       Cooling_Function_Density_Log_Slope=Cooling_Function_Density_Log_Slope*numberDensityHydrogen/coolingFunction
    else
       Cooling_Function_Density_Log_Slope=0.0d0
    end if

    return
  end function Cooling_Function_Density_Log_Slope

  double precision function Cooling_Function_Temperature_Log_Slope(temperature,numberDensityHydrogen,abundances,radiation)
    !% Return $\d\ln\Lambda/\d\ln T$ for a cooling function at the given temperature and hydrogen density for the specified set
    !% of abundances and radiation field.
    !# <include directive="coolingFunctionTemperatureSlopeCompute" type="moduleUse">
    include 'cooling.cooling_function.computeTemperatureSlope.moduleUse.inc'
    !# </include>
    implicit none
    double precision,          intent(in) :: temperature,numberDensityHydrogen
    type(abundancesStructure), intent(in) :: abundances
    type(radiationStructure),  intent(in) :: radiation
    double precision                      :: thisCoolingFunctionTemperatureSlope,coolingFunction

    ! Initialize the module.
    call Cooling_Function_Initialize
  
    Cooling_Function_Temperature_Log_Slope=0.0d0
    !# <include directive="coolingFunctionTemperatureSlopeCompute" type="code" action="subroutine">
    !#  <subroutineArgs>thisCoolingFunctionTemperatureSlope,temperature,numberDensityHydrogen,abundances,radiation</subroutineArgs>
    !#  <subroutineAction>Cooling_Function_Temperature_Log_Slope=Cooling_Function_Temperature_Log_Slope+thisCoolingFunctionTemperatureSlope</subroutineAction>
    include 'cooling.cooling_function.computeTemperatureSlope.inc'
    !# </include>
    
    ! Get cooling function.
    coolingFunction=Cooling_Function(temperature,numberDensityHydrogen,abundances,radiation)

    ! Convert to logarithmic slope.
    if (coolingFunction > 0.0d0) then
       Cooling_Function_Temperature_Log_Slope=Cooling_Function_Temperature_Log_Slope*temperature/coolingFunction
    else
       Cooling_Function_Temperature_Log_Slope=0.0d0
    end if

    return
  end function Cooling_Function_Temperature_Log_Slope

end module Cooling_Functions
