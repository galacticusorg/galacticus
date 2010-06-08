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






!% Contains a module that implements calculations of the ionization state.

module Ionization_States
  !% Implements calculations of the ionization state.
  use Abundances_Structure
  use Radiation_Structure
  use ISO_Varying_String 
  private
  public :: Electron_Density, Electron_Density_Temperature_Log_Slope, Electron_Density_Density_Log_Slope

  ! Flag to indicate if this module has been initialized.  
  logical              :: ionizationStateInitialized=.false.

  ! Name of ionization state method used.
  type(varying_string) :: ionizationStateMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Electron_Density_Get_Template), pointer :: Electron_Density_Get                       => null()
  procedure(Electron_Density_Get_Template), pointer :: Electron_Density_Temperature_Log_Slope_Get => null()
  procedure(Electron_Density_Get_Template), pointer :: Electron_Density_Density_Log_Slope_Get     => null()
  abstract interface
     double precision function Electron_Density_Get_Template(temperature,numberDensityHydrogen,abundances,radiation)
       import abundancesStructure,radiationStructure
       double precision,          intent(in) :: temperature,numberDensityHydrogen
       type(abundancesStructure), intent(in) :: abundances
       type(radiationStructure),  intent(in) :: radiation
     end function Electron_Density_Get_Template
  end interface
  
contains

  subroutine Ionization_State_Initialize
    !% Initialize the ionization state module.
    use Galacticus_Error
    use Input_Parameters
    use Memory_Management
    !# <include directive="ionizationStateMethod" type="moduleUse">
    include 'atomic.ionization_state.modules.inc'
    !# </include>
    implicit none
    
    !$omp critical(Ionization_State_Initialization) 
    ! Initialize if necessary.
    if (.not.ionizationStateInitialized) then
       ! Get the ionization state method parameter.
       !@ <inputParameter>
       !@   <name>ionizationStateMethod</name>
       !@   <defaultValue>atomic CIE Cloudy</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of the method to be used for computing the ionization state.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('ionizationStateMethod',ionizationStateMethod,defaultValue='atomic_CIE_Cloudy')

       ! Include file that makes calls to all available method initialization routines.
       !# <include directive="ionizationStateMethod" type="code" action="subroutine">
       !#  <subroutineArgs>ionizationStateMethod,Electron_Density_Get,Electron_Density_Temperature_Log_Slope_Get,Electron_Density_Density_Log_Slope_Get</subroutineArgs>
       include 'atomic.ionization_state.inc'
       !# </include>
       if (.not.(associated(Electron_Density_Get).and.associated(Electron_Density_Temperature_Log_Slope_Get) &
            & .and.associated(Electron_Density_Density_Log_Slope_Get))) call&
            & Galacticus_Error_Report('Ionization_State_Initialize','method '//char(ionizationStateMethod)//' is unrecognized')
       ionizationStateInitialized=.true.
    end if
    !$omp end critical(Ionization_State_Initialization) 
    return
  end subroutine Ionization_State_Initialize

  double precision function Electron_Density(temperature,numberDensityHydrogen,abundances,radiation)
    !% Return the electron density at the given temperature and hydrogen density for the specified set of abundances and radiation
    !% field. Units of the returned electron density are cm$^-3$.
    !# <include directive="ionizationStateCompute" type="moduleUse">
    include 'atomic.ionization_state.compute.modules.inc'
    !# </include>
    implicit none
    double precision,          intent(in) :: temperature,numberDensityHydrogen
    type(abundancesStructure), intent(in) :: abundances
    type(radiationStructure),  intent(in) :: radiation

    ! Initialize the module.
    call Ionization_State_Initialize

    ! Call the routine to do the calculation.
    Electron_Density=Electron_Density_Get(temperature,numberDensityHydrogen,abundances,radiation)
    
    return
  end function Electron_Density

  double precision function Electron_Density_Temperature_Log_Slope(temperature,numberDensityHydrogen,abundances,radiation)
    !% Return the logarithmic gradient of electron density with temperature at the given temperature and hydrogen density for the
    !% specified set of abundances and radiation field. Units of the returned electron density are cm$^-3$.
    !# <include directive="ionizationStateCompute" type="moduleUse">
    include 'atomic.ionization_state.compute.modules.inc'
    !# </include>
    implicit none
    double precision,          intent(in) :: temperature,numberDensityHydrogen
    type(abundancesStructure), intent(in) :: abundances
    type(radiationStructure),  intent(in) :: radiation

    ! Initialize the module.
    call Ionization_State_Initialize

    ! Call the routine to do the calculation.
    Electron_Density_Temperature_Log_Slope=Electron_Density_Temperature_Log_Slope_Get(temperature,numberDensityHydrogen&
         &,abundances,radiation)
    
    return
  end function Electron_Density_Temperature_Log_Slope

  double precision function Electron_Density_Density_Log_Slope(temperature,numberDensityHydrogen,abundances,radiation)
    !% Return the logarithmic gradient of electron density with respect to density at the given temperature and hydrogen density
    !% for the specified set of abundances and radiation field. Units of the returned electron density are cm$^-3$.
    !# <include directive="ionizationStateCompute" type="moduleUse">
    include 'atomic.ionization_state.compute.modules.inc'
    !# </include>
    implicit none
    double precision,          intent(in) :: temperature,numberDensityHydrogen
    type(abundancesStructure), intent(in) :: abundances
    type(radiationStructure),  intent(in) :: radiation

    ! Initialize the module.
    call Ionization_State_Initialize

    ! Call the routine to do the calculation.
    Electron_Density_Density_Log_Slope=Electron_Density_Density_Log_Slope_Get(temperature,numberDensityHydrogen&
         &,abundances,radiation)
    
    return
  end function Electron_Density_Density_Log_Slope

end module Ionization_States
