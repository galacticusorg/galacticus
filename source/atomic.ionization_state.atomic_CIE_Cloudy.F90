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








!% Contains a module which generates a tabulated atomic collisional ionization equilibrium ionization state using {\sc Cloudy}.

module Ionization_States_Atomic_CIE_Cloudy
  !% Generates a tabulated atomic collisional ionization equilibrium ionization state using {\sc Cloudy}.
  use ISO_Varying_String
  private
  public :: Ionization_State_Atomic_CIE_Cloudy_Initialize
  
  ! Flag to indicate if this module has been initialized.
  logical                     :: ionizationStateInitialized=.false.

  ! File name for the ionization state and ionization state data.
  character(len=50)           :: coolingFunctionFile='data/cooling_function_Atomic_CIE_Cloudy.xml'
  character(len=50)           :: ionizationStateFile='data/ionization_state_Atomic_CIE_Cloudy.xml'

  ! Maximum tabulated metallicity.
  double precision, parameter :: metallicityMaximumDefault=10.0d0 ! Ten times Solar.
  double precision            :: metallicityMaximum

  ! Maximum metallicity that we will ever tabulate. More than this should be physically implausible.
  double precision, parameter :: metallicityMaximumLimit  =30.0d0 ! Thirty times Solar.

contains
  
  !# <ionizationStateMethod>
  !#  <unitName>Ionization_State_Atomic_CIE_Cloudy_Initialize</unitName>
  !# </ionizationStateMethod>
  subroutine Ionization_State_Atomic_CIE_Cloudy_Initialize(ionizationStateMethod,Electron_Density_Get&
       &,Electron_Density_Temperature_Log_Slope_Get,Electron_Density_Density_Log_Slope_Get)
    !% Initializes the ``atomic CIE ionization state from {\sc Cloudy}'' module.
    implicit none
    type(varying_string),          intent(in)    :: ionizationStateMethod
    procedure(),          pointer, intent(inout) :: Electron_Density_Get,Electron_Density_Temperature_Log_Slope_Get&
         &,Electron_Density_Density_Log_Slope_Get
 
    ! Check if this ionization state has been selected.
    if (ionizationStateMethod == 'atomic_CIE_Cloudy') then
       Electron_Density_Get                       => Electron_Density_Atomic_CIE_Cloudy
       Electron_Density_Temperature_Log_Slope_Get => Electron_Density_Temperature_Log_Slope_Atomic_CIE_Cloudy
       Electron_Density_Density_Log_Slope_Get     => Electron_Density_Density_Log_Slope_Atomic_CIE_Cloudy
    end if

    return
  end subroutine Ionization_State_Atomic_CIE_Cloudy_Initialize

  subroutine Ionization_State_Atomic_CIE_Cloudy_Create(abundances)
    !% Create the ionization state.
    use Ionization_States_CIE_File
    use Abundances_Structure
    use System_Command
    implicit none
    type(abundancesStructure), intent(in) :: abundances
    logical                               :: makeFile
    character(len=32)                     :: metallicityLabel
    type(varying_string)                  :: command,ionizationStateFileVarString

    ! Generate the name of the data file and an XML input parameter file.
    !$omp critical (Ionization_State_Atomic_CIE_Cloudy_Initialize)
    ! Determine if we need to reinitialize this module.
    if (.not.ionizationStateInitialized) then
       makeFile=.true.
    else
       if (Abundances_Get_Metallicity(abundances) > 0.0d0) then
          makeFile=(min(Abundances_Get_Metallicity(abundances),metallicityMaximumLimit) > metallicityMaximum)
       else
          makeFile=.false.
       end if
       if (makeFile) then
          ! Remove the transfer function file so that a new one will be created.
          command='rm -f '//trim(ionizationStateFile)
          call System_Command_Do(command)
       end if
    end if
    ! Read the file if this module has not been initialized or if the metallicity is out of range.
    if (makeFile) then
       ! Determine maximum metallicity to which we should tabulate.
       if (Abundances_Get_Metallicity(abundances) > 0.0d0) then
          metallicityMaximum=max(max(metallicityMaximumDefault,3.0d0*Abundances_Get_Metallicity(abundances)),metallicityMaximumLimit)
       else
          metallicityMaximum=metallicityMaximumDefault
       end if
       write (metallicityLabel,'(e12.6)') dlog10(metallicityMaximum)

       ! Run Atomic_CIE_Cloudy wrapper script.
       command='./scripts/aux/Atomic_CIE_Cloudy_Driver.pl '//metallicityLabel//' '//trim(coolingFunctionFile)//' '&
            &//trim(ionizationStateFile)
       call System_Command_Do(command)

       ! Call routine to read in the tabulated data.
       ionizationStateFileVarString=trim(ionizationStateFile)
       call Ionization_State_CIE_File_Read(ionizationStateFileVarString,metallicityMaximumTabulated=metallicityMaximum)

       ! Flag that transfer function is now initialized.
       ionizationStateInitialized=.true.
    end if
    !$omp end critical (Ionization_State_Atomic_CIE_Cloudy_Initialize)
    return
  end subroutine Ionization_State_Atomic_CIE_Cloudy_Create

  double precision function Electron_Density_Atomic_CIE_Cloudy(temperature,numberDensityHydrogen,abundances,radiation)
    !% Return the electron density assuming atomic CIE as computed by {\sc Cloudy}.
    use Ionization_States_CIE_File
    use Abundances_Structure
    use Radiation_Structure
    implicit none
    double precision,          intent(in) :: temperature,numberDensityHydrogen
    type(abundancesStructure), intent(in) :: abundances
    type(radiationStructure),  intent(in) :: radiation

    ! Create the ionization state.
    call Ionization_State_Atomic_CIE_Cloudy_Create(abundances)
    
    ! Call routine to interpolate in the tabulated function.
    Electron_Density_Atomic_CIE_Cloudy=Electron_Density_CIE_File_Interpolate(temperature,numberDensityHydrogen,abundances,radiation)
    
    return
  end function Electron_Density_Atomic_CIE_Cloudy

  double precision function Electron_Density_Temperature_Log_Slope_Atomic_CIE_Cloudy(temperature,numberDensityHydrogen,abundances&
       &,radiation)
    !% Return the logarithmic slope of the electron density with respect to temperature assuming atomic CIE as computed by {\sc Cloudy}.
    use Ionization_States_CIE_File
    use System_Command
    use Abundances_Structure
    use Radiation_Structure
    implicit none
    double precision,          intent(in)  :: temperature,numberDensityHydrogen
    type(abundancesStructure), intent(in)  :: abundances
    type(radiationStructure),  intent(in)  :: radiation
    
    ! Create the ionization state.
    call Ionization_State_Atomic_CIE_Cloudy_Create(abundances)
    
    ! Call routine to interpolate in the tabulated function.
    Electron_Density_Temperature_Log_Slope_Atomic_CIE_Cloudy=Electron_Density_CIE_File_logTemperature_Interpolate(temperature&
         &,numberDensityHydrogen ,abundances,radiation)
 
    return
  end function Electron_Density_Temperature_Log_Slope_Atomic_CIE_Cloudy
      
  double precision function Electron_Density_Density_Log_Slope_Atomic_CIE_Cloudy(temperature,numberDensityHydrogen,abundances&
       &,radiation)
    !% Return the logarithmic slope of the electron density with respect to density assuming atomic CIE as computed by {\sc Cloudy}.
    use Abundances_Structure
    use Radiation_Structure
    implicit none
    double precision,          intent(in)  :: temperature,numberDensityHydrogen
    type(abundancesStructure), intent(in)  :: abundances
    type(radiationStructure),  intent(in)  :: radiation
    
    ! Electron density always scales as total density under CIE conditions.
    Electron_Density_Density_Log_Slope_Atomic_CIE_Cloudy=1.0d0
 
    return
  end function Electron_Density_Density_Log_Slope_Atomic_CIE_Cloudy
      
end module Ionization_States_Atomic_CIE_Cloudy
