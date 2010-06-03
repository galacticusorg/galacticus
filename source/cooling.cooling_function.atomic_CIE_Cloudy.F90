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








!% Contains a module which generates a tabulated atomic collisional ionization equilibrium cooling function using {\sc Cloudy}.

module Cooling_Function_Atomic_CIE_Cloudy
  !% Generates a tabulated atomic collisional ionization equilibrium cooling function using {\sc Cloudy}.
  use ISO_Varying_String
  private
  public :: Cooling_Function_Atomic_CIE_Cloudy_Initialize
  
  ! Flag to indicate if this module has been initialized.
  logical                     :: coolingFunctionInitialized=.false.

  ! File name for the cooling function data.
  character(len=50)           :: coolingFunctionFile='data/cooling_function_Atomic_CIE_Cloudy.xml'

  ! Maximum tabulated metallicity.
  double precision, parameter :: metallicityMaximumDefault=10.0d0 ! Ten times Solar.
  double precision            :: metallicityMaximum

  ! Maximum metallicity that we will ever tabulate. More than this should be physically implausible.
  double precision, parameter :: metallicityMaximumLimit  =30.0d0 ! Thirty times Solar.

contains
  
  !# <coolingFunctionMethod>
  !#  <unitName>Cooling_Function_Atomic_CIE_Cloudy_Initialize</unitName>
  !# </coolingFunctionMethod>
  subroutine Cooling_Function_Atomic_CIE_Cloudy_Initialize(coolingFunctionMethod,Cooling_Function_Get&
       &,Cooling_Function_Density_Log_Slope_Get,Cooling_Function_Temperature_Log_Slope_Get)
    !% Initializes the ``atomic CIE cooling function from {\sc Cloudy}'' module.
    implicit none
    type(varying_string),          intent(in)    :: coolingFunctionMethod
    procedure(),          pointer, intent(inout) :: Cooling_Function_Get,Cooling_Function_Density_Log_Slope_Get&
         &,Cooling_Function_Temperature_Log_Slope_Get
    
    if (coolingFunctionMethod.eq.'atomic CIE Cloudy') then
       Cooling_Function_Get                       => Cooling_Function_Atomic_CIE_Cloudy_Get
       Cooling_Function_Density_Log_Slope_Get     => Cooling_Function_Density_Log_Slope_Atomic_CIE_Cloudy_Get
       Cooling_Function_Temperature_Log_Slope_Get => Cooling_Function_Temperature_Log_Slope_Atomic_CIE_Cloudy_Get
    end if
    return
  end subroutine Cooling_Function_Atomic_CIE_Cloudy_Initialize

  subroutine Cooling_Function_Atomic_CIE_Cloudy_Get_Create(abundances)
    !% Create the cooling function.
    use Cooling_Function_CIE_File
    use Abundances_Structure
    use System_Command
    implicit none
    type(abundancesStructure), intent(in) :: abundances
    logical                               :: makeFile
    character(len=32)                     :: metallicityLabel
    type(varying_string)                  :: command,coolingFunctionFileVarString

    ! Generate the name of the data file and an XML input parameter file.
    !$omp critical (Atomic_CIE_Cloudy_Initialize)
    ! Determine if we need to reinitialize this module.
    if (.not.coolingFunctionInitialized) then
       makeFile=.true.
    else
       if (Abundances_Get_Metallicity(abundances) > 0.0d0) then
          makeFile=(min(Abundances_Get_Metallicity(abundances),metallicityMaximumLimit) > metallicityMaximum)
       else
          makeFile=.false.
       end if
       if (makeFile) then
          ! Remove the transfer function file so that a new one will be created.
          command='rm -f '//trim(coolingFunctionFile)
          call System_Command_Do(command)
       end if
    end if
    ! Read the file if this module has not been initialized or if the wavenumber is out of range.
    if (makeFile) then
       ! Determine maximum metallicity to which we should tabulate.
       if (Abundances_Get_Metallicity(abundances) > 0.0d0) then
          metallicityMaximum=max(max(metallicityMaximumDefault,3.0d0*Abundances_Get_Metallicity(abundances)),metallicityMaximumLimit)
       else
          metallicityMaximum=metallicityMaximumDefault
       end if
       write (metallicityLabel,'(e12.6)') dlog10(metallicityMaximum)

       ! Run Atomic_CIE_Cloudy wrapper script.
       command='./scripts/aux/Atomic_CIE_Cloudy_Driver.pl '//metallicityLabel//' '//trim(coolingFunctionFile)
       call System_Command_Do(command)

       ! Call routine to read in the tabulated data.
       coolingFunctionFileVarString=trim(coolingFunctionFile)
       call Cooling_Function_CIE_File_Read(coolingFunctionFileVarString,metallicityMaximumTabulated=metallicityMaximum)

       ! Flag that transfer function is now initialized.
       coolingFunctionInitialized=.true.
    end if
    !$omp end critical (Atomic_CIE_Cloudy_Initialize)
    return
  end subroutine Cooling_Function_Atomic_CIE_Cloudy_Get_Create

  double precision function Cooling_Function_Atomic_CIE_Cloudy_Get(temperature,numberDensityHydrogen,abundances,radiation)
    !% Return the cooling function assuming atomic CIE as computed by {\sc Cloudy}.
    use Cooling_Function_CIE_File
    use Abundances_Structure
    use Radiation_Structure
    implicit none
    double precision,          intent(in) :: temperature,numberDensityHydrogen
    type(abundancesStructure), intent(in) :: abundances
    type(radiationStructure),  intent(in) :: radiation

    ! Create the cooling function.
    call Cooling_Function_Atomic_CIE_Cloudy_Get_Create(abundances)

    ! Call routine to interpolate in the tabulated function.
    Cooling_Function_Atomic_CIE_Cloudy_Get=Cooling_Function_CIE_File_Interpolate(temperature,numberDensityHydrogen,abundances,radiation)

    return
  end function Cooling_Function_Atomic_CIE_Cloudy_Get
  
  double precision function Cooling_Function_Density_Log_Slope_Atomic_CIE_Cloudy_Get(temperature,numberDensityHydrogen,abundances&
       &,radiation)
    !% Return the logarithmic gradient with respect to density of cooling function assuming atomic CIE as computed by {\sc Cloudy}.
    use Abundances_Structure
    use Radiation_Structure
    implicit none
    double precision,          intent(in) :: temperature,numberDensityHydrogen
    type(abundancesStructure), intent(in) :: abundances
    type(radiationStructure),  intent(in) :: radiation

    ! Always 2 for a CIE cooling function.
    Cooling_Function_Density_Log_Slope_Atomic_CIE_Cloudy_Get=2.0d0
    return
  end function Cooling_Function_Density_Log_Slope_Atomic_CIE_Cloudy_Get
  
  double precision function Cooling_Function_Temperature_Log_Slope_Atomic_CIE_Cloudy_Get(temperature,numberDensityHydrogen,abundances,radiation)
    !% Return the cooling function assuming atomic CIE as computed by {\sc Cloudy}.
    use Cooling_Function_CIE_File
    use System_Command
    use Abundances_Structure
    use Radiation_Structure
    implicit none
    double precision,          intent(in) :: temperature,numberDensityHydrogen
    type(abundancesStructure), intent(in) :: abundances
    type(radiationStructure),  intent(in) :: radiation

    ! Create the cooling function.
    call Cooling_Function_Atomic_CIE_Cloudy_Get_Create(abundances)

    ! Call routine to interpolate in the tabulated function.
    Cooling_Function_Temperature_Log_Slope_Atomic_CIE_Cloudy_Get=Cooling_Function_CIE_File_logTemperature_Interpolate(temperature&
         &,numberDensityHydrogen,abundances,radiation)

    return
  end function Cooling_Function_Temperature_Log_Slope_Atomic_CIE_Cloudy_Get
  
end module Cooling_Function_Atomic_CIE_Cloudy
