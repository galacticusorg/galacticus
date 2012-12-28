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

!% Contains a module which generates a tabulated atomic collisional ionization equilibrium cooling function using {\sc Cloudy}.

module Cooling_Functions_Atomic_CIE_Cloudy
  !% Generates a tabulated atomic collisional ionization equilibrium cooling function using {\sc Cloudy}.
  use ISO_Varying_String
  implicit none
  private
  public :: Cooling_Function_Atomic_CIE_Cloudy_Initialize, Cooling_Function_Atomic_CIE_Cloudy,&
       & Cooling_Function_Density_Slope_Atomic_CIE_Cloudy, Cooling_Function_Temperature_Slope_Atomic_CIE_Cloudy
  
  ! Flag indicating whether or not this cooling function is selected.
  logical                     :: functionSelected=.false.

  ! Flag to indicate if this module has been initialized.
  logical                     :: coolingFunctionInitialized=.false.

  ! File name for the cooling function and ionization state data.
  character(len=55)           :: coolingFunctionFile='data/cooling/cooling_function_Atomic_CIE_Cloudy.xml'
  character(len=55)           :: chemicalStateFile  ='data/chemicalState/chemical_state_Atomic_CIE_Cloudy.xml'

  ! Maximum tabulated metallicity.
  double precision, parameter :: metallicityMaximumDefault=30.0d0 ! Thirty times Solar.
  double precision            :: metallicityMaximum

  ! Maximum metallicity that we will ever tabulate. More than this should be physically implausible.
  double precision, parameter :: metallicityMaximumLimit  =30.0d0 ! Thirty times Solar.

  ! Factor by which metallicity must exceed currently tabulated maximum before we retabulate.
  double precision, parameter :: metallicityTolerance     =0.1d0

contains
  
  !# <coolingFunctionMethods>
  !#  <unitName>Cooling_Function_Atomic_CIE_Cloudy_Initialize</unitName>
  !# </coolingFunctionMethods>
  subroutine Cooling_Function_Atomic_CIE_Cloudy_Initialize(coolingFunctionMethods,coolingFunctionsMatched)
    !% Initializes the ``atomic CIE cooling function from {\sc Cloudy}'' module.
    implicit none
    type(varying_string), intent(in   ) :: coolingFunctionMethods(:)
    integer,              intent(inout) :: coolingFunctionsMatched
 
    ! Check if this cooling function has been selected.
    if (any(coolingFunctionMethods == 'atomicCIECloudy')) then
       functionSelected=.true.
       coolingFunctionsMatched=coolingFunctionsMatched+1
    end if

    return
  end subroutine Cooling_Function_Atomic_CIE_Cloudy_Initialize

  subroutine Cooling_Function_Atomic_CIE_Cloudy_Create(gasAbundances)
    !% Create the cooling function.
    use Cooling_Functions_CIE_File
    use Abundances_Structure
    use System_Command
    use Galacticus_Input_Paths
    use String_Handling
    implicit none
    type(abundances), intent(in) :: gasAbundances
    logical                               :: makeFile
    character(len=32)                     :: metallicityLabel
    type(varying_string)                  :: command,coolingFunctionFileVarString

    ! Generate the name of the data file and an XML input parameter file.
    !$omp critical (Cooling_Function_Atomic_CIE_Cloudy_Initialize)
    ! Determine if we need to reinitialize this module.
    if (.not.coolingFunctionInitialized) then
       makeFile=.true.
    else
       if (Abundances_Get_Metallicity(gasAbundances,linearByMassSolar) > 0.0d0) then
          makeFile=(min(Abundances_Get_Metallicity(gasAbundances,linearByMassSolar),metallicityMaximumLimit) > metallicityMaximum&
               &*(1.0d0+metallicityTolerance))
       else
          makeFile=.false.
       end if
       if (makeFile) then
          ! Remove the transfer function file so that a new one will be created.
          command='rm -f '//char(Galacticus_Input_Path())//trim(coolingFunctionFile)
          call System_Command_Do(command)
       end if
    end if
    ! Read the file if this module has not been initialized or if the metallicity is out of range.
    if (makeFile) then
       ! Determine maximum metallicity to which we should tabulate.
       if (Abundances_Get_Metallicity(gasAbundances,linearByMassSolar) > 0.0d0) then
          metallicityMaximum=min(max(metallicityMaximumDefault,3.0d0*Abundances_Get_Metallicity(gasAbundances,linearByMassSolar))&
               &,metallicityMaximumLimit)
       else
          metallicityMaximum=metallicityMaximumDefault
       end if
       write (metallicityLabel,'(e12.6)') dlog10(metallicityMaximum)

       ! Run Atomic_CIE_Cloudy wrapper script.
       command=char(Galacticus_Input_Path())//'scripts/aux/Atomic_CIE_Cloudy_Driver.pl '//metallicityLabel//' '//char(Galacticus_Input_Path())//trim(coolingFunctionFile)//' '&
            &//char(Galacticus_Input_Path())//trim(chemicalStateFile)
       command=command//" "//Cooling_Function_CIE_File_Format_Version()
       call System_Command_Do(command)

       ! Call routine to read in the tabulated data.
       coolingFunctionFileVarString=char(Galacticus_Input_Path())//trim(coolingFunctionFile)
       call Cooling_Function_CIE_File_Read(coolingFunctionFileVarString,metallicityMaximumTabulated=metallicityMaximum)

       ! Flag that transfer function is now initialized.
       coolingFunctionInitialized=.true.
    end if
    !$omp end critical (Cooling_Function_Atomic_CIE_Cloudy_Initialize)
    return
  end subroutine Cooling_Function_Atomic_CIE_Cloudy_Create

  !# <coolingFunctionCompute>
  !#   <unitName>Cooling_Function_Atomic_CIE_Cloudy</unitName>
  !# </coolingFunctionCompute>
  subroutine Cooling_Function_Atomic_CIE_Cloudy(coolingFunction,temperature,numberDensityHydrogen,gasAbundances,chemicalDensities,radiation)
    !% Return the cooling function assuming atomic CIE as computed by {\sc Cloudy}.
    use Cooling_Functions_CIE_File
    use Abundances_Structure
    use Chemical_Abundances_Structure
    use Radiation_Structure
    implicit none
    double precision,                  intent(in)  :: temperature,numberDensityHydrogen
    type(abundances),         intent(in)  :: gasAbundances
    type(chemicalAbundances), intent(in)  :: chemicalDensities
    type(radiationStructure),          intent(in)  :: radiation
    double precision,                  intent(out) :: coolingFunction

    ! Check if this cooling function has been selected.
    if (functionSelected) then
       
       ! Create the cooling function.
       call Cooling_Function_Atomic_CIE_Cloudy_Create(gasAbundances)
       
       ! Call routine to interpolate in the tabulated function.
       coolingFunction=Cooling_Function_CIE_File_Interpolate(temperature,numberDensityHydrogen,gasAbundances,radiation)
       
    else

       ! Not selected, return zero.
       coolingFunction=0.0d0

    end if

    return
  end subroutine Cooling_Function_Atomic_CIE_Cloudy
  
  !# <coolingFunctionDensitySlopeCompute>
  !#   <unitName>Cooling_Function_Density_Slope_Atomic_CIE_Cloudy</unitName>
  !# </coolingFunctionDensitySlopeCompute>
  subroutine Cooling_Function_Density_Slope_Atomic_CIE_Cloudy(coolingFunctionDensitySlope,temperature,numberDensityHydrogen,gasAbundances&
       &,chemicalDensities,radiation)
    !% Return the gradient with respect to density of cooling function assuming atomic CIE as computed by {\sc Cloudy}.
    use Abundances_Structure
    use Chemical_Abundances_Structure
    use Radiation_Structure
    implicit none
    double precision,                  intent(in)  :: temperature,numberDensityHydrogen
    type(abundances),         intent(in)  :: gasAbundances
    type(chemicalAbundances), intent(in)  :: chemicalDensities
    type(radiationStructure),          intent(in)  :: radiation
    double precision,                  intent(out) :: coolingFunctionDensitySlope
    double precision                                :: coolingFunction

    ! Check if this cooling function has been selected.
    if (functionSelected) then
       
       ! Get the cooling function.
       call Cooling_Function_Atomic_CIE_Cloudy(coolingFunction,temperature,numberDensityHydrogen,gasAbundances,chemicalDensities,radiation)
       
       ! Logarithmic slope is always 2 for a CIE cooling function.
       coolingFunctionDensitySlope=2.0d0*coolingFunction/numberDensityHydrogen
       
    else
       
       ! Not selected, return zero.
       coolingFunctionDensitySlope=0.0d0
       
    end if
       
    return
  end subroutine Cooling_Function_Density_Slope_Atomic_CIE_Cloudy
  
  !# <coolingFunctionTemperatureSlopeCompute>
  !#   <unitName>Cooling_Function_Temperature_Slope_Atomic_CIE_Cloudy</unitName>
  !# </coolingFunctionTemperatureSlopeCompute>
  subroutine Cooling_Function_Temperature_Slope_Atomic_CIE_Cloudy(coolingFunctionTemperatureSlope,temperature&
       &,numberDensityHydrogen,gasAbundances,chemicalDensities,radiation)
    !% Return the cooling function assuming atomic CIE as computed by {\sc Cloudy}.
    use Cooling_Functions_CIE_File
    use Abundances_Structure
    use Chemical_Abundances_Structure
    use Radiation_Structure
    implicit none
    double precision,                  intent(in)  :: temperature,numberDensityHydrogen
    type(abundances),         intent(in)  :: gasAbundances
    type(chemicalAbundances), intent(in)  :: chemicalDensities
    type(radiationStructure),          intent(in)  :: radiation
    double precision,                  intent(out) :: coolingFunctionTemperatureSlope
    double precision                               :: coolingFunction

    ! Check if this cooling function has been selected.
    if (functionSelected) then
       
       ! Create the cooling function.
       call Cooling_Function_Atomic_CIE_Cloudy_Create(gasAbundances)
       
       ! Get the cooling function.
       call Cooling_Function_Atomic_CIE_Cloudy(coolingFunction,temperature,numberDensityHydrogen,gasAbundances,chemicalDensities,radiation)
       
       ! Call routine to interpolate in the tabulated function.
       coolingFunctionTemperatureSlope=Cooling_Function_CIE_File_logTemperature_Interpolate(temperature,numberDensityHydrogen&
            &,gasAbundances,radiation)*coolingFunction/temperature
       
    else
       
       ! Not selected, return zero.
       coolingFunctionTemperatureSlope=0.0d0
       
    end if
       
    return
  end subroutine Cooling_Function_Temperature_Slope_Atomic_CIE_Cloudy
  
end module Cooling_Functions_Atomic_CIE_Cloudy
