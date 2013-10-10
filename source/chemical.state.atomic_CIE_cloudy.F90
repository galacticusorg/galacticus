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

!% Contains a module which generates a tabulated atomic collisional ionization equilibrium ionization state using {\sc Cloudy}.

module Chemical_States_Atomic_CIE_Cloudy
  !% Generates a tabulated atomic collisional ionization equilibrium ionization state using {\sc Cloudy}.
  use ISO_Varying_String
  implicit none
  private
  public :: Chemical_State_Atomic_CIE_Cloudy_Initialize

  ! Flag to indicate if this module has been initialized.
  logical :: chemicalStateInitialized=.false.

  ! File names for the cooling function and chemical state data.
  character(len=55)           :: coolingFunctionFile='data/cooling/cooling_function_Atomic_CIE_Cloudy.xml'
  character(len=55)           :: chemicalStateFile='data/chemicalState/chemical_state_Atomic_CIE_Cloudy.xml'

  ! Maximum tabulated metallicity.
  double precision, parameter :: metallicityMaximumDefault=30.0d0 !   Thirty times Solar.
  double precision            :: metallicityMaximum

  ! Maximum metallicity that we will ever tabulate. More than this should be physically implausible.
  double precision, parameter :: metallicityMaximumLimit  =30.0d0 !   Thirty times Solar.

  ! Factor by which metallicity must exceed currently tabulated maximum before we retabulate.
  double precision, parameter :: metallicityTolerance     =0.1d0

contains

  !# <chemicalStateMethod>
  !#  <unitName>Chemical_State_Atomic_CIE_Cloudy_Initialize</unitName>
  !# </chemicalStateMethod>
  subroutine Chemical_State_Atomic_CIE_Cloudy_Initialize(chemicalStateMethod,Electron_Density_Get &
       &,Electron_Density_Temperature_Log_Slope_Get,Electron_Density_Density_Log_Slope_Get,Chemical_Densities_Get)
    !% Initializes the ``atomic CIE ionization state from {\sc Cloudy}'' module.
    implicit none
    type     (varying_string  ), intent(in   )          :: chemicalStateMethod
    procedure(Electron_Density_Atomic_CIE_Cloudy), intent(inout), pointer :: Electron_Density_Density_Log_Slope_Get
    procedure(Electron_Density_Temperature_Log_Slope_Atomic_CIE_Cloudy), intent(inout), pointer :: Electron_Density_Get
    procedure(Electron_Density_Density_Log_Slope_Atomic_CIE_Cloudy), intent(inout), pointer :: Electron_Density_Temperature_Log_Slope_Get
    procedure(Chemical_Densities_Atomic_CIE_Cloudy), intent(inout), pointer :: Chemical_Densities_Get

    ! Check if this chemical state has been selected.
    if (chemicalStateMethod == 'atomicCIECloudy') then
       Electron_Density_Get                       => Electron_Density_Atomic_CIE_Cloudy
       Electron_Density_Temperature_Log_Slope_Get => Electron_Density_Temperature_Log_Slope_Atomic_CIE_Cloudy
       Electron_Density_Density_Log_Slope_Get     => Electron_Density_Density_Log_Slope_Atomic_CIE_Cloudy
       Chemical_Densities_Get                    => Chemical_Densities_Atomic_CIE_Cloudy
   end if

    return
  end subroutine Chemical_State_Atomic_CIE_Cloudy_Initialize

  subroutine Chemical_State_Atomic_CIE_Cloudy_Create(gasAbundances)
    !% Create the chemical state.
    use Chemical_States_CIE_File
    use Abundances_Structure
    use System_Command
    use Galacticus_Input_Paths
    use String_Handling
    implicit none
    type     (abundances    ), intent(in   ) :: gasAbundances
    logical                                  :: makeFile
    character(len=32        )                :: metallicityLabel
    type     (varying_string)                :: chemicalStateFileVarString, command

    ! Generate the name of the data file and an XML input parameter file.
    !$omp critical (Chemical_State_Atomic_CIE_Cloudy_Initialize)
    ! Determine if we need to reinitialize this module.
    if (.not.chemicalStateInitialized) then
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
          command='rm -f '//char(Galacticus_Input_Path())//trim(chemicalStateFile)
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
       write (metallicityLabel,'(e12.6)') log10(metallicityMaximum)

       ! Run Atomic_CIE_Cloudy wrapper script.
       command=char(Galacticus_Input_Path())//'scripts/aux/Atomic_CIE_Cloudy_Driver.pl '//metallicityLabel//' '//char(Galacticus_Input_Path())//trim(coolingFunctionFile)//' '&
            &//char(Galacticus_Input_Path())//trim(chemicalStateFile)
       command=command//" "//Chemical_State_CIE_File_Format_Version()
       call System_Command_Do(command)

       ! Call routine to read in the tabulated data.
       chemicalStateFileVarString=char(Galacticus_Input_Path())//trim(chemicalStateFile)
       call Chemical_State_CIE_File_Read(chemicalStateFileVarString,metallicityMaximumTabulated=metallicityMaximum)

       ! Flag that transfer function is now initialized.
       chemicalStateInitialized=.true.
    end if
    !$omp end critical (Chemical_State_Atomic_CIE_Cloudy_Initialize)
    return
  end subroutine Chemical_State_Atomic_CIE_Cloudy_Create

  double precision function Electron_Density_Atomic_CIE_Cloudy(temperature,numberDensityHydrogen,gasAbundances,radiation)
    !% Return the electron density assuming atomic CIE as computed by {\sc Cloudy}.
    use Chemical_States_CIE_File
    use Abundances_Structure
    use Radiation_Structure
    implicit none
    double precision                    , intent(in   ) :: numberDensityHydrogen, temperature
    type            (abundances        ), intent(in   ) :: gasAbundances
    type            (radiationStructure), intent(in   ) :: radiation

    ! Create the chemical state.
    call Chemical_State_Atomic_CIE_Cloudy_Create(gasAbundances)

    ! Call routine to interpolate in the tabulated function.
    Electron_Density_Atomic_CIE_Cloudy=Electron_Density_CIE_File_Interpolate(temperature,numberDensityHydrogen,gasAbundances,radiation)

    return
  end function Electron_Density_Atomic_CIE_Cloudy

  double precision function Electron_Density_Temperature_Log_Slope_Atomic_CIE_Cloudy(temperature,numberDensityHydrogen,gasAbundances&
       &,radiation)
    !% Return the logarithmic slope of the electron density with respect to temperature assuming atomic CIE as computed by {\sc Cloudy}.
    use Chemical_States_CIE_File
    use Abundances_Structure
    use Radiation_Structure
    implicit none
    double precision                    , intent(in   ) :: numberDensityHydrogen, temperature
    type            (abundances        ), intent(in   ) :: gasAbundances
    type            (radiationStructure), intent(in   ) :: radiation

    ! Create the chemical state.
    call Chemical_State_Atomic_CIE_Cloudy_Create(gasAbundances)

    ! Call routine to interpolate in the tabulated function.
    Electron_Density_Temperature_Log_Slope_Atomic_CIE_Cloudy=Electron_Density_CIE_File_logTemperature_Interpolate(temperature&
         &,numberDensityHydrogen ,gasAbundances,radiation)

    return
  end function Electron_Density_Temperature_Log_Slope_Atomic_CIE_Cloudy

  double precision function Electron_Density_Density_Log_Slope_Atomic_CIE_Cloudy(temperature,numberDensityHydrogen,gasAbundances&
       &,radiation)
    !% Return the logarithmic slope of the electron density with respect to density assuming atomic CIE as computed by {\sc Cloudy}.
    use Abundances_Structure
    use Radiation_Structure
    implicit none
    double precision                    , intent(in   ) :: numberDensityHydrogen, temperature
    type            (abundances        ), intent(in   ) :: gasAbundances
    type            (radiationStructure), intent(in   ) :: radiation

    ! Electron density always scales as total density under CIE conditions.
    Electron_Density_Density_Log_Slope_Atomic_CIE_Cloudy=1.0d0

    return
  end function Electron_Density_Density_Log_Slope_Atomic_CIE_Cloudy

  subroutine Chemical_Densities_Atomic_CIE_Cloudy(theseAbundances,temperature,numberDensityHydrogen,gasAbundances,radiation)
    !% Return the densities of chemical species at the given temperature and hydrogen density for the specified set of abundances
    !% and radiation field. Units of the returned electron density are cm$^-3$.
    use Chemical_States_CIE_File
    use Abundances_Structure
    use Radiation_Structure
    use Chemical_Abundances_Structure
    implicit none
    type            (chemicalAbundances), intent(inout) :: theseAbundances
    double precision                    , intent(in   ) :: numberDensityHydrogen, temperature
    type            (abundances        ), intent(in   ) :: gasAbundances
    type            (radiationStructure), intent(in   ) :: radiation

    ! Create the chemical state.
    call Chemical_State_Atomic_CIE_Cloudy_Create(gasAbundances)

    ! Call routine to interpolate in the tabulated function.
    call Chemical_Densities_CIE_File_Interpolate(theseAbundances,temperature,numberDensityHydrogen,gasAbundances,radiation)

    return
  end subroutine Chemical_Densities_Atomic_CIE_Cloudy

end module Chemical_States_Atomic_CIE_Cloudy
