!! Copyright 2009, 2010, 2011 Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


!% Contains a module which generates a tabulated atomic collisional ionization equilibrium ionization state using {\sc Cloudy}.

module Chemical_States_Atomic_CIE_Cloudy
  !% Generates a tabulated atomic collisional ionization equilibrium ionization state using {\sc Cloudy}.
  use ISO_Varying_String
  implicit none
  private
  public :: Chemical_State_Atomic_CIE_Cloudy_Initialize
  
  ! Flag to indicate if this module has been initialized.
  logical                     :: chemicalStateInitialized=.false.

  ! File names for the cooling function and chemical state data.
  character(len=50)           :: coolingFunctionFile='data/cooling_function_Atomic_CIE_Cloudy.xml'
  character(len=50)           :: chemicalStateFile='data/chemical_state_Atomic_CIE_Cloudy.xml'

  ! Maximum tabulated metallicity.
  double precision, parameter :: metallicityMaximumDefault=30.0d0 ! Thirty times Solar.
  double precision            :: metallicityMaximum

  ! Maximum metallicity that we will ever tabulate. More than this should be physically implausible.
  double precision, parameter :: metallicityMaximumLimit  =30.0d0 ! Thirty times Solar.

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
    type(varying_string),                 intent(in)    :: chemicalStateMethod
    procedure(double precision), pointer, intent(inout) :: Electron_Density_Get,Electron_Density_Temperature_Log_Slope_Get&
         &,Electron_Density_Density_Log_Slope_Get
    procedure(),                 pointer, intent(inout) :: Chemical_Densities_Get
 
    ! Check if this chemical state has been selected.
    if (chemicalStateMethod == 'atomicCIECloudy') then
       Electron_Density_Get                       => Electron_Density_Atomic_CIE_Cloudy
       Electron_Density_Temperature_Log_Slope_Get => Electron_Density_Temperature_Log_Slope_Atomic_CIE_Cloudy
       Electron_Density_Density_Log_Slope_Get     => Electron_Density_Density_Log_Slope_Atomic_CIE_Cloudy
       Chemical_Densities_Get                    => Chemical_Densities_Atomic_CIE_Cloudy
   end if

    return
  end subroutine Chemical_State_Atomic_CIE_Cloudy_Initialize

  subroutine Chemical_State_Atomic_CIE_Cloudy_Create(abundances)
    !% Create the chemical state.
    use Chemical_States_CIE_File
    use Abundances_Structure
    use System_Command
    implicit none
    type(abundancesStructure), intent(in) :: abundances
    logical                               :: makeFile
    character(len=32)                     :: metallicityLabel
    type(varying_string)                  :: command,chemicalStateFileVarString

    ! Generate the name of the data file and an XML input parameter file.
    !$omp critical (Chemical_State_Atomic_CIE_Cloudy_Initialize)
    ! Determine if we need to reinitialize this module.
    if (.not.chemicalStateInitialized) then
       makeFile=.true.
    else
       if (Abundances_Get_Metallicity(abundances,linearByMassSolar) > 0.0d0) then
          makeFile=(min(Abundances_Get_Metallicity(abundances,linearByMassSolar),metallicityMaximumLimit) > metallicityMaximum&
               &*(1.0d0+metallicityTolerance))
       else
          makeFile=.false.
       end if
       if (makeFile) then
          ! Remove the transfer function file so that a new one will be created.
          command='rm -f '//trim(chemicalStateFile)
          call System_Command_Do(command)
       end if
    end if
    ! Read the file if this module has not been initialized or if the metallicity is out of range.
    if (makeFile) then
       ! Determine maximum metallicity to which we should tabulate.
       if (Abundances_Get_Metallicity(abundances,linearByMassSolar) > 0.0d0) then
          metallicityMaximum=min(max(metallicityMaximumDefault,3.0d0*Abundances_Get_Metallicity(abundances,linearByMassSolar))&
               &,metallicityMaximumLimit)
       else
          metallicityMaximum=metallicityMaximumDefault
       end if
       write (metallicityLabel,'(e12.6)') dlog10(metallicityMaximum)

       ! Run Atomic_CIE_Cloudy wrapper script.
       command='./scripts/aux/Atomic_CIE_Cloudy_Driver.pl '//metallicityLabel//' '//trim(coolingFunctionFile)//' '&
            &//trim(chemicalStateFile)
       call System_Command_Do(command)

       ! Call routine to read in the tabulated data.
       chemicalStateFileVarString=trim(chemicalStateFile)
       call Chemical_State_CIE_File_Read(chemicalStateFileVarString,metallicityMaximumTabulated=metallicityMaximum)

       ! Flag that transfer function is now initialized.
       chemicalStateInitialized=.true.
    end if
    !$omp end critical (Chemical_State_Atomic_CIE_Cloudy_Initialize)
    return
  end subroutine Chemical_State_Atomic_CIE_Cloudy_Create

  double precision function Electron_Density_Atomic_CIE_Cloudy(temperature,numberDensityHydrogen,abundances,radiation)
    !% Return the electron density assuming atomic CIE as computed by {\sc Cloudy}.
    use Chemical_States_CIE_File
    use Abundances_Structure
    use Radiation_Structure
    implicit none
    double precision,          intent(in) :: temperature,numberDensityHydrogen
    type(abundancesStructure), intent(in) :: abundances
    type(radiationStructure),  intent(in) :: radiation

    ! Create the chemical state.
    call Chemical_State_Atomic_CIE_Cloudy_Create(abundances)
    
    ! Call routine to interpolate in the tabulated function.
    Electron_Density_Atomic_CIE_Cloudy=Electron_Density_CIE_File_Interpolate(temperature,numberDensityHydrogen,abundances,radiation)
    
    return
  end function Electron_Density_Atomic_CIE_Cloudy

  double precision function Electron_Density_Temperature_Log_Slope_Atomic_CIE_Cloudy(temperature,numberDensityHydrogen,abundances&
       &,radiation)
    !% Return the logarithmic slope of the electron density with respect to temperature assuming atomic CIE as computed by {\sc Cloudy}.
    use Chemical_States_CIE_File
    use Abundances_Structure
    use Radiation_Structure
    implicit none
    double precision,          intent(in)  :: temperature,numberDensityHydrogen
    type(abundancesStructure), intent(in)  :: abundances
    type(radiationStructure),  intent(in)  :: radiation
    
    ! Create the chemical state.
    call Chemical_State_Atomic_CIE_Cloudy_Create(abundances)
    
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

  subroutine Chemical_Densities_Atomic_CIE_Cloudy(theseAbundances,temperature,numberDensityHydrogen,abundances,radiation)
    !% Return the densities of chemical species at the given temperature and hydrogen density for the specified set of abundances
    !% and radiation field. Units of the returned electron density are cm$^-3$.
    use Chemical_States_CIE_File
    use Abundances_Structure
    use Radiation_Structure
    use Chemical_Abundances_Structure
    implicit none
    type(chemicalAbundancesStructure), intent(inout) :: theseAbundances
    double precision,                  intent(in)    :: temperature,numberDensityHydrogen
    type(abundancesStructure),         intent(in)    :: abundances
    type(radiationStructure),          intent(in)    :: radiation

    ! Create the chemical state.
    call Chemical_State_Atomic_CIE_Cloudy_Create(abundances)
    
    ! Call routine to interpolate in the tabulated function.
    call Chemical_Densities_CIE_File_Interpolate(theseAbundances,temperature,numberDensityHydrogen,abundances,radiation)
 
    return
  end subroutine Chemical_Densities_Atomic_CIE_Cloudy
      
end module Chemical_States_Atomic_CIE_Cloudy
