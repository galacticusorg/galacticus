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

!% Contains a module that implements calculations of the chemical state.

module Chemical_States
  !% Implements calculations of the chemical state.
  use Abundances_Structure
  use Radiation_Structure
  use Chemical_Abundances_Structure
  use ISO_Varying_String
  implicit none
  private
  public :: Electron_Density, Electron_Density_Temperature_Log_Slope, Electron_Density_Density_Log_Slope, Chemical_Densities

  ! Flag to indicate if this module has been initialized.
  logical                                              :: chemicalStateInitialized                  =.false.

  ! Name of chemical state method used.
  type     (varying_string                  )          :: chemicalStateMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Electron_Density_Get_Template   ), pointer :: Electron_Density_Get                      =>null()
  procedure(Electron_Density_Get_Template   ), pointer :: Electron_Density_Temperature_Log_Slope_Get=>null()
  procedure(Electron_Density_Get_Template   ), pointer :: Electron_Density_Density_Log_Slope_Get    =>null()
  procedure(Chemical_Densities_Get_Template ), pointer :: Chemical_Densities_Get                    =>null()
  abstract interface
     double precision function Electron_Density_Get_Template(temperature,numberDensityHydrogen,gasAbundances,radiation)
       import abundances,radiationStructure
       double precision                    , intent(in   ) :: numberDensityHydrogen, temperature
       type            (abundances        ), intent(in   ) :: gasAbundances
       type            (radiationStructure), intent(in   ) :: radiation
     end function Electron_Density_Get_Template
  end interface
  abstract interface
     subroutine Chemical_Densities_Get_Template(theseAbundances,temperature,numberDensityHydrogen,gasAbundances,radiation)
       import abundances,radiationStructure,chemicalAbundances
       type            (chemicalAbundances), intent(inout) :: theseAbundances
       double precision                    , intent(in   ) :: numberDensityHydrogen, temperature
       type            (abundances        ), intent(in   ) :: gasAbundances
       type            (radiationStructure), intent(in   ) :: radiation
     end subroutine Chemical_Densities_Get_Template
  end interface

contains

  subroutine Chemical_State_Initialize
    !% Initialize the chemical state module.
    use Galacticus_Error
    use Input_Parameters
    use Memory_Management
    !# <include directive="chemicalStateMethod" type="moduleUse">
    include 'chemical.state.modules.inc'
    !# </include>
    implicit none

    ! Initialize if necessary.
    if (.not.chemicalStateInitialized) then
       !$omp critical(Chemical_State_Initialization)
       if (.not.chemicalStateInitialized) then
          ! Get the chemical state method parameter.
          !@ <inputParameter>
          !@   <name>chemicalStateMethod</name>
          !@   <defaultValue>atomicCIECloudy</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be used for computing the chemical state.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('chemicalStateMethod',chemicalStateMethod,defaultValue='atomicCIECloudy')

          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="chemicalStateMethod" type="functionCall" functionType="void">
          !#  <functionArgs>chemicalStateMethod,Electron_Density_Get,Electron_Density_Temperature_Log_Slope_Get,Electron_Density_Density_Log_Slope_Get,Chemical_Densities_Get</functionArgs>
          include 'atomic.chemical_state.inc'
          !# </include>
          if (.not.(associated(Electron_Density_Get).and.associated(Electron_Density_Temperature_Log_Slope_Get) &
               & .and.associated(Electron_Density_Density_Log_Slope_Get).and.associated(Chemical_Densities_Get))) call&
               & Galacticus_Error_Report('Chemical_State_Initialize','method '//char(chemicalStateMethod)//' is unrecognized')
          chemicalStateInitialized=.true.
       end if
       !$omp end critical(Chemical_State_Initialization)
    end if
    return
  end subroutine Chemical_State_Initialize

  double precision function Electron_Density(temperature,numberDensityHydrogen,gasAbundances,radiation)
    !% Return the electron density at the given temperature and hydrogen density for the specified set of abundances and radiation
    !% field. Units of the returned electron density are cm$^-3$.
    implicit none
    double precision                    , intent(in   ) :: numberDensityHydrogen, temperature
    type            (abundances        ), intent(in   ) :: gasAbundances
    type            (radiationStructure), intent(in   ) :: radiation

    ! Initialize the module.
    call Chemical_State_Initialize

    ! Call the routine to do the calculation.
    Electron_Density=Electron_Density_Get(temperature,numberDensityHydrogen,gasAbundances,radiation)

    return
  end function Electron_Density

  double precision function Electron_Density_Temperature_Log_Slope(temperature,numberDensityHydrogen,gasAbundances,radiation)
    !% Return the logarithmic gradient of electron density with temperature at the given temperature and hydrogen density for the
    !% specified set of abundances and radiation field. Units of the returned electron density are cm$^-3$.
    implicit none
    double precision                    , intent(in   ) :: numberDensityHydrogen, temperature
    type            (abundances        ), intent(in   ) :: gasAbundances
    type            (radiationStructure), intent(in   ) :: radiation

    ! Initialize the module.
    call Chemical_State_Initialize

    ! Call the routine to do the calculation.
    Electron_Density_Temperature_Log_Slope=Electron_Density_Temperature_Log_Slope_Get(temperature,numberDensityHydrogen&
         &,gasAbundances,radiation)

    return
  end function Electron_Density_Temperature_Log_Slope

  double precision function Electron_Density_Density_Log_Slope(temperature,numberDensityHydrogen,gasAbundances,radiation)
    !% Return the logarithmic gradient of electron density with respect to density at the given temperature and hydrogen density
    !% for the specified set of abundances and radiation field. Units of the returned electron density are cm$^-3$.
    implicit none
    double precision                    , intent(in   ) :: numberDensityHydrogen, temperature
    type            (abundances        ), intent(in   ) :: gasAbundances
    type            (radiationStructure), intent(in   ) :: radiation

    ! Initialize the module.
    call Chemical_State_Initialize

    ! Call the routine to do the calculation.
    Electron_Density_Density_Log_Slope=Electron_Density_Density_Log_Slope_Get(temperature,numberDensityHydrogen&
         &,gasAbundances,radiation)

    return
  end function Electron_Density_Density_Log_Slope

  subroutine Chemical_Densities(theseAbundances,temperature,numberDensityHydrogen,gasAbundances,radiation)
    !% Return the densities of chemical species at the given temperature and hydrogen density for the specified set of abundances
    !% and radiation field. Units of the returned electron density are cm$^-3$.
    implicit none
    type            (chemicalAbundances), intent(inout) :: theseAbundances
    double precision                    , intent(in   ) :: numberDensityHydrogen, temperature
    type            (abundances        ), intent(in   ) :: gasAbundances
    type            (radiationStructure), intent(in   ) :: radiation

    ! Initialize the module.
    call Chemical_State_Initialize

    ! Call the routine to do the calculation.
    call Chemical_Densities_Get(theseAbundances,temperature,numberDensityHydrogen,gasAbundances,radiation)

    return
  end subroutine Chemical_Densities

end module Chemical_States
