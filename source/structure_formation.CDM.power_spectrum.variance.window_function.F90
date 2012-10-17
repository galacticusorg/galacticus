!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements window functions for computing the variance of the power spectrum.

module Power_Spectrum_Window_Functions
  !% Implements window functions for computing the variance of the power spectrum.
  use, intrinsic :: ISO_C_Binding
  use ISO_Varying_String
  implicit none
  private
  public :: Power_Spectrum_Window_Function
  
  ! Flag to indicate if this module has been initialized.  
  logical              :: moduleInitialized=.false.

  ! Name of cooling rate available method used.
  type(varying_string) :: powerSpectrumWindowFunctionMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Power_Spectrum_Window_Function), pointer :: Power_Spectrum_Window_Function_Get => null()

contains

  subroutine Power_Spectrum_Window_Functions_Initialize
    !% Initialize the power spectrum window function module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="powerSpectrumWindowFunctionMethod" type="moduleUse">
    include 'structure_formation.CDM.power_spectrum.variance.window_function.modules.inc'
    !# </include>
    implicit none

    ! Initialize if necessary.
    if (.not.moduleInitialized) then
       !$omp critical(Power_Spectrum_Window_Functions_Initialization) 
       if (.not.moduleInitialized) then
          ! Get the window function method parameter.
          !@ <inputParameter>
          !@   <name>powerSpectrumWindowFunctionMethod</name>
          !@   <defaultValue>topHat</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be used for computing window functions to estimate variance from the power spectrum.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('powerSpectrumWindowFunctionMethod',powerSpectrumWindowFunctionMethod,defaultValue='topHat')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="powerSpectrumWindowFunctionMethod" type="code" action="subroutine">
          !#  <subroutineArgs>powerSpectrumWindowFunctionMethod,Power_Spectrum_Window_Function_Get</subroutineArgs>
          include 'structure_formation.CDM.power_spectrum.variance.window_function.inc'
          !# </include>
          if (.not.associated(Power_Spectrum_Window_Function_Get)) call Galacticus_Error_Report('Power_Spectrum_Window_Functions_Initialize'&
               &,'method '//char(powerSpectrumWindowFunctionMethod)//' is unrecognized')
          moduleInitialized=.true.
       end if
       !$omp end critical(Power_Spectrum_Window_Functions_Initialization) 
    end if
    return
  end subroutine Power_Spectrum_Window_Functions_Initialize

  double precision function Power_Spectrum_Window_Function(wavenumber,smoothingMass)
    !% Returns the window function for power spectrum variance computation at the specified {\tt wavenumber} (in Mpc$^{-1}$) for a
    !% given {\tt smoothingMass} (in $M_\odot$).
    implicit none
    double precision, intent(in) :: wavenumber,smoothingMass

    ! Initialize the module.
    call Power_Spectrum_Window_Functions_Initialize

    ! Call the function that does the work.
    Power_Spectrum_Window_Function=Power_Spectrum_Window_Function_Get(wavenumber,smoothingMass)
    return
  end function Power_Spectrum_Window_Function
  
end module Power_Spectrum_Window_Functions
