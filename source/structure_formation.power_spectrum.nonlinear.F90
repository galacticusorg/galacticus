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

!% Contains a module which implements the nonlinear power spectrum.

module Power_Spectra_Nonlinear
  !% Implements the nonlinear power spectrum.
  use ISO_Varying_String
  use Galacticus_Nodes
  implicit none
  private
  public :: Power_Spectrum_Nonlinear
  
  ! Flag to indicate if this module has been initialized.  
  logical              :: moduleIsInitialized=.false.

  ! Name of cooling rate available method used.
  type(varying_string) :: powerSpectrumNonlinearMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Power_Spectrum_Nonlinear), pointer :: Power_Spectrum_Nonlinear_Get => null()

contains

  subroutine Power_Spectrum_Nonlinear_Initialize
    !% Initialize the nonlinear power spectrum module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="powerSpectrumNonlinearMethod" type="moduleUse">
    include 'structure_formation.power_spectrum.nonlinear.modules.inc'
    !# </include>
    implicit none

    ! Initialize if necessary.
    if (.not.moduleIsInitialized) then
       !$omp critical(Power_Spectra_Nonlinear_Initialization) 
       if (.not.moduleIsInitialized) then
          ! Get the spheroid star formation timescale method parameter.
          !@ <inputParameter>
          !@   <name>powerSpectrumNonlinearMethod</name>
          !@   <defaultValue>CosmicEmu</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be used for computing the nonlinear power spectrum.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@   <group>starFormation</group>
          !@ </inputParameter>
          call Get_Input_Parameter('powerSpectrumNonlinearMethod',powerSpectrumNonlinearMethod,defaultValue='CosmicEmu')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="powerSpectrumNonlinearMethod" type="functionCall" functionType="void">
          !#  <functionArgs>powerSpectrumNonlinearMethod,Power_Spectrum_Nonlinear_Get</functionArgs>
          include 'structure_formation.power_spectrum.nonlinear.inc'
          !# </include>
          if (.not.associated(Power_Spectrum_Nonlinear_Get)) call Galacticus_Error_Report('Power_Spectrum_Nonlinear_Initialize'&
               &,'method '//char(powerSpectrumNonlinearMethod)//' is unrecognized')
          moduleIsInitialized=.true.
       end if
       !$omp end critical(Power_Spectra_Nonlinear_Initialization) 
    end if
    return
  end subroutine Power_Spectrum_Nonlinear_Initialize

  double precision function Power_Spectrum_Nonlinear(waveNumber,time)
    !% Return the nonlinear power spectrum for $k=${\tt wavenumber} [Mpc$^{-1}$] at the given cosmic {\tt time} [Gyr].
    implicit none
    double precision, intent(in) :: waveNumber,time

    ! Initialize the module.
    call Power_Spectrum_Nonlinear_Initialize()

    ! Get the power spectrum using the selected method.
    Power_Spectrum_Nonlinear=Power_Spectrum_Nonlinear_Get(waveNumber,time)

    return
  end function Power_Spectrum_Nonlinear
    
end module Power_Spectra_Nonlinear
