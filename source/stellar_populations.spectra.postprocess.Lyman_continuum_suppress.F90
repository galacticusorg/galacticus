!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@caltech.edu>
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


!% Contains a module which implements suppression of the Lyman continuum in galaxy spectra.

module Stellar_Population_Spectra_Postprocessing_Lyc_Suppress
  !% Implements suppression of the Lyman continuum in galaxy spectra.
  use ISO_Varying_String
  public :: Stellar_Population_Spectra_Postprocess_Lyc_Suppress_Initialize,Stellar_Population_Spectra_Postprocess_Lyc_Suppress

  ! Record of whether this method is active.
  logical :: methodIsActive

contains
  
  !# <stellarPopulationSpectraPostprocessInitialize>
  !#  <unitName>Stellar_Population_Spectra_Postprocess_Lyc_Suppress_Initialize</unitName>
  !# </stellarPopulationSpectraPostprocessInitialize>
  subroutine Stellar_Population_Spectra_Postprocess_Lyc_Suppress_Initialize(stellarPopulationSpectraPostprocessMethod,postprocessingFunction)
    !% Initializes the ``Lyman-continuum suppression'' stellar spectrum postprocessing module.
    implicit none
    type     (varying_string),          intent(in   ) :: stellarPopulationSpectraPostprocessMethod
    procedure(              ), pointer, intent(inout) :: postprocessingFunction

    if (stellarPopulationSpectraPostprocessMethod == 'lymanContinuumSuppress') postprocessingFunction => Stellar_Population_Spectra_Postprocess_Lyc_Suppress
    return
  end subroutine Stellar_Population_Spectra_Postprocess_Lyc_Suppress_Initialize

  subroutine Stellar_Population_Spectra_Postprocess_Lyc_Suppress(wavelength,age,redshift,modifier)
    !% Suppresses all starlight in the Lyman continuum.
    use Numerical_Constants_Atomic
    implicit none
    double precision, intent(in)    :: wavelength,age,redshift
    double precision, intent(inout) :: modifier

    ! Suppress all emission in the Lyman continuum.
    if (methodIsActive .and. wavelength < ionizationWavelengthHydrogen) modifier=0.0d0
    return
  end subroutine Stellar_Population_Spectra_Postprocess_Lyc_Suppress

end module Stellar_Population_Spectra_Postprocessing_Lyc_Suppress
