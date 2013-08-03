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


!% Contains a module which implements postprocessing of stellar spectra to keep only recent populations.

module Stellar_Population_Spectra_Postprocessing_Identity
  !% Implements postprocessing of stellar spectra to keep only identity populations.
  use ISO_Varying_String
  public :: Stellar_Population_Spectra_Postprocess_Identity_Init

contains

  !# <stellarPopulationSpectraPostprocessInitialize>
  !#  <unitName>Stellar_Population_Spectra_Postprocess_Identity_Init</unitName>
  !# </stellarPopulationSpectraPostprocessInitialize>
  subroutine Stellar_Population_Spectra_Postprocess_Identity_Init(stellarPopulationSpectraPostprocessMethod,postprocessingFunction)
    !% Initializes the ``identity'' stellar spectrum postprocessing module.
    implicit none
    type     (varying_string),          intent(in   ) :: stellarPopulationSpectraPostprocessMethod
    procedure(              ), pointer, intent(inout) :: postprocessingFunction

    if (stellarPopulationSpectraPostprocessMethod == 'identity') postprocessingFunction => Stellar_Population_Spectra_Postprocess_Identity       
    return
  end subroutine Stellar_Population_Spectra_Postprocess_Identity_Init
  
  subroutine Stellar_Population_Spectra_Postprocess_Identity(wavelength,age,redshift,modifier)
    !% An identity operator for postprocessing of stellar spectra (i.e. does nothing).
    use Numerical_Constants_Atomic
    implicit none
    double precision, intent(in   ) :: wavelength,age,redshift
    double precision, intent(inout) :: modifier

    return
  end subroutine Stellar_Population_Spectra_Postprocess_Identity

end module Stellar_Population_Spectra_Postprocessing_Identity
