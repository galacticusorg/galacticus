!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements a null post-processing of stellar spectra.

module Stellar_Population_Spectra_Postprocess_Null
  !% Implements a null post-processing of stellar spectra.
  implicit none
  private
  public :: Stellar_Population_Spectra_Postprocess_Null_Initialize

contains

  !# <stellarPopulationSpectraPostprocessMethod>
  !#  <unitName>Stellar_Population_Spectra_Postprocess_Null_Initialize</unitName>
  !# </stellarPopulationSpectraPostprocessMethod>
  subroutine Stellar_Population_Spectra_Postprocess_Null_Initialize(stellarPopulationSpectraPostprocessMethod,Stellar_Population_Spectra_Postprocess_Get)
    !% Initializes the ``Null'' stellar spectrum postprocessing module.
    use ISO_Varying_string
    implicit none
    type     (varying_string  ), intent(in   )          :: stellarPopulationSpectraPostprocessMethod
    procedure(double precision), intent(inout), pointer :: Stellar_Population_Spectra_Postprocess_Get

    if (stellarPopulationSpectraPostprocessMethod == 'null') Stellar_Population_Spectra_Postprocess_Get =>&
         & Stellar_Population_Spectra_Postprocess_Null_Get
    return
  end subroutine Stellar_Population_Spectra_Postprocess_Null_Initialize

  double precision function Stellar_Population_Spectra_Postprocess_Null_Get(wavelength,redshift)
    !% Null post-processing of stellar spectra: always returns unity.
    implicit none
    double precision, intent(in   ) :: redshift, wavelength

    Stellar_Population_Spectra_Postprocess_Null_Get=1.0d0
    return
  end function Stellar_Population_Spectra_Postprocess_Null_Get

end module Stellar_Population_Spectra_Postprocess_Null
