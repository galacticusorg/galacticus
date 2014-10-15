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

!% Implements a null modifier for power spectra in the halo model of clustering.

  !# <haloModelPowerSpectrumModifier name="haloModelPowerSpectrumModifierNull">
  !#  <description>A null modifier for power spectra in the halo model of clustering.</description>
  !# </haloModelPowerSpectrumModifier>

  use Tables

  type, extends(haloModelPowerSpectrumModifierClass) :: haloModelPowerSpectrumModifierNull
   contains
     procedure :: modify => nullModify
  end type haloModelPowerSpectrumModifierNull
  
  interface haloModelPowerSpectrumModifierNull
     !% Constructor for the null halo model power spectra modifier class.
     module procedure nullConstructor
  end interface haloModelPowerSpectrumModifierNull

contains

  function nullConstructor()
    !% Default constructor for the null halo model power spectra modifier class.
    use Input_Parameters
    implicit none
    type(haloModelPowerSpectrumModifierNull) :: nullConstructor

    return
  end function nullConstructor
  
  subroutine nullModify(self,wavenumber,term,powerSpectrum,powerSpectrumCovariance,mass)
    !% Applies a null modification to a halo model power spectrum.
    implicit none
    class           (haloModelPowerSpectrumModifierNull), intent(inout)                           :: self
    double precision                                    , intent(in   ), dimension(:  )           :: wavenumber
    integer                                             , intent(in   )                           :: term
    double precision                                    , intent(inout), dimension(:  )           :: powerSpectrum
    double precision                                    , intent(inout), dimension(:,:), optional :: powerSpectrumCovariance
    double precision                                    , intent(in   )                , optional :: mass
    
    ! Do nothing, except to set covariance to zero.
    if (present(powerSpectrumCovariance)) powerSpectrumCovariance=0.0d0
    return
  end subroutine nullModify
