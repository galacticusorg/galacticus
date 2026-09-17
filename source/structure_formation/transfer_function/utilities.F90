!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
!!    Andrew Benson <abenson@carnegiescience.edu>
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

!+    Contributions to this file made by: Claude.

!!{RST
Contains a module of utilities needed by ``transferFunction`` classes.
!!}

module Transfer_Function_Utilities
  !!{RST
  Provides utilities needed by ``transferFunction`` classes.
  !!}
  private
  public :: Transfer_Function_Mass_From_Wavenumber

contains

  double precision function Transfer_Function_Mass_From_Wavenumber(wavenumber,cosmologyParameters_) result(mass)
    !!{RST
    Return the mass scale corresponding to the given ``wavenumber``,

    .. math::
     M = \frac{4 \pi}{3} \Omega_\mathrm{M} \rho_\mathrm{crit} \left( \frac{\pi}{k} \right)^3.

    As a default choice, the wavenumber is converted to a length scale assuming :math:`R = \lambda/2 = \pi/k` [see Eq.(9) of
    :cite:t:`schneider_non-linear_2012`].
    !!}
    use :: Cosmology_Parameters    , only : cosmologyParametersClass
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision                          , intent(in   ) :: wavenumber
    class           (cosmologyParametersClass), intent(inout) :: cosmologyParameters_
    double precision                                          :: matterDensity

    matterDensity=+cosmologyParameters_%OmegaMatter    () &
         &        *cosmologyParameters_%densityCritical()
    mass         =+4.0d0            &
         &        *Pi               &
         &        /3.0d0            &
         &        *matterDensity    &
         &        *(                &
         &          +Pi             &
         &          /wavenumber     &
         &         )**3
    return
  end function Transfer_Function_Mass_From_Wavenumber

end module Transfer_Function_Utilities
