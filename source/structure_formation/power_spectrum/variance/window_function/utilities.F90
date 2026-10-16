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
Contains a module of utilities needed by ``powerSpectrumWindowFunction`` classes.
!!}

module Power_Spectrum_Window_Function_Utilities
  !!{RST
  Provides utilities needed by ``powerSpectrumWindowFunction`` classes.
  !!}
  private
  public :: Window_Function_Top_Hat, Window_Function_Radius_Lagrangian

contains

  double precision function Window_Function_Top_Hat(x) result(windowFunction)
    !!{RST
    Return the Fourier transform into :math:`k`-space of a top hat in real space,

    .. math::
     W(x) = \frac{3 \left( \sin x - x \cos x \right)}{x^3},

    where :math:`x = k r` for a top hat of radius :math:`r`. A series expansion,

    .. math::
     W(x) = 1 - \frac{x^2}{10} + \frac{x^4}{280} - \frac{x^6}{15120},

    is used for small :math:`x`, where the full expression suffers from cancellation, and zero is returned for
    non-positive :math:`x`.
    !!}
    implicit none
    double precision, intent(in   ) :: x
    double precision, parameter     :: xSeriesMaximum=1.0d-3
    double precision                :: xSquared

    if      (x <= 0.0d0         ) then
       windowFunction=+0.0d0
    else if (x <= xSeriesMaximum) then
       ! Use a series expansion of the window function for small x.
       xSquared      =+x**2
       windowFunction=+1.0d0                        &
            &         +xSquared*(  -1.0d0/   10.0d0 &
            &         +xSquared* ( +1.0d0/  280.0d0 &
            &         +xSquared*  (-1.0d0/15120.0d0 &
            &                     )                 &
            &                    )                  &
            &                   )
    else
       ! For larger x, use the full expression.
       windowFunction=3.0d0*(sin(x)-x*cos(x))/(x**3)
    end if
    return
  end function Window_Function_Top_Hat

  double precision function Window_Function_Radius_Lagrangian(smoothingMass,cosmologyParameters_) result(radiusLagrangian)
    !!{RST
    Return the Lagrangian radius enclosing the given ``smoothingMass`` at the mean density of the universe,

    .. math::
     r = \left( \frac{3 M}{4 \pi \Omega_\mathrm{M} \rho_\mathrm{crit}} \right)^{1/3}.
    !!}
    use :: Cosmology_Parameters    , only : cosmologyParametersClass
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision                          , intent(in   ) :: smoothingMass
    class           (cosmologyParametersClass), intent(inout) :: cosmologyParameters_

    radiusLagrangian=+(                                        &
         &             +3.0d0                                  &
         &             /4.0d0                                  &
         &             /Pi                                     &
         &             *smoothingMass                          &
         &             /cosmologyParameters_%OmegaMatter    () &
         &             /cosmologyParameters_%densityCritical() &
         &            )**(1.0d0/3.0d0)
    return
  end function Window_Function_Radius_Lagrangian

end module Power_Spectrum_Window_Function_Utilities
