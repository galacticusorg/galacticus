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

!!{
Contains a module which implements Rankine-Hugoniot jump conditions for one-dimensional (planar) shocks, returning the
density compression ratio as a function of upstream Mach number and adiabatic index.
!!}

module Shocks_1D
  !!{
  Implements calculations of one-dimensional shock jump conditions. Provides the Rankine-Hugoniot density jump ratio
  across a planar shock as a function of the upstream Mach number and adiabatic index, including the strong-shock
  (infinite Mach number) limit.
  !!}
  private
  public :: Shocks_1D_Density_Jump

  ! Effective infinite Mach number for strong shocks.
  double precision, parameter, public :: machNumberInfinite=huge(1.0d0)

contains

  double precision function Shocks_1D_Density_Jump(adiabaticIndex,machNumberPreShock)
    !!{
    Computes the density compression ratio $\rho_2/\rho_1$ across a one-dimensional planar shock using the
    Rankine-Hugoniot jump conditions. For a finite upstream Mach number $\mathcal{M}$, the ratio is
    $(\gamma+1)\mathcal{M}^2 / [(\gamma-1)\mathcal{M}^2 + 2]$; in the strong-shock limit
    ($\mathcal{M} \to \infty$) this reduces to $(\gamma+1)/(\gamma-1)$, where $\gamma$ is the adiabatic index.
    !!}
    implicit none
    double precision, intent(in   ) :: adiabaticIndex,machNumberPreShock

    if (machNumberPreShock >= machNumberInfinite) then
       Shocks_1D_Density_Jump=                               &
            &                  (adiabaticIndex       +1.0d0) &
            &                 /(adiabaticIndex       -1.0d0)
    else
       Shocks_1D_Density_Jump=                               &
            &                  (adiabaticIndex       +1.0d0) &
            &                 /(adiabaticIndex       -1.0d0) &
            &                 * machNumberPreShock**2        &
            &                 /(machNumberPreShock**2+2.0d0)
    end if
    return
  end function Shocks_1D_Density_Jump

end module Shocks_1D
