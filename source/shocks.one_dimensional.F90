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
Contains a module which implements calculations of one-dimensional shocks.
!!}

module Shocks_1D
  !!{
  iImplements calculations of one-dimensional shocks.
  !!}
  private
  public :: Shocks_1D_Density_Jump

  ! Effective infinite Mach number for strong shocks.
  double precision, parameter, public :: machNumberInfinite=huge(1.0d0)

contains

  double precision function Shocks_1D_Density_Jump(adiabaticIndex,machNumberPreShock)
    !!{
    Computes the density jump across a one-dimensional shock.
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
