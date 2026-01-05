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
Contains a module that implements various useful utility functions for calculations of chemical abundances and rates.
!!}

module Chemical_Reaction_Rates_Utilities
  !!{
  Implements various useful utility functions for calculations of chemical abundances and rates.
  !!}
  implicit none
  private
  public :: Chemicals_Mass_To_Density_Conversion, Chemicals_Mass_To_Fraction_Conversion

contains

  double precision function Chemicals_Mass_To_Density_Conversion(radius)
    !!{
    Returns the conversion factor from mass of chemicals in ($M_\odot$) to number density in cm$^{-3}$ assuming
    that the mass is distributed uniformly in a sphere of the given {\normalfont \ttfamily radius} (in Mpc).
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar     , megaParsec
    use :: Numerical_Constants_Atomic      , only : atomicMassUnit
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Prefixes    , only : hecto
    implicit none
    double precision, intent(in   ) :: radius

    Chemicals_Mass_To_Density_Conversion=3.0d0*massSolar/atomicMassUnit/4.0d0/Pi/(hecto*megaParsec*radius)**3
    return
  end function Chemicals_Mass_To_Density_Conversion

  double precision function Chemicals_Mass_To_Fraction_Conversion(massTotal)
    !!{
    Returns the conversion factor from mass of chemicals in ($M_\odot$) to a number density (in cm$^{-3}$) per unit total mass
    density given a total {\normalfont \ttfamily mass} (in $\mathrm{M}_\odot$).
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar     , megaParsec
    use :: Numerical_Constants_Atomic      , only : atomicMassUnit
    use :: Numerical_Constants_Prefixes    , only : hecto
    implicit none
    double precision, intent(in   ) :: massTotal

    Chemicals_Mass_To_Fraction_Conversion=massSolar/atomicMassUnit/(hecto*megaParsec)**3/massTotal
    return
  end function Chemicals_Mass_To_Fraction_Conversion

end module Chemical_Reaction_Rates_Utilities
