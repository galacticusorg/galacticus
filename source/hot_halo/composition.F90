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
Contains a module which provides the composition of gas in the hot atmosphere of a node.
!!}

module Hot_Halo_Composition
  !!{RST
  Provides the composition of gas in the hot atmosphere of a node, as needed to evaluate its cooling and radiative properties.
  !!}
  implicit none
  private
  public :: Hot_Halo_Gas_Composition, Hydrogen_Number_Density

contains

  subroutine Hot_Halo_Gas_Composition(node,abundancesGas,fractionsChemical)
    !!{RST
    Return the abundances of gas in the hot atmosphere of ``node`` as mass fractions and, optionally, its chemical abundances as
    number densities (in cm⁻³) per unit total mass density (in M☉ Mpc⁻³).
    !!}
    use :: Abundances_Structure             , only : abundances
    use :: Chemical_Abundances_Structure    , only : chemicalAbundances
    use :: Chemical_Reaction_Rates_Utilities, only : Chemicals_Mass_To_Fraction_Conversion
    use :: Galacticus_Nodes                 , only : nodeComponentHotHalo                 , treeNode
    implicit none
    type            (treeNode            ), intent(inout)           :: node
    type            (abundances          ), intent(  out)           :: abundancesGas
    type            (chemicalAbundances  ), intent(  out), optional :: fractionsChemical
    class           (nodeComponentHotHalo), pointer                 :: hotHalo
    type            (chemicalAbundances  )                          :: massChemical
    double precision                                                :: massToDensityConversion

    ! Get the abundances as mass fractions.
    hotHalo       => node   %hotHalo   ()
    abundancesGas =  hotHalo%abundances()
    call abundancesGas%massToMassFraction(hotHalo%mass())
    ! Get the chemicals, if requested.
    if (present(fractionsChemical)) then
       ! Scale all chemical masses by their mass in atomic mass units to get a number density.
       massChemical=hotHalo%chemicals()
       call massChemical%massToNumber(fractionsChemical)
       ! Compute factor converting mass of chemicals in (M☉) to number density in cm⁻³ per total mass density.
       if (hotHalo%mass() > 0.0d0) then
          massToDensityConversion=Chemicals_Mass_To_Fraction_Conversion(hotHalo%mass())
       else
          massToDensityConversion=0.0d0
       end if
       ! Convert to number density per unit total mass density.
       fractionsChemical=fractionsChemical*massToDensityConversion
    end if
    return
  end subroutine Hot_Halo_Gas_Composition

  double precision function Hydrogen_Number_Density(density,abundancesGas) result(numberDensityHydrogen)
    !!{RST
    Return the number density of hydrogen (in cm⁻³) in gas of the given ``density`` (in M☉ Mpc⁻³) and abundances (given as mass
    fractions).
    !!}
    use :: Abundances_Structure            , only : abundances
    use :: Numerical_Constants_Astronomical, only : massSolar       , megaParsec
    use :: Numerical_Constants_Atomic      , only : massHydrogenAtom
    use :: Numerical_Constants_Prefixes    , only : hecto
    implicit none
    double precision            , intent(in   ) :: density
    type            (abundances), intent(in   ) :: abundancesGas

    numberDensityHydrogen=+density                                    &
         &                *abundancesGas   %hydrogenMassFraction()    &
         &                *massSolar                                  &
         &                /massHydrogenAtom                           &
         &                /hecto                                  **3 &
         &                /megaParsec                             **3
    return
  end function Hydrogen_Number_Density

end module Hot_Halo_Composition
