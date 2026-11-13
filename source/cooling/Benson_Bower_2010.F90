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

!+    Contributions to this file made by: Andrew Benson, Claude.

!!{RST
Contains a module that computes the properties of the notional hot halo used by the cooling model of :cite:t:`benson_galaxy_2010-1`.
!!}

module Cooling_Benson_Bower_2010
  !!{RST
  Computes the properties of the notional hot halo used by the cooling model of :cite:t:`benson_galaxy_2010-1`.
  !!}
  private
  public :: coolingBensonBower2010Properties

contains

  logical function coolingBensonBower2010Properties(node,coolingFunction_,chemicalState_,radiation,coolingFunction,numberDensityHydrogen,countParticles) result(exists)
    !!{RST
    Compute the properties of the notional hot halo used by the cooling model of :cite:t:`benson_galaxy_2010-1`. The notional
    hot halo contains the current hot halo gas, the gas outflowed from it, and all galactic mass, spread uniformly within the outer
    radius of the hot halo, at the temperature of the hot halo gas at that radius. Returns the cooling function, :math:`\lambda(T,Z)
    = \Lambda(T,Z) n_\mathrm{H}^2` (in ergs cm⁻³ s⁻¹), the number density of hydrogen, :math:`n_\mathrm{H}` (in cm⁻³), and the number
    of particles, :math:`N`, in the notional hot halo. The function result is ``.false.`` (and these properties are not set) if the
    node has no hot halo, if the notional hot halo has no mass, or if the hot halo has no mass distribution.
    !!}
    use :: Galacticus_Nodes                 , only : treeNode                               , nodeComponentBasic         , nodeComponentHotHalo
    use :: Abundances_Structure             , only : abundances
    use :: Chemical_Abundances_Structure    , only : chemicalAbundances                     , Chemicals_Property_Count
    use :: Chemical_Reaction_Rates_Utilities, only : Chemicals_Mass_To_Density_Conversion
    use :: Chemical_States                  , only : chemicalStateClass
    use :: Cooling_Functions                , only : coolingFunctionClass
    use :: Mass_Distributions               , only : massDistributionClass                  , kinematicsDistributionClass, massDistributionZero
    use :: Coordinates                      , only : coordinateSpherical                    , assignment(=)
    use :: Galactic_Structure_Options       , only : componentTypeHotHalo                   , massTypeGaseous            , massTypeGalactic
    use :: Numerical_Constants_Astronomical , only : megaParsec
    use :: Numerical_Constants_Prefixes     , only : centi
    use :: Numerical_Constants_Math         , only : Pi
    use :: Hot_Halo_Composition             , only : Hydrogen_Number_Density
    use :: Radiation_Fields                 , only : radiationFieldCosmicMicrowaveBackground
    implicit none
    type            (treeNode                               ), intent(inout) :: node
    class           (coolingFunctionClass                   ), intent(inout) :: coolingFunction_
    class           (chemicalStateClass                     ), intent(inout) :: chemicalState_
    type            (radiationFieldCosmicMicrowaveBackground), intent(inout) :: radiation
    double precision                                         , intent(  out) :: coolingFunction        , numberDensityHydrogen  , &
         &                                                                      countParticles
    class           (nodeComponentBasic                     ), pointer       :: basic
    class           (nodeComponentHotHalo                   ), pointer       :: hotHalo
    class           (massDistributionClass                  ), pointer       :: massDistribution_
    class           (kinematicsDistributionClass            ), pointer       :: kinematicsDistribution_
    type            (coordinateSpherical                    )                :: coordinates
    double precision                                                         :: density                , temperature            , &
         &                                                                      massToDensityConversion, numberDensityAllSpecies, &
         &                                                                      massNotional
    type            (abundances                             )                :: abundances_
    type            (chemicalAbundances                     )                :: chemicalDensities_     , chemicalMasses_

    exists  =  .false.
    hotHalo => node%hotHalo()
    select type (hotHalo)
    type is (nodeComponentHotHalo)
       ! Hot halo does not exist - nothing to do here.
       return
    end select
    ! Compute the mass in the notional hot halo.
    basic             =>  node             %basic           (                         )
    massDistribution_ =>  node             %massDistribution(massType=massTypeGalactic)
    massNotional      =  +hotHalo          %mass            (                         ) &
         &               +hotHalo          %outflowedMass   (                         ) &
         &               +massDistribution_%massTotal       (                         )
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    if (massNotional <= 0.0d0) return
    ! Compute the mean density and temperature of the hot halo.
    massDistribution_       => node             %massDistribution      (componentType=componentTypeHotHalo,massType=massTypeGaseous)
    kinematicsDistribution_ => massDistribution_%kinematicsDistribution(                                                           )
    select type (massDistribution_)
    type is (massDistributionZero)
       ! No mass distribution exists for the hot halo (most likely it is in an unphysical state).
       !![
       <objectDestructor name="massDistribution_"      />
       <objectDestructor name="kinematicsDistribution_"/>
       !!]
       return
    end select
    density    =+massNotional             &
         &      *3.0d0                    &
         &      /4.0d0                    &
         &      /Pi                       &
         &      /hotHalo%outerRadius()**3
    coordinates=[hotHalo%outerRadius(),0.0d0,0.0d0]
    temperature=+kinematicsDistribution_%temperature(coordinates)
    !![
    <objectDestructor name="massDistribution_"      />
    <objectDestructor name="kinematicsDistribution_"/>
    !!]
    ! Get the abundances for this node.
    abundances_=hotHalo%abundances()
    call abundances_%massToMassFraction(hotHalo%mass())
    ! Get the chemicals for this node.
    if (Chemicals_Property_Count () > 0) then
       chemicalMasses_=hotHalo%chemicals()
       ! Scale all chemical masses by their mass in atomic mass units to get a number density.
       call chemicalMasses_%massToNumber(chemicalDensities_)
       ! Compute factor converting mass of chemicals in (M☉/mₐ) to number density in cm⁻³.
       if (hotHalo%outerRadius() > 0.0d0) then
          massToDensityConversion=Chemicals_Mass_To_Density_Conversion(hotHalo%outerRadius())
       else
          massToDensityConversion=0.0d0
       end if
       ! Convert to number density.
       chemicalDensities_=chemicalDensities_*massToDensityConversion
    end if
    ! Set epoch for radiation field.
    call radiation%timeSet(basic%time())
    ! Compute number density of hydrogen (in cm⁻³).
    numberDensityHydrogen  =  Hydrogen_Number_Density             (            density                    ,abundances_                             )
    ! Get the number density of all species, including electrons.
    numberDensityAllSpecies=+                                            numberDensityHydrogen                                                       &
         &                  /abundances_    %hydrogenNumberFraction(                                                                               ) &
         &                  +chemicalState_ %electronDensity       (     numberDensityHydrogen,temperature,abundances_                   ,radiation)
    ! Get the cooling function (in ergs cm⁻³ s⁻¹).
    coolingFunction        =+coolingFunction_%coolingFunction      (node,numberDensityHydrogen,temperature,abundances_,chemicalDensities_,radiation)
    ! Compute the number of particles in the halo.
    countParticles         =+numberDensityAllSpecies &
         &                  *4.0d0                   &
         &                  *Pi                      &
         &                  /3.0d0                   &
         &                  *(                       &
         &                    +hotHalo%outerRadius() &
         &                    *megaParsec            &
         &                    /centi                 &
         &                   )**3
    exists=.true.
    return
  end function coolingBensonBower2010Properties

end module Cooling_Benson_Bower_2010
