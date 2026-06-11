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

  !!{RST
  Implements a null mass distribution heating class.
  !!}

  !![
  <massDistributionHeating name="massDistributionHeatingNull" docformat="rst">
    <description>
    A null mass distribution heating class. The heating energy is always zero.
    </description>
  </massDistributionHeating>
  !!]
  type, extends(massDistributionHeatingClass) :: massDistributionHeatingNull
     !!{RST
     Implementation of a null mass distribution heating class.
     !!}
     private
   contains
     procedure :: specificEnergy                 => nullSpecificEnergy
     procedure :: specificEnergyGradient         => nullSpecificEnergyGradient
     procedure :: specificEnergyIsEveryWhereZero => nullSpecificEnergyIsEverywhereZero
  end type massDistributionHeatingNull

  interface massDistributionHeatingNull
     !!{RST
     Constructors for the ``massDistributionHeatingNull`` mass distribution heating class.
     !!}
     module procedure nullConstructorParameters
  end interface massDistributionHeatingNull

contains

  function nullConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the ``massDistributionHeatingNull`` mass distribution heating class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(massDistributionHeatingNull)                :: self
    type(inputParameters            ), intent(inout) :: parameters

    self=massDistributionHeatingNull()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function nullConstructorParameters
  
  double precision function nullSpecificEnergy(self,radius,massDistribution_) result(energySpecific)
    !!{RST
    Compute the specific energy in a zero-heating mass distribution.
    !!}
    implicit none
    class           (massDistributionHeatingNull), intent(inout) :: self
    double precision                             , intent(in   ) :: radius
    class           (massDistributionClass      ), intent(inout) :: massDistribution_

    energySpecific=+0.0d0
    return
  end function nullSpecificEnergy

  double precision function nullSpecificEnergyGradient(self,radius,massDistribution_) result(energySpecificGradient)
    !!{RST
    Returns the gradient of the specific energy of heating.
    !!}
    implicit none
    class           (massDistributionHeatingNull), intent(inout) :: self
    double precision                             , intent(in   ) :: radius
    class           (massDistributionClass      ), intent(inout) :: massDistribution_

    energySpecificGradient=+0.0d0
    return
  end function nullSpecificEnergyGradient

  logical function nullSpecificEnergyIsEverywhereZero(self) result(energySpecificIsEverywhereZero)
    !!{RST
    Returns true if the specific energy is everywhere zero.
    !!}
    implicit none
    class(massDistributionHeatingNull), intent(inout) :: self

    energySpecificIsEverywhereZero=.true.
    return
  end function nullSpecificEnergyIsEverywhereZero
