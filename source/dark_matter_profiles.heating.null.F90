!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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

  !% A null dark matter halo profile heating class.

  !# <darkMatterProfileHeating name="darkMatterProfileHeatingNull">
  !#  <description>A dark matter profile heating model in which the heating is always zero.</description>
  !# </darkMatterProfileHeating>

  type, extends(darkMatterProfileHeatingClass) :: darkMatterProfileHeatingNull
     !% A dark matter profile heating class with zero heating.
     private
   contains
     procedure :: specificEnergy                 => nullSpecificEnergy
     procedure :: specificEnergyGradient         => nullSpecificEnergyGradient
     procedure :: specificEnergyIsEverywhereZero => nullSpecificEnergyIsEverywhereZero
  end type darkMatterProfileHeatingNull

  interface darkMatterProfileHeatingNull
     !% Constructors for the {\normalfont \ttfamily null} dark matter profile heating class.
     module procedure nullConstructorParameters
  end interface darkMatterProfileHeatingNull

contains

  function nullConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily null} dark matter profile heating scales class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(darkMatterProfileHeatingNull), target        :: self
    type(inputParameters             ), intent(inout) :: parameters
    !$GLC attributes unused :: parameters

    self=darkMatterProfileHeatingNull()
    return
  end function nullConstructorParameters

  double precision function nullSpecificEnergy(self,node,darkMatterProfileDMO_,radius)
    !% Returns the specific energy of heating in the given {\normalfont \ttfamily node}.
    implicit none
    class           (darkMatterProfileHeatingNull), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    class           (darkMatterProfileDMOClass      ), intent(inout) :: darkMatterProfileDMO_
    double precision                              , intent(in   ) :: radius
    !$GLC attributes unused :: self, node, radius, darkMatterProfileDMO_

    nullSpecificEnergy=0.0d0
    return
  end function nullSpecificEnergy

  double precision function nullSpecificEnergyGradient(self,node,darkMatterProfileDMO_,radius)
    !% Returns the gradient of the specific energy of heating in the given {\normalfont \ttfamily node}.
    implicit none
    class           (darkMatterProfileHeatingNull), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    class           (darkMatterProfileDMOClass      ), intent(inout) :: darkMatterProfileDMO_
    double precision                              , intent(in   ) :: radius
    !$GLC attributes unused :: self, node, darkMatterProfileDMO_, radius

    nullSpecificEnergyGradient=0.0d0
    return
  end function nullSpecificEnergyGradient

  logical function nullSpecificEnergyIsEverywhereZero(self,node,darkMatterProfileDMO_)
    !% Returns true if the specific energy is everywhere zero in the given {\normalfont \ttfamily node}.
    implicit none
    class(darkMatterProfileHeatingNull), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node
    class(darkMatterProfileDMOClass      ), intent(inout) :: darkMatterProfileDMO_
    !$GLC attributes unused :: self, node, darkMatterProfileDMO_

    nullSpecificEnergyIsEverywhereZero=.true.
    return
  end function nullSpecificEnergyIsEverywhereZero
