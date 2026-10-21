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
  Implementation of a galactic component heating class which applies no heating.
  !!}

  !![
  <galacticComponentHeating name="galacticComponentHeatingNone" docformat="rst">
    <description>
    A galactic component heating class which applies no heating: the mass distributions of galactic components are unaffected
    by tidal heating.
    </description>
  </galacticComponentHeating>
  !!]
  type, extends(galacticComponentHeatingClass) :: galacticComponentHeatingNone
     !!{RST
     A galactic component heating class which applies no heating.
     !!}
     private
   contains
     procedure :: isActive         => noneIsActive
     procedure :: massDistribution => noneMassDistribution
  end type galacticComponentHeatingNone

  interface galacticComponentHeatingNone
     !!{RST
     Constructors for the :galacticus-class:`galacticComponentHeatingNone` galactic component heating class.
     !!}
     module procedure noneConstructorParameters
  end interface galacticComponentHeatingNone

contains

  function noneConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`galacticComponentHeatingNone` galactic component heating class which takes a parameter
    set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(galacticComponentHeatingNone)                :: self
    type(inputParameters             ), intent(inout) :: parameters

    self=galacticComponentHeatingNone()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function noneConstructorParameters

  logical function noneIsActive(self) result(isActive)
    !!{RST
    Return false, as this class applies no heating.
    !!}
    implicit none
    class(galacticComponentHeatingNone), intent(inout) :: self
    !$GLC attributes unused :: self

    isActive=.false.
    return
  end function noneIsActive

  function noneMassDistribution(self,massDistribution_,massDistributionDarkMatter_) result(massDistributionHeated)
    !!{RST
    Return the unheated mass distribution.
    !!}
    implicit none
    class(massDistributionSpherical   ), pointer               :: massDistributionHeated
    class(galacticComponentHeatingNone), intent(inout)         :: self
    class(massDistributionSpherical   ), intent(inout), target :: massDistribution_     , massDistributionDarkMatter_
    !$GLC attributes unused :: self, massDistributionDarkMatter_

    massDistributionHeated => massDistribution_
    !![
    <referenceCountIncrement object="massDistributionHeated"/>
    !!]
    return
  end function noneMassDistribution
