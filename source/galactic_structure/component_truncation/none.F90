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
  Implementation of a galactic component truncation class which applies no truncation.
  !!}

  !![
  <galacticComponentTruncation name="galacticComponentTruncationNone" docformat="rst">
    <description>
    A galactic component truncation class which applies no truncation: the mass distributions of galactic components are
    always their untruncated mass distributions, scaled to the current mass of the component, irrespective of any tidal mass
    loss.
    </description>
  </galacticComponentTruncation>
  !!]
  type, extends(galacticComponentTruncationClass) :: galacticComponentTruncationNone
     !!{RST
     A galactic component truncation class which applies no truncation.
     !!}
     private
   contains
     procedure :: isActive         => noneIsActive
     procedure :: massDistribution => noneMassDistribution
  end type galacticComponentTruncationNone

  interface galacticComponentTruncationNone
     !!{RST
     Constructors for the :galacticus-class:`galacticComponentTruncationNone` galactic component truncation class.
     !!}
     module procedure noneConstructorParameters
  end interface galacticComponentTruncationNone

contains

  function noneConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`galacticComponentTruncationNone` galactic component truncation class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(galacticComponentTruncationNone)                :: self
    type(inputParameters                ), intent(inout) :: parameters

    self=galacticComponentTruncationNone()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function noneConstructorParameters

  logical function noneIsActive(self) result(isActive)
    !!{RST
    Return false, as this class applies no truncation.
    !!}
    implicit none
    class(galacticComponentTruncationNone), intent(inout) :: self
    !$GLC attributes unused :: self

    isActive=.false.
    return
  end function noneIsActive

  function noneMassDistribution(self,massDistribution_,fractionMassRetained) result(massDistributionTruncated)
    !!{RST
    Return the untruncated mass distribution.
    !!}
    implicit none
    class           (massDistributionSpherical      ), pointer               :: massDistributionTruncated
    class           (galacticComponentTruncationNone), intent(inout)         :: self
    class           (massDistributionSpherical      ), intent(inout), target :: massDistribution_
    double precision                                 , intent(in   )         :: fractionMassRetained
    !$GLC attributes unused :: self, fractionMassRetained

    massDistributionTruncated => massDistribution_
    !![
    <referenceCountIncrement object="massDistributionTruncated"/>
    !!]
    return
  end function noneMassDistribution
