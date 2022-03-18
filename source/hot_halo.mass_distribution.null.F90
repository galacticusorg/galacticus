!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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
A null implementation of the hot halo mass distribution class.
!!}

  !![
  <hotHaloMassDistribution name="hotHaloMassDistributionNull">
   <description>
    A hot halo mass distribution class that assumes no hot halo mass distribution. It is useful, for example, when performing
    dark matter-only calculations.
   </description>
  </hotHaloMassDistribution>
  !!]
  type, extends(hotHaloMassDistributionClass) :: hotHaloMassDistributionNull
     !!{
     A null implementation of the hot halo mass distribution class.
     !!}
     private
   contains
     procedure :: density               => nullDensity
     procedure :: densityLogSlope       => nullDensityLogSlope
     procedure :: enclosedMass          => nullEnclosedMass
     procedure :: radialMoment          => nullRadialMoment
     procedure :: rotationNormalization => nullRotationNormalization
  end type hotHaloMassDistributionNull

  interface hotHaloMassDistributionNull
     !!{
     Constructors for the null hot halo mass distribution class.
     !!}
     module procedure nullConstructorParameters
  end interface hotHaloMassDistributionNull

contains

  function nullConstructorParameters(parameters) result(self)
    !!{
    Constructor for the null hot halo mass distribution class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(hotHaloMassDistributionNull)                :: self
    type(inputParameters            ), intent(inout) :: parameters

    self=hotHaloMassDistributionNull()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function nullConstructorParameters

  double precision function nullDensity(self,node,radius)
    !!{
    Return the density in a null hot halo mass distribution.
    !!}
    implicit none
    class           (hotHaloMassDistributionNull), intent(inout) :: self
    type            (treeNode                   ), intent(inout) :: node
    double precision                             , intent(in   ) :: radius
    !$GLC attributes unused :: self, node, radius

    nullDensity=0.0d0
    return
  end function nullDensity

  double precision function nullDensityLogSlope(self,node,radius)
    !!{
    Return the logarithmic slope of the density of the hot halo at the given {\normalfont \ttfamily radius}.
    !!}
    implicit none
    class           (hotHaloMassDistributionNull), intent(inout) :: self
    type            (treeNode                   ), intent(inout) :: node
    double precision                             , intent(in   ) :: radius
    !$GLC attributes unused :: self, node, radius

    nullDensityLogSlope=0.0d0
    return
  end function nullDensityLogSlope

  double precision function nullEnclosedMass(self,node,radius)
    !!{
    Return the mass enclosed in the hot halo at the given {\normalfont \ttfamily radius}.
    !!}
    implicit none
    class           (hotHaloMassDistributionNull), intent(inout)         :: self
    type            (treeNode                   ), intent(inout), target :: node
    double precision                             , intent(in   )         :: radius
    !$GLC attributes unused :: self, node, radius

    nullEnclosedMass=0.0d0
    return
  end function nullEnclosedMass

  double precision function nullRadialMoment(self,node,moment,radius)
    !!{
    Return the density of the hot halo at the given {\normalfont \ttfamily radius}.
    !!}
    implicit none
    class           (hotHaloMassDistributionNull), intent(inout) :: self
    type            (treeNode                   ), intent(inout) :: node
    double precision                             , intent(in   ) :: moment, radius
    !$GLC attributes unused :: self, node, radius, moment

    nullRadialMoment=0.0d0
    return
  end function nullRadialMoment

  double precision function nullRotationNormalization(self,node)
    !!{
    Returns the relation between specific angular momentum and rotation velocity (assuming a rotation velocity that is constant
    in radius) for {\normalfont \ttfamily node}. Specifically, the normalization, $A$, returned is such that $V_\mathrm{rot} =
    A J/M$.
    !!}
    implicit none
    class(hotHaloMassDistributionNull), intent(inout) :: self
    type (treeNode                   ), intent(inout) :: node
    !$GLC attributes unused :: self, node

    nullRotationNormalization=0.0d0
  return
  end function nullRotationNormalization

