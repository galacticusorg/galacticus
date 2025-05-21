!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
  Implementation of a zero mass distribution class.
  !!}

  !![
  <massDistribution name="massDistributionZero">
   <description>A zero mass distribution class.</description>
  </massDistribution>
  !!]
  type, public, extends(massDistributionSpherical) :: massDistributionZero
     !!{
     A zero mass distribution.
     !!}
   contains
     procedure :: massTotal             => zeroMassTotal
     procedure :: density               => zeroDensity
     procedure :: densityGradientRadial => zeroDensityGradientRadial
     procedure :: densityRadialMoment   => zeroDensityRadialMoment
     procedure :: massEnclosedBySphere  => zeroMassEnclosedBySphere
     procedure :: rotationCurve         => zeroRotationCurve
     procedure :: rotationCurveGradient => zeroRotationCurveGradient
     procedure :: potentialIsAnalytic   => zeroPotentialIsAnalytic
     procedure :: potential             => zeroPotential
  end type massDistributionZero

  interface massDistributionZero
     !!{
     Constructors for the \refClass{massDistributionZero} mass distribution class.
     !!}
     module procedure zeroConstructorParameters
     module procedure zeroConstructorInternal
  end interface massDistributionZero

contains

  function zeroConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{massDistributionZero} mass distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (massDistributionZero)                :: self
    type   (inputParameters     ), intent(inout) :: parameters
    logical                                      :: dimensionless

    !![
    <inputParameter>
     <name>dimensionless</name>
     <defaultValue>.true.</defaultValue>
     <description>If true the null profile is considered to be dimensionless.</description>
     <source>parameters</source>
    </inputParameter>
    !!]
    self=massDistributionZero(dimensionless)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function zeroConstructorParameters
  
  function zeroConstructorInternal(dimensionless) result(self)
    !!{
    Constructor for the \refClass{massDistributionZero} mass distribution class.
    !!}
    use :: Galactic_Structure_Options, only : componentTypeUnknown, massTypeUnknown
    implicit none
    type   (massDistributionZero)                :: self
    logical                      , intent(in   ) :: dimensionless
    !![
    <constructorAssign variables="dimensionless"/>
    !!]
    
    self%componentType=componentTypeUnknown
    self%massType     =massTypeUnknown
    return
  end function zeroConstructorInternal

  double precision function zeroMassTotal(self)
    !!{
    Return the total mass in the zero distribution.
    !!}
    implicit none
    class(massDistributionZero), intent(inout) :: self
    !$GLC attributes unused :: self

    zeroMassTotal=0.0d0
    return
  end function zeroMassTotal

  double precision function zeroDensity(self,coordinates)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in a zero distribution.
    !!}
    implicit none
    class(massDistributionZero), intent(inout) :: self
    class(coordinate          ), intent(in   ) :: coordinates
    !$GLC attributes unused :: self, coordinates

    zeroDensity=0.0d0
    return
  end function zeroDensity

  double precision function zeroDensityGradientRadial(self,coordinates,logarithmic)
    !!{
    Return the density gradient in the radial direction for a point mass.
    !!}
    implicit none
    class  (massDistributionZero), intent(inout), target   :: self
    class  (coordinate          ), intent(in   )           :: coordinates
    logical                      , intent(in   ), optional :: logarithmic
    !$GLC attributes unused :: self, coordinates, logarithmic
    
    zeroDensityGradientRadial=0.0d0
    return
  end function zeroDensityGradientRadial

  double precision function zeroMassEnclosedBySphere(self,radius)
    !!{
    Computes the mass enclosed within a sphere of given {\normalfont \ttfamily radius} for a zero distribution.
    !!}
    implicit none
    class           (massDistributionZero), intent(inout), target :: self
    double precision                      , intent(in   )         :: radius
     !$GLC attributes unused :: self, radius

    zeroMassEnclosedBySphere=0.0d0
    return
  end function zeroMassEnclosedBySphere

  double precision function zeroRotationCurve(self,radius)
    !!{
    Return the rotation curve for a zero mass distribution.
    !!}
    implicit none
    class           (massDistributionZero), intent(inout) :: self
    double precision                      , intent(in   ) :: radius
    !$GLC attributes unused :: self, radius

    zeroRotationCurve=0.0d0
    return
  end function zeroRotationCurve

  double precision function zeroRotationCurveGradient(self,radius)
    !!{
    Return the rotation curve gradient for a spherical mass distribution.
    !!}
    implicit none
    class           (massDistributionZero), intent(inout) :: self
    double precision                      , intent(in   ) :: radius
    !$GLC attributes unused :: self, radius

    zeroRotationCurveGradient=0.0d0
    return
  end function zeroRotationCurveGradient

  logical function zeroPotentialIsAnalytic(self) result(isAnalytic)
    !!{
    Return that the potential has an analytic form.
    !!}
    implicit none
    class(massDistributionZero), intent(inout) :: self

    isAnalytic=.true.
    return
  end function zeroPotentialIsAnalytic

  double precision function zeroPotential(self,coordinates,status)
    !!{
    Return the potential at the specified {\normalfont \ttfamily coordinates} for a point mass.
    !!}
    use :: Galactic_Structure_Options, only : structureErrorCodeSuccess
    implicit none
    class(massDistributionZero             ), intent(inout), target   :: self
    class(coordinate                       ), intent(in   )           :: coordinates
    type (enumerationStructureErrorCodeType), intent(  out), optional :: status
    !$GLC attributes unused :: self, coordinates

    if (present(status)) status=structureErrorCodeSuccess
    zeroPotential=0.0d0
    return
  end function zeroPotential

  double precision function zeroDensityRadialMoment(self,moment,radiusMinimum,radiusMaximum,isInfinite)
    !!{
    Computes radial moments of the density for a point mass.
    !!}
    implicit none
    class           (massDistributionZero), intent(inout)           :: self
    double precision                      , intent(in   )           :: moment
    double precision                      , intent(in   ), optional :: radiusMinimum, radiusMaximum
    logical                               , intent(  out), optional :: isInfinite
    !$GLC attributes unused :: self, moment, radiusMinimum, radiusMaximum

    zeroDensityRadialMoment=0.0d0
    if (present(isInfinite)) isInfinite=.false.
    return    
  end function zeroDensityRadialMoment
