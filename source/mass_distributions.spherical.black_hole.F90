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
  Implementation of a black hole distribution class.
  !!}

  !![
  <massDistribution name="massDistributionBlackHole">
   <description>A mass distribution class for point masses.</description>
  </massDistribution>
  !!]
  type, public, extends(massDistributionSpherical) :: massDistributionBlackHole
     !!{
     A black hole distribution.
     !!}
     double precision :: mass, radiusGravitational
   contains
     procedure :: massTotal             => blackHoleMassTotal
     procedure :: density               => blackHoleDensity
     procedure :: densityGradientRadial => blackHoleDensityGradientRadial
     procedure :: densityRadialMoment   => blackHoleDensityRadialMoment
     procedure :: massEnclosedBySphere  => blackHoleMassEnclosedBySphere
     procedure :: rotationCurve         => blackHoleRotationCurve
     procedure :: rotationCurveGradient => blackHoleRotationCurveGradient
     procedure :: potentialIsAnalytic   => blackHolePotentialIsAnalytic
     procedure :: potential             => blackHolePotential
  end type massDistributionBlackHole

  interface massDistributionBlackHole
     !!{
     Constructors for the \refClass{massDistributionBlackHole} mass distribution class.
     !!}
     module procedure blackHoleConstructorParameters
     module procedure blackHoleConstructorInternal
  end interface massDistributionBlackHole

contains

  function blackHoleConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{massDistributionBlackHole} mass distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters          , only : inputParameters
    use :: Galactic_Structure_Options, only : enumerationComponentTypeEncode, enumerationMassTypeEncode
    implicit none
    type            (massDistributionBlackHole)                :: self
    type            (inputParameters          ), intent(inout) :: parameters
    double precision                                           :: mass
    logical                                                    :: dimensionless
    type            (varying_string           )                :: componentType
    type            (varying_string           )                :: massType

    !![
    <inputParameter>
      <name>mass</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The mass of the black hole.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>dimensionless</name>
      <defaultValue>.true.</defaultValue>
      <description>If true the point mass distribution is considered to be dimensionless.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>componentType</name>
      <defaultValue>var_str('unknown')</defaultValue>
      <description>The component type that this mass distribution represents.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massType</name>
      <defaultValue>var_str('unknown')</defaultValue>
      <description>The mass type that this mass distribution represents.</description>
      <source>parameters</source>
    </inputParameter>
    <conditionalCall>
     <call>self=massDistributionBlackHole(componentType=enumerationComponentTypeEncode(componentType,includesPrefix=.false.),massType=enumerationMassTypeEncode(massType,includesPrefix=.false.){conditions})</call>
     <argument name="mass"          value="mass"          parameterPresent="parameters"/>
     <argument name="dimensionless" value="dimensionless" parameterPresent="parameters"/>
    </conditionalCall>
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function blackHoleConstructorParameters
  
  function blackHoleConstructorInternal(mass,dimensionless,componentType,massType) result(self)
    !!{
    Constructor for the \refClass{massDistributionBlackHole} mass distribution class.
    !!}
    use :: Error                           , only : Error_Report
    use :: Numerical_Comparison            , only : Values_Differ
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    use :: Numerical_Constants_Physical    , only : speedLight
    use :: Numerical_Constants_Prefixes    , only : milli
    implicit none
    type            (massDistributionBlackHole   )                          :: self
    double precision                              , intent(in   ), optional :: mass
    logical                                       , intent(in   ), optional :: dimensionless
    type            (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional :: massType
    !![
    <constructorAssign variables="componentType, massType"/>
    !!]
    ! Determine if profile is dimensionless.
    self%dimensionless=.false.
    if (present(dimensionless)) self%dimensionless=dimensionless
    ! If dimensionless, then set scale length and mass to unity.
    if (self%dimensionless) then
       if (present(mass)) then
          if (Values_Differ(mass,1.0d0,absTol=1.0d-6)) call Error_Report('mass should be unity for a dimensionless profile (or simply do not specify a mass)'//{introspection:location})
       end if
       self%mass         =1.0d0
    else if (present(mass)) then
       self%mass         =mass
       self%dimensionless=.false.
    else
       call Error_Report('either specify a mass, or declare the distribution to be dimensionless'//{introspection:location})
    end if
    ! Compute the gravitational radius for the black hole.
    if (self%dimensionless) then
       self%radiusGravitational=+1.0d0
    else
       self%radiusGravitational=+     gravitationalConstant_internal     &
            &                   *self%mass                               &
            &                   /     milli                          **2 &
            &                   /     speedLight                     **2
    end if
    return
  end function blackHoleConstructorInternal

  double precision function blackHoleMassTotal(self)
    !!{
    Return the total mass in the black hole.
    !!}
    implicit none
    class(massDistributionBlackHole), intent(inout) :: self
 
    blackHoleMassTotal=self%mass
    return
  end function blackHoleMassTotal

  double precision function blackHoleDensity(self,coordinates)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in a $\beta$-profile mass distribution.
    !!}
    implicit none
    class(massDistributionBlackHole), intent(inout) :: self
    class(coordinate               ), intent(in   ) :: coordinates
    
    blackHoleDensity=     0.0d0
    if (coordinates%rSphericalSquared() > 0.0d0) return
    blackHoleDensity=huge(0.0d0)
    return
  end function blackHoleDensity

  double precision function blackHoleDensityGradientRadial(self,coordinates,logarithmic)
    !!{
    Return the density gradient in the radial direction for a point mass.
    !!}
    implicit none
    class  (massDistributionBlackHole), intent(inout), target   :: self
    class  (coordinate               ), intent(in   )           :: coordinates
    logical                           , intent(in   ), optional :: logarithmic
    !$GLC attributes unused :: logarithmic
    
    blackHoleDensityGradientRadial=      0.0d0
    if (coordinates%rSphericalSquared() > 0.0d0) return
    blackHoleDensityGradientRadial=-huge(0.0d0)
    return
  end function blackHoleDensityGradientRadial

  double precision function blackHoleMassEnclosedBySphere(self,radius)
    !!{
    Computes the mass enclosed within a sphere of given {\normalfont \ttfamily radius} for a black hole.
    !!}
    implicit none
    class           (massDistributionBlackHole), intent(inout), target :: self
    double precision                           , intent(in   )         :: radius
 
    blackHoleMassEnclosedBySphere=+self%mass
    return
  end function blackHoleMassEnclosedBySphere

  double precision function blackHoleRotationCurve(self,radius)
    !!{
    Return the rotation curve for a blackHole mass distribution.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    use :: Numerical_Constants_Physical    , only : speedLight
    use :: Numerical_Constants_Prefixes    , only : milli
    implicit none
    class           (massDistributionBlackHole), intent(inout) :: self
    double precision                           , intent(in   ) :: radius

    if (self%mass <= 0.0d0) then
       blackHoleRotationCurve=0.0d0
    else if (radius <= self%radiusGravitational) then
       if (self%dimensionless) then
          blackHoleRotationCurve=+1.0d0
       else
          blackHoleRotationCurve=+milli      &
               &                 *speedLight
       end if
    else
       blackHoleRotationCurve=+sqrt(                                   &
            &                       +self%massEnclosedBySphere(radius) &
            &                       /                          radius  &
            &                      )
       ! Make dimensionful if necessary.
       if (.not.self%dimensionless) blackHoleRotationCurve=+sqrt(gravitationalConstant_internal) &
            &                                              *blackHoleRotationCurve
    end if
    return
  end function blackHoleRotationCurve

  double precision function blackHoleRotationCurveGradient(self,radius)
    !!{
    Return the rotation curve gradient for a spherical mass distribution.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    use :: Numerical_Constants_Math        , only : Pi
    implicit none
    class           (massDistributionBlackHole), intent(inout) :: self
    double precision                           , intent(in   ) :: radius
 
    if     (                                         &
         &        self%mass                <= 0.0d0  &
         &  .or.                                     &
         &        self%radiusGravitational >  radius &
         &) then
       blackHoleRotationCurveGradient=+0.0d0
    else
       blackHoleRotationCurveGradient=-self%mass      &
            &                         /     radius**2
       ! Make dimensionful if necessary.
       if (.not.self%dimensionless) blackHoleRotationCurveGradient=+gravitationalConstant_internal &
            &                                                      *blackHoleRotationCurveGradient
    end if
    return
  end function blackHoleRotationCurveGradient

  logical function blackHolePotentialIsAnalytic(self) result(isAnalytic)
    !!{
    Return that the potential has an analytic form.
    !!}
    implicit none
    class(massDistributionBlackHole), intent(inout) :: self

    isAnalytic=.true.
    return
  end function blackHolePotentialIsAnalytic

  double precision function blackHolePotential(self,coordinates,status)
    !!{
    Return the potential at the specified {\normalfont \ttfamily coordinates} for a point mass.
    !!}
    use :: Galactic_Structure_Options      , only : structureErrorCodeSuccess
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    class(massDistributionBlackHole        ), intent(inout), target   :: self
    class(coordinate                       ), intent(in   )           :: coordinates
    type (enumerationStructureErrorCodeType), intent(  out), optional :: status

    if (present(status)) status=structureErrorCodeSuccess
    blackHolePotential=-     self       %mass                   &
         &             /max(                                    &
         &                   coordinates%rSpherical         (), &
         &                   self       %radiusGravitational    &
         &                  )
    if (.not.self%dimensionless)                               &
         & blackHolePotential=+blackHolePotential            &
         &                    *gravitationalConstant_internal
    return
  end function blackHolePotential

  double precision function blackHoleDensityRadialMoment(self,moment,radiusMinimum,radiusMaximum,isInfinite)
    !!{
    Computes radial moments of the density for a point mass.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (massDistributionBlackHole), intent(inout)           :: self
    double precision                           , intent(in   )           :: moment
    double precision                           , intent(in   ), optional :: radiusMinimum , radiusMaximum
    logical                                    , intent(  out), optional :: isInfinite
     
    blackHoleDensityRadialMoment=0.0d0
    ! Moment is zero if:
    if    (moment        > 0.0d0) return ! The moment is positive (in which case ∫ δ(r) rᵐ dr = 0).
    if (present(radiusMinimum)) then
       if (radiusMinimum > 0.0d0) return ! The lower limit of the integral does not extend to zero.
    end if
    if (present(isInfinite)) then
       isInfinite=.true.
    else
       call Error_Report('radial moment is infinite'//{introspection:location})
    end if
    return    
  end function blackHoleDensityRadialMoment
