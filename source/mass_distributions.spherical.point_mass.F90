!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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
  Implementation of a point mass distribution class.
  !!}

  !![
  <massDistribution name="massDistributionPointMass">
   <description>A mass distribution class for point masses.</description>
  </massDistribution>
  !!]
  type, public, extends(massDistributionSpherical) :: massDistributionPointMass
     !!{
     A point mass distribution.
     !!}
     double precision :: mass
   contains
     procedure :: density               => pointMassDensity
     procedure :: densityGradientRadial => pointMassDensityGradientRadial
     procedure :: densityRadialMoment   => pointMassDensityRadialMoment
     procedure :: massEnclosedBySphere  => pointMassMassEnclosedBySphere
     procedure :: potential             => pointMassPotential
  end type massDistributionPointMass

  interface massDistributionPointMass
     !!{
     Constructors for the {\normalfont \ttfamily pointMass} mass distribution class.
     !!}
     module procedure pointMassConstructorParameters
     module procedure pointMassConstructorInternal
  end interface massDistributionPointMass

contains

  function pointMassConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily pointMass} mass distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters          , only : inputParameters
    use :: Galactic_Structure_Options, only : enumerationComponentTypeEncode, enumerationMassTypeEncode
    implicit none
    type            (massDistributionPointMass)                :: self
    type            (inputParameters          ), intent(inout) :: parameters
    double precision                                           :: mass
    logical                                                    :: dimensionless
    type            (varying_string           )                :: componentType
    type            (varying_string           )                :: massType

    !![
    <inputParameter>
      <name>mass</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The mass of the point.</description>
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
     <call>self=massDistributionPointMass(componentType=enumerationComponentTypeEncode(componentType,includesPrefix=.false.),massType=enumerationMassTypeEncode(massType,includesPrefix=.false.){conditions})</call>
     <argument name="mass"          value="mass"          parameterPresent="parameters"/>
     <argument name="dimensionless" value="dimensionless" parameterPresent="parameters"/>
    </conditionalCall>
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function pointMassConstructorParameters
  
  function pointMassConstructorInternal(mass,dimensionless,componentType,massType) result(self)
    !!{
    Constructor for {\normalfont \ttfamily pointMass} mass distribution class.
    !!}
    use :: Error               , only : Error_Report
    use :: Numerical_Comparison, only : Values_Differ
    implicit none
    type            (massDistributionPointMass   )                          :: self
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
    return
  end function pointMassConstructorInternal

  double precision function pointMassDensity(self,coordinates,componentType,massType)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in a $\beta$-profile mass distribution.
    !!}
    implicit none
    class(massDistributionPointMass   ), intent(inout)           :: self
    class(coordinate                  ), intent(in   )           :: coordinates
    type (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type (enumerationMassTypeType     ), intent(in   ), optional :: massType
    
    pointMassDensity=     0.0d0
    if (.not.self       %matches          (componentType,massType)        ) return
    if (     coordinates%rSphericalSquared(                      ) > 0.0d0) return
    pointMassDensity=huge(0.0d0)
    return
  end function pointMassDensity

  double precision function pointMassDensityGradientRadial(self,coordinates,logarithmic,componentType,massType)
    !!{
    Return the density gradient in the radial direction for a point mass.
    !!}
    implicit none
    class  (massDistributionPointMass   ), intent(inout)           :: self
    class  (coordinate                  ), intent(in   )           :: coordinates
    logical                              , intent(in   ), optional :: logarithmic
    type   (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type   (enumerationMassTypeType     ), intent(in   ), optional :: massType
    !$GLC attributes unused :: logarithmic
    
    pointMassDensityGradientRadial=      0.0d0
    if (.not.self       %matches          (componentType,massType)        ) return
    if (     coordinates%rSphericalSquared(                      ) > 0.0d0) return
    pointMassDensityGradientRadial=-huge(0.0d0)
    return
  end function pointMassDensityGradientRadial

  double precision function pointMassMassEnclosedBySphere(self,radius,componentType,massType)
    !!{
    Computes the mass enclosed within a sphere of given {\normalfont \ttfamily radius} for a constant density cloud.
    !!}
    implicit none
    class           (massDistributionPointMass   ), intent(inout), target   :: self
    double precision                              , intent(in   )           :: radius
    type            (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional :: massType

    if (.not.self%matches(componentType,massType)) then
       pointMassMassEnclosedBySphere=0.0d0
       return
    end if
    pointMassMassEnclosedBySphere=+self%mass
    return
  end function pointMassMassEnclosedBySphere

  double precision function pointMassPotential(self,coordinates,componentType,massType)
    !!{
    Return the potential at the specified {\normalfont \ttfamily coordinates} for a point mass.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    implicit none
    class(massDistributionPointMass   ), intent(inout)           :: self
    class(coordinate                  ), intent(in   )           :: coordinates
    type (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type (enumerationMassTypeType     ), intent(in   ), optional :: massType

    if (.not.self%matches(componentType,massType)) then
       pointMassPotential=0.0d0
       return
    end if
    pointMassPotential=-self%mass                              &
         &             /coordinates%rSpherical()
    if (.not.self%dimensionless)                               &
         & pointMassPotential=+pointMassPotential              &
         &                    *gravitationalConstantGalacticus
    return
  end function pointMassPotential

  double precision function pointMassDensityRadialMoment(self,moment,radiusMinimum,radiusMaximum,isInfinite,componentType,massType)
    !!{
    Computes radial moments of the density for a point mass.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (massDistributionPointMass   ), intent(inout)           :: self
    double precision                              , intent(in   )           :: moment
    double precision                              , intent(in   ), optional :: radiusMinimum , radiusMaximum
    logical                                       , intent(  out), optional :: isInfinite
    type            (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional :: massType
    
    pointMassDensityRadialMoment=0.0d0
    ! Moment is zero if:
    if    (.not.self%matches(componentType,massType)        ) return ! Mass or component types do not match.
    if    (     moment                               > 0.0d0) return ! The moment is positive (in which case ∫ δ(r) rᵐ dr = 0).
    if (present(radiusMinimum)) then
       if (     radiusMinimum                        > 0.0d0) return ! The lower limit of the integral does not extend to zero.
    end if
    if (present(isInfinite)) then
       isInfinite=.true.
    else
       call Error_Report('radial moment is infinite'//{introspection:location})
    end if
    return    
  end function pointMassDensityRadialMoment
