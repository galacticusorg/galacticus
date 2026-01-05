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
  Implementation of a mass distribution class for an infinite radial extent gas slab with a Gaussian vertical distribution.
  !!}

  !![
  <massDistribution name="massDistributionGaussianSlab">
   <description>An infinite radial extent gas slab with a Gaussian vertical distribution: $\rho(r,z)=\rho_0 \exp(z^2/2 z^2_\mathrm{s})$.</description>
  </massDistribution>
  !!]
  type, public, extends(massDistributionCylindrical) :: massDistributionGaussianSlab
     !!{
     An infinite radial extent gas slab with a Gaussian vertical distribution: $\rho(r,z)=\rho_0 \exp(z^2/2 z^2_\mathrm{s})$.
     !!}
     private
     double precision :: scaleHeight, densityCentral
   contains
     procedure :: density                    => gaussianSlabDensity
     procedure :: densitySphericalAverage    => gaussianSlabDensitySphericalAverage
     procedure :: rotationCurve              => gaussianSlabRotationCurve
     procedure :: rotationCurveGradient      => gaussianSlabRotationCurveGradient
     procedure :: surfaceDensity             => gaussianSlabSurfaceDensity
     procedure :: surfaceDensityRadialMoment => gaussianSlabSurfaceDensityRadialMoment
     procedure :: radiusHalfMass             => gaussianSlabRadiusHalfMass
   end type massDistributionGaussianSlab

  interface massDistributionGaussianSlab
     !!{
     Constructors for the \refClass{massDistributionGaussianSlab} mass distribution class.
     !!}
     module procedure gaussianSlabConstructorParameters
     module procedure gaussianSlabConstructorInternal
  end interface massDistributionGaussianSlab

contains

  function gaussianSlabConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{massDistributionGaussianSlab} mass distribution class which builds the object from a parameter
    set.
    !!}
     use :: Galactic_Structure_Options, only : enumerationComponentTypeEncode, enumerationMassTypeEncode
     use :: Input_Parameters          , only : inputParameter                , inputParameters
    implicit none
    type            (massDistributionGaussianSlab)                :: self
    type            (inputParameters             ), intent(inout) :: parameters
    double precision                                              :: scaleHeight   , densityCentral
    logical                                                       :: dimensionless
    type            (varying_string              )                :: componentType
    type            (varying_string              )                :: massType

    !![
    <inputParameter>
      <name>scaleHeight</name>
      <defaultSource>\citep{kregel_flattening_2002}</defaultSource>
      <defaultValue>0.137d0</defaultValue>
      <description>The scale height of the Gaussian slab profile.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>densityCentral</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The density in the slab mid-plane.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>dimensionless</name>
      <defaultValue>.true.</defaultValue>
      <description>If true the Gaussian slab profile is considered to be dimensionless.</description>
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
    !!]
    self=massDistributionGaussianSlab(scaleHeight,densityCentral,dimensionless,enumerationComponentTypeEncode(componentType,includesPrefix=.false.),enumerationMassTypeEncode(massType,includesPrefix=.false.))
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function gaussianSlabConstructorParameters

  function gaussianSlabConstructorInternal(scaleHeight,densityCentral,dimensionless,componentType,massType) result(self)
    !!{
    Internal constructor for ``gaussianSlab'' mass distribution class.
    !!}
    use :: Error                   , only : Error_Report
    use :: Numerical_Comparison    , only : Values_Differ
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    type            (massDistributionGaussianSlab)                          :: self
    double precision                              , intent(in   ), optional :: scaleHeight  , densityCentral
    logical                                       , intent(in   ), optional :: dimensionless
    type            (enumerationComponentTypeType), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType     ), intent(in   ), optional :: massType
    !![
    <constructorAssign variables="componentType, massType"/>
    !!]
    
    ! Determine if profile is dimensionless.
    if (present(dimensionless)) then
       self%dimensionless=dimensionless
    else
       self%dimensionless=.false.
    end if
    ! If dimensionless, then set scale length and mass to unity.
    if (self%dimensionless) then
       if (present(scaleHeight   )) then
          if (Values_Differ(scaleHeight   ,1.0d0,absTol=1.0d-6)) call Error_Report('scale height should be unity for a dimensionless profile (or simply do not specify a scale height)'//{introspection:location})
       end if
       if (present(densityCentral)) then
          if (Values_Differ(densityCentral,1.0d0,absTol=1.0d-6)) call Error_Report('central density should be unity for a dimensionless profile (or simply do not specify a density)'  //{introspection:location})
       end if
       self%scaleHeight   =1.0d0
       self%densityCentral=1.0d0
    else
       ! Set properties radius.
       if (.not.present(scaleHeight   )) call Error_Report('scale height must be specified for dimensionful profiles'   //{introspection:location})
       if (.not.present(densityCentral)) call Error_Report('central density must be specified for dimensionful profiles'//{introspection:location})
       self%scaleHeight   =scaleHeight
       self%densityCentral=densityCentral
    end if
    return
  end function gaussianSlabConstructorInternal

  double precision function gaussianSlabDensity(self,coordinates)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in a Gaussian slab mass distribution.
    !!}
    use :: Coordinates, only : assignment(=), coordinateCylindrical
    use :: Error      , only : Error_Report
    implicit none
    class           (massDistributionGaussianSlab), intent(inout) :: self
    class           (coordinate                  ), intent(in   ) :: coordinates
    type            (coordinateCylindrical       )                :: position
    double precision                                              :: z

    ! If disk is razor thin, density is undefined.
    if (self%scaleHeight <= 0.0d0) call Error_Report('density undefined for razor-thin slab'//{introspection:location})
    ! Get position in cylindrical coordinate system.
    position=coordinates
    ! Compute density.
    z=abs(position%z())/self%scaleHeight
    gaussianSlabDensity=self%densityCentral*exp(-0.5d0*z**2)
    return
  end function gaussianSlabDensity

  double precision function gaussianSlabDensitySphericalAverage(self,radius)
    !!{
    Return the spherically-averaged density at the specified {\normalfont \ttfamily radius} in a Gaussian slab mass distribution.
    !!}
    implicit none
    class           (massDistributionGaussianSlab), intent(inout) :: self
    double precision                              , intent(in   ) :: radius
    !$GLC attributes unused :: self, radius

    gaussianSlabDensitySphericalAverage=0.0d0
    call Error_Report('spherically-averaged density profile is not implemented'//{introspection:location})
    return
  end function gaussianSlabDensitySphericalAverage

  double precision function gaussianSlabRotationCurve(self,radius)
    !!{
    Rotation curve for a infinite extent Gaussian slab.
    !!}
    implicit none
    class           (massDistributionGaussianSlab), intent(inout) :: self
    double precision                              , intent(in   ) :: radius
    !$GLC attributes unused :: self, radius
    
    gaussianSlabRotationCurve=0.0d0
    return
  end function gaussianSlabRotationCurve

  double precision function gaussianSlabRotationCurveGradient(self,radius)
    !!{
    Rotation curve gradient for a infinite extent Gaussian slab.
    !!}
    implicit none
    class           (massDistributionGaussianSlab), intent(inout) :: self
    double precision                              , intent(in   ) :: radius
    !$GLC attributes unused :: self, radius
    
    gaussianSlabRotationCurveGradient=0.0d0
    return
  end function gaussianSlabRotationCurveGradient

  double precision function gaussianSlabSurfaceDensity(self,coordinates)
    !!{
    Return the surface density at the specified {\normalfont \ttfamily coordinates} in a Gaussian slab mass distribution.
    !!}
    use :: Coordinates             , only : coordinate
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class(massDistributionGaussianSlab), intent(inout) :: self
    class(coordinate                  ), intent(in   ) :: coordinates
    !$GLC attributes unused :: coordinates
    
    gaussianSlabSurfaceDensity=+sqrt(2.0d0*Pi)      &
         &                     *self%scaleHeight    &
         &                     *self%densityCentral
    return
  end function gaussianSlabSurfaceDensity

  double precision function gaussianSlabSurfaceDensityRadialMoment(self,moment,radiusMinimum,radiusMaximum,isInfinite)
    !!{
    Compute radial moments of the Gaussian slab mass distribution surface density profile.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (massDistributionGaussianSlab), intent(inout)           :: self
    double precision                              , intent(in   )           :: moment
    double precision                              , intent(in   ), optional :: radiusMinimum, radiusMaximum
    logical                                       , intent(  out), optional :: isInfinite
    !$GLC attributes unused :: self, moment, radiusMinimum, radiusMaximum
    
    ! All moments are infinite.
    gaussianSlabSurfaceDensityRadialMoment=huge(0.0d0)
    if (present(isInfinite)) then
       isInfinite=.true.
    else
       call Error_Report('moment is infinite'//{introspection:location})
    end if
    return
  end function gaussianSlabSurfaceDensityRadialMoment

  double precision function gaussianSlabRadiusHalfMass(self)
    !!{
    Return the half-mass radius for an infinite extent Gaussian slab mass distribution.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(massDistributionGaussianSlab), intent(inout) :: self
    !$GLC attributes unused :: self

    gaussianSlabRadiusHalfMass=0.0d0
    call Error_Report('half mass radius is undefined'//{introspection:location})
    return
  end function gaussianSlabRadiusHalfMass
