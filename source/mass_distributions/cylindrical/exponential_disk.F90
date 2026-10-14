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
  Implementation of an exponential disk mass distribution class.
  !!}
  
  use :: Tables, only : table1DLogarithmicLinear

  !![
  <massDistribution name="massDistributionExponentialDisk" docformat="rst">
   <description>
   The exponential disk mass distribution: :math:`\rho(r,z)=\rho_0 \exp(-r/r_\mathrm{s}) \hbox{sech}^2(z/z_\mathrm{s})`.
   </description>
  </massDistribution>
  !!]
  type, public, extends(massDistributionCylindrical) :: massDistributionExponentialDisk
     !!{RST
     The exponential disk mass distribution: :math:`\rho(r,z)=\rho_0 \exp(-r/r_\mathrm{s}) \hbox{sech}^2(z/z_\mathrm{s})`.
     !!}
     private
     double precision                                                        :: radiusScale                              , scaleHeight                              , &
          &                                                                     densityNormalization                     , surfaceDensityNormalization              , &
          &                                                                     mass
     ! Tables used for rotation curves and potential.
     logical                                                                 :: scaleLengthFactorSet                     , rotationCurveInitialized         =.false., &
          &                                                                     rotationCurveGradientInitialized =.false., potentialInitialized             =.false., &
          &                                                                     accelerationInitialized          =.false.
     double precision                                                        :: scaleLengthFactor
     type            (table1DLogarithmicLinear)                              :: rotationCurveTable                       , rotationCurveGradientTable               , &
          &                                                                     potentialTable
     double precision                          , allocatable, dimension(:  ) :: accelerationRadii                        , accelerationHeights
     double precision                          , allocatable, dimension(:,:) :: accelerationRadial                       , accelerationVertical                     , &
          &                                                                     tidalTensorRadialRadial                  , tidalTensorVerticalVertical              , &
          &                                                                     tidalTensorCross
     double precision                                                        :: accelerationRadiusMinimumLog             , accelerationRadiusMaximumLog             , &
          &                                                                     accelerationHeightMinimumLog             , accelerationHeightMaximumLog             , &
          &                                                                     accelerationRadiusInverseInterval        , accelerationHeightInverseInterval
   contains
     !![
     <methods docformat="rst">
       <method description="Tabulates the potential for an exponential disk mass distribution."                            method="tabulate"                         />
       <method description="Compute the Bessel function factor appearing in the exponential disk rotation curve."          method="besselFactorRotationCurve"        />
       <method description="Compute the Bessel function factor appearing in the exponential disk rotation curve gradient." method="besselFactorRotationCurveGradient"/>
       <method description="Compute the Bessel function factor appearing in the exponential disk potential."               method="besselFactorPotential"            />
       <method description="Tabulate the gravitational acceleration and tidal tensor due to the disk."                     method="accelerationTabulate"             />
       <method description="Interpolate in the tabulated gravitational acceleration and/or tidal tensor due to the disk."  method="accelerationInterpolate"          />
     </methods>
     !!]
     final     ::                                            exponentialDiskDestructor
     procedure :: tabulate                                => exponentialDiskTabulate
     procedure :: besselFactorRotationCurve               => exponentialDiskBesselFactorRotationCurve
     procedure :: besselFactorRotationCurveGradient       => exponentialDiskBesselFactorRotationCurveGradient
     procedure :: besselFactorPotential                   => exponentialDiskBesselFactorPotential
     procedure :: assumeMonotonicDecreasingSurfaceDensity => exponentialDiskAssumeMonotonicDecreasingSurfaceDensity
     procedure :: massTotal                               => exponentialDiskMassTotal
     procedure :: density                                 => exponentialDiskDensity
     procedure :: densityGradientRadial                   => exponentialDensityGradientRadial
     procedure :: densitySlopeLogarithmicCentral          => exponentialDiskDensitySlopeLogarithmicCentral
     procedure :: densitySphericalAverage                 => exponentialDiskDensitySphericalAverage
     procedure :: surfaceDensity                          => exponentialDiskSurfaceDensity
     procedure :: massEnclosedBySphere                    => exponentialDiskMassEnclosedBySphere
     procedure :: massEnclosedByCylinder                  => exponentialDiskMassEnclosedByCylinder
     procedure :: radiusEnclosingSurfaceDensity           => exponentialDiskRadiusEnclosingSurfaceDensity
     procedure :: potentialIsAnalytic                     => exponentialDiskPotentialIsAnalytic
     procedure :: potential                               => exponentialDiskPotential
     procedure :: rotationCurve                           => exponentialDiskRotationCurve
     procedure :: rotationCurveGradient                   => exponentialDiskRotationCurveGradient
     procedure :: radiusHalfMass                          => exponentialDiskRadiusHalfMass
     procedure :: surfaceDensityRadialMoment              => exponentialDiskSurfaceDensityRadialMoment
     procedure :: acceleration                            => exponentialDiskAcceleration
     procedure :: tidalTensor                             => exponentialDiskTidalTensor
     procedure :: accelerationTabulate                    => exponentialDiskAccelerationTabulate
     procedure :: accelerationInterpolate                 => exponentialDiskAccelerationInterpolate
     procedure :: positionSample                          => exponentialDiskPositionSample
  end type massDistributionExponentialDisk

  interface massDistributionExponentialDisk
     !!{RST
     Constructors for the :galacticus-class:`massDistributionExponentialDisk` mass distribution class.
     !!}
     module procedure exponentialDiskConstructorParameters
     module procedure exponentialDiskConstructorInternal
  end interface massDistributionExponentialDisk

  ! The radius (in units of the disk scale length) beyond which the disk is treated as a point mass for the purposes of computing
  ! rotation curves.
  double precision, parameter :: radiusMaximum           =3.0d+1

  ! Potential tabulation.
  integer         , parameter :: potentialPointsPerDecade=10
  double precision, parameter :: potentialRadiusMinimum  =1.0d-3
  double precision, parameter :: potentialRadiusMaximum  =5.0d+2

  ! Seed range for the rotation curve tabulations. Every range begins from this one, so that any two tabulations - built in
  ! different runs, at different requested half-radii - overlap and can be merged onto the common lattice.
  double precision, parameter :: rotationCurveHalfRadiusMinimumDefault=1.0d-6, rotationCurveHalfRadiusMaximumDefault=1.0d+1

contains

  function exponentialDiskConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`massDistributionExponentialDisk` mass distribution class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters          , only : inputParameter                , inputParameters
    use :: Galactic_Structure_Options, only : enumerationComponentTypeEncode, enumerationMassTypeEncode
    implicit none
    type            (massDistributionExponentialDisk)                :: self
    type            (inputParameters                ), intent(inout) :: parameters
    double precision                                                 :: mass         , radiusScale, &
         &                                                              scaleHeight
    logical                                                          :: dimensionless
    type            (varying_string                 )                :: componentType
    type            (varying_string                 )                :: massType

    !![
    <inputParameter docformat="rst">
      <name>scaleHeight</name>
      <defaultSource>
      :cite:p:`kregel_flattening_2002`
      </defaultSource>
      <defaultValue>0.137d0</defaultValue>
      <description>
      The scale height of the exponential disk profile.
      </description>
      <source>parameters</source>
      <minimum inclusive="false">0.0</minimum>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>radiusScale</name>
      <defaultValue>1.0d0</defaultValue>
      <description>
      The scale radius of the exponential disk profile.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>mass</name>
      <defaultValue>1.0d0</defaultValue>
      <description>
      The mass of the exponential disk profile.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>dimensionless</name>
      <defaultValue>.true.</defaultValue>
      <description>
      If true the exponential disk profile is considered to be dimensionless.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>componentType</name>
      <defaultValue>var_str('unknown')</defaultValue>
      <description>
      The component type that this mass distribution represents.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>massType</name>
      <defaultValue>var_str('unknown')</defaultValue>
      <description>
      The mass type that this mass distribution represents.
      </description>
      <source>parameters</source>
    </inputParameter>
    <conditionalCall>
     <call>self=massDistributionExponentialDisk(scaleHeight=scaleHeight,componentType=enumerationComponentTypeEncode(componentType,includesPrefix=.false.),massType=enumerationMassTypeEncode(massType,includesPrefix=.false.){conditions})</call>
     <argument name="mass"          value="mass"          parameterPresent="parameters"/>
     <argument name="radiusScale"   value="radiusScale"   parameterPresent="parameters"/>
     <argument name="dimensionless" value="dimensionless" parameterPresent="parameters"/>
    </conditionalCall>
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function exponentialDiskConstructorParameters

  function exponentialDiskConstructorInternal(radiusScale,scaleHeight,mass,dimensionless,componentType,massType) result(self)
    !!{RST
    Internal constructor for "exponentialDisk" mass distribution class.
    !!}
    use :: Error                   , only : Error_Report
    use :: Numerical_Comparison    , only : Values_Differ
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    type            (massDistributionExponentialDisk)                          :: self
    double precision                                 , intent(in   ), optional :: radiusScale  , scaleHeight, &
         &                                                                        mass
    logical                                          , intent(in   ), optional :: dimensionless
    type            (enumerationComponentTypeType   ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType        ), intent(in   ), optional :: massType
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
       if (present(radiusScale)) then
          if (Values_Differ(radiusScale,1.0d0,absTol=1.0d-6)) call Error_Report('radiusScale should be unity for a dimensionless profile (or simply do not specify a scale length)'//{introspection:location})
       end if
       if (present(mass       )) then
          if (Values_Differ(mass       ,1.0d0,absTol=1.0d-6)) call Error_Report('mass should be unity for a dimensionless profile (or simply do not specify a mass)'               //{introspection:location})
       end if
       self%radiusScale                =1.0d0
       self%mass                       =1.0d0
       self%surfaceDensityNormalization=1.0d0/2.0d0/Pi
    else
       ! Set scale radius.
       if (.not.present(radiusScale)) call Error_Report('scale radius must be specified for dimensionful profiles'//{introspection:location})
       if (.not.present(mass       )) call Error_Report('mass must be specified for dimensionful profiles'        //{introspection:location})
       self%radiusScale                =radiusScale
       self%mass                       =mass
       self%surfaceDensityNormalization=self%mass/2.0d0/Pi/self%radiusScale**2
    end if
    ! Set the scale height.
    if (present(scaleHeight)) then
       self%scaleHeight         = scaleHeight
       self%densityNormalization= self%surfaceDensityNormalization &
            &                    /self%scaleHeight                 &
            &                    /2.0d0
    else
       ! No scale height given, assume a razor-thin disk.
       self%scaleHeight         =0.0d0
       self%densityNormalization=0.0d0
    end if
    ! Initialize rotation curve tables.
    self%scaleLengthFactorSet            =.false.
    self%rotationCurveInitialized        =.false.
    self%rotationCurveGradientInitialized=.false.
    self%potentialInitialized            =.false.
    self%accelerationInitialized         =.false.
    self%scaleLengthFactor               =0.0d0
    return
  end function exponentialDiskConstructorInternal

  subroutine exponentialDiskDestructor(self)
    !!{RST
    Destructor for exponential disk mass distributions.
    !!}
    implicit none
    type(massDistributionExponentialDisk), intent(inout) :: self

    if (self%rotationCurveInitialized        ) call self%rotationCurveTable        %destroy()
    if (self%rotationCurveGradientInitialized) call self%rotationCurveGradientTable%destroy()
    if (self%potentialInitialized            ) call self%potentialTable            %destroy()
    return
  end subroutine exponentialDiskDestructor

  subroutine exponentialDiskTabulate(self)
    !!{RST
    Build tables used for exponential disk mass distributions.
    !!}
    use :: Bessel_Functions, only : Bessel_Function_I0    , Bessel_Function_I1, Bessel_Function_K0, Bessel_Function_K1
    use :: Table_Labels    , only : extrapolationTypeAbort
    implicit none
    class           (massDistributionExponentialDisk), intent(inout) :: self
    integer                                                          :: i   , potentialPointsCount
    double precision                                                 :: x

    ! Build table if necessary.
    if (.not.self%potentialInitialized) then
       ! Determine how many points to tabulate.
       potentialPointsCount=int(log10(potentialRadiusMaximum/potentialRadiusMinimum)*dble(potentialPointsPerDecade))+1
       ! Create the table.
       call self%potentialTable%destroy()
       call self%potentialTable%create(potentialRadiusMinimum,potentialRadiusMaximum,potentialPointsCount,extrapolationType=spread(extrapolationTypeAbort,1,2))
       ! Compute Bessel factors.
       do i=1,potentialPointsCount
          x=self%potentialTable%x(i)
          call self%potentialTable%populate(                                               &
               &                            +x                                             &
               &                            *(                                             &
               &                              +Bessel_Function_I0(x)*Bessel_Function_K1(x) &
               &                              -Bessel_Function_I1(x)*Bessel_Function_K0(x) &
               &                             ),                                            &
               &                            i                                              &
               &                           )
       end do
       ! Record that table is initialized.
       self%potentialInitialized=.true.
    end if
    return
  end subroutine exponentialDiskTabulate

  logical function exponentialDiskAssumeMonotonicDecreasingSurfaceDensity(self) result(assumeMonotonicDecreasingSurfaceDensity)
    !!{RST
    Return true indicating that this distribution has a monotonically-decreasing surface density.
    !!}
    implicit none
    class(massDistributionExponentialDisk), intent(inout) :: self

    assumeMonotonicDecreasingSurfaceDensity=.true.
    return
  end function exponentialDiskAssumeMonotonicDecreasingSurfaceDensity

  double precision function exponentialDiskMassTotal(self)
    !!{RST
    Return the total mass in an exponential disk distribution.
    !!}
    implicit none
    class(massDistributionExponentialDisk), intent(inout) :: self

    exponentialDiskMassTotal=self%mass
    return
  end function exponentialDiskMassTotal

  double precision function exponentialDiskRadiusHalfMass(self)
    !!{RST
    Return the half-mass radius in an exponential disk mass distribution.
    !!}
    implicit none
    class           (massDistributionExponentialDisk), intent(inout) :: self
    double precision                                 , parameter     :: radiusHalfMassToScaleRadius=1.678346990d0

    exponentialDiskRadiusHalfMass=+radiusHalfMassToScaleRadius &
         &                        *self%radiusScale
    return
  end function exponentialDiskRadiusHalfMass

  double precision function exponentialDiskDensity(self,coordinates)
    !!{RST
    Return the density at the specified ``coordinates`` in an exponential disk mass distribution.
    !!}
    use :: Coordinates, only : assignment(=), coordinateCylindrical
    use :: Error      , only : Error_Report
    implicit none
    class           (massDistributionExponentialDisk), intent(inout) :: self
    class           (coordinate                     ), intent(in   ) :: coordinates
    type            (coordinateCylindrical          )                :: position
    double precision                                 , parameter     :: coshArgumentMaximum=50.0d0
    double precision                                                 :: r                         , z, &
         &                                                              coshTerm

    ! If disk is razor thin, density is undefined.
    if (self%scaleHeight <= 0.0d0) call Error_Report('density undefined for razor-thin disk'//{introspection:location})
    ! Get position in cylindrical coordinate system.
    position=coordinates
    ! Compute density.
    r=    position%r() /self%radiusScale
    z=abs(position%z())/self%scaleHeight
    if (z > coshArgumentMaximum) then
       coshTerm=(2.0d0*exp(-z)/(1.0d0+exp(-2.0d0*z)))**2
    else
       coshTerm=1.0d0/cosh(z)**2
    end if
    exponentialDiskDensity=self%densityNormalization*exp(-r)*coshTerm
    return
  end function exponentialDiskDensity

  double precision function exponentialDensityGradientRadial(self,coordinates,logarithmic)
    !!{RST
    Return the density gradient in the radial direction in a scaled spherical mass distribution.
    !!}
    use :: Coordinates, only : assignment(=), coordinateCylindrical
    use :: Error      , only : Error_Report
    implicit none
    class           (massDistributionExponentialDisk), intent(inout), target   :: self
    class           (coordinate                     ), intent(in   )           :: coordinates
    logical                                          , intent(in   ), optional :: logarithmic
    double precision                                 , parameter               :: coshArgumentMaximum=50.0d0
    type            (coordinateCylindrical          )                          :: position
    double precision                                                           :: r                         , z, &
         &                                                                        coshTerm
    !![
    <optionalArgument name="logarithmic" defaultsTo=".false."/>
    !!]
    
    ! If disk is razor thin, density is undefined.
    if (self%scaleHeight <= 0.0d0) call Error_Report('density undefined for razor-thin disk'//{introspection:location})
    ! Get position in cylindrical coordinate system.
    position=coordinates
    ! Compute density.
    r=    position%r() /self%radiusScale
    z=abs(position%z())/self%scaleHeight
    if (z > coshArgumentMaximum) then
       coshTerm=(2.0d0*exp(-z)/(1.0d0+exp(-2.0d0*z)))**2
    else
       coshTerm=1.0d0/cosh(z)**2
    end if
    exponentialDensityGradientRadial=-(r+2.0d0*z*tanh(z))
    if (.not.logarithmic_)                                                                          &
         & exponentialDensityGradientRadial=+         exponentialDensityGradientRadial              &
         &                                  *self    %density                         (coordinates) &
         &                                  /position%rSpherical                      (           )
    return
  end function exponentialDensityGradientRadial

  double precision function exponentialDiskDensitySlopeLogarithmicCentral(self) result(slope)
    !!{RST
    Return the central logarithmic slope of the density profile in an exponential disk mass distribution. The density is
    finite at the center of the disk, and so is evaluable there.
    !!}
    implicit none
    class(massDistributionExponentialDisk), intent(inout) :: self
    !$GLC attributes unused :: self

    slope=0.0d0
    return
  end function exponentialDiskDensitySlopeLogarithmicCentral
  
  double precision function exponentialDiskDensitySphericalAverage(self,radius)
    !!{RST
    Return the spherically-averaged density at the specified ``coordinates`` in an exponential disk mass distribution. Note that this assumes the thin-disk approximation.
    !!}
    implicit none
    class           (massDistributionExponentialDisk), intent(inout) :: self
    double precision                                 , intent(in   ) :: radius

    exponentialDiskDensitySphericalAverage=+0.5d0                                 &
         &                                 *     self%surfaceDensityNormalization &
         &                                 /                               radius &
         &                                 *exp(                                  &
         &                                      -                          radius &
         &                                      /self                %radiusScale &
         &                                     )
    return
  end function exponentialDiskDensitySphericalAverage

  double precision function exponentialDiskMassEnclosedBySphere(self,radius) result(mass)
    !!{RST
    Computes the mass enclosed within a sphere of given ``radius`` for exponential disk mass distributions. Note that this assumes the thin-disk approximation.
    !!}
    implicit none
    class           (massDistributionExponentialDisk), intent(inout), target :: self
    double precision                                 , intent(in   )         :: radius

    mass=self%massEnclosedByCylinder(radius)
    return
  end function exponentialDiskMassEnclosedBySphere

  double precision function exponentialDiskMassEnclosedByCylinder(self,radius) result(mass)
    !!{RST
    Computes the mass enclosed within a cylinder of given ``radius`` for exponential disk mass distributions.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (massDistributionExponentialDisk), intent(inout), target :: self
    double precision                                 , intent(in   )         :: radius
    double precision                                                         :: fractionalRadius

    fractionalRadius=+radius                              &
         &           /self%radiusScale
    mass            =+2.0d0                               &
         &           *Pi                                  &
         &           *self%radiusScale                **2 &
         &           *self%surfaceDensityNormalization    &
         &           *(                                   &
         &             +1.0d0                             &
         &             -(                                 &
         &               +1.0d0                           &
         &               +fractionalRadius                &
         &              )                                 &
         &             *exp(-fractionalRadius)            &
         &           )
    return
  end function exponentialDiskMassEnclosedByCylinder

  double precision function exponentialDiskSurfaceDensity(self,coordinates)
    !!{RST
    Return the surface density at the specified ``coordinates`` in an exponential disk mass distribution.
    !!}
    use :: Coordinates, only : coordinate
    implicit none
    class           (massDistributionExponentialDisk), intent(inout) :: self
    class           (coordinate                     ), intent(in   ) :: coordinates
    double precision                                                 :: r

    ! Get the radial coordinate.
    r=coordinates%rCylindrical()/self%radiusScale
    ! Compute the density.
    exponentialDiskSurfaceDensity=self%surfaceDensityNormalization*exp(-r)
    return
  end function exponentialDiskSurfaceDensity

  double precision function exponentialDiskRadiusEnclosingSurfaceDensity(self,densitySurface,radiusGuess) result(radius)
    !!{RST
    Computes the radius enclosing a given surface density for exponential disk mass distributions.
    !!}    
    implicit none
    class           (massDistributionExponentialDisk), intent(inout), target   :: self
    double precision                                 , intent(in   )           :: densitySurface
    double precision                                 , intent(in   ), optional :: radiusGuess
    
    radius=-     self%radiusScale                 &
    &      *log(                                  &
    &           +     densitySurface              &
    &           /self%surfaceDensityNormalization &
    &          )
    return
  end function exponentialDiskRadiusEnclosingSurfaceDensity
  
  double precision function exponentialDiskRotationCurve(self,radius)
    !!{RST
    Return the mid-plane rotation curve for an exponential disk.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    class           (massDistributionExponentialDisk), intent(inout) :: self
    double precision                                 , intent(in   ) :: radius
    double precision                                                 :: r           , halfRadius, &
         &                                                              radiusFactor

    ! Get scale-free radius.
    r=radius/self%radiusScale
    ! Compute rotation curve.
    if (r > radiusMaximum) then
       ! Beyond some maximum radius, approximate the disk as a spherical distribution to avoid evaluating Bessel functions for
       ! very large arguments.
       exponentialDiskRotationCurve=sqrt(self%mass/radius)
    else
       ! We are often called at precisely one scale length. Use pre-computed factors in that case.
       if (r == 1.0d0) then
          if (.not.self%scaleLengthFactorSet) then
             halfRadius               =0.5d0
             self%scaleLengthFactor   =self%besselFactorRotationCurve(halfRadius)
             self%scaleLengthFactorSet=.true.
          end if
          radiusFactor=self%scaleLengthFactor
       else
          halfRadius       =0.5d0*r
          radiusFactor=self%besselFactorRotationCurve(halfRadius)
       end if
       exponentialDiskRotationCurve=sqrt(2.0d0*(self%mass/self%radiusScale)*radiusFactor)
    end if
    ! Make dimensionful if necessary.
    if (.not.self%dimensionless) exponentialDiskRotationCurve= &
         & +sqrt(gravitationalConstant_internal)               &
         & *exponentialDiskRotationCurve
    return
  end function exponentialDiskRotationCurve

  double precision function exponentialDiskRotationCurveGradient(self,radius)
    !!{RST
    Return the mid-plane rotation curve gradient for an exponential disk.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    class           (massDistributionExponentialDisk), intent(inout) :: self
    double precision                                 , intent(in   ) :: radius
    double precision                                 , parameter     :: fractionalRadiusMaximum=30.0d0
    double precision                                                 :: besselArgument                , besselFactor

    ! Compute Bessel functions argument.
    besselArgument=+radius           &
         &         /2.0d0            &
         &         /self%radiusScale
    if (2.0d0*besselArgument > fractionalRadiusMaximum) then
       ! Beyond some maximum radius, approximate the disk as a point mass to avoid evaluating Bessel functions for
       ! very large arguments.
       exponentialDiskRotationCurveGradient=-self%mass/radius**2
    else
       ! Compute the gradient.
       besselFactor=self%besselFactorRotationCurveGradient(besselArgument)
       exponentialDiskRotationCurveGradient=+self%mass           &
            &                               *besselFactor        &
            &                               /self%radiusScale**2
    end if
    ! Make dimensionful if necessary.
    if (.not.self%dimensionless) exponentialDiskRotationCurveGradient= &
         &  +gravitationalConstant_internal                            &
         &  *exponentialDiskRotationCurveGradient
    return
  end function exponentialDiskRotationCurveGradient

  logical function exponentialDiskPotentialIsAnalytic(self) result(isAnalytic)
    !!{RST
    Return that the potential has an analytic form.
    !!}
    implicit none
    class(massDistributionExponentialDisk), intent(inout) :: self

    isAnalytic=.true.
    return
  end function exponentialDiskPotentialIsAnalytic

  double precision function exponentialDiskPotential(self,coordinates,status)
    !!{RST
    Return the gravitational potential for an exponential disk.
    !!}
    use :: Coordinates                     , only : assignment(=)                 , coordinateCylindrical
    use :: Galactic_Structure_Options      , only : structureErrorCodeSuccess
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    class           (massDistributionExponentialDisk  ), intent(inout), target   :: self
    class           (coordinate                       ), intent(in   )           :: coordinates
    type            (enumerationStructureErrorCodeType), intent(  out), optional :: status
    type            (coordinateCylindrical            )                          :: position
    double precision                                                             :: correctionSmallRadius, halfRadius, &
         &                                                                          radius

    if (present(status)) status=structureErrorCodeSuccess
    ! Get position in cylindrical coordinate system.
    position=coordinates
    ! Compute density.
    radius=position%r()
    ! If the radius is sufficiently large, treat the disk as a point mass.
    if (radius > potentialRadiusMaximum*self%radiusScale) then
       exponentialDiskPotential=-self%mass/radius
    else
       ! Radius is sufficiently small to use the full calculation.
       ! Compute the potential. If the radius is lower than the height then approximate the disk
       ! mass as being spherically distributed.
       if (radius > self%scaleHeight) then
          halfRadius           =radius/2.0d0/self%radiusScale
          correctionSmallRadius=0.0d0
       else
          halfRadius           =+self%scaleHeight/self%radiusScale/2.0d0
          correctionSmallRadius=+self%massEnclosedBySphere(self%scaleHeight) &
               &                /self%scaleHeight
          if (radius > 0.0d0) correctionSmallRadius=+correctionSmallRadius             &
               &                                    -self%massEnclosedBySphere(radius) &
               &                                    /radius
       end if
       ! Compute the potential including the correction to small radii.
       exponentialDiskPotential=                                  &
            &             -self%mass                              &
            &             /self%radiusScale                       &
            &             *self%besselFactorPotential(halfRadius) &
            &             +correctionSmallRadius
    end if
    ! Make dimensionful if necessary.
    if (.not.self%dimensionless) exponentialDiskPotential= &
         &  +gravitationalConstant_internal                &
         &  *exponentialDiskPotential
    return
  end function exponentialDiskPotential

  double precision function exponentialDiskBesselFactorPotential(self,halfRadius)
    !!{RST
    Compute Bessel function factors appearing in the expression for an razor-thin exponential disk gravitational potential.
    !!}
    use :: Numerical_Constants_Math, only : eulersConstant, ln2
    implicit none
    class           (massDistributionExponentialDisk), intent(inout) :: self
    double precision                                 , intent(in   ) :: halfRadius

    ! For small half-radii, use a series expansion for a more accurate result.
    if      (halfRadius <=                 0.0d0) then
       exponentialDiskBesselFactorPotential=1.0d0
    else if (halfRadius < potentialRadiusMinimum) then
       exponentialDiskBesselFactorPotential=1.0d0+(eulersConstant-ln2+log(halfRadius))*halfRadius**2
    else
       call self%tabulate()
       exponentialDiskBesselFactorPotential=self%potentialTable%interpolate(halfRadius)
    end if
    return
  end function exponentialDiskBesselFactorPotential

  double precision function exponentialDiskBesselFactorRotationCurve(self,halfRadius)
    !!{RST
    Compute Bessel function factors appearing in the expression for an razor-thin exponential disk rotation curve.
    !!}
    use :: Bessel_Functions        , only : Bessel_Function_I0, Bessel_Function_I1, Bessel_Function_K0, Bessel_Function_K1
    use :: Numerical_Constants_Math, only : eulersConstant    , ln2
    use :: Numerical_Ranges        , only : Range_Pinned      , rangeLattice      , gridSchemePerDecade
    implicit none
    class           (massDistributionExponentialDisk), intent(inout)               :: self
    double precision                                 , intent(in   )               :: halfRadius
    double precision                                 , parameter                   :: halfRadiusSmall             =1.0d-3
    integer                                          , parameter                   :: rotationCurvePointsPerDecade=100
    integer                                                                        :: iPoint
    double precision                                                               :: x
    type            (rangeLattice                   )                              :: lattice
    logical                                          , allocatable  , dimension(:) :: isComputed

    ! For small half-radii, use a series expansion for a more accurate result.
    if (halfRadius <= 0.0d0) then
       exponentialDiskBesselFactorRotationCurve=0.0d0
       return
    else if (halfRadius < halfRadiusSmall) then
       exponentialDiskBesselFactorRotationCurve=(ln2-eulersConstant-0.5d0-log(halfRadius))*halfRadius**2
       return
    end if
    ! Find the range of half-radii to tabulate, pinned to an absolute lattice so that the tabulation - and hence every value
    ! interpolated from it - is independent of the half-radius at which it happened to be first requested, and so that it can be
    ! extended without recomputing any value already found. The safety margin which `Range_Pinned` applies by default is a factor
    ! of two at each end, which is the margin this tabulation always used.
    lattice=Range_Pinned(                                                                                              &
         &                              halfRadius                                                                   , &
         &                              rotationCurvePointsPerDecade                                                 , &
         &                              gridSchemePerDecade                                                          , &
         &               rangeCurrent  =[rotationCurveHalfRadiusMinimumDefault,rotationCurveHalfRadiusMaximumDefault], &
         &               latticeCurrent=self%rotationCurveTable%lattice                                                &
         &              )
    if (.not.self%rotationCurveInitialized.or..not.self%rotationCurveTable%lattice%covers(lattice)) then
       ! Extend the tabulation onto the new lattice, preserving every value already computed.
       call self%rotationCurveTable%extend(lattice,isComputed)
       ! Compute Bessel factors, skipping any point whose value was preserved by the extension.
       do iPoint=1,lattice%count
          if (isComputed(iPoint)) cycle
          x=self%rotationCurveTable%x(iPoint)
          call self%rotationCurveTable%populate(                                               &
               &                                +x**2                                          &
               &                                *(                                             &
               &                                  +Bessel_Function_I0(x)*Bessel_Function_K0(x) &
               &                                  -Bessel_Function_I1(x)*Bessel_Function_K1(x) &
               &                                 ),                                            &
               &                                iPoint                                         &
               &                               )
       end do
       ! Flag that the rotation curve is now initialized.
       self%rotationCurveInitialized=.true.
    end if
    ! Interpolate in the tabulated function.
    exponentialDiskBesselFactorRotationCurve=self%rotationCurveTable%interpolate(halfRadius)
    return
  end function exponentialDiskBesselFactorRotationCurve

  double precision function exponentialDiskBesselFactorRotationCurveGradient(self,halfRadius)
    !!{RST
    Compute Bessel function factors appearing in the expression for a razor-thin exponential disk rotation curve gradient.
    !!}
    use :: Bessel_Functions        , only : Bessel_Function_I0, Bessel_Function_I1, Bessel_Function_K0 , Bessel_Function_K1
    use :: Numerical_Constants_Math, only : eulersConstant    , ln2
    use :: Numerical_Ranges        , only : Range_Pinned      , rangeLattice      , gridSchemePerDecade
    implicit none
    class           (massDistributionExponentialDisk), intent(inout)               :: self
    double precision                                 , intent(in   )               :: halfRadius
    double precision                                 , parameter                   :: halfRadiusSmall                     =1.0d-3
    double precision                                 , parameter                   :: halfRadiusLarge                     =1.0d+2
    integer                                          , parameter                   :: rotationCurveGradientPointsPerDecade=100
    integer                                                                        :: iPoint
    double precision                                                               :: x
    type            (rangeLattice                   )                              :: lattice
    logical                                          , allocatable  , dimension(:) :: isComputed

    ! For small and large half-radii, use a series expansion for a more accurate result.
    if (halfRadius == 0.0d0) then
       exponentialDiskBesselFactorRotationCurveGradient=0.0d0
       return
    else if (halfRadius < halfRadiusSmall) then
       exponentialDiskBesselFactorRotationCurveGradient=(ln2-log(halfRadius)-eulersConstant-1.0d0)*halfRadius &
            & +(1.5d0*ln2-1.5d0*log(halfRadius)+0.25d0-1.5d0*eulersConstant)*halfRadius**3
       return
    else if (halfRadius > halfRadiusLarge) then
       exponentialDiskBesselFactorRotationCurveGradient=-0.125d0-27.0d0/64.0d0/halfRadius**2
       return
    end if
    ! As for the rotation curve itself, pin the range to an absolute lattice so that the tabulation does not depend on the
    ! half-radius at which it was first requested, and can be extended without recomputing any value already found.
    lattice=Range_Pinned(                                                                                              &
         &                              halfRadius                                                                   , &
         &                              rotationCurveGradientPointsPerDecade                                         , &
         &                              gridSchemePerDecade                                                          , &
         &               rangeCurrent  =[rotationCurveHalfRadiusMinimumDefault,rotationCurveHalfRadiusMaximumDefault], &
         &               latticeCurrent=self%rotationCurveGradientTable%lattice                                        &
         &              )
    if (.not.self%rotationCurveGradientInitialized.or..not.self%rotationCurveGradientTable%lattice%covers(lattice)) then
       ! Extend the tabulation onto the new lattice, preserving every value already computed.
       call self%rotationCurveGradientTable%extend(lattice,isComputed)
       ! Compute Bessel factors, skipping any point whose value was preserved by the extension.
       do iPoint=1,lattice%count
          if (isComputed(iPoint)) cycle
          x=self%rotationCurveGradientTable%x(iPoint)
          call self%rotationCurveGradientTable%populate                                      &
               &  (                                                                          &
               &   +2.0d0                                                                    &
               &   *x                                                                        &
               &   *(                                                                        &
               &     +  Bessel_Function_I0(x)                         *Bessel_Function_K0(x) &
               &     -  Bessel_Function_I1(x)                         *Bessel_Function_K1(x) &
               &    )                                                                        &
               &   +x**2                                                                     &
               &   *(                                                                        &
               &     +  Bessel_Function_I1(x)                         *Bessel_Function_K0(x) &
               &     -  Bessel_Function_K1(x)                         *Bessel_Function_I0(x) &
               &     -(+Bessel_Function_I0(x)-Bessel_Function_I1(x)/x)*Bessel_Function_K1(x) &
               &     -(-Bessel_Function_K0(x)-Bessel_Function_K1(x)/x)*Bessel_Function_I1(x) &
               &    ),                                                                       &
               &   iPoint                                                                    &
               &  )
       end do
       ! Flag that the rotation curve is now initialized.
       self%rotationCurveGradientInitialized=.true.
    end if
    ! Interpolate in the tabulated function.
    exponentialDiskBesselFactorRotationCurveGradient=self%rotationCurveGradientTable%interpolate(halfRadius)
    return
  end function exponentialDiskBesselFactorRotationCurveGradient

  double precision function exponentialDiskSurfaceDensityRadialMoment(self,moment,radiusMinimum,radiusMaximum,isInfinite)
    !!{RST
    Compute radial moments of the exponential disk mass distribution surface density profile.
    !!}
    use :: Error          , only : Error_Report
    use :: Gamma_Functions, only : Gamma_Function, Gamma_Function_Incomplete
    implicit none
    class           (massDistributionExponentialDisk), intent(inout)           :: self
    double precision                                 , intent(in   )           :: moment
    double precision                                 , intent(in   ), optional :: radiusMinimum, radiusMaximum
    logical                                          , intent(  out), optional :: isInfinite
    double precision                                                           :: integralLow  , integralHigh

    ! All moments n>-1 are finite.
    if (present(isInfinite)) isInfinite=(moment <= -1.0d0)
    if (moment <= -1.0d0) then
       exponentialDiskSurfaceDensityRadialMoment=0.0d0
       if (present(isInfinite)) return
       call Error_Report('moment is infinite'//{introspection:location})
    end if
    ! Compute the moment.
    if (present(radiusMinimum)) then
       integralLow =-Gamma_Function_Incomplete(moment+1.0d0,radiusMinimum)
    else
       integralLow =-Gamma_Function_Incomplete(moment+1.0d0,0.0d0        )
    end if
    if (present(radiusMaximum)) then
       integralHigh=-Gamma_Function_Incomplete(moment+1.0d0,radiusMaximum)
    else
       integralHigh=+0.0d0
    end if
    exponentialDiskSurfaceDensityRadialMoment=(integralHigh-integralLow)*Gamma_Function(moment+1.0d0)*self%radiusScale**(moment+1.0d0)
    return
  end function exponentialDiskSurfaceDensityRadialMoment

  function exponentialDiskAcceleration(self,coordinates)
    !!{RST
    Computes the gravitational acceleration at ``coordinates`` for exponential disk mass distributions.
    !!}
    use :: Coordinates                     , only : assignment(=), coordinateCartesian           , coordinateCylindrical
    use :: Numerical_Constants_Astronomical, only : gigaYear     , gravitationalConstant_internal, megaParsec
    use :: Numerical_Constants_Prefixes    , only : kilo
    implicit none
    type            (coordinateCartesian            )                :: exponentialDiskAcceleration
    class           (massDistributionExponentialDisk), intent(inout) :: self
    class           (coordinate                     ), intent(in   ) :: coordinates
    double precision                                 , dimension(3)  :: positionCartesian         , accelerationVector
    type            (coordinateCylindrical          )                :: coordinatesCylindrical
    type            (coordinateCartesian            )                :: coordinatesCartesian
    double precision                                                 :: accelerationRadial        , accelerationVertical

    ! Get position in cylindrical and Cartesian coordinate systems.
    coordinatesCylindrical=coordinates
    coordinatesCartesian  =coordinates
    positionCartesian     =coordinatesCartesian
    ! Ensure that acceleration is tabulated.
    call self%accelerationTabulate()
    ! Interpolate in the tables.
    call self%accelerationInterpolate(coordinatesCylindrical,accelerationRadial,accelerationVertical)
    ! Convert components of the acceleration vector to Cartesian coordinate system.
    if (abs(accelerationRadial) > 0.0d0) then
       accelerationVector(1:2)=+           accelerationRadial             &
            &                  *           positionCartesian       (1:2)  &
            &                  /           coordinatesCylindrical%r(   )
    else
       accelerationVector(1:2)=+0.0d0
    end if
    accelerationVector   (3  )=+           accelerationVertical           &
         &                     *sign(1.0d0,coordinatesCylindrical%z(   ))
    ! For dimensionful profiles apply the unit conversions and scalings.
    if (.not.self%isDimensionless())                     &
         & accelerationVector=+accelerationVector        &
         &                    *kilo                           &
         &                    *gigaYear                       &
         &                    /megaParsec                     &
         &                    *gravitationalConstant_internal &
         &                    *self%mass                      &
         &                    /self%radiusScale**2
    exponentialDiskAcceleration=accelerationVector
    return
  end function exponentialDiskAcceleration

  function exponentialDiskTidalTensor(self,coordinates)
    !!{RST
    Computes the gravitational tidal tensor at ``coordinates`` for exponential disk mass distributions.
    !!}
    use :: Coordinates                     , only : assignment(=)                 , coordinateCartesian, coordinateCylindrical
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    type            (tensorRank2Dimension3Symmetric )                :: exponentialDiskTidalTensor
    class           (massDistributionExponentialDisk), intent(inout) :: self
    class           (coordinate                     ), intent(in   ) :: coordinates
    double precision                                 , dimension(3)  :: positionCartesian
    double precision                                 , parameter     :: radiusCylindricalSmall    =1.0d-6
    type            (coordinateCartesian            )                :: coordinatesCartesian
    type            (coordinateCylindrical          )                :: coordinatesCylindrical
    double precision                                                 :: accelerationRadial               , accelerationVertical       , &
         &                                                              tidalTensorRadialRadial          , tidalTensorVerticalVertical, &
         &                                                              tidalTensorCross                 , radiusCylindrical

    ! Get position in cylindrical and Cartesian coordinate systems.
    coordinatesCylindrical=coordinates
    coordinatesCartesian  =coordinatesCylindrical
    positionCartesian     =coordinatesCartesian
    radiusCylindrical     =coordinatesCylindrical%r()/self%radiusScale
    positionCartesian     =positionCartesian         /self%radiusScale
    ! Ensure that acceleration is tabulated.
    call self%accelerationTabulate()
    ! Interpolate in the tables.
    call self%accelerationInterpolate(coordinatesCylindrical,accelerationRadial,accelerationVertical,tidalTensorRadialRadial,tidalTensorVerticalVertical,tidalTensorCross)
    if (coordinatesCylindrical%z() < 0.0d0) then
       accelerationVertical=-accelerationVertical
       tidalTensorCross    =-tidalTensorCross
    end if
    ! Convert tensor components to Cartesian coordinate system.
    if (radiusCylindrical < radiusCylindricalSmall) then
       exponentialDiskTidalTensor=tensorRank2Dimension3Symmetric(                                                                                                                                                                                                                              &
            &                                                    x00=+tidalTensorRadialRadial                                                                                                                                                                                                , &
            &                                                    x01=+0.0d0                                                                                                                                                                                                                  , &
            &                                                    x11=+tidalTensorRadialRadial                                                                                                                                                                                                , &
            &                                                    x02=+0.0d0                                                                                                                                                                                                                  , &
            &                                                    x12=+0.0d0                                                                                                                                                                                                                  , &
            &                                                    x22=                                                                                                                     +tidalTensorVerticalVertical                                                                         &
            &                                                   )
    else
       exponentialDiskTidalTensor=tensorRank2Dimension3Symmetric(                                                                                                                                                                                                                              &
            &                                                    x00=+accelerationRadial/radiusCylindrical*(+1.0d0-(positionCartesian(1)**2                        /radiusCylindrical**2))+tidalTensorRadialRadial    *(positionCartesian(1)**2                        /radiusCylindrical**2), &
            &                                                    x01=+accelerationRadial/radiusCylindrical*(      -(positionCartesian(1)   *positionCartesian(2)   /radiusCylindrical**2))+tidalTensorRadialRadial    *(positionCartesian(1)   *positionCartesian(2)   /radiusCylindrical**2), &
            &                                                    x11=+accelerationRadial/radiusCylindrical*(+1.0d0-(                        positionCartesian(2)**2/radiusCylindrical**2))+tidalTensorRadialRadial    *(                        positionCartesian(2)**2/radiusCylindrical**2), &
            &                                                    x02=                                                                                                                     +tidalTensorCross           * positionCartesian(1)                           /radiusCylindrical    , &
            &                                                    x12=                                                                                                                     +tidalTensorCross           *                         positionCartesian(2)   /radiusCylindrical    , &
            &                                                    x22=                                                                                                                     +tidalTensorVerticalVertical                                                                         &
            &                                                   )
    end if
    ! For dimensionful profiles apply the unit conversions and scalings.
    if (.not.self%isDimensionless())                                   &
         & exponentialDiskTidalTensor=+exponentialDiskTidalTensor      &
         &                             *gravitationalConstant_internal &
         &                             *self%mass                      &
         &                             /self%radiusScale**3
    return
  end function exponentialDiskTidalTensor
  
  subroutine exponentialDiskAccelerationInterpolate(self,coordinatesCylindrical,accelerationRadial,accelerationVertical,tidalTensorRadialRadial,tidalTensorVerticalVertical,tidalTensorCross)
    !!{RST
    Interpolate gravitational accelerations and tidal tensors in the tabulated solutions for an exponential disk.
    !!}
    use :: Coordinates, only : assignment(=), coordinateCartesian, coordinateCylindrical
    implicit none
    class           (massDistributionExponentialDisk), intent(inout)            :: self
    type            (coordinateCylindrical          ), intent(in   )            :: coordinatesCylindrical
    double precision                                 , intent(  out) , optional :: accelerationRadial     , accelerationVertical       , &
         &                                                                         tidalTensorRadialRadial, tidalTensorVerticalVertical, &
         &                                                                         tidalTensorCross
    type            (coordinateCartesian            )                           :: coordinatesCartesian
    double precision                                 , dimension(3  )           :: positionCartesian
    double precision                                 , dimension(0:1)           :: hRadius                , hHeight
    double precision                                                            :: radiusCylindrical      , radiusCylindricalLog       , &
         &                                                                         heightCylindrical      , heightCylindricalLog       , &
         &                                                                         radiusSpherical
    integer                                                                     :: iRadius                , iHeight                    , &
         &                                                                         jRadius                , jHeight

    ! Find interpolating factors.
    radiusCylindrical=    coordinatesCylindrical%r()/self%radiusScale
    heightCylindrical=abs(coordinatesCylindrical%z()/self%radiusScale)
    if     (                                                                              &
         &   radiusCylindrical > self%accelerationRadii  (size(self%accelerationRadii  )) &
         &  .or.                                                                          &
         &   heightCylindrical > self%accelerationHeights(size(self%accelerationHeights)) &
         & ) then
       ! Use spherical approximation.
       coordinatesCartesian= coordinatesCylindrical
       positionCartesian   = coordinatesCartesian
       radiusSpherical     =+sqrt(sum(positionCartesian**2))     
       if (present(accelerationVertical))                                                                             &
            & accelerationVertical       =-                                                      heightCylindrical    &
            &                             /             radiusSpherical**3                                            &
            &                             *        self%radiusScale    **2
       if (present(accelerationRadial  ))                                                                             &
            & accelerationRadial         =-                                 radiusCylindrical                         &
            &                             /             radiusSpherical**3                                            &
            &                             *        self%radiusScale    **2
       if (present(tidalTensorVerticalVertical))                                                                      &
            & tidalTensorVerticalVertical=+(                                                                          &
            &                               -(1.0d0/    radiusSpherical**3)                                           &
            &                               +(3.0d0/    radiusSpherical**5)*                     heightCylindrical**2 &
            &                              )                                                                          &
            &                             *        self%radiusScale    **3
       if (present(tidalTensorRadialRadial    ))                                                                      &
            & tidalTensorRadialRadial    =+(                                                                          &
            &                               -(1.0d0/    radiusSpherical**3)                                           &
            &                               +(3.0d0/    radiusSpherical**5)*radiusCylindrical**2                      &
            &                              )                                                                          &
            &                             *        self%radiusScale    **3
       if (present(tidalTensorCross           ))                                                                      &
            & tidalTensorCross           =+(                                                                          &
            &                               +(3.0d0/    radiusSpherical**5)*radiusCylindrical   *heightCylindrical    &
            &                              )                                                                          &
            &                             *        self%radiusScale    **3
    else
       ! Interpolate in tabulated solution.
       if (radiusCylindrical < self%accelerationRadii  (1)) then
          ! The tidal tensor must become independent of R as R→0 (since the radial acceleration goes to zero at R=0).
          iRadius   = 1
          hRadius   =[0.0d0,1.0d0]
       else
          radiusCylindricalLog   =log(radiusCylindrical)
          hRadius             (0)=+(                                        &
               &                    +     radiusCylindricalLog              &
               &                    -self%accelerationRadiusMinimumLog      &
               &                   )                                        &
               &                  *  self%accelerationRadiusInverseInterval
          iRadius                =+int(hRadius(0))+             1
          hRadius             (0)=+    hRadius(0) -dble(iRadius-1    )
          hRadius             (1)=-    hRadius(0) +             1.0d0
       end if
       if (heightCylindrical < self%accelerationHeights(1)) then
          ! The tidal tensor must become independent of z as z→0 (since the vertical acceleration goes to zero at z=0).
          iHeight                = 1
          hHeight                =[0.0d0,1.0d0]
       else
          heightCylindricalLog   =log(heightCylindrical)
          hHeight             (0)=+(                                        &
               &                    +     heightCylindricalLog              &
               &                    -self%accelerationHeightMinimumLog      &
               &                   )                                        &
               &                  *  self%accelerationHeightInverseInterval
          iHeight                =+int(hHeight(0))+             1
          hHeight             (0)=+    hHeight(0) -dble(iHeight-1    )
          hHeight             (1)=-    hHeight(0) +             1.0d0
       end if
       if (present(accelerationRadial         )) accelerationRadial         =0.0d0
       if (present(accelerationVertical       )) accelerationVertical       =0.0d0
       if (present(tidalTensorRadialRadial    )) tidalTensorRadialRadial    =0.0d0
       if (present(tidalTensorVerticalVertical)) tidalTensorVerticalVertical=0.0d0
       if (present(tidalTensorCross           )) tidalTensorCross           =0.0d0
       do jRadius=0,1
          do jHeight=0,1
             if (present(accelerationRadial         )) accelerationRadial         =accelerationRadial         +self%accelerationRadial         (iRadius+jRadius,iHeight+jHeight)*(1.0d0-hRadius(jRadius))*(1.0d0-hHeight(jHeight))
             if (present(accelerationVertical       )) accelerationVertical       =accelerationVertical       +self%accelerationVertical       (iRadius+jRadius,iHeight+jHeight)*(1.0d0-hRadius(jRadius))*(1.0d0-hHeight(jHeight))
             if (present(tidalTensorRadialRadial    )) tidalTensorRadialRadial    =tidalTensorRadialRadial    +self%tidalTensorRadialRadial    (iRadius+jRadius,iHeight+jHeight)*(1.0d0-hRadius(jRadius))*(1.0d0-hHeight(jHeight))
             if (present(tidalTensorVerticalVertical)) tidalTensorVerticalVertical=tidalTensorVerticalVertical+self%tidalTensorVerticalVertical(iRadius+jRadius,iHeight+jHeight)*(1.0d0-hRadius(jRadius))*(1.0d0-hHeight(jHeight))
             if (present(tidalTensorCross           )) tidalTensorCross           =tidalTensorCross           +self%tidalTensorCross           (iRadius+jRadius,iHeight+jHeight)*(1.0d0-hRadius(jRadius))*(1.0d0-hHeight(jHeight))
          end do
       end do
       if (present(accelerationRadial  ) .and. radiusCylindrical < self%accelerationRadii  (1)) accelerationRadial  =accelerationRadial  *radiusCylindrical/self%accelerationRadii  (1)
       if (present(accelerationVertical) .and. heightCylindrical < self%accelerationHeights(1)) accelerationVertical=accelerationVertical*heightCylindrical/self%accelerationHeights(1)
    end if
    return
  end subroutine exponentialDiskAccelerationInterpolate

  subroutine exponentialDiskAccelerationTabulate(self)
    !!{RST
    Tabulate the acceleration (and tidal tensor) due to the exponential disk mass distribution. Uses the approach of
    :cite:t:`kuijken_mass_1989`, expressing each quantity as a Hankel transform of the surface density multiplied by a vertical
    kernel, :math:`K(k,z)`, which accounts for the :math:`\mathrm{sech}^2` vertical structure of the disk. The vertical kernel is
    evaluated by direct numerical integration (see ``kernelVertical``) rather than by the series expansion of
    :cite:t:`kuijken_mass_1989`, the terms of which do not decay with order. The tabulation is built for a dimensionless disk.
    !!}
    use :: Bessel_Functions        , only : Bessel_Function_J0_Zero, Bessel_Function_J1_Zero, Bessel_Function_Jn_Zero
    use :: Display                 , only : displayCounter         , displayCounterClear    , displayIndent          , displayUnindent, &
          &                                 verbosityLevelWorking
    use :: File_Utilities          , only : Directory_Make         , File_Exists            , File_Lock              , File_Path      , &
          &                                 File_Unlock            , lockDescriptor
    use :: Error                   , only : Error_Report
    use :: Input_Paths             , only : inputPath              , pathTypeDataDynamic
    use :: HDF5_Access             , only : hdf5Access
    use :: IO_HDF5                 , only : hdf5File
    use :: ISO_Varying_String      , only : char                   , operator(//)           , varying_string
    use :: Numerical_Integration   , only : integrator
    use :: Numerical_Ranges        , only : Make_Range             , rangeTypeLogarithmic
    use :: String_Handling         , only : operator(//)
    implicit none
    class           (massDistributionExponentialDisk), intent(inout)              :: self
    double precision                                 , parameter                  :: radiusMinimum                    = 1.0d-2, radiusMaximum    =5.0d1
    double precision                                 , parameter                  :: radiiPerDecade                   =30.0d+0
    double precision                                 , parameter                  :: wavenumberMaximumFactor          =10.0d+0
    integer                                          , parameter                  :: fileFormatVersionCurrent         = 2
    double precision                                 , parameter                  :: wavenumberKernelMinimum          = 1.0d-3, wavenumberKernelMaximum=1.0d5
    double precision                                 , parameter                  :: wavenumberKernelPerDecade        =32.0d+0
    type            (integrator                     ), save         , allocatable :: integratorAccelerationRadial             , integratorAccelerationVertical       , &
         &                                                                           integratorTidalTensorRadialRadial        , integratorTidalTensorVerticalVertical, &
         &                                                                           integratorTidalTensorCross               , integratorKernel                     , &
         &                                                                           integratorKernelGradient                 , integratorKernelGradient2
    logical                                          , save                       :: converged
    integer                                          , save                       :: iBesselZero                              , besselOrder                          , &
         &                                                                           indexHeight
    double precision                                 , save                       :: height                                   , accelerationDelta                    , &
         &                                                                           wavenumberLow                            , wavenumberHigh                       , &
         &                                                                           tidalTensorDelta                         , tidalTensorRadialRadial              , &
         &                                                                           wavenumberScaled                         , heightScaled
    !$omp threadprivate(height,indexHeight,accelerationDelta,tidalTensorDelta,tidalTensorRadialRadial,wavenumberLow,wavenumberHigh,wavenumberScaled,heightScaled,iBesselZero,besselOrder,converged,integratorAccelerationRadial,integratorAccelerationVertical,integratorTidalTensorRadialRadial,integratorTidalTensorVerticalVertical,integratorTidalTensorCross,integratorKernel,integratorKernelGradient,integratorKernelGradient2)
    integer                                                                       :: countRadii                               , iRadius                              , &
         &                                                                           iHeight                                  , countWork                            , &
         &                                                                           countWavenumbers                         , iWavenumber
    double precision                                                              :: radius                                   , scaleHeight                          , &
         &                                                                           wavenumberKernelLogarithmicMinimum       , wavenumberKernelLogarithmicInterval
    double precision                                 , allocatable , dimension(:,:,:) :: kernelTable

    ! Return if acceleration is initialized.
    if (self%accelerationInitialized) return
    accelerationTabulate: block
      type     (varying_string) :: fileName
      character(len=8         ) :: label
      type     (lockDescriptor) :: fileLock
      
      ! Construct a file name for the table. The file format version is included in the name so that tables built by earlier
      ! versions (which used an incorrect vertical kernel - see below) are never reused, while leaving those tables in place for
      ! any older executables which share the same dynamic data path.
      write (label,'(f8.6)') self%scaleHeight/self%radiusScale
      fileName=inputPath(pathTypeDataDynamic)// &
           &   'galacticStructure/'          // &
           &   self%objectType()             // &
           &   '_h'                          // &
           &   trim(adjustl(label))          // &
           &   '_v'                          // &
           &   fileFormatVersionCurrent      // &
           &   '.hdf5'
      call Directory_Make(File_Path(fileName))
      ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
      call File_Lock(fileName,fileLock,lockIsShared=.true.)
      if (File_Exists(fileName)) then
         !$ call hdf5Access%set()
         hdf5ReadScope: block
           type(hdf5File  ) :: file
           file=hdf5File(fileName,readOnly=.true.)
           call file%readDataset('radii'                      ,self%accelerationRadii          )
           call file%readDataset('heights'                    ,self%accelerationHeights        )
           call file%readDataset('accelerationRadial'         ,self%accelerationRadial         )
           call file%readDataset('accelerationVertical'       ,self%accelerationVertical       )
           call file%readDataset('tidalTensorRadialRadial'    ,self%tidalTensorRadialRadial    )
           call file%readDataset('tidalTensorVerticalVertical',self%tidalTensorVerticalVertical)
           call file%readDataset('tidalTensorCross'           ,self%tidalTensorCross           )
         end block hdf5ReadScope
         !$ call hdf5Access%unset()
      else
         ! Generate grid in radius and height.
         countRadii=int(log10(radiusMaximum/radiusMinimum)*radiiPerDecade)+1
         allocate(self%accelerationRadii          (countRadii           ))
         allocate(self%accelerationHeights        (countRadii           ))
         allocate(self%accelerationRadial         (countRadii,countRadii))
         allocate(self%accelerationVertical       (countRadii,countRadii))
         allocate(self%tidalTensorRadialRadial    (countRadii,countRadii))
         allocate(self%tidalTensorVerticalVertical(countRadii,countRadii))
         allocate(self%tidalTensorCross           (countRadii,countRadii))
         self%accelerationRadii  =Make_Range(radiusMinimum,radiusMaximum,countRadii,rangeTypeLogarithmic)
         self%accelerationHeights=Make_Range(radiusMinimum,radiusMaximum,countRadii,rangeTypeLogarithmic)
         ! Compute the dimensionless scale height (in units of the radial scale length), as the tabulation is built for a
         ! dimensionless disk.
         scaleHeight=+self%scaleHeight &
              &      /self%radiusScale
         ! Tabulate the vertical kernel, and its derivatives, as a function of wavenumber at each height. The kernel depends only
         ! on wavenumber and height, so tabulating it once here avoids recomputing it (each evaluation of which requires a
         ! numerical integration) for every radius, and for every component of the acceleration and tidal tensor.
         countWavenumbers                   =int(log10(wavenumberKernelMaximum/wavenumberKernelMinimum)*wavenumberKernelPerDecade)+1
         wavenumberKernelLogarithmicMinimum =log(wavenumberKernelMinimum)
         wavenumberKernelLogarithmicInterval=log(wavenumberKernelMaximum/wavenumberKernelMinimum)/dble(countWavenumbers-1)
         allocate(kernelTable(countWavenumbers,countRadii,3))
         !$omp parallel private(iWavenumber)
         allocate(integratorKernel         )
         allocate(integratorKernelGradient )
         allocate(integratorKernelGradient2)
         integratorKernel         =integrator(kernelIntegrand         ,toleranceAbsolute=1.0d-12,toleranceRelative=1.0d-6)
         integratorKernelGradient =integrator(kernelGradientIntegrand ,toleranceAbsolute=1.0d-12,toleranceRelative=1.0d-6)
         integratorKernelGradient2=integrator(kernelGradient2Integrand,toleranceAbsolute=1.0d-12,toleranceRelative=1.0d-6)
         !$omp do schedule(dynamic)
         do iHeight=1,countRadii
            height=self%accelerationHeights(iHeight)
            do iWavenumber=1,countWavenumbers
               wavenumberLow=exp(wavenumberKernelLogarithmicMinimum+dble(iWavenumber-1)*wavenumberKernelLogarithmicInterval)
               kernelTable(iWavenumber,iHeight,1)=kernelVerticalDirect         (wavenumberLow)
               kernelTable(iWavenumber,iHeight,2)=kernelVerticalGradientDirect (wavenumberLow)
               kernelTable(iWavenumber,iHeight,3)=kernelVerticalGradient2Direct(wavenumberLow)
            end do
         end do
         !$omp end do
         deallocate(integratorKernel         )
         deallocate(integratorKernelGradient )
         deallocate(integratorKernelGradient2)
         !$omp end parallel
         ! Iterate over radii and heights.
         call displayIndent("tabulating gravitational accelerations for exponential disk",verbosityLevelWorking)
         countWork=0
         do iRadius=1,countRadii
            radius=self%accelerationRadii(iRadius)
            !$omp parallel
            allocate(integratorAccelerationRadial         )
            allocate(integratorAccelerationVertical       )
            allocate(integratorTidalTensorRadialRadial    )
            allocate(integratorTidalTensorVerticalVertical)
            allocate(integratorTidalTensorCross           )
            allocate(integratorKernel                     )
            allocate(integratorKernelGradient             )
            allocate(integratorKernelGradient2            )
            integratorAccelerationRadial         =integrator(accelerationRadialIntegrand         ,toleranceAbsolute=1.0d-06,toleranceRelative=1.0d-3)
            integratorAccelerationVertical       =integrator(accelerationVerticalIntegrand       ,toleranceAbsolute=1.0d-06,toleranceRelative=1.0d-3)
            integratorTidalTensorRadialRadial    =integrator(tidalTensorRadialRadialIntegrand    ,toleranceAbsolute=1.0d-06,toleranceRelative=1.0d-3)
            integratorTidalTensorVerticalVertical=integrator(tidalTensorVerticalVerticalIntegrand,toleranceAbsolute=1.0d-06,toleranceRelative=1.0d-3)
            integratorTidalTensorCross           =integrator(tidalTensorCrossIntegrand           ,toleranceAbsolute=1.0d-06,toleranceRelative=1.0d-3)
            integratorKernel                     =integrator(kernelIntegrand                     ,toleranceAbsolute=1.0d-12,toleranceRelative=1.0d-6)
            integratorKernelGradient             =integrator(kernelGradientIntegrand             ,toleranceAbsolute=1.0d-12,toleranceRelative=1.0d-6)
            integratorKernelGradient2            =integrator(kernelGradient2Integrand            ,toleranceAbsolute=1.0d-12,toleranceRelative=1.0d-6)
            !$omp do schedule(dynamic)
            do iHeight=1,countRadii
               !$omp atomic
               countWork=countWork+1
               call displayCounter(int(100.0d0*dble(countWork)/dble(countRadii**2)),iRadius == 1 .and. iHeight == 1,verbosityLevelWorking)
               height     =self%accelerationHeights(iHeight)
               indexHeight=                         iHeight
               ! Evaluate the integral for the radial component of acceleration.
               self%accelerationRadial(iRadius,iHeight)=0.0d0
               wavenumberHigh                          =0.0d0
               iBesselZero                             =0
               converged                               =.false.
               do while (.not.converged)
                  iBesselZero      =+                        iBesselZero  &
                       &            +                        1
                  wavenumberLow    =+wavenumberHigh
                  wavenumberHigh   =+Bessel_Function_J1_Zero(iBesselZero) &
                       &            /radius        
                  accelerationDelta=+integratorAccelerationRadial%integrate(wavenumberLow,wavenumberHigh)
                  converged=abs(accelerationDelta) < 1.0d-6*abs(self%accelerationRadial(iRadius,iHeight))
                  self%accelerationRadial(iRadius,iHeight)=+self%accelerationRadial(iRadius,iHeight) &
                       &                                   +     accelerationDelta
               end do
               ! Evaluate the integral for the vertical component of acceleration.
               self%accelerationVertical(iRadius,iHeight)=0.0d0
               wavenumberHigh                            =0.0d0
               iBesselZero                               =0
               converged                                 =.false.
               do while (.not.converged)
                  iBesselZero      =+                       iBesselZero   &
                       &            +                       1
                  wavenumberLow    =+wavenumberHigh
                  wavenumberHigh   =+Bessel_Function_J0_Zero(iBesselZero) &
                       &            /radius
                  accelerationDelta=+integratorAccelerationVertical%integrate(wavenumberLow,wavenumberHigh)
                  converged=abs(accelerationDelta) < 1.0d-6*abs(self%accelerationVertical(iRadius,iHeight))
                  self%accelerationVertical(iRadius,iHeight)=+self%accelerationVertical(iRadius,iHeight) &
                       &                                     +     accelerationDelta
               end do
               ! Evaluate the integral for the radial component of the tidal tensor.
               self%tidalTensorRadialRadial(iRadius,iHeight)=0.0d0
               do besselOrder=0,2,2                
                  tidalTensorRadialRadial=0.0d0
                  wavenumberHigh         =0.0d0
                  iBesselZero            =0
                  converged              =.false.
                  do while (.not.converged)
                     iBesselZero  =+iBesselZero &
                          &        +1
                     wavenumberLow=+wavenumberHigh
                     select case (besselOrder)
                     case (0)
                        wavenumberHigh=+Bessel_Function_J0_Zero(      iBesselZero) &
                             &         /radius
                     case (2)
                        wavenumberHigh=+Bessel_Function_Jn_Zero(2.0d0,iBesselZero) &
                             &         /radius
                     case default
                        call Error_Report('incorrect Bessel function order'//{introspection:location})
                     end select
                     tidalTensorDelta=+integratorTidalTensorRadialRadial%integrate(wavenumberLow,wavenumberHigh)
                     converged=abs(tidalTensorDelta) < 1.0d-6*abs(tidalTensorRadialRadial)                   
                     tidalTensorRadialRadial=+tidalTensorRadialRadial &
                          &                  +     tidalTensorDelta
                  end do
                  self%tidalTensorRadialRadial(iRadius,iHeight)=+self%tidalTensorRadialRadial(iRadius,iHeight) &
                       &                                        +     tidalTensorRadialRadial
               end do
               ! Evaluate the integral for the vertical-vertical component of the tidal tensor.
               self%tidalTensorVerticalVertical(iRadius,iHeight)=0.0d0
               wavenumberHigh                                   =0.0d0
               iBesselZero                                      =0
               converged                                        =.false.
               do while (.not.converged)
                  iBesselZero     =+                        iBesselZero  &
                       &           +                        1
                  wavenumberLow   =+wavenumberHigh
                  wavenumberHigh  =+Bessel_Function_J0_Zero(iBesselZero) &
                       &           /radius
                  tidalTensorDelta=+integratorTidalTensorVerticalVertical%integrate(wavenumberLow,wavenumberHigh)
                  converged=abs(tidalTensorDelta) < 1.0d-6*abs(self%tidalTensorVerticalVertical(iRadius,iHeight))
                  self%tidalTensorVerticalVertical(iRadius,iHeight)=+self%tidalTensorVerticalVertical(iRadius,iHeight) &
                       &                                            +     tidalTensorDelta
               end do
               ! Evaluate the integral for the cross component of the tidal tensor.
               self%tidalTensorCross(iRadius,iHeight)=0.0d0
               wavenumberHigh                        =0.0d0
               iBesselZero                           =0
               converged                             =.false.
               do while (.not.converged)
                  iBesselZero     =+                        iBesselZero  &
                       &           +                        1
                  wavenumberLow   =+wavenumberHigh
                  wavenumberHigh  =+Bessel_Function_J1_Zero(iBesselZero) &
                       &           /radius
                  tidalTensorDelta=+integratorTidalTensorCross%integrate(wavenumberLow,wavenumberHigh)
                  converged=abs(tidalTensorDelta) < 1.0d-6*abs(self%tidalTensorCross(iRadius,iHeight))
                  self%tidalTensorCross(iRadius,iHeight)=+self%tidalTensorCross(iRadius,iHeight) &
                       &                                 +     tidalTensorDelta
               end do
            end do
            !$omp end do
            deallocate(integratorAccelerationRadial         )
            deallocate(integratorAccelerationVertical       )
            deallocate(integratorTidalTensorRadialRadial    )
            deallocate(integratorTidalTensorVerticalVertical)
            deallocate(integratorTidalTensorCross           )
            deallocate(integratorKernel                     )
            deallocate(integratorKernelGradient             )
            deallocate(integratorKernelGradient2            )
            !$omp end parallel
         end do
         call displayCounterClear(       verbosityLevelWorking)
         call displayUnindent     ("done",verbosityLevelWorking)
         !$ call hdf5Access%set()
         hdf5WriteScope: block
           type(hdf5File  ) :: file
           file=hdf5File(fileName,overWrite=.true.,readOnly=.false.)
           call file%writeDataset  (self%accelerationRadii          ,'radii'                      )
           call file%writeDataset  (self%accelerationHeights        ,'heights'                    )
           call file%writeDataset  (self%accelerationRadial         ,'accelerationRadial'         )
           call file%writeDataset  (self%accelerationVertical       ,'accelerationVertical'       )
           call file%writeDataset  (self%tidalTensorRadialRadial    ,'tidalTensorRadialRadial'    )
           call file%writeDataset  (self%tidalTensorVerticalVertical,'tidalTensorVerticalVertical')
           call file%writeDataset  (self%tidalTensorCross           ,'tidalTensorCross'           )
           call file%writeAttribute(     fileFormatVersionCurrent   ,'fileFormat'                 )
         end block hdf5WriteScope
         !$ call hdf5Access%unset()
      end if
      call File_Unlock(fileLock)
      ! Compute factors needed for interpolation.
      self%accelerationRadiusMinimumLog     =      log(self%accelerationRadii  (                            1 ))
      self%accelerationRadiusMaximumLog     =      log(self%accelerationRadii  (size(self%accelerationRadii  )))
      self%accelerationHeightMinimumLog     =      log(self%accelerationHeights(                            1 ))
      self%accelerationHeightMaximumLog     =      log(self%accelerationHeights(size(self%accelerationHeights)))
      self%accelerationRadiusInverseInterval=1.0d0/log(self%accelerationRadii  (2)/self%accelerationRadii  (1))
      self%accelerationHeightInverseInterval=1.0d0/log(self%accelerationHeights(2)/self%accelerationHeights(1))
      ! Record that the acceleration table is initialized.
      self%accelerationInitialized=.true.
    end block accelerationTabulate
    return

  contains

    double precision function accelerationRadialIntegrand(wavenumber)
      !!{RST
      Integrand for the radial component of the acceleration.
      !!}
      use :: Bessel_Functions, only : Bessel_Function_J1
      implicit none
      double precision, intent(in   ) :: wavenumber

      accelerationRadialIntegrand=-                   wavenumber            &
           &                      *Bessel_Function_J1(wavenumber*radius)    &
           &                      *kernelVertical    (wavenumber       )    &
           &                      /(                                        &
           &                        +                 1.0d0                 &
           &                        +                 wavenumber        **2 &
           &                       )**1.5d0
      return
    end function accelerationRadialIntegrand

    double precision function accelerationVerticalIntegrand(wavenumber)
      !!{RST
      Integrand for the vertical component of the acceleration.
      !!}
      use :: Bessel_Functions, only : Bessel_Function_J0
      implicit none
      double precision, intent(in   ) :: wavenumber

      accelerationVerticalIntegrand=+Bessel_Function_J0    (wavenumber*radius)    &
           &                        *kernelVerticalGradient(wavenumber       )    &
           &                        /(                                            &
           &                          +                     1.0d0                 &
           &                          +                     wavenumber        **2 &
           &                         )**1.5d0
      return
    end function accelerationVerticalIntegrand

    double precision function tidalTensorRadialRadialIntegrand(wavenumber)
      !!{RST
      Integrand for the of the :math:`\partial^2 \Phi \over \partial R^2` component of the tidal tensor.
      !!}
      use :: Bessel_Functions, only : Bessel_Function_J0, Bessel_Function_Jn
      use :: Error           , only : Error_Report
      implicit none
      double precision, intent(in   ) :: wavenumber
      double precision                :: besselFunction

      select case (besselOrder)
      case (0)
         besselFunction=+Bessel_Function_J0(  wavenumber*radius)
      case (2)
         besselFunction=-Bessel_Function_Jn(2,wavenumber*radius)
      case default
         besselFunction=+0.0d0
         call Error_Report('incorrect Bessel function order'//{introspection:location})
      end select
      tidalTensorRadialRadialIntegrand=-0.5d0                         &
           &                           *              wavenumber  **2 &
           &                           *besselFunction                &
           &                           *kernelVertical(wavenumber)    &
           &                           /(                             &
           &                             +            1.0d0           &
           &                             +            wavenumber  **2 &
           &                            )**1.5d0
      return
    end function tidalTensorRadialRadialIntegrand

    double precision function tidalTensorCrossIntegrand(wavenumber)
      !!{RST
      Integrand for the of the :math:`\partial^2 \Phi \over \partial R \partial z` component of the tidal tensor.
      !!}
      use :: Bessel_Functions, only : Bessel_Function_J1
      implicit none
      double precision, intent(in   ) :: wavenumber

      tidalTensorCrossIntegrand=-wavenumber                                   &
           &                    *Bessel_Function_J1    (wavenumber*radius)    &
           &                    *kernelVerticalGradient(wavenumber       )    &
           &                    /(                                            &
           &                      +                     1.0d0                 &
           &                      +                     wavenumber        **2 &
           &                     )**1.5d0
      return
    end function tidalTensorCrossIntegrand

    double precision function tidalTensorVerticalVerticalIntegrand(wavenumber)
      !!{RST
      Integrand for the of the :math:`\partial^2 \Phi \over \partial z^2` component of the tidal tensor.
      !!}
      use :: Bessel_Functions, only : Bessel_Function_J0
      implicit none
      double precision, intent(in   ) :: wavenumber

      tidalTensorVerticalVerticalIntegrand=+Bessel_Function_J0     (wavenumber*radius)    &
           &                               *kernelVerticalGradient2(wavenumber       )    &
           &                               /(                                             &
           &                                 +                      1.0d0                 &
           &                                 +                      wavenumber        **2 &
           &                                )**1.5d0
      return
    end function tidalTensorVerticalVerticalIntegrand

    double precision function kernelVerticalDirect(wavenumber)
      !!{RST
      The vertical kernel, :math:`K(k,z)`, appearing in the expression for the potential of the disk,
      :math:`\Phi(R,z) = -2\pi G \int_0^\infty \mathrm{d}k \, J_0(kR) \tilde{\Sigma}(k) K(k,z)`, where :math:`\tilde{\Sigma}(k)`
      is the Hankel transform of the surface density. :math:`K(k,z)` is the convolution of the normalized vertical density
      profile, :math:`\rho_z(z)=\mathrm{sech}^2(z/h)/2h`, with the razor-thin kernel, :math:`\exp(-k|z|)`. Writing
      :math:`u=z/h` and :math:`\kappa=kh` this can be expressed as a single Laplace-type integral:

      .. math::

         K(k,z) = {1 \over 2} \int_0^\infty \mathrm{d}x \, \mathrm{e}^{-\kappa x} \left[ \mathrm{sech}^2(u+x) + \mathrm{sech}^2(u-x) \right],

      which is normalized such that :math:`K(0,z)=1`. This function evaluates the kernel directly, by numerical integration. During
      tabulation the kernel is instead interpolated from a pre-computed table (see ``kernelVertical``).
      !!}
      implicit none
      double precision, intent(in   ) :: wavenumber

      kernelVerticalDirect=+0.5d0                                        &
           &               *kernelIntegrate(wavenumber,integratorKernel)
      return
    end function kernelVerticalDirect

    double precision function kernelVerticalGradientDirect(wavenumber)
      !!{RST
      The :math:`z` derivative of the vertical kernel, :math:`K(k,z)`:

      .. math::

         {\partial K \over \partial z} = {k \over 2} \int_0^\infty \mathrm{d}x \, \mathrm{e}^{-\kappa x} \left[ \mathrm{sech}^2(u+x) - \mathrm{sech}^2(u-x) \right].
      !!}
      implicit none
      double precision, intent(in   ) :: wavenumber

      kernelVerticalGradientDirect=+0.5d0                                                &
           &                       *wavenumber                                           &
           &                       *kernelIntegrate(wavenumber,integratorKernelGradient)
      return
    end function kernelVerticalGradientDirect

    double precision function kernelVerticalGradient2Direct(wavenumber)
      !!{RST
      The second :math:`z` derivative of the vertical kernel, :math:`K(k,z)`, found by differentiating the integrand of
      :math:`\partial K/\partial z` with respect to :math:`u`:

      .. math::

         {\partial^2 K \over \partial z^2} = {k \over 2h} \int_0^\infty \mathrm{d}x \, \mathrm{e}^{-\kappa x} \left[ S^\prime(u+x) - S^\prime(u-x) \right],

      where :math:`S(y)=\mathrm{sech}^2 y`. This is equivalent to :math:`k^2 K - 2 k \rho_z(z)` (which follows from
      :math:`(\partial^2/\partial z^2 - k^2)K = -2k\rho_z`), but avoids the catastrophic cancellation between those two terms
      which occurs for :math:`\kappa \gg 1`.
      !!}
      implicit none
      double precision, intent(in   ) :: wavenumber

      kernelVerticalGradient2Direct=+0.5d0                                                 &
           &                        *wavenumber                                            &
           &                        /scaleHeight                                           &
           &                        *kernelIntegrate(wavenumber,integratorKernelGradient2)
      return
    end function kernelVerticalGradient2Direct

    double precision function kernelVertical(wavenumber)
      !!{RST
      The vertical kernel, :math:`K(k,z)`, at the current height (see ``kernelVerticalDirect``), interpolated from the table
      computed at the start of the tabulation.
      !!}
      implicit none
      double precision, intent(in   ) :: wavenumber

      kernelVertical=kernelInterpolate(wavenumber,1)
      return
    end function kernelVertical

    double precision function kernelVerticalGradient(wavenumber)
      !!{RST
      The :math:`z` derivative of the vertical kernel, :math:`K(k,z)`, at the current height (see
      ``kernelVerticalGradientDirect``), interpolated from the table computed at the start of the tabulation.
      !!}
      implicit none
      double precision, intent(in   ) :: wavenumber

      kernelVerticalGradient=kernelInterpolate(wavenumber,2)
      return
    end function kernelVerticalGradient

    double precision function kernelVerticalGradient2(wavenumber)
      !!{RST
      The second :math:`z` derivative of the vertical kernel, :math:`K(k,z)`, at the current height (see
      ``kernelVerticalGradient2Direct``), interpolated from the table computed at the start of the tabulation.
      !!}
      implicit none
      double precision, intent(in   ) :: wavenumber

      kernelVerticalGradient2=kernelInterpolate(wavenumber,3)
      return
    end function kernelVerticalGradient2

    double precision function kernelInterpolate(wavenumber,indexKernel)
      !!{RST
      Interpolate the vertical kernel (``indexKernel=1``), or its first (``indexKernel=2``) or second (``indexKernel=3``)
      :math:`z` derivative, at the current height in the table computed at the start of the tabulation. Four-point Lagrange
      interpolation in :math:`\log k` is used, as this requires no state (and so is safe for concurrent use by multiple OpenMP
      threads). Outside of the tabulated range of wavenumber the kernel is evaluated directly.
      !!}
      use :: Error, only : Error_Report
      implicit none
      double precision              , intent(in   ) :: wavenumber
      integer                       , intent(in   ) :: indexKernel
      double precision, dimension(4)                :: weight
      double precision                              :: position
      integer                                       :: i

      if (wavenumber < wavenumberKernelMinimum .or. wavenumber > wavenumberKernelMaximum) then
         select case (indexKernel)
         case (1)
            kernelInterpolate=kernelVerticalDirect         (wavenumber)
         case (2)
            kernelInterpolate=kernelVerticalGradientDirect (wavenumber)
         case (3)
            kernelInterpolate=kernelVerticalGradient2Direct(wavenumber)
         case default
            kernelInterpolate=0.0d0
            call Error_Report('incorrect kernel index'//{introspection:location})
         end select
      else
         ! Find the position in the table, and the index of the second of the four points used in the interpolation.
         position=+(                                         &
              &     +log(wavenumber)                         &
              &     -    wavenumberKernelLogarithmicMinimum  &
              &    )                                         &
              &   /      wavenumberKernelLogarithmicInterval
         i       =min(max(int(position)+1,2),countWavenumbers-2)
         position=position-dble(i-1)
         ! Evaluate the Lagrange interpolation weights.
         weight  =[                                                           &
              &    -(position      )*(position-1.0d0)*(position-2.0d0)/6.0d0, &
              &    +(position+1.0d0)*(position-1.0d0)*(position-2.0d0)/2.0d0, &
              &    -(position+1.0d0)*(position      )*(position-2.0d0)/2.0d0, &
              &    +(position+1.0d0)*(position      )*(position-1.0d0)/6.0d0  &
              &   ]
         kernelInterpolate=sum(weight*kernelTable(i-1:i+2,indexHeight,indexKernel))
      end if
      return
    end function kernelInterpolate

    double precision function kernelIntegrate(wavenumber,integrator_)
      !!{RST
      Evaluate the integral over :math:`x` appearing in the vertical kernel, :math:`K(k,z)`, or its derivatives, using the
      provided integrator.

      The integrand has structure on three scales: near :math:`x=0` (where it decays as :math:`\mathrm{e}^{-(\kappa+2)x}`), and
      near :math:`x=u` (where :math:`\mathrm{sech}^2(u-x)` peaks, with width of order unity). For large :math:`\kappa` the
      structure near :math:`x=0` (and the onset of the peak at :math:`x=u`) has width :math:`1/\kappa`. Since an adaptive
      integrator can entirely miss a narrow feature lying between its initial sample points (and so report a confident, but
      incorrect, result) the range of integration is broken at each of these scales, so that every feature lies close to the
      boundary of some sub-interval.
      !!}
      implicit none
      double precision            , intent(in   ) :: wavenumber
      type            (integrator), intent(inout) :: integrator_
      ! Range (in x) beyond which sech²(u±x) has decayed by a factor e⁻⁵⁰ from its peak.
      double precision            , parameter     :: rangeSech       =25.0d0
      ! Range (in κx) beyond which exp(-κx) has decayed by a factor e⁻⁵⁰.
      double precision            , parameter     :: rangeExponential=50.0d0
      double precision            , dimension(7)  :: limits
      double precision                            :: rangeDecay
      integer                                     :: i

      wavenumberScaled=+wavenumber  &
           &           *scaleHeight
      heightScaled    =+height      &
           &           /scaleHeight
      if (wavenumberScaled*rangeSech > rangeExponential) then
         rangeDecay=+rangeExponential &
              &     /wavenumberScaled
      else
         rangeDecay=+rangeSech
      end if
      limits=[                                                         &
           &  0.0d0                                                  , &
           &  min(heightScaled,rangeSech,rangeDecay)                 , &
           &  min(heightScaled,rangeSech           )                 , &
           &  max(min(heightScaled,rangeSech),heightScaled-rangeSech), &
           &  heightScaled                                           , &
           &  heightScaled+rangeDecay                                , &
           &  heightScaled+rangeSech                                   &
           & ]
      kernelIntegrate=0.0d0
      do i=1,size(limits)-1
         if (limits(i+1) > limits(i)) kernelIntegrate=+kernelIntegrate                              &
              &                                       +integrator_%integrate(limits(i),limits(i+1))
      end do
      return
    end function kernelIntegrate

    double precision function kernelIntegrand(x)
      !!{RST
      Integrand for the vertical kernel, :math:`K(k,z)`.
      !!}
      implicit none
      double precision, intent(in   ) :: x

      kernelIntegrand=+exp(-wavenumberScaled*x)      &
           &          *(                             &
           &            +sechSquared(heightScaled+x) &
           &            +sechSquared(heightScaled-x) &
           &           )
      return
    end function kernelIntegrand

    double precision function kernelGradientIntegrand(x)
      !!{RST
      Integrand for the :math:`z` derivative of the vertical kernel, :math:`K(k,z)`.
      !!}
      implicit none
      double precision, intent(in   ) :: x

      kernelGradientIntegrand=+exp(-wavenumberScaled*x)      &
           &                  *(                             &
           &                    +sechSquared(heightScaled+x) &
           &                    -sechSquared(heightScaled-x) &
           &                   )
      return
    end function kernelGradientIntegrand

    double precision function kernelGradient2Integrand(x)
      !!{RST
      Integrand for the second :math:`z` derivative of the vertical kernel, :math:`K(k,z)`.
      !!}
      implicit none
      double precision, intent(in   ) :: x

      kernelGradient2Integrand=+exp(-wavenumberScaled*x)                &
           &                   *(                                       &
           &                     +sechSquaredDerivative(heightScaled+x) &
           &                     -sechSquaredDerivative(heightScaled-x) &
           &                    )
      return
    end function kernelGradient2Integrand

    double precision function sechSquared(y)
      !!{RST
      Evaluate :math:`\mathrm{sech}^2 y`, written in terms of :math:`\mathrm{e}^{-2|y|}` so that it can not overflow at large
      :math:`|y|`.
      !!}
      implicit none
      double precision, intent(in   ) :: y
      double precision                :: expTerm

      expTerm    =+exp(-2.0d0*abs(y))
      sechSquared=+4.0d0           &
           &      *expTerm         &
           &      /(               &
           &        +1.0d0         &
           &        +expTerm       &
           &       )**2
      return
    end function sechSquared

    double precision function sechSquaredDerivative(y)
      !!{RST
      Evaluate :math:`\mathrm{d}\,\mathrm{sech}^2 y/\mathrm{d}y = -2\,\mathrm{sech}^2 y \tanh y`, written in terms of
      :math:`\mathrm{e}^{-2|y|}` so that it can not overflow at large :math:`|y|`.
      !!}
      implicit none
      double precision, intent(in   ) :: y
      double precision                :: expTerm

      expTerm              =+exp(-2.0d0*abs(y))
      sechSquaredDerivative=-2.0d0          &
           &                *sechSquared(y) &
           &                *sign(1.0d0,y)  &
           &                *(              &
           &                  +1.0d0        &
           &                  -expTerm      &
           &                 )              &
           &                /(              &
           &                  +1.0d0        &
           &                  +expTerm      &
           &                 )
      return
    end function sechSquaredDerivative
    
  end subroutine exponentialDiskAccelerationTabulate
  
  function exponentialDiskPositionSample(self,randomNumberGenerator_)
    !!{RST
    Sample a position from an exponential disk distribution.
    !!}
    use :: Coordinates             , only : coordinateCartesian, assignment(=)
    use :: Lambert_Ws              , only : Lambert_Wm1
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    type            (coordinateCartesian            )                :: exponentialDiskPositionSample
    class           (massDistributionExponentialDisk), intent(inout) :: self
    class           (randomNumberGeneratorClass     ), intent(inout) :: randomNumberGenerator_
    double precision                                                 :: radius                       , height, &
         &                                                              phi

    ! Select a radial coordinate.
    radius=(-1.0d0-Lambert_Wm1((-1.0d0+      randomNumberGenerator_%uniformSample())/exp(1.0d0)))*self%radiusScale
    ! Select a vertical coordinate.
    height=(      -atanh      ( +1.0d0-2.0d0*randomNumberGenerator_%uniformSample()            ))*self%scaleHeight
    ! Angular coordinate is uniformly distributed between 0 and 2π.
    phi   =+                                 randomNumberGenerator_%uniformSample()              *2.0d0*Pi
    ! Return Cartesian coordinates.
    exponentialDiskPositionSample=[radius*cos(phi),radius*sin(phi),height]
    return
  end function exponentialDiskPositionSample
