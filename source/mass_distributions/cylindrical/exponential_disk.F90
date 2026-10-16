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
  
  use :: Tables, only : table1DLinearCSpline, table1DLogarithmicLinear

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
          &                                                                     accelerationInitialized          =.false., rotationCurveMidplaneInitialized =.false.
     double precision                                                        :: scaleLengthFactor
     type            (table1DLogarithmicLinear)                              :: rotationCurveTable                       , rotationCurveGradientTable               , &
          &                                                                     potentialTable
     double precision                          , allocatable, dimension(:  ) :: accelerationRadii                        , accelerationHeights
     double precision                          , allocatable, dimension(:,:) :: accelerationRadial                       , accelerationVertical                     , &
          &                                                                     tidalTensorRadialRadial                  , tidalTensorVerticalVertical              , &
          &                                                                     tidalTensorCross                         , potentialTabulated
     ! Table of the logarithm of the square of the mid-plane rotation curve for a thick disk, as a function of the logarithm of
     ! radius.
     type            (table1DLinearCSpline    )                              :: rotationCurveMidplaneTable
     double precision                                                        :: rotationCurveMidplaneRadiusMinimum       , rotationCurveMidplaneSquaredMinimum
     double precision                                                        :: accelerationRadiusMinimumLog             , accelerationRadiusMaximumLog             , &
          &                                                                     accelerationHeightMinimumLog             , accelerationHeightMaximumLog             , &
          &                                                                     accelerationRadiusInverseInterval        , accelerationHeightInverseInterval
   contains
     !![
     <methods docformat="rst">
       <method description="Tabulates the potential for an exponential disk mass distribution."                                       method="tabulate"                         />
       <method description="Compute the Bessel function factor appearing in the exponential disk rotation curve."                     method="besselFactorRotationCurve"        />
       <method description="Compute the Bessel function factor appearing in the exponential disk rotation curve gradient."            method="besselFactorRotationCurveGradient"/>
       <method description="Compute the Bessel function factor appearing in the exponential disk potential."                          method="besselFactorPotential"            />
       <method description="Tabulate the gravitational acceleration, tidal tensor, and potential due to the disk."                    method="accelerationTabulate"             />
       <method description="Tabulate the mid-plane rotation curve of a thick disk."                                                   method="rotationCurveMidplaneTabulate"    />
       <method description="Interpolate in the tabulated gravitational acceleration, tidal tensor, and/or potential due to the disk." method="accelerationInterpolate"          />
       <method description="Interpolate in the tabulated mid-plane rotation curve (and its gradient) of a thick disk."                method="rotationCurveSquaredMidplane"     />
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
     procedure :: rotationCurveMidplaneTabulate           => exponentialDiskRotationCurveMidplaneTabulate
     procedure :: accelerationInterpolate                 => exponentialDiskAccelerationInterpolate
     procedure :: rotationCurveSquaredMidplane            => exponentialDiskRotationCurveSquaredMidplane
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

  ! Tabulation of the mid-plane rotation curve of a thick disk (radii in units of the disk scale length). The table extends to
  ! `radiusMaximum`, beyond which the disk is treated as a point mass.
  double precision, parameter :: radiusRotationCurveMidplaneMinimum   =1.0d-4
  integer         , parameter :: countRotationCurveMidplanePerDecade  =100

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
    self%rotationCurveMidplaneInitialized=.false.
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
    if (self%rotationCurveMidplaneInitialized) call self%rotationCurveMidplaneTable%destroy()
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
       exponentialDiskRotationCurve=sqrt( self%mass/     radius                                           )
    else if (self%scaleHeight > 0.0d0) then
       ! For a thick disk, interpolate in the tabulated mid-plane rotation curve.
       exponentialDiskRotationCurve=sqrt((self%mass/self%radiusScale)*self%rotationCurveSquaredMidplane(r))
    else
       ! For a razor-thin disk, use the analytic solution.
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
       exponentialDiskRotationCurveGradient=-self%mass                                                               &
            &                               /     radius     **2
    else if (self%scaleHeight > 0.0d0) then
       ! For a thick disk, use the gradient of the interpolant in the tabulated mid-plane rotation curve, so that the gradient
       ! is exactly consistent with the rotation curve.
       exponentialDiskRotationCurveGradient=+self%mass                                                               &
            &                               *self%rotationCurveSquaredMidplane(2.0d0*besselArgument,gradient=.true.) &
            &                               /self%radiusScale**2
    else
       ! For a razor-thin disk, use the analytic solution.
       ! Compute the gradient.
       besselFactor                        =self%besselFactorRotationCurveGradient(besselArgument)
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

  double precision function exponentialDiskRotationCurveSquaredMidplane(self,radius,gradient) result(rotationCurveSquared)
    !!{RST
    Return the square of the mid-plane rotation curve, :math:`V_\mathrm{c}^2(R)`, of a dimensionless thick disk (or, if
    ``gradient`` is true, its gradient :math:`\mathrm{d}V_\mathrm{c}^2/\mathrm{d}R`), interpolated from the tabulated
    solution. The tabulation is of :math:`\log V_\mathrm{c}^2` as a function of :math:`\log R`, interpolated with a cubic
    spline, so that the rotation curve and its gradient are continuous, and the gradient is exactly the derivative of the
    interpolated rotation curve. Below the smallest tabulated radius the disk density is effectively uniform within the
    cylinder, such that :math:`V_\mathrm{c}^2 \propto R^2`.
    !!}
    implicit none
    class           (massDistributionExponentialDisk), intent(inout)           :: self
    double precision                                 , intent(in   )           :: radius
    logical                                          , intent(in   ), optional :: gradient
    double precision                                                           :: radiusLogarithmic
    !![
    <optionalArgument name="gradient" defaultsTo=".false."/>
    !!]

    call self%rotationCurveMidplaneTabulate()
    if (radius < self%rotationCurveMidplaneRadiusMinimum) then
       if (gradient_) then
          rotationCurveSquared=+2.0d0                                         &
               &               *self%rotationCurveMidplaneSquaredMinimum      &
               &               *     radius                                   &
               &               /self%rotationCurveMidplaneRadiusMinimum   **2
       else
          rotationCurveSquared=+self%rotationCurveMidplaneSquaredMinimum      &
               &               *(                                             &
               &                 +     radius                                 &
               &                 /self%rotationCurveMidplaneRadiusMinimum     &
               &                )**2
       end if
    else
       radiusLogarithmic                  =+log(                                                    radius               )
       rotationCurveSquared               =+exp(self%rotationCurveMidplaneTable%interpolate        (radiusLogarithmic   ))
       if (gradient_) rotationCurveSquared=+                                                        rotationCurveSquared   &
            &                              *    self%rotationCurveMidplaneTable%interpolateGradient(radiusLogarithmic   )  &
            &                              /                                                        radius
    end if
    return
  end function exponentialDiskRotationCurveSquaredMidplane

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
    Return the gravitational potential for an exponential disk. For a disk of finite thickness the potential is interpolated
    from the tabulated solution (see ``accelerationTabulate``). For a razor-thin disk the analytic solution is used.
    !!}
    use :: Coordinates                     , only : assignment(=)                 , coordinateCylindrical
    use :: Galactic_Structure_Options      , only : structureErrorCodeSuccess
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    class           (massDistributionExponentialDisk  ), intent(inout), target   :: self
    class           (coordinate                       ), intent(in   )           :: coordinates
    type            (enumerationStructureErrorCodeType), intent(  out), optional :: status
    type            (coordinateCylindrical            )                          :: position
    double precision                                                             :: halfRadius  , radius, &
         &                                                                          potential

    if (present(status)) status=structureErrorCodeSuccess
    ! Get position in cylindrical coordinate system.
    position=coordinates
    ! Compute density.
    radius=position%r()
    if (self%scaleHeight > 0.0d0) then
       ! For a thick disk, interpolate in the tabulated potential.
       call self%accelerationTabulate   (                            )
       call self%accelerationInterpolate(position,potential=potential)
       exponentialDiskPotential=+self%mass        &
            &                   /self%radiusScale &
            &                   *     potential
    else if (radius > potentialRadiusMaximum*self%radiusScale) then
       ! If the radius is sufficiently large, treat the disk as a point mass.
       exponentialDiskPotential=-self%mass        &
            &                   /     radius
    else
       ! Radius is sufficiently small to use the full calculation for a razor-thin disk. Note that this potential is evaluated
       ! in the plane of the disk - it is independent of height.
       halfRadius              =+radius           &
            &                   /2.0d0            &
            &                   /self%radiusScale
       exponentialDiskPotential=-self%mass                              &
            &                   /self%radiusScale                       &
            &                   *self%besselFactorPotential(halfRadius)
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
  
  subroutine exponentialDiskAccelerationInterpolate(self,coordinatesCylindrical,accelerationRadial,accelerationVertical,tidalTensorRadialRadial,tidalTensorVerticalVertical,tidalTensorCross,potential)
    !!{RST
    Interpolate gravitational accelerations, tidal tensors, and the potential in the tabulated solutions for an exponential disk.
    Results are for a disk of unit mass and unit scale length (with :math:`G=1`), at positive height. Beyond the tabulated range the
    disk is treated as a point mass. Below the smallest tabulated radius the radial acceleration is taken to decrease linearly to
    zero, and the tidal tensor to be constant. The potential there is found by integrating that acceleration from the edge of the
    table, so that it remains consistent with the acceleration. The table includes the mid-plane. Between the mid-plane and the
    smallest non-zero height tabulated, quantities which are odd functions of :math:`z` (the vertical acceleration, and the
    cross term of the tidal tensor) are interpolated linearly in :math:`z`, while those which are even functions of :math:`z` are
    interpolated linearly in :math:`z^2`, such that their gradients vanish in the mid-plane.
    !!}
    use :: Coordinates, only : coordinateCylindrical
    implicit none
    class           (massDistributionExponentialDisk), intent(inout)            :: self
    type            (coordinateCylindrical          ), intent(in   )            :: coordinatesCylindrical
    double precision                                 , intent(  out) , optional :: accelerationRadial     , accelerationVertical       , &
         &                                                                         tidalTensorRadialRadial, tidalTensorVerticalVertical, &
         &                                                                         tidalTensorCross       , potential
    double precision                                 , dimension(0:1)           :: hRadius                , hHeightEven                , &
         &                                                                         hHeightOdd
    double precision                                                            :: radiusCylindrical      , heightCylindrical          , &
         &                                                                         radiusSpherical
    integer                                                                     :: iRadius                , iHeight

    ! Find the scale-free position.
    radiusCylindrical=    coordinatesCylindrical%r()/self%radiusScale
    heightCylindrical=abs(coordinatesCylindrical%z()/self%radiusScale)
    if     (                                                                              &
         &   radiusCylindrical > self%accelerationRadii  (size(self%accelerationRadii  )) &
         &  .or.                                                                          &
         &   heightCylindrical > self%accelerationHeights(size(self%accelerationHeights)) &
         & ) then
       ! Use point mass approximation.
       radiusSpherical=sqrt(radiusCylindrical**2+heightCylindrical**2)
       if (present(accelerationVertical       )) accelerationVertical       =-                           heightCylindrical   /radiusSpherical**3
       if (present(accelerationRadial         )) accelerationRadial         =-       radiusCylindrical                       /radiusSpherical**3
       if (present(tidalTensorVerticalVertical)) tidalTensorVerticalVertical=-1.0d0                                          /radiusSpherical**3 &
            &                                                                +3.0d0*                     heightCylindrical**2/radiusSpherical**5
       if (present(tidalTensorRadialRadial    )) tidalTensorRadialRadial    =-1.0d0                                          /radiusSpherical**3 &
            &                                                                +3.0d0*radiusCylindrical**2                     /radiusSpherical**5
       if (present(tidalTensorCross           )) tidalTensorCross           =+3.0d0*radiusCylindrical   *heightCylindrical   /radiusSpherical**5
       if (present(potential                  )) potential                  =-1.0d0                                          /radiusSpherical
    else
       ! Find interpolating factors.
       if (radiusCylindrical < self%accelerationRadii  (1)) then
          ! The tidal tensor must become independent of R as R→0 (since the radial acceleration goes to zero at R=0).
          iRadius   = 1
          hRadius   =[0.0d0,1.0d0]
       else
          hRadius(0)=+(                                        &
               &       +log(radiusCylindrical)                 &
               &       -self%accelerationRadiusMinimumLog      &
               &      )                                        &
               &     *  self%accelerationRadiusInverseInterval
          iRadius   =min(int(hRadius(0))+1              ,size(self%accelerationRadii)-1)
          hRadius(0)=+       hRadius(0) -dble(iRadius-1)
          hRadius(1)=-       hRadius(0) +1.0d0
       end if
       if (heightCylindrical < self%accelerationHeights(2)) then
          ! Between the mid-plane (the first height tabulated) and the smallest non-zero height tabulated, interpolate linearly in
          ! z for quantities that are odd functions of z, and linearly in z² for those that are even functions of z.
          iHeight       = 1
          hHeightOdd (0)=+heightCylindrical/self%accelerationHeights(2)
          hHeightEven(0)=+hHeightOdd(0)**2
       else
          ! Other heights are tabulated logarithmically.
          hHeightOdd (0)=+(                                        &
               &           +log(heightCylindrical)                 &
               &           -self%accelerationHeightMinimumLog      &
               &          )                                        &
               &         *  self%accelerationHeightInverseInterval
          iHeight       =min(int(hHeightOdd(0))+1,size(self%accelerationHeights)-2)
          hHeightOdd (0)=+    hHeightOdd(0) -dble(iHeight-1)
          iHeight       =+iHeight+1
          hHeightEven(0)=+hHeightOdd(0)
       end if
       hHeightOdd (1)=-hHeightOdd (0)+1.0d0
       hHeightEven(1)=-hHeightEven(0)+1.0d0
       ! Interpolate in tabulated solution.
       if (present(accelerationRadial         )) accelerationRadial         =interpolate(self%accelerationRadial         ,hHeightEven)
       if (present(accelerationVertical       )) accelerationVertical       =interpolate(self%accelerationVertical       ,hHeightOdd )
       if (present(tidalTensorRadialRadial    )) tidalTensorRadialRadial    =interpolate(self%tidalTensorRadialRadial    ,hHeightEven)
       if (present(tidalTensorVerticalVertical)) tidalTensorVerticalVertical=interpolate(self%tidalTensorVerticalVertical,hHeightEven)
       if (present(tidalTensorCross           )) tidalTensorCross           =interpolate(self%tidalTensorCross           ,hHeightOdd )
       if (present(potential                  )) then
          potential=interpolate(self%potentialTabulated,hHeightEven)
          ! Below the tabulated range the acceleration decreases linearly to zero. Integrate this acceleration to find the
          ! potential.
          if (radiusCylindrical < self%accelerationRadii(1))                   &
               & potential=+potential                                          &
               &           +interpolate(self%accelerationRadial,hHeightEven)   &
               &           *(                                                  &
               &             +self%accelerationRadii(1)**2                     &
               &             -     radiusCylindrical   **2                     &
               &            )                                                  &
               &           /2.0d0                                              &
               &           /  self%accelerationRadii(1)
       end if
       if (present(accelerationRadial) .and. radiusCylindrical < self%accelerationRadii(1)) accelerationRadial=accelerationRadial*radiusCylindrical/self%accelerationRadii(1)
    end if
    return

  contains

    double precision function interpolate(table,hHeight)
      !!{RST
      Perform bilinear interpolation in the given table, using the given interpolating factors in height.
      !!}
      implicit none
      double precision, intent(in   ), dimension(:,:) :: table
      double precision, intent(in   ), dimension(0:1) :: hHeight
      integer                                         :: jRadius, jHeight

      interpolate=0.0d0
      do jRadius=0,1
         do jHeight=0,1
            interpolate=+interpolate                            &
                 &      +table(iRadius+jRadius,iHeight+jHeight) &
                 &      *(1.0d0-hRadius(jRadius))               &
                 &      *(1.0d0-hHeight(jHeight))
         end do
      end do
      return
    end function interpolate

  end subroutine exponentialDiskAccelerationInterpolate

  subroutine exponentialDiskAccelerationTabulate(self)
    !!{RST
    Tabulate the acceleration, tidal tensor, and potential due to the exponential disk mass distribution (see
    ``exponentialDiskSolutionTabulate``).
    !!}
    implicit none
    class(massDistributionExponentialDisk), intent(inout) :: self

    if (.not.self%accelerationInitialized) call exponentialDiskSolutionTabulate(self,midplane=.false.)
    return
  end subroutine exponentialDiskAccelerationTabulate

  subroutine exponentialDiskRotationCurveMidplaneTabulate(self)
    !!{RST
    Tabulate the mid-plane rotation curve of the exponential disk mass distribution (see ``exponentialDiskSolutionTabulate``).
    This is tabulated separately from the acceleration, tidal tensor, and potential as it requires far fewer integrals to be
    evaluated, and is required much more frequently.
    !!}
    implicit none
    class(massDistributionExponentialDisk), intent(inout) :: self

    if (.not.self%rotationCurveMidplaneInitialized) call exponentialDiskSolutionTabulate(self,midplane=.true. )
    return
  end subroutine exponentialDiskRotationCurveMidplaneTabulate

  subroutine exponentialDiskSolutionTabulate(self,midplane)
    !!{RST
    Tabulate either the acceleration, tidal tensor, and potential due to the exponential disk mass distribution on a grid of
    radius and height, or (if ``midplane`` is true) its mid-plane rotation curve. Uses the approach of :cite:t:`kuijken_mass_1989`, expressing each quantity as a Hankel transform of
    the surface density multiplied by a vertical kernel, :math:`K(k,z)`, which accounts for the :math:`\mathrm{sech}^2` vertical
    structure of the disk. For example, the potential is

    .. math::

       \Phi(R,z) = -G M \int_0^\infty \mathrm{d}k \, {J_0(kR) K(k,z) \over (1+k^2)^{3/2}},

    in units where the disk scale length is unity, and the mid-plane rotation curve is :math:`V_\mathrm{c}^2(R) = -R
    a_R(R,0)`. The vertical kernel is evaluated by direct numerical integration (see ``kernelVertical``) rather than by the series
    expansion of :cite:t:`kuijken_mass_1989`, the terms of which do not decay with order. The tabulation is built for a
    dimensionless disk.
    !!}
    use    :: Display              , only : displayCounter              , displayCounterClear , displayIndent , displayUnindent, &
          &                                 verbosityLevelWorking
    use    :: File_Utilities       , only : Directory_Make              , File_Exists         , File_Lock     , File_Path      , &
          &                                 File_Unlock                 , lockDescriptor
    use    :: Input_Paths          , only : inputPath                   , pathTypeDataDynamic
    use    :: HDF5_Access          , only : hdf5Access
    use    :: IO_HDF5              , only : hdf5File
    use    :: ISO_Varying_String   , only : char                        , operator(//)        , varying_string
    use    :: Numerical_Integration, only : integrator
    use    :: Numerical_Ranges     , only : Make_Range                  , rangeTypeLogarithmic
    use    :: String_Handling      , only : operator(//)
    use    :: Table_Labels         , only : extrapolationTypeExtrapolate
    !$ use :: OMP_Lib              , only : OMP_In_Parallel
    implicit none
    class           (massDistributionExponentialDisk), intent(inout)                   :: self
    logical                                          , intent(in   )                   :: midplane
    double precision                                 , parameter                       :: radiusTableMinimum                 = 1.0d-2, radiusTableMaximum                   =5.0d1
    double precision                                 , parameter                       :: radiiPerDecade                     =30.0d+0
    integer                                          , parameter                       :: fileFormatVersionCurrent           = 5
    double precision                                 , parameter                       :: wavenumberKernelMinimum            = 1.0d-3, wavenumberKernelMaximum              =1.0d7
    double precision                                 , parameter                       :: wavenumberKernelPerDecade          =32.0d+0
    type            (integrator                     ), save         , allocatable      :: integratorAccelerationRadial               , integratorAccelerationVertical             , &
         &                                                                                integratorTidalTensorRadialRadial          , integratorTidalTensorVerticalVertical      , &
         &                                                                                integratorTidalTensorCross                 , integratorPotential                        , &
         &                                                                                integratorKernel                           , integratorKernelGradient                   , &
         &                                                                                integratorKernelGradient2
    integer                                          , save                            :: besselOrder                                , indexHeight
    double precision                                 , save                            :: radius                                     , height                                     , &
         &                                                                                wavenumberScaled                           , heightScaled
    !$omp threadprivate(radius,height,indexHeight,wavenumberScaled,heightScaled,besselOrder,integratorAccelerationRadial,integratorAccelerationVertical,integratorTidalTensorRadialRadial,integratorTidalTensorVerticalVertical,integratorTidalTensorCross,integratorPotential,integratorKernel,integratorKernelGradient,integratorKernelGradient2)
    integer                                                                            :: countRadii                                 , iRadius                                    , &
         &                                                                                iHeight                                    , countWork                                  , &
         &                                                                                countWavenumbers                           , iWavenumber                                , &
         &                                                                                countRotationCurve                         , iRotationCurve                             , &
         &                                                                                countHeights                               , countHeightsKernel
    double precision                                                                   :: scaleHeight                                , wavenumber                                 , &
         &                                                                                wavenumberKernelLogarithmicMinimum         , wavenumberKernelLogarithmicInterval
    double precision                                 , allocatable  , dimension(:,:,:) :: kernelTable
    double precision                                 , allocatable  , dimension(:    ) :: rotationCurveRadii                         , rotationCurveSquared                       , &
         &                                                                                heightsKernel
    logical                                                                            :: tableRead                                  , tabulateInParallel

    solutionTabulate: block
      type     (varying_string) :: fileName
      character(len=8         ) :: label
      type     (lockDescriptor) :: fileLock
      
      ! Construct a file name for the table. The file format version is included in the name so that tables built by earlier
      ! versions (which used an incorrect vertical kernel, or lack the potential and mid-plane row of the table) are never reused,
      ! while leaving those tables in place for any older executables which share the same dynamic data path.
      write (label,'(f8.6)') self%scaleHeight/self%radiusScale
      fileName=inputPath(pathTypeDataDynamic)// &
           &   'galacticStructure/'          // &
           &   self%objectType()             // &
           &   '_h'                          // &
           &   trim(adjustl(label))
      if (midplane) fileName=fileName//'_midplane'
      fileName=fileName                      // &
           &   '_v'                          // &
           &   fileFormatVersionCurrent      // &
           &   '.hdf5'
      call Directory_Make(File_Path(fileName))
      ! Read the table from file if it exists, holding a shared lock while doing so. Otherwise, obtain an exclusive lock before
      ! building and writing the table. Once that lock is obtained the file must be checked for again, as another process (e.g. a
      ! different MPI process) may have built it in the meantime. (A shared lock must not be held while writing, as concurrent
      ! writers then collide when creating the file.) Always obtain the file lock before the hdf5Access lock to avoid deadlocks
      ! between OpenMP threads.
      tableRead=.false.
      if (File_Exists(fileName)) then
         call File_Lock  (fileName,fileLock,lockIsShared=.true.)
         call tableReadFile(fileName)
         call File_Unlock(         fileLock                    )
         tableRead=.true.
      end if
      if (.not.tableRead) call File_Lock(fileName,fileLock,lockIsShared=.false.)
      if (tableRead) then
         ! Table was read from file - nothing more to do.
      else if (File_Exists(fileName)) then
         call tableReadFile(fileName)
      else
         if (midplane) then
            ! Generate grid in radius for the mid-plane rotation curve. The vertical kernel is needed in the mid-plane only.
            countRotationCurve=int(log10(radiusMaximum/radiusRotationCurveMidplaneMinimum)*dble(countRotationCurveMidplanePerDecade))+1
            allocate(rotationCurveRadii  (countRotationCurve))
            allocate(rotationCurveSquared(countRotationCurve))
            rotationCurveRadii=Make_Range(radiusRotationCurveMidplaneMinimum,radiusMaximum,countRotationCurve,rangeTypeLogarithmic)
            heightsKernel     =[0.0d0]
         else
            ! Generate grid in radius and height. The grid in height consists of the mid-plane, followed by logarithmically-spaced
            ! heights.
            countRadii  =int(log10(radiusTableMaximum/radiusTableMinimum)*radiiPerDecade)+1
            countHeights=countRadii+1
            allocate(self%accelerationRadii          (countRadii             ))
            allocate(self%accelerationHeights        (countHeights           ))
            allocate(self%accelerationRadial         (countRadii,countHeights))
            allocate(self%accelerationVertical       (countRadii,countHeights))
            allocate(self%tidalTensorRadialRadial    (countRadii,countHeights))
            allocate(self%tidalTensorVerticalVertical(countRadii,countHeights))
            allocate(self%tidalTensorCross           (countRadii,countHeights))
            allocate(self%potentialTabulated         (countRadii,countHeights))
            self%accelerationRadii                  =Make_Range(radiusTableMinimum,radiusTableMaximum,countRadii,rangeTypeLogarithmic)
            self%accelerationHeights(1             )=0.0d0
            self%accelerationHeights(2:countHeights)=Make_Range(radiusTableMinimum,radiusTableMaximum,countRadii,rangeTypeLogarithmic)
            heightsKernel                           =self%accelerationHeights
         end if
         countHeightsKernel=size(heightsKernel)
         ! Compute the dimensionless scale height (in units of the radial scale length), as the tabulation is built for a
         ! dimensionless disk.
         scaleHeight=+self%scaleHeight &
              &      /self%radiusScale
         ! Tabulate the vertical kernel, and its derivatives, as a function of wavenumber at each height. The kernel depends only on wavenumber and height, so tabulating it once
         ! here avoids recomputing it (each evaluation of which requires a numerical integration) for every radius, and for every
         ! component of the acceleration, tidal tensor, and potential.
         countWavenumbers                   =int(log10(wavenumberKernelMaximum/wavenumberKernelMinimum)*wavenumberKernelPerDecade)+1
         wavenumberKernelLogarithmicMinimum =log(wavenumberKernelMinimum)
         wavenumberKernelLogarithmicInterval=log(wavenumberKernelMaximum/wavenumberKernelMinimum)/dble(countWavenumbers-1)
         allocate(kernelTable(countWavenumbers,countHeightsKernel,3))
         ! Parallelize the tabulation only if not already within a parallel region. The tables are usually first needed during
         ! tree evolution, within the parallel region over trees, where all threads are already occupied. A nested parallel region
         ! there would start new threads, whose threadprivate data (e.g. the copies of event hooks made on entry to each parallel
         ! region) are never freed when those threads exit.
         tabulateInParallel=.true.
         !$ tabulateInParallel=.not.OMP_In_Parallel()
         !$omp parallel private(iWavenumber,wavenumber) if (tabulateInParallel)
         allocate(integratorKernel         )
         allocate(integratorKernelGradient )
         allocate(integratorKernelGradient2)
         integratorKernel         =integrator(kernelIntegrand         ,toleranceAbsolute=1.0d-12,toleranceRelative=1.0d-6)
         integratorKernelGradient =integrator(kernelGradientIntegrand ,toleranceAbsolute=1.0d-12,toleranceRelative=1.0d-6)
         integratorKernelGradient2=integrator(kernelGradient2Integrand,toleranceAbsolute=1.0d-12,toleranceRelative=1.0d-6)
         !$omp do schedule(dynamic)
         do iHeight=1,countHeightsKernel
            height=heightsKernel(iHeight)
            do iWavenumber=1,countWavenumbers
               wavenumber=exp(wavenumberKernelLogarithmicMinimum+dble(iWavenumber-1)*wavenumberKernelLogarithmicInterval)
               kernelTable(iWavenumber,iHeight,1)=kernelVerticalDirect         (wavenumber)
               kernelTable(iWavenumber,iHeight,2)=kernelVerticalGradientDirect (wavenumber)
               kernelTable(iWavenumber,iHeight,3)=kernelVerticalGradient2Direct(wavenumber)
            end do
         end do
         !$omp end do
         deallocate(integratorKernel         )
         deallocate(integratorKernelGradient )
         deallocate(integratorKernelGradient2)
         !$omp end parallel
         if (midplane) then
            ! Tabulate the mid-plane rotation curve, V_c²(R) = -R a_R(R,z=0).
            !$omp parallel if (tabulateInParallel)
            call integratorsAllocate()
            height     =heightsKernel(1)
            indexHeight=              1
            !$omp do schedule(dynamic)
            do iRotationCurve=1,countRotationCurve
               radius                              =+rotationCurveRadii(iRotationCurve)
               rotationCurveSquared(iRotationCurve)=-radius                                          &
                    &                               *hankelIntegrate(integratorAccelerationRadial,1)
            end do
            !$omp end do
            call integratorsDeallocate()
            !$omp end parallel
         else
            ! Iterate over radii and heights.
            call displayIndent("tabulating gravitational accelerations for exponential disk",verbosityLevelWorking)
            countWork=0
            do iRadius=1,countRadii
               !$omp parallel if (tabulateInParallel)
               radius=self%accelerationRadii(iRadius)
               call integratorsAllocate()
               !$omp do schedule(dynamic)
               do iHeight=1,countHeights
                  !$omp atomic
                  countWork=countWork+1
                  call displayCounter(int(100.0d0*dble(countWork)/dble(countRadii*countHeights)),iRadius == 1 .and. iHeight == 1,verbosityLevelWorking)
                  height     =self%accelerationHeights(iHeight)
                  indexHeight=                         iHeight
                  ! Evaluate the integrals for the components of the acceleration. In the mid-plane the vertical acceleration (and the
                  ! cross term of the tidal tensor) vanish by symmetry - their integrands are identically zero there, so are not
                  ! integrated.
                  self%accelerationRadial         (iRadius,iHeight)=+hankelIntegrate(integratorAccelerationRadial         ,1          )
                  if (height > 0.0d0) then
                  self%accelerationVertical       (iRadius,iHeight)=+hankelIntegrate(integratorAccelerationVertical       ,0          )
                  else
                     self%accelerationVertical    (iRadius,iHeight)=+0.0d0
                  end if
                  ! Evaluate the integrals for the components of the tidal tensor. The radial component is the sum of terms involving
                  ! Bessel functions of order 0 and 2.
                  self%tidalTensorRadialRadial    (iRadius,iHeight)=+0.0d0
                  do besselOrder=0,2,2
                     self%tidalTensorRadialRadial (iRadius,iHeight)=+self%tidalTensorRadialRadial(iRadius,iHeight)                      &
                          &                                         +hankelIntegrate(integratorTidalTensorRadialRadial    ,besselOrder)
                  end do
                  self%tidalTensorVerticalVertical(iRadius,iHeight)=+hankelIntegrate(integratorTidalTensorVerticalVertical,0          )
                  if (height > 0.0d0) then
                     self%tidalTensorCross        (iRadius,iHeight)=+hankelIntegrate(integratorTidalTensorCross           ,1          )
                  else
                     self%tidalTensorCross        (iRadius,iHeight)=+0.0d0
                  end if
                  ! Evaluate the integral for the potential.
                  self%potentialTabulated         (iRadius,iHeight)=hankelIntegrate(integratorPotential                  ,0          )
               end do
               !$omp end do
               call integratorsDeallocate()
               !$omp end parallel
            end do
            call displayCounterClear(       verbosityLevelWorking)
            call displayUnindent     ("done",verbosityLevelWorking)
         end if
         !$ call hdf5Access%set()
         hdf5WriteScope: block
           type(hdf5File) :: file
           file=hdf5File(fileName,overWrite=.true.,readOnly=.false.)
           if (midplane) then
              call file%writeDataset  (     rotationCurveRadii         ,'rotationCurveRadii'         )
              call file%writeDataset  (     rotationCurveSquared       ,'rotationCurveSquared'       )
           else
              call file%writeDataset  (self%accelerationRadii          ,'radii'                      )
              call file%writeDataset  (self%accelerationHeights        ,'heights'                    )
              call file%writeDataset  (self%accelerationRadial         ,'accelerationRadial'         )
              call file%writeDataset  (self%accelerationVertical       ,'accelerationVertical'       )
              call file%writeDataset  (self%tidalTensorRadialRadial    ,'tidalTensorRadialRadial'    )
              call file%writeDataset  (self%tidalTensorVerticalVertical,'tidalTensorVerticalVertical')
              call file%writeDataset  (self%tidalTensorCross           ,'tidalTensorCross'           )
              call file%writeDataset  (self%potentialTabulated         ,'potential'                  )
           end if
           call    file%writeAttribute(     fileFormatVersionCurrent   ,'fileFormat'                 )
         end block hdf5WriteScope
         !$ call hdf5Access%unset()
      end if
      if (.not.tableRead) call File_Unlock(fileLock)
      if (midplane) then
         ! Build the spline interpolation of the mid-plane rotation curve.
         call self%rotationCurveMidplaneTable%create  (                                                            &
              &                                        log(rotationCurveRadii(                      1 ))         , &
              &                                        log(rotationCurveRadii(size(rotationCurveRadii)))         , &
              &                                                              size(rotationCurveRadii)            , &
              &                                        extrapolationType=spread(extrapolationTypeExtrapolate,1,2)  &
              &                                       )
         call self%rotationCurveMidplaneTable%populate(log(rotationCurveSquared))
         self%rotationCurveMidplaneRadiusMinimum =rotationCurveRadii  (1)
         self%rotationCurveMidplaneSquaredMinimum=rotationCurveSquared(1)
         ! Record that the mid-plane rotation curve table is initialized.
         self%rotationCurveMidplaneInitialized=.true.
      else
         ! Compute factors needed for interpolation.
         self%accelerationRadiusMinimumLog     =      log(self%accelerationRadii  (                            1 ))
         self%accelerationRadiusMaximumLog     =      log(self%accelerationRadii  (size(self%accelerationRadii  )))
         self%accelerationHeightMinimumLog     =      log(self%accelerationHeights(                            2 ))
         self%accelerationHeightMaximumLog     =      log(self%accelerationHeights(size(self%accelerationHeights)))
         self%accelerationRadiusInverseInterval=1.0d0/log(self%accelerationRadii  (2)/self%accelerationRadii  (1))
         self%accelerationHeightInverseInterval=1.0d0/log(self%accelerationHeights(3)/self%accelerationHeights(2))
         ! Record that the acceleration table is initialized.
         self%accelerationInitialized=.true.
      end if
    end block solutionTabulate
    return

  contains

    subroutine tableReadFile(fileName)
      !!{RST
      Read the tabulated solution from file. The caller must hold a lock on the file.
      !!}
      implicit none
      type(varying_string), intent(in   ) :: fileName

      !$ call hdf5Access%set()
      hdf5ReadScope: block
        type(hdf5File) :: file
        file=hdf5File(fileName,readOnly=.true.)
        if (midplane) then
           call file%readDataset('rotationCurveRadii'         ,     rotationCurveRadii         )
           call file%readDataset('rotationCurveSquared'       ,     rotationCurveSquared       )
        else
           call file%readDataset('radii'                      ,self%accelerationRadii          )
           call file%readDataset('heights'                    ,self%accelerationHeights        )
           call file%readDataset('accelerationRadial'         ,self%accelerationRadial         )
           call file%readDataset('accelerationVertical'       ,self%accelerationVertical       )
           call file%readDataset('tidalTensorRadialRadial'    ,self%tidalTensorRadialRadial    )
           call file%readDataset('tidalTensorVerticalVertical',self%tidalTensorVerticalVertical)
           call file%readDataset('tidalTensorCross'           ,self%tidalTensorCross           )
           call file%readDataset('potential'                  ,self%potentialTabulated         )
        end if
      end block hdf5ReadScope
      !$ call hdf5Access%unset()
      return
    end subroutine tableReadFile

    subroutine integratorsAllocate()
      !!{RST
      Allocate and construct the (thread private) integrators used in the tabulation.
      !!}
      implicit none

      allocate(integratorAccelerationRadial         )
      allocate(integratorAccelerationVertical       )
      allocate(integratorTidalTensorRadialRadial    )
      allocate(integratorTidalTensorVerticalVertical)
      allocate(integratorTidalTensorCross           )
      allocate(integratorPotential                  )
      allocate(integratorKernel                     )
      allocate(integratorKernelGradient             )
      allocate(integratorKernelGradient2            )
      integratorAccelerationRadial         =integrator(accelerationRadialIntegrand         ,toleranceAbsolute=1.0d-06,toleranceRelative=1.0d-3)
      integratorAccelerationVertical       =integrator(accelerationVerticalIntegrand       ,toleranceAbsolute=1.0d-06,toleranceRelative=1.0d-3)
      integratorTidalTensorRadialRadial    =integrator(tidalTensorRadialRadialIntegrand    ,toleranceAbsolute=1.0d-06,toleranceRelative=1.0d-3)
      integratorTidalTensorVerticalVertical=integrator(tidalTensorVerticalVerticalIntegrand,toleranceAbsolute=1.0d-06,toleranceRelative=1.0d-3)
      integratorTidalTensorCross           =integrator(tidalTensorCrossIntegrand           ,toleranceAbsolute=1.0d-06,toleranceRelative=1.0d-3)
      ! The potential is integrated to a tighter tolerance than other quantities. Close to the center of the disk the potential
      ! varies by only a small fraction of its value between grid points, and so must be computed precisely to give an accurate
      ! (and monotonic) radial gradient.
      integratorPotential                  =integrator(potentialIntegrand                  ,toleranceAbsolute=1.0d-10,toleranceRelative=1.0d-7)
      integratorKernel                     =integrator(kernelIntegrand                     ,toleranceAbsolute=1.0d-12,toleranceRelative=1.0d-6)
      integratorKernelGradient             =integrator(kernelGradientIntegrand             ,toleranceAbsolute=1.0d-12,toleranceRelative=1.0d-6)
      integratorKernelGradient2            =integrator(kernelGradient2Integrand            ,toleranceAbsolute=1.0d-12,toleranceRelative=1.0d-6)
      return
    end subroutine integratorsAllocate

    subroutine integratorsDeallocate()
      !!{RST
      Destroy the (thread private) integrators used in the tabulation.
      !!}
      implicit none

      deallocate(integratorAccelerationRadial         )
      deallocate(integratorAccelerationVertical       )
      deallocate(integratorTidalTensorRadialRadial    )
      deallocate(integratorTidalTensorVerticalVertical)
      deallocate(integratorTidalTensorCross           )
      deallocate(integratorPotential                  )
      deallocate(integratorKernel                     )
      deallocate(integratorKernelGradient             )
      deallocate(integratorKernelGradient2            )
      return
    end subroutine integratorsDeallocate

    double precision function hankelIntegrate(integrator_,order) result(integral)
      !!{RST
      Evaluate the integral over wavenumber, from zero to infinity, of the integrand of the given integrator at the current
      radius. The integral is accumulated over the intervals between successive zeros of the Bessel function of the given
      ``order`` (as a function of :math:`kR`). The partial sums so obtained form a (typically alternating) sequence which
      converges slowly - for the components of the tidal tensor, whose integrands carry additional powers of :math:`k`, many
      thousands of intervals would be needed for the contribution from an interval to fall below :math:`10^{-6}` of the total.
      The limit of the sequence is therefore estimated using the :math:`\epsilon`-algorithm of :cite:t:`wynn_device_1956`
      applied to the most recent partial sums, and accumulation stops once successive estimates agree to a relative precision
      of :math:`10^{-8}`. Accumulation also stops (and the partial sum is returned) if the contribution from an interval falls
      below :math:`10^{-6}` of the total, which occurs first if the integrand decays rapidly.
      !!}
      use :: Bessel_Functions, only : Bessel_Function_J0_Zero, Bessel_Function_J1_Zero, Bessel_Function_Jn_Zero
      use :: Error           , only : Error_Report
      implicit none
      type            (integrator), intent(inout)               :: integrator_
      integer                     , intent(in   )               :: order
      integer                     , parameter                   :: countPartialSums=16     , countIntervalsMinimum=8
      double precision            , parameter                   :: toleranceDirect = 1.0d-6, toleranceAccelerated =1.0d-8
      double precision            , dimension(countPartialSums) :: partialSums
      double precision            , dimension(3               ) :: estimates
      double precision                                          :: wavenumberLow           , wavenumberHigh              , &
           &                                                       integralDelta
      integer                                                   :: iBesselZero             , countStored
      logical                                                   :: converged

      integral      =0.0d0
      wavenumberHigh=0.0d0
      iBesselZero   =0
      countStored   =0
      estimates     =huge(0.0d0)
      converged     =.false.
      do while (.not.converged)
         iBesselZero  =+iBesselZero    &
              &        +1
         wavenumberLow=+wavenumberHigh
         select case (order)
         case (0)
            wavenumberHigh=+Bessel_Function_J0_Zero(      iBesselZero) &
                 &         /radius
         case (1)
            wavenumberHigh=+Bessel_Function_J1_Zero(      iBesselZero) &
                 &         /radius
         case (2)
            wavenumberHigh=+Bessel_Function_Jn_Zero(2.0d0,iBesselZero) &
                 &         /radius
         case default
            wavenumberHigh=+0.0d0
            call Error_Report('incorrect Bessel function order'//{introspection:location})
         end select
         integralDelta=+integrator_%integrate(wavenumberLow,wavenumberHigh)
         converged    = abs(integralDelta) < toleranceDirect*abs(integral)
         integral     =+integral      &
              &        +integralDelta
         ! Store the partial sum, retaining only the most recent.
         if (countStored < countPartialSums) then
            countStored             =countStored+1
         else
            partialSums(1:countStored-1)=partialSums(2:countStored)
         end if
         partialSums(countStored)=integral
         ! Estimate the limit of the sequence of partial sums, and test for convergence of successive estimates.
         if (.not.converged .and. iBesselZero >= countIntervalsMinimum) then
            estimates=[estimates(2:3),wynnEpsilon(partialSums(1:countStored))]
            if     (                                                                         &
                 &   abs(estimates(3)-estimates(2)) < toleranceAccelerated*abs(estimates(3)) &
                 &  .and.                                                                    &
                 &   abs(estimates(3)-estimates(1)) < toleranceAccelerated*abs(estimates(3)) &
                 & ) then
               integral =estimates(3)
               converged=.true.
            end if
         end if
      end do
      return
    end function hankelIntegrate

    double precision function wynnEpsilon(sequence) result(estimate)
      !!{RST
      Estimate the limit of the given sequence using the :math:`\epsilon`-algorithm of :cite:t:`wynn_device_1956`. The estimate is
      the final element of the highest even-order column of the :math:`\epsilon` table which can be constructed. Construction
      stops early if two successive elements of a column are equal (to within rounding), as the sequence has then converged.
      !!}
      implicit none
      double precision, intent(in   ), dimension(:             ) :: sequence
      double precision               , dimension(size(sequence)) :: columnPrevious            , columnCurrent, &
           &                                                        columnNext
      double precision, parameter                                :: differenceMinimum=1.0d-280
      double precision                                           :: difference
      integer                                                    :: countSequence             , iColumn      , &
           &                                                        i
      logical                                                    :: degenerate

      countSequence =size(sequence)
      estimate      =sequence(countSequence)
      columnPrevious=0.0d0
      columnCurrent =sequence
      degenerate    =.false.
      do iColumn=1,countSequence-1
         do i=1,countSequence-iColumn
            difference=columnCurrent(i+1)-columnCurrent(i)
            if (abs(difference) <= max(epsilon(0.0d0)*max(abs(columnCurrent(i+1)),abs(columnCurrent(i))),differenceMinimum)) then
               degenerate=.true.
               exit
            end if
            columnNext(i)=columnPrevious(i+1)+1.0d0/difference
         end do
         if (degenerate) then
            ! Two successive elements of this column are equal. If this is an even-order column, this is the limit.
            if (mod(iColumn-1,2) == 0) estimate=columnCurrent(i+1)
            exit
         end if
         if (mod(iColumn,2) == 0) estimate=columnNext(countSequence-iColumn)
         columnPrevious(1:countSequence-iColumn+1)=columnCurrent(1:countSequence-iColumn+1)
         columnCurrent (1:countSequence-iColumn  )=columnNext   (1:countSequence-iColumn  )
      end do
      return
    end function wynnEpsilon

    double precision function potentialIntegrand(wavenumber)
      !!{RST
      Integrand for the gravitational potential.
      !!}
      use :: Bessel_Functions, only : Bessel_Function_J0
      implicit none
      double precision, intent(in   ) :: wavenumber

      potentialIntegrand=-Bessel_Function_J0(wavenumber*radius)    &
           &             *kernelVertical    (wavenumber       )    &
           &             /(                                        &
           &               +                 1.0d0                 &
           &               +                 wavenumber        **2 &
           &              )**1.5d0
      return
    end function potentialIntegrand

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
    
  end subroutine exponentialDiskSolutionTabulate
  
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
