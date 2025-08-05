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
  Implementation of an exponential disk mass distribution class.
  !!}
  
  !$ use :: OMP_Lib, only : omp_lock_kind
  use    :: Tables , only : table1DLogarithmicLinear

  !![
  <massDistribution name="massDistributionExponentialDisk">
   <description>The exponential disk mass distribution: $\rho(r,z)=\rho_0 \exp(-r/r_\mathrm{s}) \hbox{sech}^2(z/z_\mathrm{s})$.</description>
  </massDistribution>
  !!]
  type, public, extends(massDistributionCylindrical) :: massDistributionExponentialDisk
     !!{
     The exponential disk mass distribution: $\rho(r,z)=\rho_0 \exp(-r/r_\mathrm{s}) \hbox{sech}^2(z/z_\mathrm{s})$.
     !!}
     private
     double precision                                                        :: scaleRadius                                   , scaleHeight                                   , &
          &                                                                     densityNormalization                          , surfaceDensityNormalization                   , &
          &                                                                     mass
     ! Tables used for rotation curves and potential.
     logical                                                                 :: scaleLengthFactorSet                          , rotationCurveInitialized              =.false., &
          &                                                                     rotationCurveGradientInitialized      =.false., potentialInitialized                  =.false., &
          &                                                                     accelerationInitialized               =.false.
     double precision                                                        :: scaleLengthFactor
     double precision                                                        :: rotationCurveHalfRadiusMinimum                , rotationCurveHalfRadiusMaximum
     double precision                                                        :: rotationCurveGradientHalfRadiusMinimum        , rotationCurveGradientHalfRadiusMaximum
     type            (table1DLogarithmicLinear)                              :: rotationCurveTable                            , rotationCurveGradientTable                    , &
          &                                                                     potentialTable
     double precision                          , allocatable, dimension(:  ) :: accelerationRadii                             , accelerationHeights
     double precision                          , allocatable, dimension(:,:) :: accelerationRadial                            , accelerationVertical                          , &
          &                                                                     tidalTensorRadialRadial                       , tidalTensorVerticalVertical                   , &
          &                                                                     tidalTensorCross
     double precision                                                        :: accelerationRadiusMinimumLog                  , accelerationRadiusMaximumLog                  , &
          &                                                                     accelerationHeightMinimumLog                  , accelerationHeightMaximumLog                  , &
          &                                                                     accelerationRadiusInverseInterval             , accelerationHeightInverseInterval
     ! Locks.
     !$ integer      (omp_lock_kind           )                              :: factorComputeLock                             , rotationCurveLock                             , &
     !$   &                                                                     rotationCurveGradientLock                     , potentialLock
   contains
     !![
     <methods>
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
     procedure :: densitySphericalAverage                 => exponentialDiskDensitySphericalAverage
     procedure :: surfaceDensity                          => exponentialDiskSurfaceDensity
     procedure :: massEnclosedBySphere                    => exponentialDiskMassEnclosedBySphere
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
     !!{
     Constructors for the \refClass{massDistributionExponentialDisk} mass distribution class.
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

contains

  function exponentialDiskConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{massDistributionExponentialDisk} mass distribution class which builds the object from a parameter
    set.
    !!}
    use :: Input_Parameters          , only : inputParameter                , inputParameters
    use :: Galactic_Structure_Options, only : enumerationComponentTypeEncode, enumerationMassTypeEncode
    implicit none
    type            (massDistributionExponentialDisk)                :: self
    type            (inputParameters                ), intent(inout) :: parameters
    double precision                                                 :: mass         , scaleRadius, &
         &                                                              scaleHeight
    logical                                                          :: dimensionless
    type            (varying_string                 )                :: componentType
    type            (varying_string                 )                :: massType

    !![
    <inputParameter>
      <name>scaleHeight</name>
      <defaultSource>\citep{kregel_flattening_2002}</defaultSource>
      <defaultValue>0.137d0</defaultValue>
      <description>The scale height of the exponential disk profile.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>scaleRadius</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The scale radius of the exponential disk profile.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>mass</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The mass of the exponential disk profile.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>dimensionless</name>
      <defaultValue>.true.</defaultValue>
      <description>If true the exponential disk profile is considered to be dimensionless.</description>
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
     <call>self=massDistributionExponentialDisk(scaleHeight=scaleHeight,componentType=enumerationComponentTypeEncode(componentType,includesPrefix=.false.),massType=enumerationMassTypeEncode(massType,includesPrefix=.false.){conditions})</call>
     <argument name="mass"          value="mass"          parameterPresent="parameters"/>
     <argument name="scaleRadius"   value="scaleRadius"   parameterPresent="parameters"/>
     <argument name="dimensionless" value="dimensionless" parameterPresent="parameters"/>
    </conditionalCall>
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function exponentialDiskConstructorParameters

  function exponentialDiskConstructorInternal(scaleRadius,scaleHeight,mass,dimensionless,componentType,massType) result(self)
    !!{
    Internal constructor for ``exponentialDisk'' mass distribution class.
    !!}
    use :: Error                   , only : Error_Report
    use :: Numerical_Comparison    , only : Values_Differ
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    type            (massDistributionExponentialDisk)                          :: self
    double precision                                 , intent(in   ), optional :: scaleRadius                                 , scaleHeight                                 , &
         &                                                                        mass
    logical                                          , intent(in   ), optional :: dimensionless
    type            (enumerationComponentTypeType   ), intent(in   ), optional :: componentType
    type            (enumerationMassTypeType        ), intent(in   ), optional :: massType
    double precision                                 , parameter               :: rotationCurveHalfRadiusMaximumDefault=1.0d+1, rotationCurveHalfRadiusMinimumDefault=1.0d-6
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
       if (present(scaleRadius)) then
          if (Values_Differ(scaleRadius,1.0d0,absTol=1.0d-6)) call Error_Report('scaleRadius should be unity for a dimensionless profile (or simply do not specify a scale length)'//{introspection:location})
       end if
       if (present(mass       )) then
          if (Values_Differ(mass       ,1.0d0,absTol=1.0d-6)) call Error_Report('mass should be unity for a dimensionless profile (or simply do not specify a mass)'               //{introspection:location})
       end if
       self%scaleRadius                =1.0d0
       self%mass                       =1.0d0
       self%surfaceDensityNormalization=1.0d0/2.0d0/Pi
    else
       ! Set scale radius.
       if (.not.present(scaleRadius)) call Error_Report('scale radius must be specified for dimensionful profiles'//{introspection:location})
       if (.not.present(mass       )) call Error_Report('mass must be specified for dimensionful profiles'        //{introspection:location})
       self%scaleRadius                =scaleRadius
       self%mass                       =mass
       self%surfaceDensityNormalization=self%mass/2.0d0/Pi/self%scaleRadius**2
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
    self%rotationCurveHalfRadiusMinimum        =rotationCurveHalfRadiusMinimumDefault
    self%rotationCurveHalfRadiusMaximum        =rotationCurveHalfRadiusMaximumDefault
    self%rotationCurveGradientHalfRadiusMinimum=rotationCurveHalfRadiusMinimumDefault
    self%rotationCurveGradientHalfRadiusMaximum=rotationCurveHalfRadiusMaximumDefault
    self%scaleLengthFactorSet                  =.false.
    self%rotationCurveInitialized              =.false.
    self%rotationCurveGradientInitialized      =.false.
    self%potentialInitialized                  =.false.
    self%accelerationInitialized               =.false.
    self%scaleLengthFactor                     =0.0d0
    ! Initialize locks.
    !$ call OMP_Init_Lock(self%factorComputeLock        )
    !$ call OMP_Init_Lock(self%rotationCurveLock        )
    !$ call OMP_Init_Lock(self%rotationCurveGradientLock)
    !$ call OMP_Init_Lock(self%potentialLock            )
    return
  end function exponentialDiskConstructorInternal

  subroutine exponentialDiskDestructor(self)
    !!{
    Destructor for exponential disk mass distributions.
    !!}
    implicit none
    type(massDistributionExponentialDisk), intent(inout) :: self

    if (self%rotationCurveInitialized        ) call self%rotationCurveTable        %destroy()
    if (self%rotationCurveGradientInitialized) call self%rotationCurveGradientTable%destroy()
    if (self%potentialInitialized            ) call self%potentialTable            %destroy()
    !$ call OMP_Destroy_Lock(self%factorComputeLock        )
    !$ call OMP_Destroy_Lock(self%rotationCurveLock        )
    !$ call OMP_Destroy_Lock(self%rotationCurveGradientLock)
    !$ call OMP_Destroy_Lock(self%potentialLock            )
    return
  end subroutine exponentialDiskDestructor

  subroutine exponentialDiskTabulate(self)
    !!{
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
    !!{
    Return true indicating that this distribution has a monotonically-decreasing surface density.
    !!}
    implicit none
    class(massDistributionExponentialDisk), intent(inout) :: self

    assumeMonotonicDecreasingSurfaceDensity=.true.
    return
  end function exponentialDiskAssumeMonotonicDecreasingSurfaceDensity

  double precision function exponentialDiskMassTotal(self)
    !!{
    Return the total mass in an exponential disk distribution.
    !!}
    implicit none
    class(massDistributionExponentialDisk), intent(inout) :: self

    exponentialDiskMassTotal=self%mass
    return
  end function exponentialDiskMassTotal

  double precision function exponentialDiskRadiusHalfMass(self)
    !!{
    Return the half-mass radius in an exponential disk mass distribution.
    !!}
    implicit none
    class           (massDistributionExponentialDisk), intent(inout) :: self
    double precision                                 , parameter     :: radiusHalfMassToScaleRadius=1.678346990d0

    exponentialDiskRadiusHalfMass=+radiusHalfMassToScaleRadius &
         &                        *self%scaleRadius
    return
  end function exponentialDiskRadiusHalfMass

  double precision function exponentialDiskDensity(self,coordinates)
    !!{
    Return the density at the specified {\normalfont \ttfamily coordinates} in an exponential disk mass distribution.
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
    r=    position%r() /self%scaleRadius
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
    !!{
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
    r=    position%r() /self%scaleRadius
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
  
  double precision function exponentialDiskDensitySphericalAverage(self,radius)
    !!{
    Return the spherically-averaged density at the specified {\normalfont \ttfamily coordinates} in an exponential disk mass
    distribution. Note that this assumes the thin-disk approximation.
    !!}
    implicit none
    class           (massDistributionExponentialDisk), intent(inout) :: self
    double precision                                 , intent(in   ) :: radius

    exponentialDiskDensitySphericalAverage=+0.5d0                                 &
         &                                 *     self%surfaceDensityNormalization &
         &                                 /                               radius &
         &                                 *exp(                                  &
         &                                      -                          radius &
         &                                      /self                %scaleRadius &
         &                                     )
    return
  end function exponentialDiskDensitySphericalAverage

  double precision function exponentialDiskMassEnclosedBySphere(self,radius)
    !!{
    Computes the mass enclosed within a sphere of given {\normalfont \ttfamily radius} for exponential disk mass
    distributions. Note that this assumes the thin-disk approximation.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (massDistributionExponentialDisk), intent(inout), target :: self
    double precision                                 , intent(in   )         :: radius
    double precision                                                         :: fractionalRadius

    fractionalRadius                   =+radius                              &
         &                              /self%scaleRadius
    exponentialDiskMassEnclosedBySphere=+2.0d0                               &
         &                              *Pi                                  &
         &                              *self%scaleRadius                **2 &
         &                              *self%surfaceDensityNormalization    &
         &                              *(                                   &
         &                                +1.0d0                             &
         &                                -(                                 &
         &                                  +1.0d0                           &
         &                                  +fractionalRadius                &
         &                                 )                                 &
         &                                *exp(-fractionalRadius)            &
         &                              )
    return
  end function exponentialDiskMassEnclosedBySphere

  double precision function exponentialDiskSurfaceDensity(self,coordinates)
    !!{
    Return the surface density at the specified {\normalfont \ttfamily coordinates} in an exponential disk mass distribution.
    !!}
    use :: Coordinates, only : coordinate
    implicit none
    class           (massDistributionExponentialDisk), intent(inout) :: self
    class           (coordinate                     ), intent(in   ) :: coordinates
    double precision                                                 :: r

    ! Get the radial coordinate.
    r=coordinates%rCylindrical()/self%scaleRadius
    ! Compute the density.
    exponentialDiskSurfaceDensity=self%surfaceDensityNormalization*exp(-r)
    return
  end function exponentialDiskSurfaceDensity

  double precision function exponentialDiskRadiusEnclosingSurfaceDensity(self,densitySurface,radiusGuess) result(radius)
    !!{
    Computes the radius enclosing a given surface density for exponential disk mass distributions.
    !!}    
    implicit none
    class           (massDistributionExponentialDisk), intent(inout), target   :: self
    double precision                                 , intent(in   )           :: densitySurface
    double precision                                 , intent(in   ), optional :: radiusGuess
    
    radius=-     self%scaleRadius                 &
    &      *log(                                  &
    &           +     densitySurface              &
    &           /self%surfaceDensityNormalization &
    &          )
    return
  end function exponentialDiskRadiusEnclosingSurfaceDensity
  
  double precision function exponentialDiskRotationCurve(self,radius)
    !!{
    Return the mid-plane rotation curve for an exponential disk.
    !!}
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    implicit none
    class           (massDistributionExponentialDisk), intent(inout) :: self
    double precision                                 , intent(in   ) :: radius
    double precision                                                 :: r           , halfRadius, &
         &                                                              radiusFactor

    ! Get scale-free radius.
    r=radius/self%scaleRadius
    ! Compute rotation curve.
    if (r > radiusMaximum) then
       ! Beyond some maximum radius, approximate the disk as a spherical distribution to avoid evaluating Bessel functions for
       ! very large arguments.
       exponentialDiskRotationCurve=sqrt(self%mass/radius)
    else
       ! We are often called at precisely one scale length. Use pre-computed factors in that case.
       if (r == 1.0d0) then
          if (.not.self%scaleLengthFactorSet) then
             !$ call OMP_Set_Lock(self%factorComputeLock)
             !$ if (.not.self%scaleLengthFactorSet) then
                halfRadius               =0.5d0
                self%scaleLengthFactor   =self%besselFactorRotationCurve(halfRadius)
                self%scaleLengthFactorSet=.true.
             !$ end if
             !$ call OMP_Unset_Lock(self%factorComputeLock)
          end if
          radiusFactor=self%scaleLengthFactor
       else
          halfRadius       =0.5d0*r
          radiusFactor=self%besselFactorRotationCurve(halfRadius)
       end if
       exponentialDiskRotationCurve=sqrt(2.0d0*(self%mass/self%scaleRadius)*radiusFactor)
    end if
    ! Make dimensionful if necessary.
    if (.not.self%dimensionless) exponentialDiskRotationCurve= &
         & +sqrt(gravitationalConstant_internal)               &
         & *exponentialDiskRotationCurve
    return
  end function exponentialDiskRotationCurve

  double precision function exponentialDiskRotationCurveGradient(self,radius)
    !!{
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
         &         /self%scaleRadius
    if (2.0d0*besselArgument > fractionalRadiusMaximum) then
       ! Beyond some maximum radius, approximate the disk as a point mass to avoid evaluating Bessel functions for
       ! very large arguments.
       exponentialDiskRotationCurveGradient=-self%mass/radius**2
    else
       ! Compute the gradient.
       besselFactor=self%besselFactorRotationCurveGradient(besselArgument)
       exponentialDiskRotationCurveGradient=+self%mass           &
            &                               *besselFactor        &
            &                               /self%scaleRadius**2
    end if
    ! Make dimensionful if necessary.
    if (.not.self%dimensionless) exponentialDiskRotationCurveGradient= &
         &  +gravitationalConstant_internal                            &
         &  *exponentialDiskRotationCurveGradient
    return
  end function exponentialDiskRotationCurveGradient

  logical function exponentialDiskPotentialIsAnalytic(self) result(isAnalytic)
    !!{
    Return that the potential has an analytic form.
    !!}
    implicit none
    class(massDistributionExponentialDisk), intent(inout) :: self

    isAnalytic=.true.
    return
  end function exponentialDiskPotentialIsAnalytic

  double precision function exponentialDiskPotential(self,coordinates,status)
    !!{
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
    if (radius > potentialRadiusMaximum*self%scaleRadius) then
       exponentialDiskPotential=-self%mass/radius
    else
       ! Radius is sufficiently small to use the full calculation.
       ! Compute the potential. If the radius is lower than the height then approximate the disk
       ! mass as being spherically distributed.
       if (radius > self%scaleHeight) then
          halfRadius           =radius/2.0d0/self%scaleRadius
          correctionSmallRadius=0.0d0
       else
          halfRadius           =+self%scaleHeight/self%scaleRadius/2.0d0
          correctionSmallRadius=+self%massEnclosedBySphere(self%scaleHeight) &
               &                /self%scaleHeight
          if (radius > 0.0d0) correctionSmallRadius=+correctionSmallRadius             &
               &                                    -self%massEnclosedBySphere(radius) &
               &                                    /radius
       end if
       ! Compute the potential including the correction to small radii.
       exponentialDiskPotential=                                  &
            &             -self%mass                              &
            &             /self%scaleRadius                       &
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
    !!{
    Compute Bessel function factors appearing in the expression for an razor-thin exponential
    disk gravitational potential.
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
       !$ call OMP_Set_Lock(self%potentialLock)
       call self%tabulate()
       exponentialDiskBesselFactorPotential=self%potentialTable%interpolate(halfRadius)
       !$ call OMP_Unset_Lock(self%potentialLock)
    end if
    return
  end function exponentialDiskBesselFactorPotential

  double precision function exponentialDiskBesselFactorRotationCurve(self,halfRadius)
    !!{
    Compute Bessel function factors appearing in the expression for an razor-thin exponential disk rotation curve.
    !!}
    use :: Bessel_Functions        , only : Bessel_Function_I0, Bessel_Function_I1, Bessel_Function_K0, Bessel_Function_K1
    use :: Numerical_Constants_Math, only : eulersConstant    , ln2
    implicit none
    class           (massDistributionExponentialDisk), intent(inout) :: self
    double precision                                 , intent(in   ) :: halfRadius
    double precision                                 , parameter     :: halfRadiusSmall             =1.0d-3
    integer                                          , parameter     :: rotationCurvePointsPerDecade=100
    integer                                                          :: iPoint                             , rotationCurvePointsCount
    double precision                                                 :: x
    logical                                                          :: makeTable

    ! For small half-radii, use a series expansion for a more accurate result.
    if (halfRadius <= 0.0d0) then
       exponentialDiskBesselFactorRotationCurve=0.0d0
       return
    else if (halfRadius < halfRadiusSmall) then
       exponentialDiskBesselFactorRotationCurve=(ln2-eulersConstant-0.5d0-log(halfRadius))*halfRadius**2
       return
    end if
    !$ call OMP_Set_Lock(self%rotationCurveLock)
    if (.not.self%rotationCurveInitialized) then
       makeTable=.true.
    else
       makeTable= halfRadius < self%rotationCurveTable%x(+1) &
         &       .or.                                        &
         &        halfRadius > self%rotationCurveTable%x(-1)
    end if
    if (makeTable) then
       ! Find the minimum and maximum half-radii to tabulate.
       self%rotationCurveHalfRadiusMinimum=min(self%rotationCurveHalfRadiusMinimum,0.5d0*halfRadius)
       self%rotationCurveHalfRadiusMaximum=max(self%rotationCurveHalfRadiusMaximum,2.0d0*halfRadius)
       ! Determine how many points to tabulate.
       rotationCurvePointsCount=int(log10(self%rotationCurveHalfRadiusMaximum/self%rotationCurveHalfRadiusMinimum)*dble(rotationCurvePointsPerDecade))+1
       ! Allocate table arrays.
       call self%rotationCurveTable%destroy()
       call self%rotationCurveTable%create(self%rotationCurveHalfRadiusMinimum,self%rotationCurveHalfRadiusMaximum,rotationCurvePointsCount)
       ! Compute Bessel factors.
       do iPoint=1,rotationCurvePointsCount
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
    !$ call OMP_Unset_Lock(self%rotationCurveLock)
    return
  end function exponentialDiskBesselFactorRotationCurve

  double precision function exponentialDiskBesselFactorRotationCurveGradient(self,halfRadius)
    !!{
    Compute Bessel function factors appearing in the expression for a razor-thin exponential disk rotation curve gradient.
    !!}
    use :: Bessel_Functions        , only : Bessel_Function_I0, Bessel_Function_I1, Bessel_Function_K0, Bessel_Function_K1
    use :: Numerical_Constants_Math, only : eulersConstant    , ln2
    implicit none
    class           (massDistributionExponentialDisk), intent(inout) :: self
    double precision                                 , intent(in   ) :: halfRadius
    double precision                                 , parameter     :: halfRadiusSmall                     =1.0d-3
    double precision                                 , parameter     :: halfRadiusLarge                     =1.0d+2
    integer                                          , parameter     :: rotationCurveGradientPointsPerDecade=100
    integer                                                          :: iPoint                                      , rotationCurveGradientPointsCount
    double precision                                                 :: x

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
    !$ call OMP_Set_Lock(self%rotationCurveGradientLock)
    if     (                                                    &
         &   .not.self%rotationCurveGradientInitialized         &
         &  .or.                                                &
         &   halfRadius < self%rotationCurveGradientTable%x(+1) &
         &  .or.                                                &
         &   halfRadius > self%rotationCurveGradientTable%x(-1) &
         & ) then
       ! Find the minimum and maximum half-radii to tabulate.
       self%rotationCurveGradientHalfRadiusMinimum=min(self%rotationCurveGradientHalfRadiusMinimum,0.5d0*halfRadius)
       self%rotationCurveGradientHalfRadiusMaximum=max(self%rotationCurveGradientHalfRadiusMaximum,2.0d0*halfRadius)
       ! Determine how many points to tabulate.
       rotationCurveGradientPointsCount=int(log10(self%rotationCurveGradientHalfRadiusMaximum/self%rotationCurveGradientHalfRadiusMinimum)*dble(rotationCurveGradientPointsPerDecade))+1
       ! Allocate table arrays.
       call self%rotationCurveGradientTable%destroy()
       call self%rotationCurveGradientTable%create(self%rotationCurveGradientHalfRadiusMinimum,self%rotationCurveGradientHalfRadiusMaximum,rotationCurveGradientPointsCount)
       ! Compute Bessel factors.
       do iPoint=1,rotationCurveGradientPointsCount
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
    !$ call OMP_Unset_Lock(self%rotationCurveGradientLock)
    return
  end function exponentialDiskBesselFactorRotationCurveGradient

  double precision function exponentialDiskSurfaceDensityRadialMoment(self,moment,radiusMinimum,radiusMaximum,isInfinite)
    !!{
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
    exponentialDiskSurfaceDensityRadialMoment=(integralHigh-integralLow)*Gamma_Function(moment+1.0d0)*self%scaleRadius**(moment+1.0d0)
    return
  end function exponentialDiskSurfaceDensityRadialMoment

  function exponentialDiskAcceleration(self,coordinates)
    !!{
    Computes the gravitational acceleration at {\normalfont \ttfamily coordinates} for exponential disk mass distributions.
    !!}
    use :: Coordinates                     , only : assignment(=), coordinateCartesian           , coordinateCylindrical
    use :: Numerical_Constants_Astronomical, only : gigaYear     , gravitationalConstant_internal, megaParsec
    use :: Numerical_Constants_Prefixes    , only : kilo
    implicit none
    double precision                                 , dimension(3)  :: exponentialDiskAcceleration
    class           (massDistributionExponentialDisk), intent(inout) :: self
    class           (coordinate                     ), intent(in   ) :: coordinates
    double precision                                 , dimension(3)  :: positionCartesian
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
       exponentialDiskAcceleration(1:2)=+           accelerationRadial             &
            &                           *           positionCartesian       (1:2)  &
            &                           /           coordinatesCylindrical%r(   )
    else
       exponentialDiskAcceleration(1:2)=+0.0d0
    end if
    exponentialDiskAcceleration   (3  )=+           accelerationVertical           &
         &                              *sign(1.0d0,coordinatesCylindrical%z(   ))
    ! For dimensionful profiles apply the unit conversions and scalings.
    if (.not.self%isDimensionless())                                   &
         & exponentialDiskAcceleration=+exponentialDiskAcceleration    &
         &                             *kilo                           &
         &                             *gigaYear                       &
         &                             /megaParsec                     &
         &                             *gravitationalConstant_internal &
         &                             *self%mass                      &
         &                             /self%scaleRadius**2
    return
  end function exponentialDiskAcceleration

  function exponentialDiskTidalTensor(self,coordinates)
    !!{
    Computes the gravitational tidal tensor at {\normalfont \ttfamily coordinates} for exponential disk mass distributions.
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
    radiusCylindrical     =coordinatesCylindrical%r()/self%scaleRadius
    positionCartesian     =positionCartesian         /self%scaleRadius
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
         &                             /self%scaleRadius**3
    return
  end function exponentialDiskTidalTensor
  
  subroutine exponentialDiskAccelerationInterpolate(self,coordinatesCylindrical,accelerationRadial,accelerationVertical,tidalTensorRadialRadial,tidalTensorVerticalVertical,tidalTensorCross)
    !!{
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
    radiusCylindrical=    coordinatesCylindrical%r()/self%scaleRadius
    heightCylindrical=abs(coordinatesCylindrical%z()/self%scaleRadius)
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
            &                             *        self%scaleRadius    **2
       if (present(accelerationRadial  ))                                                                             &
            & accelerationRadial         =-                                 radiusCylindrical                         &
            &                             /             radiusSpherical**3                                            &
            &                             *        self%scaleRadius    **2
       if (present(tidalTensorVerticalVertical))                                                                      &
            & tidalTensorVerticalVertical=+(                                                                          &
            &                               -(1.0d0/    radiusSpherical**3)                                           &
            &                               +(3.0d0/    radiusSpherical**5)*                     heightCylindrical**2 &
            &                              )                                                                          &
            &                             *        self%scaleRadius    **3
       if (present(tidalTensorRadialRadial    ))                                                                      &
            & tidalTensorRadialRadial    =+(                                                                          &
            &                               -(1.0d0/    radiusSpherical**3)                                           &
            &                               +(3.0d0/    radiusSpherical**5)*radiusCylindrical**2                      &
            &                              )                                                                          &
            &                             *        self%scaleRadius    **3
       if (present(tidalTensorCross           ))                                                                      &
            & tidalTensorCross           =+(                                                                          &
            &                               +(3.0d0/    radiusSpherical**5)*radiusCylindrical   *heightCylindrical    &
            &                              )                                                                          &
            &                             *        self%scaleRadius    **3
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
    !!{
    Tabulate the acceleration (and tidal tensor) due to the exponential disk mass distribution. Uses the approach of
    \cite{kuijken_mass_1989}. The tabulation is built for a dimensionless disk.
    !!}
    use :: Bessel_Functions        , only : Bessel_Function_J0_Zero, Bessel_Function_J1_Zero, Bessel_Function_Jn_Zero
    use :: Display                 , only : displayCounter         , displayCounterClear    , displayIndent          , displayUnindent, &
          &                                 verbosityLevelWorking
    use :: File_Utilities          , only : Directory_Make         , File_Exists            , File_Lock              , File_Path      , &
          &                                 File_Unlock            , lockDescriptor
    use :: Error                   , only : Error_Report
    use :: Input_Paths             , only : inputPath              , pathTypeDataDynamic
    use :: HDF5_Access             , only : hdf5Access
    use :: IO_HDF5                 , only : hdf5Object
    use :: ISO_Varying_String      , only : char                   , operator(//)           , varying_string
    use :: Numerical_Constants_Math, only : Pi
    use :: Numerical_Integration   , only : integrator
    use :: Numerical_Ranges        , only : Make_Range             , rangeTypeLogarithmic
    implicit none
    class           (massDistributionExponentialDisk), intent(inout) :: self
    double precision                                 , parameter     :: radiusMinimum                    = 1.0d-2, radiusMaximum    =5.0d1
    double precision                                 , parameter     :: radiiPerDecade                   =30.0d+0
    double precision                                 , parameter     :: wavenumberMaximumFactor          =10.0d+0
    integer                                          , parameter     :: xi                               = 2
    type            (integrator                     ), save          :: integratorAccelerationRadial             , integratorAccelerationVertical       , &
         &                                                              integratorTidalTensorRadialRadial        , integratorTidalTensorVerticalVertical, &
         &                                                              integratorTidalTensorCross
    logical                                          , save          :: converged
    integer                                          , save          :: iBesselZero                              , besselOrder
    double precision                                 , save          :: height                                   , accelerationDelta                    , &
         &                                                              wavenumberLow                            , wavenumberHigh                       , &
         &                                                              tidalTensorDelta                         , tidalTensorRadialRadial
    !$omp threadprivate(height,accelerationDelta,tidalTensorDelta,tidalTensorRadialRadial,wavenumberLow,wavenumberHigh,iBesselZero,besselOrder,converged,integratorAccelerationRadial,integratorAccelerationVertical,integratorTidalTensorRadialRadial,integratorTidalTensorVerticalVertical,integratorTidalTensorCross)
    integer                                                          :: countRadii                               , iRadius                              , &
         &                                                              iHeight                                  , countWork
    double precision                                                 :: radius                                   , beta

    ! Return if acceleration is initialized.
    if (self%accelerationInitialized) return
    block
      type     (varying_string) :: fileName
      character(len=8         ) :: label
      type     (hdf5Object    ) :: file
      type     (lockDescriptor) :: fileLock
      
      ! Construct a file name for the table.
      write (label,'(f8.6)') self%scaleHeight/self%scaleRadius
      fileName=inputPath(pathTypeDataDynamic)// &
           &   'galacticStructure/'          // &
           &   self%objectType()             // &
           &   '_h'                          // &
           &   trim(adjustl(label))          // &
           &   '.hdf5'
      call Directory_Make(char(File_Path(char(fileName))))
      ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
      call File_Lock(char(fileName),fileLock,lockIsShared=.true.)
      if (File_Exists(fileName)) then
         !$ call hdf5Access%set()
         file=hdf5Object      (char(fileName                     ),readOnly=.true.                 )
         call file%readDataset(     'radii'                       ,self%accelerationRadii          )
         call file%readDataset(     'heights'                     ,self%accelerationHeights        )
         call file%readDataset(     'accelerationRadial'          ,self%accelerationRadial         )
         call file%readDataset(     'accelerationVertical'        ,self%accelerationVertical       )
         call file%readDataset(     'tidalTensorRadialRadial'     ,self%tidalTensorRadialRadial    )
         call file%readDataset(     'tidalTensorVerticalVertical' ,self%tidalTensorVerticalVertical)
         call file%readDataset(     'tidalTensorCross'            ,self%tidalTensorCross           )
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
         ! Compute the vertical inverse scale-height. Note that our definition of β differs slightly from that of Kuijken & Gilmore
         ! (1989). They assume a density profile in the vertical direction of the form:
         !
         !  ρ(z) = sech²(βz/ξ)
         !
         ! while we use:
         !
         !  ρ(z) = sech²(z/h)
         !
         ! where h is the scale-height. Therefore:
         !
         !  β=ξ/h
         !
         ! and we then make it dimensionless by multiplying by the radial scale length.
         beta   =+dble(xi)         &
              &  *self%scaleRadius &
              &  /self%scaleHeight
         ! Iterate over radii and heights.
         call displayIndent("tabulating gravitational accelerations for exponential disk",verbosityLevelWorking)
         countWork=0
         do iRadius=1,countRadii
            radius=self%accelerationRadii(iRadius)
            !$omp parallel
            integratorAccelerationRadial         =integrator(accelerationRadialIntegrand         ,toleranceAbsolute=1.0d-6,toleranceRelative=1.0d-3)
            integratorAccelerationVertical       =integrator(accelerationVerticalIntegrand       ,toleranceAbsolute=1.0d-6,toleranceRelative=1.0d-3)
            integratorTidalTensorRadialRadial    =integrator(tidalTensorRadialRadialIntegrand    ,toleranceAbsolute=1.0d-6,toleranceRelative=1.0d-3)
            integratorTidalTensorVerticalVertical=integrator(tidalTensorVerticalVerticalIntegrand,toleranceAbsolute=1.0d-6,toleranceRelative=1.0d-3)
            integratorTidalTensorCross           =integrator(tidalTensorCrossIntegrand           ,toleranceAbsolute=1.0d-6,toleranceRelative=1.0d-3)
            !$omp do
            do iHeight=1,countRadii
               !$omp atomic
               countWork=countWork+1
               call displayCounter(int(100.0d0*dble(countWork)/dble(countRadii**2)),iRadius == 1 .and. iHeight == 1,verbosityLevelWorking)
               height=self%accelerationHeights(iHeight)          
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
            !$omp end parallel
         end do
         call displayCounterClear(       verbosityLevelWorking)
         call displayUnindent     ("done",verbosityLevelWorking)
         !$ call hdf5Access%set()
         file=hdf5Object       (char   (fileName                        )                              ,overWrite=.true.,readOnly=.false.)
         call file%writeDataset(        self%accelerationRadii           ,'radii'                                                        )
         call file%writeDataset(        self%accelerationHeights         ,'heights'                                                      )
         call file%writeDataset(        self%accelerationRadial          ,'accelerationRadial'                                           )
         call file%writeDataset(        self%accelerationVertical        ,'accelerationVertical'                                         )
         call file%writeDataset(        self%tidalTensorRadialRadial     ,'tidalTensorRadialRadial'                                      )
         call file%writeDataset(        self%tidalTensorVerticalVertical ,'tidalTensorVerticalVertical'                                  )
         call file%writeDataset(        self%tidalTensorCross            ,'tidalTensorCross'                                             )
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
    end block
    return

  contains

    double precision function accelerationRadialIntegrand(wavenumber)
      !!{
      Integrand for the radial component of the acceleration.
      !!}
      use :: Bessel_Functions, only : Bessel_Function_J1
      implicit none
      double precision, intent(in   ) :: wavenumber

      accelerationRadialIntegrand=-                   wavenumber            &
           &                      *Bessel_Function_J1(wavenumber*radius)    &
           &                      *Iz                (wavenumber       )    &
           &                      /(                                        &
           &                        +                 1.0d0                 &
           &                        +                 wavenumber        **2 &
           &                       )**1.5d0                                 &
           &                      /self%scaleHeight                         &
           &                      /4.0d0
      return
    end function accelerationRadialIntegrand

    double precision function accelerationVerticalIntegrand(wavenumber)
      !!{
      Integrand for the radial component of the acceleration.
      !!}
      use :: Bessel_Functions, only : Bessel_Function_J0
      implicit none
      double precision, intent(in   ) :: wavenumber

      accelerationVerticalIntegrand=+Bessel_Function_J0(wavenumber*radius)    &
           &                        *dIzdz             (wavenumber       )    &
           &                        /(                                        &
           &                          +                 1.0d0                 &
           &                          +                 wavenumber        **2 &
           &                         )**1.5d0                                 &
           &                        /self%scaleHeight                         &
           &                        /4.0d0
      return
    end function accelerationVerticalIntegrand

    double precision function tidalTensorRadialRadialIntegrand(wavenumber)
      !!{
      Integrand for the of the $\partial^2 \Phi \over \partial R^2$ component of the tidal tensor.
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
      tidalTensorRadialRadialIntegrand=-0.5d0             &
           &                           *   wavenumber **2 &
           &                           *besselFunction    &
           &                           *Iz(wavenumber)    &
           &                           /(                 &
           &                             + 1.0d0          &
           &                             + wavenumber **2 &
           &                            )**1.5d0          &
           &                           /self%scaleHeight  &
           &                           /4.0d0
      return
    end function tidalTensorRadialRadialIntegrand

    double precision function tidalTensorCrossIntegrand(wavenumber)
      !!{
      Integrand for the of the $\partial^2 \Phi \over \partial R \partial z$ component of the tidal tensor.
      !!}
      use :: Bessel_Functions, only : Bessel_Function_J1
      implicit none
      double precision, intent(in   ) :: wavenumber

      tidalTensorCrossIntegrand=-wavenumber                               &
           &                    *Bessel_Function_J1(wavenumber*radius)    &
           &                    *dIzdz             (wavenumber       )    &
           &                    /(                                        &
           &                      +                 1.0d0                 &
           &                      +                 wavenumber        **2 &
           &                     )**1.5d0                                 &
           &                    /self%scaleHeight                         &
           &                    /4.0d0
      return
    end function tidalTensorCrossIntegrand

    double precision function tidalTensorVerticalVerticalIntegrand(wavenumber)
      !!{
      Integrand for the of the $\partial^2 \Phi \over \partial z^2$ component of the tidal tensor.
      !!}
      use :: Bessel_Functions, only : Bessel_Function_J0
      implicit none
      double precision, intent(in   ) :: wavenumber

      tidalTensorVerticalVerticalIntegrand=+Bessel_Function_J0(wavenumber*radius)    &
           &                               *d2Izdz2           (wavenumber       )    &
           &                               /(                                        &
           &                                 +                 1.0d0                 &
           &                                 +                 wavenumber        **2 &
           &                                )**1.5d0                                 &
           &                               /self%scaleHeight                         &
           &                               /4.0d0
      return
    end function tidalTensorVerticalVerticalIntegrand
    
    double precision function Iz(wavenumber)
      !!{
      $z$-dependent term appearing in the expression for the potential of the disk.
      !!}
      implicit none
      double precision, intent(in   ) :: wavenumber
      double precision                :: IzmOdd    , IzmEven, &
           &                             IzDelta
      integer                         :: m
      logical                         :: converged

      Iz       =+0.0d0
      m        =-2
      converged=.false.
      do while (.not.converged)
         m      =+m                   &
              &  +2
         IzmOdd =+Izm(wavenumber,m-1)
         IzmEven=+Izm(wavenumber,m  )
         IzDelta=+IzmOdd              &
              &  +IzmEven
         Iz     =+Iz                  &
              &  +IzDelta
         if (abs(IzDelta) < 1.0d-3*abs(Iz) .or. (m > 100 .and. (abs(Iz) == 0.0d0 .or. abs(IzDelta) == 0.0d0))) converged=.true.
      end do
      Iz=Iz*2.0d0**(1+xi)
      return
    end function Iz

    double precision function Izm(wavenumber,m)
      !!{
      Evaluate the $m$-dependent part of the $I(z)$ integral.
      !!}
      use :: Binomial_Coefficients, only : Binomial_Coefficient
      implicit none
      double precision, intent(in   ) :: wavenumber
      integer         , intent(in   ) :: m
      double precision                :: mFactor   , divisor

      if (m < 0) then
         Izm=0.0d0
      else
         mFactor=+1.0d0         &
              &  +2.0d0         &
              &  *dble(m )      &
              &  /dble(xi)
         divisor=+mFactor   **2 &
              &  *beta      **2 &
              &  -wavenumber**2
         if (divisor == 0.0d0) then
            ! Divisor is zero, but the term is still regular - use l'Hopital's rule to derive the correct value.
            Izm     =+Binomial_Coefficient(-xi,m)                         &
                 &  *(                                                    &
                 &    +beta*mFactor*height                                &
                 &    +1.0d0                                              &
                 &   )                                                    &
                 &  *exp(-wavenumber*height)                              &
                 &  /     wavenumber&
                 &/2.0d0
         else
            Izm    =+Binomial_Coefficient(-xi,m)                          &
                 &  *(                                                    &
                 &    +beta      *mFactor*exp(-wavenumber        *height) &
                 &    -wavenumber        *exp(-beta      *mFactor*height) &
                 &   )                                                    &
                 &  /divisor
         end if
      end if
      return
    end function Izm
    
    double precision function dIzdz(wavenumber)
      !!{
      $z$ derivative of the $z$-dependent term appearing in the expression for the potential of the disk.
      !!}
      use :: Binomial_Coefficients, only : Binomial_Coefficient
      implicit none
      double precision, intent(in   ) :: wavenumber
      double precision                :: dIzdzmOdd , dIzdzmEven, &
           &                             dIzdzDelta
      integer                         :: m
      logical                         :: converged

      dIzdz    =+0.0d0
      m        =-2
      converged=.false.
      do while (.not.converged)
         m         =+m                      &
              &     +2
         dIzdzmOdd =+dIzdzm(wavenumber,m-1)
         dIzdzmEven=+dIzdzm(wavenumber,m  )
         dIzdzDelta=+dIzdzmOdd              &
              &     +dIzdzmEven
         dIzdz     =+dIzdz                  &
              &     +dIzdzDelta
         if (abs(dIzdzDelta) < 1.0d-3*abs(dIzdz) .or. (m > 100 .and. (abs(dIzdz) == 0.0d0 .or. abs(dIzdzDelta) == 0.0d0))) converged=.true.
      end do
      dIzdz=dIzdz*2.0d0**(1+xi)
      return
    end function dIzdz

    double precision function dIzdzm(wavenumber,m)
      !!{
      Evaluate the $m$-dependent part of the $\mathrm{d}I(z)/\mathrm{d}z$ integral.
      !!}
      use :: Binomial_Coefficients, only : Binomial_Coefficient
      implicit none
      double precision, intent(in   ) :: wavenumber
      integer         , intent(in   ) :: m
      double precision                :: mFactor   , divisor

      if (m < 0) then
         dIzdzm=0.0d0
      else
         mFactor=+1.0d0         &
              &  +2.0d0         &
              &  *dble(m )      &
              &  /dble(xi)
         divisor=+mFactor   **2 &
              &  *beta      **2 &
              &  -wavenumber**2
         if (divisor == 0.0d0) then
            ! Divisor is zero, but the term is still regular - use l'Hopital's rule to derive the correct value.
            dIzdzm=-Binomial_Coefficient(-xi,m)       &
                 & *beta                              &
                 & *mFactor                           &
                 & *                height            &
                 & *exp(-wavenumber*height)           &
                 & /2.0d0
         else
            dIzdzm=-Binomial_Coefficient(-xi,m)       &
                 & *       beta                       &
                 & *                 mFactor          &
                 & *       wavenumber                 &
                 & *(                                 &
                 &   +exp(-wavenumber        *height) &
                 &   -exp(-beta      *mFactor*height) &
                 &  )                                 &
                 & /divisor
         end if
      end if
      return
    end function dIzdzm

    double precision function d2Izdz2(wavenumber)
      !!{
      $z$ $2^\mathrm{nd}$ derivative of the $z$-dependent term appearing in the expression for the potential of the disk.
      !!}
      use :: Binomial_Coefficients, only : Binomial_Coefficient
      implicit none
      double precision, intent(in   ) :: wavenumber
      double precision                :: d2Izdz2mOdd , d2Izdz2mEven, &
           &                             d2Izdz2Delta
      integer                         :: m
      logical                         :: converged

      d2Izdz2  =+0.0d0
      m        =-2
      converged=.false.
      do while (.not.converged)
         m           =+m                        &
              &       +2
         d2Izdz2mOdd =+d2Izdz2m(wavenumber,m-1)
         d2Izdz2mEven=+d2Izdz2m(wavenumber,m  )
         d2Izdz2Delta=+d2Izdz2mOdd              &
              &       +d2Izdz2mEven
         d2Izdz2     =+d2Izdz2                  &
              &       +d2Izdz2Delta
         if (abs(d2Izdz2Delta) < 1.0d-4*abs(d2Izdz2) .or. (m > 100 .and. (abs(d2Izdz2) == 0.0d0 .or. abs(d2Izdz2Delta) == 0.0d0))) converged=.true.
      end do
      d2Izdz2=d2Izdz2*2.0d0**(1+xi)
      return
    end function d2Izdz2

    double precision function d2Izdz2m(wavenumber,m)
      !!{
      Evaluate the $m$-dependent part of the $\mathrm{d}^2I(z)/\mathrm{d}z^2$ integral.
      !!}
      use :: Binomial_Coefficients, only : Binomial_Coefficient
      implicit none
      double precision, intent(in   ) :: wavenumber
      integer         , intent(in   ) :: m
      double precision                :: mFactor   , divisor

      if (m < 0) then
         d2Izdz2m=0.0d0
      else
         mFactor=+1.0d0         &
              &  +2.0d0         &
              &  *dble(m )      &
              &  /dble(xi)
         divisor=+mFactor   **2 &
              &  *beta      **2 &
              &  -wavenumber**2
         if (divisor == 0.0d0) then
            ! Divisor is zero, but the term is still regular - use l'Hopital's rule to derive the correct value.
            d2Izdz2m=+Binomial_Coefficient(-xi,m)                          &
                 &   *  beta      *mFactor                                 &
                 &   *  wavenumber**2                                      &
                 &   *exp(-wavenumber*height)                              &
                 &   /2.0d0
         else
            d2Izdz2m=+Binomial_Coefficient(-xi,m)                          &
                 &   *  beta      *mFactor                                 &
                 &   *  wavenumber                                         &
                 &   *(                                                    &
                 &     +wavenumber        *exp(-wavenumber        *height) &
                 &     -beta      *mFactor*exp(-beta      *mFactor*height) &
                 &    )                                                    &
                 &   /divisor
         end if
      end if
      return
    end function d2Izdz2m
    
  end subroutine exponentialDiskAccelerationTabulate
  
  function exponentialDiskPositionSample(self,randomNumberGenerator_)
    !!{
    Sample a position from an exponential disk distribution.
    !!}
    use :: Lambert_Ws              , only : Lambert_Wm1
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision                                 , dimension(3)  :: exponentialDiskPositionSample
    class           (massDistributionExponentialDisk), intent(inout) :: self
    class           (randomNumberGeneratorClass     ), intent(inout) :: randomNumberGenerator_
    double precision                                                 :: radius                       , height, &
         &                                                              phi

    ! Select a radial coordinate.
    radius=(-1.0d0-Lambert_Wm1((-1.0d0+      randomNumberGenerator_%uniformSample())/exp(1.0d0)))*self%scaleRadius
    ! Select a vertical coordinate.
    height=(      -atanh      ( +1.0d0-2.0d0*randomNumberGenerator_%uniformSample()            ))*self%scaleHeight
    ! Angular coordinate is uniformly distributed between 0 and 2π.
    phi   =+                                 randomNumberGenerator_%uniformSample()              *2.0d0*Pi
    ! Return Cartesian coordinates.
    exponentialDiskPositionSample=[radius*cos(phi),radius*sin(phi),height]
    return
  end function exponentialDiskPositionSample
