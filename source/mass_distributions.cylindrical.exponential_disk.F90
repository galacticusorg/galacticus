!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

  !% Implementation of an exponential disk mass distribution class.

  !$ use :: OMP_Lib, only : omp_lock_kind
  use    :: Tables , only : table1DLogarithmicLinear

  !# <massDistribution name="massDistributionExponentialDisk">
  !#  <description>The exponential disk mass distribution: $\rho(r,z)=\rho_0 \exp(-r/r_\mathrm{s}) \hbox{sech}^2(z/z_\mathrm{s})$.</description>
  !# </massDistribution>
  type, public, extends(massDistributionCylindrical) :: massDistributionExponentialDisk
     !% The exponential disk mass distribution: $\rho(r,z)=\rho_0 \exp(-r/r_\mathrm{s}) \hbox{sech}^2(z/z_\mathrm{s})$.
     private
     double precision                           :: scaleRadius                           , scaleHeight                           , &
          &                                        densityNormalization                  , surfaceDensityNormalization           , &
          &                                        mass
     ! Tables used for rotation curves and potential.
     logical                                    :: scaleLengthFactorSet                  , rotationCurveInitialized              , &
          &                                        rotationCurveGradientInitialized      , potentialInitialized
     double precision                           :: scaleLengthFactor
     double precision                           :: rotationCurveHalfRadiusMinimum        , rotationCurveHalfRadiusMaximum
     double precision                           :: rotationCurveGradientHalfRadiusMinimum, rotationCurveGradientHalfRadiusMaximum
     type            (table1DLogarithmicLinear) :: rotationCurveTable                    , rotationCurveGradientTable            , &
          &                                        potentialTable
     ! Locks.
     !$ integer      (omp_lock_kind           ) :: factorComputeLock                    , rotationCurveLock                     , &
     !$   &                                        rotationCurveGradientLock            , potentialLock
   contains
     !@ <objectMethods>
     !@   <object>massDistributionExponentialDisk</object>
     !@   <objectMethod>
     !@     <method>tabulate</method>
     !@     <description>Tabulates the potential for an exponential disk mass distribution.</description>
     !@     <type>\void</type>
     !@     <arguments></arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>besselFactorRotationCurve</method>
     !@     <description>Compute the Bessel function factor appearing in the exponential disk rotation curve.</description>
     !@     <type>\doublezero</type>
     !@     <arguments>\doublezero\ halfRadius\argin</arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>besselFactorRotationCurveGradient</method>
     !@     <description>Compute the Bessel function factor appearing in the exponential disk rotation curve gradient.</description>
     !@     <type>\doublezero</type>
     !@     <arguments>\doublezero\ halfRadius\argin</arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>besselFactorPotential</method>
     !@     <description>Compute the Bessel function factor appearing in the exponential disk potential.</description>
     !@     <type>\doublezero</type>
     !@     <arguments>\doublezero\ halfRadius\argin</arguments>
     !@   </objectMethod>
     !@ </objectMethods>
     final     ::                                      exponentialDiskDestructor
     procedure :: tabulate                          => exponentialDiskTabulate
     procedure :: besselFactorRotationCurve         => exponentialDiskBesselFactorRotationCurve
     procedure :: besselFactorRotationCurveGradient => exponentialDiskBesselFactorRotationCurveGradient
     procedure :: besselFactorPotential             => exponentialDiskBesselFactorPotential
     procedure :: density                           => exponentialDiskDensity
     procedure :: surfaceDensity                    => exponentialDiskSurfaceDensity
     procedure :: massEnclosedBySphere              => exponentialDiskMassEnclosedBySphere
     procedure :: potential                         => exponentialDiskPotential
     procedure :: rotationCurve                     => exponentialDiskRotationCurve
     procedure :: rotationCurveGradient             => exponentialDiskRotationCurveGradient
     procedure :: radiusHalfMass                    => exponentialDiskRadiusHalfMass
     procedure :: surfaceDensityRadialMoment        => exponentialDiskSurfaceDensityRadialMoment
  end type massDistributionExponentialDisk

  interface massDistributionExponentialDisk
     !% Constructors for the {\normalfont \ttfamily exponentialDisk} mass distribution class.
     module procedure exponentialDiskConstructorParameters
     module procedure exponentialDiskConstructorInternal
  end interface massDistributionExponentialDisk

  ! The radius (in units of the disk scale length) beyond which the disk is treated as a point mass for the purposes of computing
  ! rotation curves.
  double precision, parameter :: exponentialDiskRadiusMaximum           =3.0d+1

  ! Potential tabulation.
  integer         , parameter :: exponentialDiskPotentialPointsPerDecade=10
  double precision, parameter :: exponentialDiskPotentialRadiusMinimum  =1.0d-3
  double precision, parameter :: exponentialDiskPotentialRadiusMaximum  =5.0d+2

contains

  function exponentialDiskConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily exponentialDisk} mass distribution class which builds the object from a parameter
    !% set.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (massDistributionExponentialDisk)                :: self
    type            (inputParameters                ), intent(inout) :: parameters
    double precision                                                 :: mass         , scaleRadius, &
         &                                                              scaleHeight
    logical                                                          :: dimensionless

    !# <inputParameter>
    !#   <name>scaleHeight</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultSource>\citep{kregel_flattening_2002}</defaultSource>
    !#   <defaultValue>0.137d0</defaultValue>
    !#   <description>The scale height of the exponential disk profile.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>scaleRadius</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>1.0d0</defaultValue>
    !#   <description>The scale radius of the exponential disk profile.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>mass</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>1.0d0</defaultValue>
    !#   <description>The mass of the exponential disk profile.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>dimensionless</name>
    !#   <defaultValue>.true.</defaultValue>
    !#   <cardinality>1</cardinality>
    !#   <description>If true the exponential disk profile is considered to be dimensionless.</description>
    !#   <source>parameters</source>
    !#   <type>boolean</type>
    !# </inputParameter>
    !# <conditionalCall>
    !#  <call>self=massDistributionExponentialDisk(scaleHeight=scaleHeight{conditions})</call>
    !#  <argument name="mass"          value="mass"          parameterPresent="parameters"/>
    !#  <argument name="scaleRadius"   value="scaleRadius"   parameterPresent="parameters"/>
    !#  <argument name="dimensionless" value="dimensionless" parameterPresent="parameters"/>
    !# </conditionalCall>
    !# <inputParametersValidate source="parameters"/>
    return
  end function exponentialDiskConstructorParameters

  function exponentialDiskConstructorInternal(scaleRadius,scaleHeight,mass,dimensionless) result(self)
    !% Internal constructor for ``exponentialDisk'' mass distribution class.
    use :: Galacticus_Error        , only : Galacticus_Error_Report
    use :: Numerical_Comparison    , only : Values_Differ
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    type            (massDistributionExponentialDisk)                          :: self
    double precision                                 , intent(in   ), optional :: scaleRadius                                 , scaleHeight                                 , &
         &                                                                        mass
    logical                                          , intent(in   ), optional :: dimensionless
    double precision                                 , parameter               :: rotationCurveHalfRadiusMaximumDefault=1.0d+1, rotationCurveHalfRadiusMinimumDefault=1.0d-6

    ! Determine if profile is dimensionless.
    if (present(dimensionless)) then
       self%dimensionless=dimensionless
    else
       self%dimensionless=.false.
    end if
    ! If dimensionless, then set scale length and mass to unity.
    if (self%dimensionless) then
       if (present(scaleRadius)) then
          if (Values_Differ(scaleRadius,1.0d0,absTol=1.0d-6)) call Galacticus_Error_Report('scaleRadius should be unity for a dimensionless profile (or simply do not specify a scale length)'//{introspection:location})
       end if
       if (present(mass       )) then
          if (Values_Differ(mass       ,1.0d0,absTol=1.0d-6)) call Galacticus_Error_Report('mass should be unity for a dimensionless profile (or simply do not specify a mass)'               //{introspection:location})
       end if
       self%scaleRadius                =1.0d0
       self%mass                       =1.0d0
       self%surfaceDensityNormalization=1.0d0/2.0d0/Pi
    else
       ! Set core radius.
       if (.not.present(scaleRadius)) call Galacticus_Error_Report('scale radius must be specified for dimensionful profiles'//{introspection:location})
       if (.not.present(mass       )) call Galacticus_Error_Report('mass must be specified for dimensionful profiles'        //{introspection:location})
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
    self%scaleLengthFactor                     =0.0d0
    ! Initialize locks.
    !$ call OMP_Init_Lock(self%factorComputeLock        )
    !$ call OMP_Init_Lock(self%rotationCurveLock        )
    !$ call OMP_Init_Lock(self%rotationCurveGradientLock)
    !$ call OMP_Init_Lock(self%potentialLock            )
    return
  end function exponentialDiskConstructorInternal

  subroutine exponentialDiskDestructor(self)
    !% Destructor for exponential disk mass distributions.
    implicit none
    type(massDistributionExponentialDisk), intent(inout) :: self

    call self%rotationCurveTable        %destroy()
    call self%rotationCurveGradientTable%destroy()
    call self%potentialTable            %destroy()
    !$ call OMP_Destroy_Lock(self%factorComputeLock        )
    !$ call OMP_Destroy_Lock(self%rotationCurveLock        )
    !$ call OMP_Destroy_Lock(self%rotationCurveGradientLock)
    !$ call OMP_Destroy_Lock(self%potentialLock            )
    return
  end subroutine exponentialDiskDestructor

  subroutine exponentialDiskTabulate(self)
    !% Build tables used for exponential disk mass distributions.
    use :: Bessel_Functions, only : Bessel_Function_I0    , Bessel_Function_I1, Bessel_Function_K0, Bessel_Function_K1
    use :: Table_Labels    , only : extrapolationTypeAbort
    implicit none
    class           (massDistributionExponentialDisk), intent(inout) :: self
    integer                                                          :: i   , potentialPointsCount
    double precision                                                 :: x

    ! Build table if necessary.
    if (.not.self%potentialInitialized) then
       ! Determine how many points to tabulate.
       potentialPointsCount=int(log10(exponentialDiskPotentialRadiusMaximum/exponentialDiskPotentialRadiusMinimum)*dble(exponentialDiskPotentialPointsPerDecade))+1
       ! Create the table.
       call self%potentialTable%destroy()
       call self%potentialTable%create(exponentialDiskPotentialRadiusMinimum,exponentialDiskPotentialRadiusMaximum,potentialPointsCount,extrapolationType=spread(extrapolationTypeAbort,1,2))
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

  double precision function exponentialDiskRadiusHalfMass(self)
    !% Return the half-mass radius in an exponential disk mass distribution.
    implicit none
    class           (massDistributionExponentialDisk), intent(inout) :: self
    double precision                                 , parameter     :: radiusHalfMassToScaleRadius=1.678346990d0

    exponentialDiskRadiusHalfMass=+radiusHalfMassToScaleRadius &
         &                        *self%scaleRadius
    return
  end function exponentialDiskRadiusHalfMass

  double precision function exponentialDiskDensity(self,coordinates)
    !% Return the density at the specified {\normalfont \ttfamily coordinates} in an exponential disk mass distribution.
    use :: Coordinates     , only : assignment(=)          , coordinateCylindrical
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class           (massDistributionExponentialDisk), intent(inout) :: self
    class           (coordinate                     ), intent(in   ) :: coordinates
    type            (coordinateCylindrical          )                :: position
    double precision                                 , parameter     :: coshArgumentMaximum=50.0d0
    double precision                                                 :: r                         , z, &
         &                                                              coshTerm

    ! If disk is razor thin, density is undefined.
    if (self%scaleHeight <= 0.0d0) call Galacticus_Error_Report('density undefined for razor-thin disk'//{introspection:location})
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

  double precision function exponentialDiskMassEnclosedBySphere(self,radius)
    !% Computes the mass enclosed within a sphere of given {\normalfont \ttfamily radius} for exponential disk mass distributions.
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
    !% Return the surface density at the specified {\normalfont \ttfamily coordinates} in an exponential disk mass distribution.
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

  double precision function exponentialDiskRotationCurve(self,radius)
    !% Return the mid-plane rotation curve for an exponential disk.
    use :: Numerical_Constants_Physical, only : gravitationalConstantGalacticus
    implicit none
    class           (massDistributionExponentialDisk), intent(inout) :: self
    double precision                                 , intent(in   ) :: radius
    double precision                                                 :: r           , halfRadius, &
         &                                                              radiusFactor

    ! Get scale-free radius.
    r=radius/self%scaleRadius
    ! Compute rotation curve.
    if (r > exponentialDiskRadiusMaximum) then
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
         & +sqrt(gravitationalConstantGalacticus)              &
         & *exponentialDiskRotationCurve
    return
  end function exponentialDiskRotationCurve

  double precision function exponentialDiskRotationCurveGradient(self,radius)
    !% Return the mid-plane rotation curve gradient for an exponential disk.
    use :: Numerical_Constants_Physical, only : gravitationalConstantGalacticus
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
         &  +sqrt(gravitationalConstantGalacticus)                     &
         &  *exponentialDiskRotationCurveGradient
    return
  end function exponentialDiskRotationCurveGradient

  double precision function exponentialDiskPotential(self,coordinates)
    !% Return the gravitational potential for an exponential disk.
    use :: Coordinates                 , only : assignment(=)                  , coordinateCylindrical
    use :: Numerical_Constants_Physical, only : gravitationalConstantGalacticus
    implicit none
    class           (massDistributionExponentialDisk), intent(inout) :: self
    class           (coordinate                     ), intent(in   ) :: coordinates
    type            (coordinateCylindrical          )                :: position
    double precision                                                 :: correctionSmallRadius, halfRadius, &
         &                                                              radius

    ! Get position in cylindrical coordinate system.
    position=coordinates
    ! Compute density.
    radius=position%r()
    ! If the radius is sufficiently large, treat the disk as a point mass.
    if (radius > exponentialDiskPotentialRadiusMaximum*self%scaleRadius) then
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
         &  +gravitationalConstantGalacticus               &
         &  *exponentialDiskPotential
    return
  end function exponentialDiskPotential

  double precision function exponentialDiskBesselFactorPotential(self,halfRadius)
    !% Compute Bessel function factors appearing in the expression for an razor-thin exponential
    !% disk gravitational potential.
    use :: Numerical_Constants_Math, only : eulersConstant, ln2
    implicit none
    class           (massDistributionExponentialDisk), intent(inout) :: self
    double precision                                 , intent(in   ) :: halfRadius

    ! For small half-radii, use a series expansion for a more accurate result.
    if      (halfRadius <=                                0.0d0) then
       exponentialDiskBesselFactorPotential=1.0d0
    else if (halfRadius < exponentialDiskPotentialRadiusMinimum) then
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
    !% Compute Bessel function factors appearing in the expression for an razor-thin exponential disk rotation curve.
    use :: Bessel_Functions        , only : Bessel_Function_I0, Bessel_Function_I1, Bessel_Function_K0, Bessel_Function_K1
    use :: Numerical_Constants_Math, only : eulersConstant    , ln2
    implicit none
    class           (massDistributionExponentialDisk), intent(inout) :: self
    double precision                                 , intent(in   ) :: halfRadius
    double precision                                 , parameter     :: halfRadiusSmall                  =1.0d-3
    integer                                          , parameter     :: rotationCurvePointsPerDecade=10
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
    !% Compute Bessel function factors appearing in the expression for a razor-thin exponential disk rotation curve gradient.
    use :: Bessel_Functions        , only : Bessel_Function_I0, Bessel_Function_I1, Bessel_Function_K0, Bessel_Function_K1
    use :: Numerical_Constants_Math, only : eulersConstant    , ln2
    implicit none
    class           (massDistributionExponentialDisk), intent(inout) :: self
    double precision                                 , intent(in   ) :: halfRadius
    double precision                                 , parameter     :: halfRadiusSmall                     =1.0d-3
    double precision                                 , parameter     :: halfRadiusLarge                     =1.0d+2
    integer                                          , parameter     :: rotationCurveGradientPointsPerDecade=10
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
               &    2.0d0                                                                    &
               &   *x                                                                        &
               &   *(                                                                        &
               &      Bessel_Function_I0(x)*Bessel_Function_K0(x)                            &
               &     -Bessel_Function_I1(x)*Bessel_Function_K1(x)                            &
               &    )                                                                        &
               &   +x**2                                                                     &
               &   *(                                                                        &
               &        Bessel_Function_I1(x)                         *Bessel_Function_K0(x) &
               &     -  Bessel_Function_K1(x)                         *Bessel_Function_I0(x) &
               &     -( Bessel_Function_I0(x)-Bessel_Function_I1(x)/x)*Bessel_Function_K1(x) &
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
    !% Compute radial moments of the exponential disk mass distribution surface density profile.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    use :: Gamma_Functions , only : Gamma_Function         , Gamma_Function_Incomplete
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
       call Galacticus_Error_Report('moment is infinite'//{introspection:location})
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
