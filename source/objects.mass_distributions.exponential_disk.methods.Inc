!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
!!    Andrew Benson <abenson@obs.carnegiescience.edu>
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

! Exponential disk functions.
subroutine Mass_Distribution_Exponential_Disk_Initialize(self,scaleRadius,scaleHeight,mass,isDimensionless)
  !% Initialize a $\beta$-profile mass distribution.
  use Numerical_Constants_Math
  use Numerical_Comparison
  use Galacticus_Error
  implicit none
  class           (massDistributionExponentialDisk), intent(  out)           :: self
  double precision                                 , intent(in   ), optional :: scaleRadius    , scaleHeight, &
       &                                                                        mass
  logical                                          , intent(in   ), optional :: isDimensionless
  double precision                                 , parameter               :: rotationCurveHalfRadiusMaximumDefault=1.0d+1, &
       &                                                                        rotationCurveHalfRadiusMinimumDefault=1.0d-6

  ! Determine if profile is dimensionless.
  if (present(isDimensionless)) then
     self%dimensionless=isDimensionless
  else
     self%dimensionless=.false.
  end if
  ! If dimensionless, then set scale length and mass to unity.
  if (self%dimensionless) then
     if (present(scaleRadius)) then
        if (Values_Differ(scaleRadius,1.0d0,absTol=1.0d-6)) call Galacticus_Error_Report('Mass_Distribution_Exponential_Disk_Initialize','scaleRadius should be unity for a dimensionless profile (or simply do not specify a scale length)')
     end if
     if (present(mass       )) then
        if (Values_Differ(mass       ,1.0d0,absTol=1.0d-6)) call Galacticus_Error_Report('Mass_Distribution_Exponential_Disk_Initialize','mass should be unity for a dimensionless profile (or simply do not specify a mass)')
     end if
     self%scaleRadius                =1.0d0
     self%mass                       =1.0d0
     self%surfaceDensityNormalization=1.0d0/2.0d0/Pi
  else
     ! Set core radius.
     if (.not.present(scaleRadius)) call Galacticus_Error_Report('Mass_Distribution_Exponential_Disk_Initialize','scale radius must be specified for dimensionful profiles')
     if (.not.present(mass       )) call Galacticus_Error_Report('Mass_Distribution_Exponential_Disk_Initialize','mass must be specified for dimensionful profiles')
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
     self%scaleHeight=0.0d0
  end if
  ! Initialize rotation curve tables.
  self%rotationCurveHalfRadiusMinimum        =rotationCurveHalfRadiusMinimumDefault
  self%rotationCurveHalfRadiusMaximum        =rotationCurveHalfRadiusMaximumDefault
  self%rotationCurveGradientHalfRadiusMinimum=rotationCurveHalfRadiusMinimumDefault
  self%rotationCurveGradientHalfRadiusMaximum=rotationCurveHalfRadiusMaximumDefault
  return
end subroutine Mass_Distribution_Exponential_Disk_Initialize

subroutine Mass_Distribution_Exponential_Disk_Build_Tables(self)
  !% Build tables used for exponential disk mass distributions.
  use Bessel_Functions
  use Table_Labels
  implicit none
  class           (massDistributionExponentialDisk), intent(inout) :: self
  integer                                                          :: i   , potentialPointsCount
  double precision                                                 :: x
  !GCC$ attributes unused :: self
  
  ! Build table if necessary.
  if (.not.exponentialDiskPotentialInitialized) then
     ! Determine how many points to tabulate.
     potentialPointsCount=int(log10(exponentialDiskPotentialRadiusMaximum/exponentialDiskPotentialRadiusMinimum)*dble(exponentialDiskPotentialPointsPerDecade))+1
     ! Create the table.
     call exponentialDiskPotentialTable%destroy()
     call exponentialDiskPotentialTable%create(exponentialDiskPotentialRadiusMinimum,exponentialDiskPotentialRadiusMaximum,potentialPointsCount,extrapolationType=spread(extrapolationTypeAbort,1,2))
     ! Compute Bessel factors.
     do i=1,potentialPointsCount
        x=exponentialDiskPotentialTable%x(i)
        call exponentialDiskPotentialTable%populate(                                               &
             &                                      +x                                             &
             &                                      *(                                             &
             &                                        +Bessel_Function_I0(x)*Bessel_Function_K1(x) &
             &                                        -Bessel_Function_I1(x)*Bessel_Function_K0(x) &
             &                                       ),                                            &
             &                                      i                                              &
             &                                     )
     end do
     ! Record that table is initialized.
     exponentialDiskPotentialInitialized=.true.
  end if
  return
end subroutine Mass_Distribution_Exponential_Disk_Build_Tables

double precision function Mass_Distribution_Exponential_Disk_Half_Mass_Radius(self)
  !% Return the half-mass radius in an exponential disk mass distribution.
  implicit none
  class           (massDistributionExponentialDisk), intent(inout) :: self
  double precision                                 , parameter     :: halfMassRadiusToScaleRadius=1.678346990d0

  Mass_Distribution_Exponential_Disk_Half_Mass_Radius=halfMassRadiusToScaleRadius*self%scaleRadius
  return
end function Mass_Distribution_Exponential_Disk_Half_Mass_Radius

double precision function Mass_Distribution_Exponential_Disk_Density(self,coordinates)
  !% Return the density at the specified {\normalfont \ttfamily coordinates} in an exponential disk mass distribution.
  use Galacticus_Error
  use Coordinates
  implicit none
  class           (massDistributionExponentialDisk), intent(in   ) :: self
  class           (coordinate                     ), intent(in   ) :: coordinates
  type            (coordinateCylindrical          )                :: position
  double precision                                 , parameter     :: coshArgumentMaximum=50.0d0
  double precision                                                 :: r                         , z, &
       &                                                              coshTerm

  ! If disk is razor thin, density is undefined.
  if (self%scaleHeight <= 0.0d0) call Galacticus_Error_Report('Mass_Distribution_Exponential_Disk_Density','density undefined for razor-thin disk')
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
  Mass_Distribution_Exponential_Disk_Density=self%densityNormalization*exp(-r)*coshTerm
  return
end function Mass_Distribution_Exponential_Disk_Density

double precision function Mass_Distribution_Exponential_Disk_Mass_Enc_By_Sphere(self,radius)
  !% Computes the mass enclosed within a sphere of given {\normalfont \ttfamily radius} for exponential disk mass distributions.
  use Numerical_Constants_Math
  implicit none
  class           (massDistributionExponentialDisk), intent(inout), target :: self
  double precision                                 , intent(in   )         :: radius
  double precision                                                         :: fractionalRadius

  fractionalRadius                                     =+radius                              &
       &                                                /self%scaleRadius
  Mass_Distribution_Exponential_Disk_Mass_Enc_By_Sphere=+2.0d0                               &
       &                                                *Pi                                  &
       &                                                *self%scaleRadius                **2 &
       &                                                *self%surfaceDensityNormalization    &
       &                                                *(                                   &
       &                                                  +1.0d0                             &
       &                                                  -(                                 &
       &                                                    +1.0d0                           &
       &                                                    +fractionalRadius                &
       &                                                   )                                 &
       &                                                  *exp(-fractionalRadius)            &
       &                                                )
  return
end function Mass_Distribution_Exponential_Disk_Mass_Enc_By_Sphere

double precision function Mass_Distribution_Exponential_Disk_Surface_Density(self,coordinates)
  !% Return the surface density at the specified {\normalfont \ttfamily coordinates} in an exponential disk mass distribution.
  use Coordinates
  implicit none
  class           (massDistributionExponentialDisk), intent(inout) :: self
  class           (coordinate                     ), intent(in   ) :: coordinates
  type            (coordinateCylindrical          )                :: position
  double precision                                                 :: r

  ! Get position in cylindrical coordinate system.
  position=coordinates
  ! Compute density.
  r=position%r()/self%scaleRadius
  Mass_Distribution_Exponential_Disk_Surface_Density=self%surfaceDensityNormalization*exp(-r)
  return
end function Mass_Distribution_Exponential_Disk_Surface_Density

double precision function Mass_Distribution_Exponential_Disk_Rotation_Curve(self,radius)
  !% Return the mid-plane rotation curve for an exponential disk.
  use Coordinates
  use Numerical_Constants_Physical
  implicit none
  class           (massDistributionExponentialDisk), intent(inout) :: self
  double precision                                 , intent(in   ) :: radius
  double precision                                                 :: r           , halfR, &
       &                                                              radiusFactor

  ! Get scale-free radius.
  r=radius/self%scaleRadius
  ! Compute rotation curve.
  if (r > radiusMaximum) then
     ! Beyond some maximum radius, approximate the disk as a spherical distribution to avoid evaluating Bessel functions for
     ! very large arguments.
     Mass_Distribution_Exponential_Disk_Rotation_Curve=sqrt(self%mass/radius)
  else
     ! We are often called at precisely one scale length. Use pre-computed factors in that case.
     if (r == 1.0d0) then
        !$omp critical (ExponentialDiskFactorCompute)
        if (.not.self%scaleLengthFactorSet) then
           halfR                    =0.5d0
           self%scaleLengthFactor   =Mass_Distribution_Exponential_Disk_Bessel_Factors(self,halfR)
           self%scaleLengthFactorSet=.true.
        end if
        !$omp end critical (ExponentialDiskFactorCompute)
        radiusFactor=self%scaleLengthFactor
     else
        halfR       =0.5d0*r
        radiusFactor=Mass_Distribution_Exponential_Disk_Bessel_Factors(self,halfR)
     end if
     Mass_Distribution_Exponential_Disk_Rotation_Curve=sqrt(2.0d0*(self%mass/self%scaleRadius)*radiusFactor)
  end if
  ! Make dimensionful if necessary.
  if (.not.self%dimensionless) Mass_Distribution_Exponential_Disk_Rotation_Curve= &
       & +sqrt(gravitationalConstantGalacticus)                                   &
       & *Mass_Distribution_Exponential_Disk_Rotation_Curve
  return
end function Mass_Distribution_Exponential_Disk_Rotation_Curve

double precision function Mass_Distribution_Exponential_Disk_Rotation_Curve_Gradient(self,radius)
  !% Return the mid-plane rotation curve gradient for an exponential disk.
  use Numerical_Constants_Physical
  use Coordinates
  implicit none
  class           (massDistributionExponentialDisk), intent(inout) :: self
  double precision                                 , intent(in   ) :: radius
  double precision                                 , parameter     :: fractionalRadiusMaximum=30.0d0
  double precision                                                 :: besselArgument                , besselFactor
  
  ! Compute Bessel functions argument.
  besselArgument= radius           &
       &         /2.0d0            &
       &         /self%scaleRadius
  if (2.0d0*besselArgument > fractionalRadiusMaximum) then
     ! Beyond some maximum radius, approximate the disk as a point mass to avoid evaluating Bessel functions for
     ! very large arguments.
     Mass_Distribution_Exponential_Disk_Rotation_Curve_Gradient=-self%mass/radius**2
  else
     ! Compute the gradient.
     besselFactor=Mass_Distribution_Exponential_Disk_Grdnt_Bssl_Fctrs(self,besselArgument)
     Mass_Distribution_Exponential_Disk_Rotation_Curve_Gradient= &
          &                         +self%mass                   &
          &                         *besselFactor                &
          &                         /self%scaleRadius**2
  end if
  ! Make dimensionful if necessary.
  if (.not.self%dimensionless) Mass_Distribution_Exponential_Disk_Rotation_Curve_Gradient= &
       &  +sqrt(gravitationalConstantGalacticus)                                           &
       &  *Mass_Distribution_Exponential_Disk_Rotation_Curve_Gradient
  return
end function Mass_Distribution_Exponential_Disk_Rotation_Curve_Gradient

double precision function Mass_Distribution_Exponential_Disk_Potential(self,coordinates)
  !% Return the gravitational potential for an exponential disk.
  use Bessel_Functions
  use Numerical_Constants_Physical
  use Coordinates
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
     Mass_Distribution_Exponential_Disk_Potential=-self%mass/radius
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
     Mass_Distribution_Exponential_Disk_Potential=                                                    &
          &             -self%mass                                                                    &
          &             /self%scaleRadius                                                             &
          &             *Mass_Distribution_Exponential_Disk_Potential_Bessel_Factors(self,halfRadius) &
          &             +correctionSmallRadius     
  end if
  ! Make dimensionful if necessary.
  if (.not.self%dimensionless) Mass_Distribution_Exponential_Disk_Potential= &
       &  +gravitationalConstantGalacticus                                   &
       &  *Mass_Distribution_Exponential_Disk_Potential
  return
end function Mass_Distribution_Exponential_Disk_Potential

double precision function Mass_Distribution_Exponential_Disk_Potential_Bessel_Factors(self,halfR)
  !% Compute Bessel function factors appearing in the expression for an razor-thin exponential
  !% disk gravitational potential.
  use Numerical_Constants_Math
  use Bessel_Functions
  implicit none
  class           (massDistributionExponentialDisk), intent(inout) :: self
  double precision                                 , intent(in   ) :: halfR

  ! For small half-radii, use a series expansion for a more accurate result.
  if      (halfR <=                                0.0d0) then
     Mass_Distribution_Exponential_Disk_Potential_Bessel_Factors=1.0d0
  else if (halfR < exponentialDiskPotentialRadiusMinimum) then
     Mass_Distribution_Exponential_Disk_Potential_Bessel_Factors=1.0d0+(eulersConstant-ln2+log(halfR))*halfR**2
  else
     call Mass_Distribution_Exponential_Disk_Build_Tables(self)
     Mass_Distribution_Exponential_Disk_Potential_Bessel_Factors=exponentialDiskPotentialTable%interpolate(halfR)
  end if
  return
end function Mass_Distribution_Exponential_Disk_Potential_Bessel_Factors

double precision function Mass_Distribution_Exponential_Disk_Bessel_Factors(self,halfR)
  !% Compute Bessel function factors appearing in the expression for an razor-thin exponential disk rotation curve.
  use Numerical_Constants_Math
  use Bessel_Functions
  implicit none
  class           (massDistributionExponentialDisk), intent(inout) :: self
  double precision                                 , intent(in   ) :: halfR
  double precision                                 , parameter     :: halfRSmall                  =1.0d-3
  integer                                          , parameter     :: rotationCurvePointsPerDecade=10
  integer                                                          :: iPoint                             , rotationCurvePointsCount
  double precision                                                 :: x

  ! For small half-radii, use a series expansion for a more accurate result.
  if (halfR <= 0.0d0) then
     Mass_Distribution_Exponential_Disk_Bessel_Factors=0.0d0
     return
  else if (halfR < halfRSmall) then
     Mass_Distribution_Exponential_Disk_Bessel_Factors=(ln2-eulersConstant-0.5d0-log(halfR))*halfR**2
     return
  end if
  !$omp critical(Exponential_Disk_Rotation_Curve_Tabulate)
  if     (                                       &
       &   .not.self%rotationCurveInitialized    &
       &  .or.                                   &
       &   halfR < self%rotationCurveTable%x(+1) &
       &  .or.                                   &
       &   halfR > self%rotationCurveTable%x(-1) &
       & ) then
     ! Find the minimum and maximum half-radii to tabulate.
     self%rotationCurveHalfRadiusMinimum=min(self%rotationCurveHalfRadiusMinimum,0.5d0*halfR)
     self%rotationCurveHalfRadiusMaximum=max(self%rotationCurveHalfRadiusMaximum,2.0d0*halfR)
     ! Determine how many points to tabulate.
     rotationCurvePointsCount=int(log10(self%rotationCurveHalfRadiusMaximum/self%rotationCurveHalfRadiusMinimum)*dble(rotationCurvePointsPerDecade))+1
     ! Allocate table arrays.
     call self%rotationCurveTable%destroy()
     call self%rotationCurveTable%create(self%rotationCurveHalfRadiusMinimum,self%rotationCurveHalfRadiusMaximum,rotationCurvePointsCount)
     ! Compute Bessel factors.
     do iPoint=1,rotationCurvePointsCount
        x=self%rotationCurveTable%x(iPoint)
        call self%rotationCurveTable%populate(                                          &
             &                            x**2                                          &
             &                           *(                                             &
             &                              Bessel_Function_I0(x)*Bessel_Function_K0(x) &
             &                             -Bessel_Function_I1(x)*Bessel_Function_K1(x) &
             &                            ),                                            &
             &                           iPoint                                         &
             &                          )
     end do
     ! Flag that the rotation curve is now initialized.
     self%rotationCurveInitialized=.true.
  end if
  ! Interpolate in the tabulated function.
  Mass_Distribution_Exponential_Disk_Bessel_Factors=self%rotationCurveTable%interpolate(halfR)
  !$omp end critical(Exponential_Disk_Rotation_Curve_Tabulate)
  return
end function Mass_Distribution_Exponential_Disk_Bessel_Factors

double precision function Mass_Distribution_Exponential_Disk_Grdnt_Bssl_Fctrs(self,halfR)
  !% Compute Bessel function factors appearing in the expression for a razor-thin exponential disk rotation curve gradient.
  use Numerical_Constants_Math
  use Bessel_Functions
  implicit none
  class           (massDistributionExponentialDisk), intent(inout) :: self
  double precision                                 , intent(in   ) :: halfR
  double precision                                 , parameter     :: halfRSmall                          =1.0d-3
  double precision                                 , parameter     :: halfRLarge                          =1.0d+2
  integer                                          , parameter     :: rotationCurveGradientPointsPerDecade=10
  integer                                                          :: iPoint                                      , rotationCurveGradientPointsCount
  double precision                                                 :: x

  ! For small and large half-radii, use a series expansion for a more accurate result.
  if (halfR == 0.0d0) then
     Mass_Distribution_Exponential_Disk_Grdnt_Bssl_Fctrs=0.0d0
     return
  else if (halfR < halfRSmall) then
     Mass_Distribution_Exponential_Disk_Grdnt_Bssl_Fctrs=(ln2-log(halfR)-eulersConstant-1.0d0)*halfR &
          & +(1.5d0*ln2-1.5d0*log(halfR)+0.25d0-1.5d0*eulersConstant)*halfR**3
     return
  else if (halfR > halfRLarge) then
     Mass_Distribution_Exponential_Disk_Grdnt_Bssl_Fctrs=-0.125d0-27.0d0/64.0d0/halfR**2
     return
  end if
  !$omp critical(Exponential_Disk_Rotation_Curve_Gradient_Tabulate)
  if     (                                               &
       &   .not.self%rotationCurveGradientInitialized    &
       &  .or.                                           &
       &   halfR < self%rotationCurveGradientTable%x(+1) &
       &  .or.                                           &
       &   halfR > self%rotationCurveGradientTable%x(-1) &
       & ) then
     ! Find the minimum and maximum half-radii to tabulate.
     self%rotationCurveGradientHalfRadiusMinimum=min(self%rotationCurveGradientHalfRadiusMinimum,0.5d0*halfR)
     self%rotationCurveGradientHalfRadiusMaximum=max(self%rotationCurveGradientHalfRadiusMaximum,2.0d0*halfR)
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
  Mass_Distribution_Exponential_Disk_Grdnt_Bssl_Fctrs=self%rotationCurveGradientTable%interpolate(halfR)
  !$omp end critical(Exponential_Disk_Rotation_Curve_Gradient_Tabulate)
  return
end function Mass_Distribution_Exponential_Disk_Grdnt_Bssl_Fctrs

double precision function Mass_Distribution_Exponential_Disk_Srfc_Dnsty_Rdl_Mmnt(self,moment,radiusMinimum,radiusMaximum,isInfinite)
  !% Compute radial moments of the exponential disk mass distribution surface density profile.
  use Gamma_Functions
  use Galacticus_Error
  implicit none
  class           (massDistributionExponentialDisk), intent(inout)           :: self
  double precision                                 , intent(in   )           :: moment 
  double precision                                 , intent(in   ), optional :: radiusMinimum, radiusMaximum
  logical                                          , intent(  out), optional :: isInfinite
  double precision                                                           :: integralLow  , integralHigh

  ! All moments n>-1 are finite.
  if (present(isInfinite)) isInfinite=(moment > -1.0d0)
  if (moment <= -1.0d0) then
     if (present(isInfinite)) return
     call Galacticus_Error_Report('Mass_Distribution_Exponential_Disk_Surface_Density_Radial_Moment','moment is infinite')
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
  Mass_Distribution_Exponential_Disk_Srfc_Dnsty_Rdl_Mmnt=(integralHigh-integralLow)*Gamma_Function(moment+1.0d0)*self%scaleRadius**(moment+1.0d0)
  return
end function Mass_Distribution_Exponential_Disk_Srfc_Dnsty_Rdl_Mmnt

subroutine Mass_Distribution_Exponential_Disk_State_Store(self,stateFile,fgslStateFile)
  !% Store the mass distribution to file.
  use FGSL
  implicit none
  class  (massDistributionExponentialDisk), intent(inout) :: self
  integer                                 , intent(in   ) :: stateFile
  type   (fgsl_file                      ), intent(in   ) :: fgslStateFile
  !GCC$ attributes unused :: fgslStateFile
  
  write (stateFile) self%scaleLengthFactor                     , self%scaleLengthFactorSet                  , &
       &            self%rotationCurveHalfRadiusMinimum        , self%rotationCurveHalfRadiusMaximum        , &
       &            self%rotationCurveGradientHalfRadiusMinimum, self%rotationCurveGradientHalfRadiusMaximum
  return
end subroutine Mass_Distribution_Exponential_Disk_State_Store

subroutine Mass_Distribution_Exponential_Disk_State_Restore(self,stateFile,fgslStateFile)
  !% Restore the mass distribution from file.
  use FGSL
  implicit none
  class  (massDistributionExponentialDisk), intent(inout) :: self
  integer                                 , intent(in   ) :: stateFile
  type   (fgsl_file                      ), intent(in   ) :: fgslStateFile
  !GCC$ attributes unused :: fgslStateFile
  
  read (stateFile) self%scaleLengthFactor                     , self%scaleLengthFactorSet                  , &
       &           self%rotationCurveHalfRadiusMinimum        , self%rotationCurveHalfRadiusMaximum        , &
       &           self%rotationCurveGradientHalfRadiusMinimum, self%rotationCurveGradientHalfRadiusMaximum
  ! Mark tables as uninitialized to trigger rebuilds.
  self%rotationCurveInitialized        =.false.
  self%rotationCurveGradientInitialized=.false.
  return
end subroutine Mass_Distribution_Exponential_Disk_State_Restore
