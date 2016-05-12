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

! Miyamoto-Nagai disk functions.
subroutine Mass_Distribution_Miyamoto_Nagai_Initialize(self,a,b,mass,isDimensionless)
  !% Initialize a $\beta$-profile mass distribution.
  use Numerical_Constants_Math
  use Numerical_Comparison
  use Galacticus_Error
  implicit none
  class           (massDistributionMiyamotoNagai), intent(  out)           :: self
  double precision                               , intent(in   ), optional :: a              ,b, &
       &                                                                      mass
  logical                                        , intent(in   ), optional :: isDimensionless

  ! Determine if profile is dimensionless.
  if (present(isDimensionless)) then
     self%dimensionless=isDimensionless
  else
     self%dimensionless=.false.
  end if
  ! If dimensionless, then set scale length and mass to unity.
  if (self%dimensionless) then
     if (.not.present(b)) call Galacticus_Error_Report('Mass_Distribution_Miyamoto_Nagai_Initialize','"b" must be specified for dimensionless profiles')
     if (present(a   )) then
        if (Values_Differ(a   ,1.0d0,absTol=1.0d-6)) call Galacticus_Error_Report('Mass_Distribution_Miyamoto_Nagai_Initialize','"a" should be unity for a dimensionless profile (or simply do not specify it)')
     end if
     if (present(mass)) then
        if (Values_Differ(mass,1.0d0,absTol=1.0d-6)) call Galacticus_Error_Report('Mass_Distribution_Miyamoto_Nagai_Initialize','mass should be unity for a dimensionless profile (or simply do not specify a mass)')
     end if
     self%a                   =1.0d0
     self%b                   =b
     self%mass                =1.0d0
  else
     if (.not.present(a   )) call Galacticus_Error_Report('Mass_Distribution_Miyamoto_Nagai_Initialize','"a" must be specified for dimensionful profiles')
     if (.not.present(b   )) call Galacticus_Error_Report('Mass_Distribution_Miyamoto_Nagai_Initialize','"b" must be specified for dimensionful profiles')
     if (.not.present(mass)) call Galacticus_Error_Report('Mass_Distribution_Miyamoto_Nagai_Initialize','mass must be specified for dimensionful profiles')
     self%a                   =a
     self%b                   =b
     self%mass                =mass
  end if
  ! Compute shape parameter and density normalization.
  self%shape               =+self%b       &
       &                    /self%a
  self%densityNormalization=              &
       &                    +self%mass    &
       &                    /self%a   **2 &
       &                    /self%b       &
       &                    /4.0d0        &
       &                    /Pi
  ! Set table initialization states.
  self%surfaceDensityInitialized=.false.
  self%  massEnclosedInitialized=.false.
  return
end subroutine Mass_Distribution_Miyamoto_Nagai_Initialize

double precision function Mass_Distribution_Miyamoto_Nagai_Density(self,coordinates)
  !% Return the density at the specified {\normalfont \ttfamily coordinates} in an exponential disk mass distribution.
  use Coordinates
  implicit none
  class           (massDistributionMiyamotoNagai), intent(in   ) :: self
  class           (coordinate                   ), intent(in   ) :: coordinates
  type            (coordinateCylindrical        )                :: position
  double precision                                               :: r          , z

  ! Get position in cylindrical coordinate system.
  position=coordinates
  ! Compute density.
  r=    position%r() /self%a
  z=abs(position%z())/self%b
  ! Evaluate the density.
  Mass_Distribution_Miyamoto_Nagai_Density=+self%densityNormalization &
       &                                   *(                         &
       &                                     +r**2                    &
       &                                     +(                       &
       &                                       +1.0d0                 &
       &                                       +self%shape            &
       &                                       *sqrt(                 &
       &                                             +1.0d0           &
       &                                             +z**2            &
       &                                            )                 &
       &                                      )**2                    &
       &                                     *(                       &
       &                                       +1.0d0                 &
       &                                       +3.0d0                 &
       &                                       *self%shape            &
       &                                       *sqrt(                 &
       &                                             +1.0d0           &
       &                                             +z**2            &
       &                                            )                 &
       &                                      )                       &
       &                                    )                         &
       &                                   /(                         &
       &                                     +r**2                    &
       &                                     +(                       &
       &                                       +1.0d0                 &
       &                                       +self%shape            &
       &                                       *sqrt(                 &
       &                                             +1.0d0           &
       &                                             +z**2            &
       &                                            )                 &
       &                                      )**2                    &
       &                                    )**2.5d0                  &
       &                                   /(                         &
       &                                     +1.0d0                   &
       &                                     +z**2                    &
       &                                    )**1.5d0
  return
end function Mass_Distribution_Miyamoto_Nagai_Density

double precision function Mass_Distribution_Miyamoto_Nagai_Mass_Enc_By_Sphere(self,radius)
  !% Computes the mass enclosed within a sphere of given {\normalfont \ttfamily radius} for exponential disk mass distributions.
  use Numerical_Constants_Math
  implicit none
  class           (massDistributionMiyamotoNagai), intent(inout), target :: self
  double precision                               , intent(in   )         :: radius

  ! Ensure mass enclosed profile is tabulated.
  call self%massEnclosedInitialize()
  ! Evaluate the mass enclosed.
  Mass_Distribution_Miyamoto_Nagai_Mass_Enc_By_Sphere=self%massEnclosedTable%interpolate(radius)
  return
end function Mass_Distribution_Miyamoto_Nagai_Mass_Enc_By_Sphere

double precision function Mass_Distribution_Miyamoto_Nagai_Half_Mass_Radius(self)
  !% Return the half-mass radius in a Miyamoto-Nagai mass distribution.
  implicit none
  class           (massDistributionMiyamotoNagai), intent(inout) :: self

  ! Ensure mass enclosed profile is tabulated.
  call self%massEnclosedInitialize()
  ! Return the half-mass radius.
  Mass_Distribution_Miyamoto_Nagai_Half_Mass_Radius=self%halfMassRadiusValue
  return
end function Mass_Distribution_Miyamoto_Nagai_Half_Mass_Radius

subroutine Mass_Distribution_Miyamoto_Nagai_Mass_Enclosed_Initialize(self)
  !% Construct a tabulation of the mass enclosed by a sphere in a Miyamoto-Nagai mass distribution.
  use Numerical_Integration
  use Numerical_Constants_Math
  use FGSL
  use Table_Labels
  implicit none
  class           (massDistributionMiyamotoNagai), intent(inout) :: self
  double precision                               , parameter     :: radiusMinimum       =1.0d-6, radiusMaximum =1.0d3
  double precision                               , parameter     :: heightStart         =0.0d+0, heightInfinity=1.0d3
  integer                                        , parameter     :: radiusCount         =300
  type            (fgsl_function                )                :: integrandFunction
  type            (fgsl_integration_workspace   )                :: integrationWorkspace
  integer                                                        :: i
  double precision                                               :: radius

  ! Mass enclosed is not analytically tractable. Tabulate the results of numerical integration.
  if (.not.self%massEnclosedInitialized) then
     ! Initialize a table of radii.
     call self%massEnclosedTable%create(radiusMinimum*self%a,radiusMaximum*self%a,radiusCount,1,extrapolationType=spread(extrapolationTypeExtrapolate,1,2))
     ! Compute surface density at each radius.
     do i=1,radiusCount
        call self                                                      &
             & %massEnclosedTable                                      &
             &  %populate(                                             &
             &            +4.0d0                                       &
             &            *Pi                                          &
             &            *Integrate(                                  &
             &                       0.0d0                           , &
             &                       self%massEnclosedTable%x(i)     , &
             &                       Mass_Enclosed_Integrand_R       , &
             &                       integrandFunction               , &
             &                       integrationWorkspace            , &
             &                       toleranceAbsolute        =0.0d+0, &
             &                       toleranceRelative        =1.0d-3  &
             &                )                                      , &
             &            i                                            &
             &           )
        call Integrate_Done(integrandFunction,integrationWorkspace)
     end do
     ! Find the half-mass radius.
     do i=2,radiusCount
        if     (                                                                     &
             &   self%massEnclosedTable%y(i  ) >= 0.5d0*self%massEnclosedTable%y(-1) &
             &  .and.                                                                &
             &   self%massEnclosedTable%y(i-1) <  0.5d0*self%massEnclosedTable%y(-1) &
             & ) self%halfMassRadiusValue=+  self%massEnclosedTable%x(i-1)           &
             &                            +(                                         &
             &                              +self%massEnclosedTable%x(i  )           &
             &                              -self%massEnclosedTable%x(i-1)           &
             &                             )                                         &
             &                            *(                                         &
             &                              +0.5d0                                   &
             &                              *self%massEnclosedTable%y( -1)           &
             &                              -self%massEnclosedTable%y(i-1)           &
             &                             )                                         &
             &                            /(                                         &
             &                              +self%massEnclosedTable%y(i  )           &
             &                              -self%massEnclosedTable%y(i-1)           &
             &                            )
     end do
     ! Record that surface density is now initialized.
     self%massEnclosedInitialized=.true.
  end if
  return

contains
  
  double precision function Mass_Enclosed_Integrand_R(r)
    !% Integrand function used for finding the mass enclosed by a sphere in Miyamoto-Nagai disks.
    use Coordinates
    implicit none
    double precision                            , intent(in   ) :: r
    type            (fgsl_function             )                :: integrandFunction1
    type            (fgsl_integration_workspace)                :: integrationWorkspace1

    radius=r
    Mass_Enclosed_Integrand_R=+radius                                                 &
         &                    *Integrate(                                             &
         &                               0.0d0                                      , &
         &                               sqrt(                                        &
         &                                    +self%massEnclosedTable%x(i)**2         &
         &                                    -radius                     **2         &
         &                                   )                                      , &
         &                               Mass_Enclosed_Integrand_Z                  , &
         &                               integrandFunction1                         , &
         &                               integrationWorkspace1                      , &
         &                               toleranceAbsolute                   =0.0d+0, &
         &                               toleranceRelative                   =1.0d-3  &
         &                              )
    call Integrate_Done(integrandFunction1,integrationWorkspace1)
    return
  end function Mass_Enclosed_Integrand_R

  double precision function Mass_Enclosed_Integrand_Z(z)
    !% Integrand function used for finding the mass enclosed by a sphere in Miyamoto-Nagai disks.
    use Coordinates
    implicit none
    double precision                       , intent(in   ) :: z
    type            (coordinateCylindrical)                :: position

    position                 =[radius,0.0d0,z]
    Mass_Enclosed_Integrand_Z=self%density(position)
    return
  end function Mass_Enclosed_Integrand_Z

end subroutine Mass_Distribution_Miyamoto_Nagai_Mass_Enclosed_Initialize

subroutine Mass_Distribution_Miyamoto_Nagai_Surface_Density_Initialize(self)
  !% Construct a tabulation of the surface density profile in a Miyamoto-Nagai mass distribution.
  use Numerical_Integration
  use FGSL
  use Table_Labels
  implicit none
  class           (massDistributionMiyamotoNagai), intent(inout) :: self
  double precision                               , parameter     :: radiusMinimum       =1.0d-6, radiusMaximum =1.0d3
  double precision                               , parameter     :: heightStart         =0.0d+0, heightInfinity=1.0d3
  integer                                        , parameter     :: radiusCount         =300
  type            (fgsl_function                )                :: integrandFunction
  type            (fgsl_integration_workspace   )                :: integrationWorkspace
  integer                                                        :: i

  ! Surface density is not analytically tractable. Tabulate the results of numerical integration.
  if (.not.self%surfaceDensityInitialized) then
     ! Initialize a table of radii.
     call self%surfaceDensityTable%create(radiusMinimum*self%a,radiusMaximum*self%a,radiusCount,1,extrapolationType=spread(extrapolationTypeExtrapolate,1,2))
     ! Compute surface density at each radius.
     do i=1,radiusCount
        call self                                                      &
             & %surfaceDensityTable                                    &
             &  %populate(                                             &
             &            +2.0d0                                       &
             &            *Integrate(                                  &
             &                       self%b*heightStart              , &
             &                       self%b*heightInfinity           , &
             &                       Surface_Density_Integrand       , &
             &                       integrandFunction               , &
             &                       integrationWorkspace            , &
             &                       toleranceAbsolute        =0.0d+0, &
             &                       toleranceRelative        =1.0d-3  &
             &                )                                      , &
             &            i                                            &
             &           )
        call Integrate_Done(integrandFunction,integrationWorkspace)
     end do
     ! Record that surface density is now initialized.
     self%surfaceDensityInitialized=.true.
  end if
  return

contains
  
  double precision function Surface_Density_Integrand(height)
    !% Integrand function used for finding the surface density of Miyamoto-Nagai disks.
    use Coordinates
    implicit none
    double precision                       , intent(in   ) :: height
    type            (coordinateCylindrical)                :: position

    position                 =[self%surfaceDensityTable%x(i),0.0d0,height]
    Surface_Density_Integrand=self%density(position)
    return
  end function Surface_Density_Integrand

end subroutine Mass_Distribution_Miyamoto_Nagai_Surface_Density_Initialize

double precision function Mass_Distribution_Miyamoto_Nagai_Surface_Density(self,coordinates)
  !% Return the surface density at the specified {\normalfont \ttfamily coordinates} in a Miyamoto-Nagai mass distribution.
  use Coordinates
  implicit none
  class           (massDistributionMiyamotoNagai), intent(inout) :: self
  class           (coordinate                   ), intent(in   ) :: coordinates
  type            (coordinateCylindrical        )                :: position

  ! Ensure surface density profile is tabulated.
  call self%surfaceDensityInitialize()
  ! Evaluate the surface density.
  position=coordinates
  Mass_Distribution_Miyamoto_Nagai_Surface_Density=self%surfaceDensityTable%interpolate(position%r())
  return
end function Mass_Distribution_Miyamoto_Nagai_Surface_Density

double precision function Mass_Distribution_Miyamoto_Nagai_Srfc_Dnsty_Rdl_Mmnt(self,moment,radiusMinimum,radiusMaximum,isInfinite)
  !% Compute radial moments of the Miyamoto-Nagai mass distribution surface density profile.
  use Galacticus_Error
  implicit none
  class           (massDistributionMiyamotoNagai), intent(inout)           :: self
  double precision                               , intent(in   )           :: moment 
  double precision                               , intent(in   ), optional :: radiusMinimum          , radiusMaximum
  logical                                        , intent(  out), optional :: isInfinite
  procedure       (integrandWeight              ), pointer                 :: integrandWeightFunction
  double precision                                                         :: radiusMinimumActual    , radiusMaximumActual

  ! Set infinity status.
  if (present(isInfinite)) isInfinite=.false.
  ! Ensure surface density profile is tabulated.
  call self%surfaceDensityInitialize()
  ! Determine the radii to use.
  if (present(radiusMinimum)) then
     radiusMinimumActual=radiusMinimum
  else
     radiusMinimumActual=self%surfaceDensityTable%x(+1)
  end if
  if (present(radiusMaximum)) then
     radiusMaximumActual=radiusMaximum
  else
     radiusMaximumActual=self%surfaceDensityTable%x(-1)
  end if
  ! Check for finite moments.
  if     (                                                                     &
       &   (moment <= -1.0d0 .and.              radiusMinimumActual  <= 0.0d0) &
       &  .or.                                                                 &
       &   (moment >=  2.0d0 .and. .not.present(radiusMaximum      )         ) &
       & ) then
     Mass_Distribution_Miyamoto_Nagai_Srfc_Dnsty_Rdl_Mmnt=0.0d0
     if (present(isInfinite)) then
        isInfinite=.true.
        return
     else
        call Galacticus_Error_Report('Mass_Distribution_Miyamoto_Nagai_Srfc_Dnsty_Rdl_Mmnt','-1 < momemnt < 2 is required for finite moment')
     end if
  end if
  ! Perform the integration.
  integrandWeightFunction => momentWeight
  Mass_Distribution_Miyamoto_Nagai_Srfc_Dnsty_Rdl_Mmnt=                                      &
       & +sum(                                                                               &
       &      +         self%surfaceDensityTable%integrationWeights(                         &
       &                                                            radiusMinimumActual    , &
       &                                                            radiusMaximumActual    , &
       &                                                            integrandWeightFunction  &
       &                                                           )                         &
       &      *reshape(                                                                      &
       &                self%surfaceDensityTable%ys                (                         &
       &                                                           )                       , &
       &               [                                                                     &
       &                self%surfaceDensityTable%size              (                         &
       &                                                           )                         &
       &               ]                                                                     &
       &              )                                                                      &
       &     )
  return

contains  
  
  double precision function momentWeight(radius)
    !% The weight function used in computing radial moments of the Miyamoto-Nagai mass distribution surface density.
    implicit none
    double precision, intent(in   ) :: radius

    momentWeight=radius**moment
    return
  end function momentWeight
  
end function Mass_Distribution_Miyamoto_Nagai_Srfc_Dnsty_Rdl_Mmnt

double precision function Mass_Distribution_Miyamoto_Nagai_Rotation_Curve(self,radius)
  !% Return the mid-plane rotation curve for an exponential disk.
  use Coordinates
  use Numerical_Constants_Physical
  implicit none
  class           (massDistributionMiyamotoNagai), intent(inout) :: self
  double precision                               , intent(in   ) :: radius
  double precision                                               :: r
  
  ! Get dimensionless radius.
  r=radius/self%a
  ! Evaluate the rotation curve.
  Mass_Distribution_Miyamoto_Nagai_Rotation_Curve=+sqrt(           &
       &                                                +self%mass &
       &                                                /self%a    &
       &                                               )           &
       &                                          *r               &
       &                                          /(               &
       &                                            +r**2          &
       &                                            +(             &
       &                                              +1.0d0       &
       &                                              +self%shape  &
       &                                             )**2          &
       &                                           )**0.75d0
  ! Make dimensionfull if necessary.
  if (.not.self%dimensionless)                             &
       & Mass_Distribution_Miyamoto_Nagai_Rotation_Curve=  &
       &  +gravitationalConstantGalacticus                 &
       &  *Mass_Distribution_Miyamoto_Nagai_Rotation_Curve
  return
end function Mass_Distribution_Miyamoto_Nagai_Rotation_Curve

double precision function Mass_Distribution_Miyamoto_Nagai_Rotation_Curve_Gradient(self,radius)
  !% Return the mid-plane rotation curve gradient for an exponential disk.
  use Numerical_Constants_Physical
  use Coordinates
  implicit none
  class           (massDistributionMiyamotoNagai), intent(inout) :: self
  double precision                               , intent(in   ) :: radius
  double precision                                               :: r
  
  ! Get dimensionless radius.
  r=radius/self%a
  ! Evaluate the rotation curve.
  Mass_Distribution_Miyamoto_Nagai_Rotation_Curve_Gradient=+self%mass        &
       &                                                   /self%a   **2     &
       &                                                   *(                &
       &                                                     +2.0d0          &
       &                                                     *r              &
       &                                                     -3.0d0          &
       &                                                     *r      **3     &
       &                                                     /(              &
       &                                                       +r    **2     &
       &                                                       +(            &
       &                                                         +1.0d0      &
       &                                                         +self%shape &
       &                                                        )**2         &
       &                                                      )              &
       &                                                    )                &
       &                                                   /(                &
       &                                                     +r**2           &
       &                                                     +(              &
       &                                                       +1.0d0        &
       &                                                       +self%shape   &
       &                                                      )**2           &
       &                                                    )**1.5d0
  ! Make dimensionfull if necessary.
  if (.not.self%dimensionless)                                      &
       & Mass_Distribution_Miyamoto_Nagai_Rotation_Curve_Gradient=  &
       &  +gravitationalConstantGalacticus                          &
       &  *Mass_Distribution_Miyamoto_Nagai_Rotation_Curve_Gradient
  return
end function Mass_Distribution_Miyamoto_Nagai_Rotation_Curve_Gradient

double precision function Mass_Distribution_Miyamoto_Nagai_Potential(self,coordinates)
  !% Return the gravitational potential for an exponential disk.
  use Numerical_Constants_Physical
  use Coordinates
  implicit none
  class           (massDistributionMiyamotoNagai), intent(inout) :: self
  class           (coordinate                   ), intent(in   ) :: coordinates
  type            (coordinateCylindrical        )                :: position
  double precision                                               :: r          , z

  ! Get position in cylindrical coordinate system.
  position=coordinates
  ! Compute density.
  r=    position%r() /self%a
  z=abs(position%z())/self%b
  ! Evaluate the potential.
  Mass_Distribution_Miyamoto_Nagai_Potential=-self%mass                       &
       &                                     /self%a                          &
       &                                     /sqrt(                           &
       &                                           +r**2                      &
       &                                           +(                         &
       &                                             +1.0d0                   &
       &                                             +self%shape              &
       &                                             *sqrt(                   &
       &                                                   +1.0d0             &
       &                                                   +z**2              &
       &                                                  )                   &
       &                                            )**2                      &
       &                                          )
  ! Make dimensionfull if necessary.
  if (.not.self%dimensionless)                        &
       & Mass_Distribution_Miyamoto_Nagai_Potential=  &
       &  +gravitationalConstantGalacticus            &
       &  *Mass_Distribution_Miyamoto_Nagai_Potential
  return
end function Mass_Distribution_Miyamoto_Nagai_Potential
