!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

! Beta-profile functions.
subroutine Mass_Distribution_Beta_Profile_Initialize(self,beta,coreRadius,densityNormalization,mass,outerRadius,isDimensionless)
  !% Initialize a $\beta$-profile mass distribution.
  use Numerical_Constants_Math
  use Hypergeometric_Functions
  use Numerical_Comparison
  use Galacticus_Error
  use Galacticus_Display
  implicit none
  class           (massDistributionBetaProfile), intent(  out)           :: self
  double precision                             , intent(in   )           :: beta
  double precision                             , intent(in   ), optional :: densityNormalization, mass      , &
       &                                                                    outerRadius         , coreRadius
  logical                                      , intent(in   ), optional :: isDimensionless
  double precision                                                       :: r
  character       (len=64                     )                          :: message

  ! Store beta parameter.
  self%beta=beta
  ! Determine if profile is dimensionless.
  if (present(isDimensionless)) then
     self%dimensionless=isDimensionless
  else
     self%dimensionless=.false.
  end if
  ! If dimensionless, then set scale length and mass to unity.
  if (self%dimensionless) then
     if (present(coreRadius         )) then
        if (Values_Differ(coreRadius         ,1.0d0,absTol=1.0d-6)) call Galacticus_Error_Report('Mass_Distribution_Beta_Profile_Initialize','coreRadius should be unity for a dimensionless profile (or simply do not specify a scale length)')
     end if
     if (present(densityNormalization)) then
        if (Values_Differ(densityNormalization,1.0d0,absTol=1.0d-6)) call Galacticus_Error_Report('Mass_Distribution_Beta_Profile_Initialize','densityNormalization should be unity for a dimensionless profile (or simply do not specify a densityNormalization)')
     end if
     if (present(mass                )) call Galacticus_Error_Report('Mass_Distribution_Beta_Profile_Initialize','mass cannot be specified for a dimensionless profile')
     if (present(outerRadius         )) call Galacticus_Error_Report('Mass_Distribution_Beta_Profile_Initialize','mass cannot be specified for a dimensionless profile')
     self%coreRadius          =1.0d0
     self%densityNormalization=1.0d0
  else
     ! Set core radius.
     if (.not.present(coreRadius)) call Galacticus_Error_Report('Mass_Distribution_Beta_Profile_Initialize','core radius must be specified for dimensionful profiles')
     self%coreRadius=coreRadius
     ! Determine density normalization.
     if      (                                   &
          &   present(densityNormalization)      &
          &  ) then
        self%densityNormalization=densityNormalization
     else if (&
          &   present(mass                ).and. &
          &   present(outerRadius         )      &
          &  ) then
        r=outerRadius/coreRadius
        if (outerRadius > 0.0d0) then
           self%densityNormalization=3.0d0*mass/4.0d0/Pi/outerRadius**3/Hypergeometric_2F1([1.5d0,1.5d0*beta],[2.5d0],-r**2)
        else
           call Galacticus_Error_Report('Mass_Distribution_Beta_Profile_Initialize','unphysical outer radius')
        end if
        ! Assert that the mass within the outer radius equals that specified.
        if (.not.Values_Agree(self%massEnclosedBySphere(outerRadius),mass,relTol=1.0d-6,absTol=tiny(0.0d0))) then
           call Galacticus_Display_Indent('beta-profile parameters:')
           write (message,'(a,e12.6)') '    coreRadius: ',coreRadius
           call Galacticus_Display_Message(message)
           write (message,'(a,e12.6)') '   outerRadius: ',outerRadius
           call Galacticus_Display_Message(message)
           write (message,'(a,e12.6)') '          mass: ',mass
           call Galacticus_Display_Message(message)
           write (message,'(a,e12.6)') '          beta: ',beta
           call Galacticus_Display_Message(message)
           write (message,'(a,e12.6)') 'mass(<r_outer): ',self%massEnclosedBySphere(outerRadius)
           call Galacticus_Display_Message(message)
           call Galacticus_Display_Unindent('done')
           call Galacticus_Error_Report('Mass_Distribution_Beta_Profile_Initialize','profile normalization failed')
        end if
     end if
  end if
  return
end subroutine Mass_Distribution_Beta_Profile_Initialize

double precision function Mass_Distribution_Beta_Profile_Density(self,coordinates)
  !% Return the density at the specified {\tt coordinates} in a $\beta$-profile mass distribution.
  use Coordinates
  implicit none
  class           (massDistributionBetaProfile), intent(in   ) :: self
  class           (coordinate                 ), intent(in   ) :: coordinates
  type            (coordinateSpherical        )                :: position
  double precision                                             :: r

  ! Get position in spherical coordinate system.
  position=coordinates
  ! Compute density.
  r=position%r()/self%coreRadius
  Mass_Distribution_Beta_Profile_Density=self%densityNormalization/(1.0d0+r**2)**(1.5d0*self%beta)
  return
end function Mass_Distribution_Beta_Profile_Density

double precision function Mass_Distribution_Beta_Profile_Density_Gradient_Radial(self,coordinates,logarithmic)
  !% Return the density at the specified {\tt coordinates} in a $\beta$-profile mass distribution.
  use Coordinates
  implicit none
  class           (massDistributionBetaProfile), intent(in   )           :: self
  class           (coordinate                 ), intent(in   )           :: coordinates
  logical                                      , intent(in   ), optional :: logarithmic
  type            (coordinateSpherical        )                          :: position
  double precision                                                       :: r
  logical                                                                :: logarithmicActual

  ! Set default options.
  logarithmicActual=.false.
  if (present(logarithmic)) logarithmicActual=logarithmic
  ! Get position in spherical coordinate system.
  position=coordinates
  r       =position%r()/self%coreRadius
  ! Compute density gradient.
  if (logarithmicActual) then
     Mass_Distribution_Beta_Profile_Density_Gradient_Radial= &
          & -3.0d0                                           &
          & *self%beta                                       &
          & * r**2                                           &
          & /(r**2+1.0d0)
  else
     Mass_Distribution_Beta_Profile_Density_Gradient_Radial= &
          & -3.0d0                                           &
          & *self%beta                                       &
          & *self%densityNormalization                       &
          & /self%coreRadius                                 &
          & * r**2                                           &
          & /(r**2+1.0d0)**(1.5d0*self%beta+1.0d0)
  end if
  return
end function Mass_Distribution_Beta_Profile_Density_Gradient_Radial

double precision function Mass_Distribution_Beta_Profile_Mass_Enc_By_Sphere(self,radius)
  !% Computes the mass enclosed within a sphere of given {\tt radius} for $\beta$-profile mass distributions. Result computed
  !% using \href{http://www.wolframalpha.com/input/?i=integrate+4*pi*r^2*rho%2F%281%2Br^2%29^%283*beta%2F2%29}{Wolfram Alpha}.
  use Numerical_Constants_Math
  use Hypergeometric_Functions
  implicit none
  class           (massDistributionBetaProfile), intent(in   ), target :: self
  double precision                             , intent(in   )         :: radius
  double precision                                                     :: fractionalRadius

  fractionalRadius=radius/self%coreRadius
  Mass_Distribution_Beta_Profile_Mass_Enc_By_Sphere=  &
       &  4.0d0                                       &
       & /3.0d0                                       &
       & *Pi                                          &
       & *self%densityNormalization                   &
       & *radius**3                                   &
       & *Hypergeometric_2F1(                         &
       &                     [1.5d0,1.5d0*self%beta], &
       &                     [2.5d0                ], &
       &                     -fractionalRadius**2     &
       &                    )
  return
end function Mass_Distribution_Beta_Profile_Mass_Enc_By_Sphere

double precision function Mass_Distribution_Beta_Profile_Potential(self,coordinates)
  !% Return the potential at the specified {\tt coordinates} in a $\beta$-profile mass distribution. Calculated using
  !% \href{http://www.wolframalpha.com/input/?i=integrate+4%2F3+%CF%80+r+%CF%81+2F1%283%2F2%2C+%283+%CE%B2%29%2F2%2C+5%2F2%2C+-r^2%29}{Wolfram
  !% Alpha}.
  use Numerical_Constants_Physical
  use Numerical_Constants_Math
  use Hypergeometric_Functions
  use Coordinates
  use Numerical_Comparison
  implicit none
  class           (massDistributionBetaProfile), intent(in   ) :: self
  class           (coordinate                 ), intent(in   ) :: coordinates
  type            (coordinateSpherical        )                :: position
  double precision                             , parameter     :: fractionalRadiusMinimum=1.0d-3
  double precision                                             :: fractionalRadius

  ! Get position in spherical coordinate system.
  position=coordinates
  ! Compute the potential at this position.
  
  fractionalRadius=position%r()/self%coreRadius
  if (Values_Agree(self%beta,2.0d0/3.0d0,absTol=1.0d-6)) then
     if (fractionalRadius < fractionalRadiusMinimum) then
        Mass_Distribution_Beta_Profile_Potential= &
             &  Pi                                &
             & *self%densityNormalization         &
             & *(                                 &
             &    4.0d0                           &
             &   +2.0d0                           &
             &   /3.0d0                           &
             &   *fractionalRadius**2             &
             &   -fractionalRadius**4             &
             &   /5.0d0                           &
             &  )
    else
        Mass_Distribution_Beta_Profile_Potential= &
             &  2.0d0                             &
             & *Pi                                &
             & *self%densityNormalization         &
             & *(                                 &
             &    log (                           &
             &          1.0d0                     &
             &         +fractionalRadius**2       &
             &        )                           &
             &   +2.0d0                           &
             &   *atan( fractionalRadius   )      &
             &   /      fractionalRadius          &
             &  )
     end if
  else
     if (fractionalRadius < fractionalRadiusMinimum) then
        Mass_Distribution_Beta_Profile_Potential= &
             &  Pi                                &
             & *self%densityNormalization         &
             & *fractionalRadius**3               &
             & *(                                 &
             &    4.0d0/3.0d0                     &
             &   -6.0d0                           &
             &   /5.0d0                           &
             &   *self%beta                       &
             &   *fractionalRadius**2             &
             &  )
     else
        Mass_Distribution_Beta_Profile_Potential=                   &
             &  2.0d0                                               &
             & /3.0d0                                               &
             & *Pi                                                  &
             & *self%densityNormalization                           &
             & *     fractionalRadius**2                            &
             & *(                                                   &
             &    6.0d0                                             &
             &   *  (fractionalRadius**2+1.0d0)**(-1.5d0*self%beta) &
             &   *(                                                 &
             &      (fractionalRadius**2+1.0d0)**(+1.5d0*self%beta) &
             &     -(fractionalRadius**2+1.0d0)                     &
             &    )                                                 &
             &   /(3.0d0*self%beta-2.0d0)                           &
             &   /   fractionalRadius**2                            &
             &   -2.0d0                                             &
             &   *Hypergeometric_2F1(                               &
             &                       [1.5d0,1.5d0*self%beta],       &
             &                       [2.5d0                ],       &
             &                       -fractionalRadius**2           &
             &                      )                               &
             &  )
     end if
  end if
  if (.not.self%isDimensionless())                  &
       & Mass_Distribution_Beta_Profile_Potential=  &
       &   Mass_Distribution_Beta_Profile_Potential &
       &   *gravitationalConstantGalacticus
  return
end function Mass_Distribution_Beta_Profile_Potential

double precision function Mass_Distribution_Beta_Profile_Density_Radial_Moment(self,moment,radiusMinimum,radiusMaximum,isInfinite)
  !% Computes radial moments of the density in a $\beta$-profile mass distribution.
  use Hypergeometric_Functions
  implicit none
  class           (massDistributionBetaProfile), intent(in   )           :: self
  double precision                             , intent(in   )           :: moment
  double precision                             , intent(in   ), optional :: radiusMinimum, radiusMaximum
  logical                                      , intent(  out), optional :: isInfinite
  double precision                                                       :: fractionalRadiusMinimum, fractionalRadiusMaximum
  
  if (present(isInfinite)) isInfinite=.false.
  if (present(radiusMaximum)) then
     fractionalRadiusMaximum=radiusMaximum/self%coreRadius
     Mass_Distribution_Beta_Profile_Density_Radial_Moment=                   &
          & +fractionalRadiusMaximum**(moment+1.0d0)                         &
          & *Hypergeometric_2F1     (                                        &
          &                          [(moment+1.0d0)/2.0d0,1.5d0*self%beta], &
          &                          [(moment+3.0d0)/2.0d0                ], &
          &                          -fractionalRadiusMaximum**2             &
          &                         )                                        &
          & /                         (moment+1.0d0)
  else
     Mass_Distribution_Beta_Profile_Density_Radial_Moment=0.0d0
  end if
  if (present(radiusMinimum)) then
     fractionalRadiusMinimum=radiusMinimum/self%coreRadius
  else
     fractionalRadiusMinimum=0.0d0
  end if
  Mass_Distribution_Beta_Profile_Density_Radial_Moment=                   &
       & +Mass_Distribution_Beta_Profile_Density_Radial_Moment            &
       & -fractionalRadiusMinimum**(moment+1.0d0)                         &
       & *Hypergeometric_2F1     (                                        &
       &                          [(moment+1.0d0)/2.0d0,1.5d0*self%beta], &
       &                          [(moment+3.0d0)/2.0d0                ], &
       &                          -fractionalRadiusMinimum**2             &
       &                         )                                        &
       & /                         (moment+1.0d0)
  ! Convert to dimensionful units.
  Mass_Distribution_Beta_Profile_Density_Radial_Moment         &
       & =Mass_Distribution_Beta_Profile_Density_Radial_Moment &
       & *self%densityNormalization                            &
       & *self%coreRadius          **moment
  return
end function Mass_Distribution_Beta_Profile_Density_Radial_Moment
