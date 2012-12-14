!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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
  implicit none
  class           (massDistributionBetaProfile), intent(  out)           :: self
  double precision                             , intent(in   )           :: beta,coreRadius
  double precision                             , intent(in   ), optional :: densityNormalization,mass,outerRadius
  logical                                      , intent(in   ), optional :: isDimensionless
  double precision                                                       :: r

  ! Store beta parameter.
  self%beta      =beta
  ! Store core radius.
  self%coreRadius=coreRadius
  ! Determine density normalization.
  if      (&
       &   present(densityNormalization)      &
       &  ) then
     self%densityNormalization=densityNormalization
  else if (&
       &   present(mass                ).and. &
       &   present(outerRadius         )      &
       &  ) then
     r=outerRadius/coreRadius
     self%densityNormalization=3.0d0*mass/4.0d0/Pi/outerRadius**3/Hypergeometric_2F1([1.5d0,1.5d0*beta],[2.5d0],-r**2)
  end if
  ! Determine if profile is dimensionless.
  if      (present(isDimensionless     )) then
     self%dimensionless=isDimensionless
  else
     self%dimensionless=.false.
  end if
  return
end subroutine Mass_Distribution_Beta_Profile_Initialize

double precision function Mass_Distribution_Beta_Profile_Density(self,coordinates)
  !% Return the density at the specified {\tt coordinates} in a $\beta$-profile mass distribution.
  use Coordinates
  implicit none
  class(massDistributionBetaProfile), intent(in) :: self
  class(coordinate                 ), intent(in) :: coordinates
  type (coordinateSpherical        )             :: position
  double precision                               :: r

  ! Get position in spherical coordinate system.
  position=coordinates
  ! Compute density.
  r=position%r()/self%coreRadius
  Mass_Distribution_Beta_Profile_Density=self%densityNormalization/(1.0d0+r**2)**(1.5d0*self%beta)
  return
end function Mass_Distribution_Beta_Profile_Density
