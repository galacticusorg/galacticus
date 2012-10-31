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

! NFW functions.
subroutine Mass_Distribution_NFW_Initialize(self,scaleLength,concentration,densityNormalization,mass,virialRadius,isDimensionless)
  !% Initialize an NFW mass distribution.
  use Numerical_Constants_Math
  use Galacticus_Error
  implicit none
  class           (massDistributionNFW), intent(  out)           :: self
  double precision                     , intent(in   ), optional :: scaleLength,concentration,densityNormalization,mass &
       &,virialRadius
  logical                              , intent(in   ), optional :: isDimensionless
  double precision                                               :: r

  ! Determine scale length
  if      (                            &
       &   present(scaleLength  )      &
       &  ) then
     self%scaleLength         =scaleLength
  else if (                            &
       &   present(concentration).and. &
       &   present(virialRadius )      &
       &  ) then
     self%scaleLength=virialRadius/concentration
  else
     call Galacticus_Error_Report('Mass_Distribution_NFW_Initialize','no means to determine scale length')
  end if
  ! Determine density normalization.
  if      (                                   &
       &   present(densityNormalization)      &
       &  ) then
     self%densityNormalization=densityNormalization
  else if (                                   &
       &   present(mass                ).and. &
       &   present(virialRadius        )      &
       &  ) then
     r=virialRadius/self%scaleLength
     self%densityNormalization=mass/4.0d0/Pi/self%scaleLength**3/(log(1.0d0+r)-r/(1.0d0+r))
  else
     call Galacticus_Error_Report('Mass_Distribution_NFW_Initialize','either "densityNormalization", or "mass" and "virialRadius" must be specified')
  end if
  ! Determine if profile is dimensionless.
  if      (present(isDimensionless     )) then
     self%dimensionless=isDimensionless
  else
     self%dimensionless=.false.
  end if
  return
end subroutine Mass_Distribution_NFW_Initialize

double precision function Mass_Distribution_NFW_Density(self,coordinates)
  !% Return the density at the specified {\tt coordinates} in an NFW mass distribution.
  use Coordinates
  implicit none
  class           (massDistributionNFW), intent(in   ) :: self
  class           (coordinate         ), intent(in   ) :: coordinates
  type            (coordinateSpherical)                :: position
  double precision                                     :: r

  ! Get position in spherical coordinate system.
  position=coordinates
  ! Compute the density at this position.
  r=position%r()/self%scaleLength
  Mass_Distribution_NFW_Density=self%densityNormalization/r/(1.0d0+r)**2
  return
end function Mass_Distribution_NFW_Density
