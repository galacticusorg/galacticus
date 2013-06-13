!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

! Hernquist functions.
subroutine Mass_Distribution_Hernquist_Initialize(self,scaleLength,densityNormalization,mass,isDimensionless)
  !% Initialize a Hernquist mass distribution.
  use Numerical_Constants_Math
  use Numerical_Comparison
  use Galacticus_Error
  implicit none
  class           (massDistributionHernquist), intent(  out)           :: self
  double precision                           , intent(in   ), optional :: scaleLength,densityNormalization,mass
  logical                                    , intent(in   ), optional :: isDimensionless

  ! Determine if profile is dimensionless.
  if      (present(isDimensionless     )) then
     self%dimensionless=isDimensionless
  else
     self%dimensionless=.false.
  end if
  ! If dimensionless, then set scale length and mass to unity.
  if (self%dimensionless) then
     if (present(scaleLength         )) then
        if (Values_Differ(scaleLength         ,1.0d0,absTol=1.0d-6)) call Galacticus_Error_Report('Mass_Distribution_Hernquist_Initialize','scaleLength should be unity for a dimensionless profile (or simply do not specify a scale length)')
     end if
     if (present(mass                )) then
        if (Values_Differ(mass                ,1.0d0,absTol=1.0d-6)) call Galacticus_Error_Report('Mass_Distribution_Hernquist_Initialize','mass should be unity for a dimensionless profile (or simply do not specify a mass)')
     end if
     if (present(densityNormalization)) then
        if (Values_Differ(densityNormalization,1.0d0,absTol=1.0d-6)) call Galacticus_Error_Report('Mass_Distribution_Hernquist_Initialize','densityNormalization should be 1/2π for a dimensionless profile (or simply do not specify a densityNormalization)')
     end if
     self%scaleLength         =1.0d0
     self%mass                =1.0d0
     self%densityNormalization=0.5d0/Pi
  else
     if      (present(scaleLength         )) then
        self%scaleLength=scaleLength
     else
        call Galacticus_Error_Report('Mass_Distribution_Hernquist_Initialize','"scaleLength" must be specified')
     end if
     if      (present(densityNormalization)) then
        self%densityNormalization=densityNormalization
        self%mass                =densityNormalization*2.0d0*Pi*scaleLength**3
     else if (present(mass                )) then
        self%densityNormalization=mass/scaleLength**3/2.0d0/Pi
        self%mass                =mass
     else
        call Galacticus_Error_Report('Mass_Distribution_Hernquist_Initialize','one of "densityNormalization" or "mass" must be specified')
     end if
  end if
  return
end subroutine Mass_Distribution_Hernquist_Initialize

double precision function Mass_Distribution_Hernquist_Density(self,coordinates)
  !% Return the density at the specified {\tt coordinates} in a Hernquist mass distribution.
  use Coordinates
  implicit none
  class           (massDistributionHernquist), intent(in   ) :: self
  class           (coordinate               ), intent(in   ) :: coordinates
  type            (coordinateSpherical      )                :: position
  double precision                                           :: r

  ! Get position in spherical coordinate system.
  position=coordinates
  ! Compute the density at this position.
  r=position%r()/self%scaleLength
  Mass_Distribution_Hernquist_Density=self%densityNormalization/r/(1.0d0+r)**3
  return
end function Mass_Distribution_Hernquist_Density

double precision function Mass_Distribution_Hernquist_Density_Radial_Moment(self,moment,isInfinite)
  !% Returns a radial density moment for the Hernquist mass distribution.
  use Numerical_Constants_Math
  use Numerical_Comparison
  use Galacticus_Error
  implicit none
  class(massDistributionHernquist), intent(in   )           :: self
  double precision                , intent(in   )           :: moment
  logical                         , intent(  out), optional :: isInfinite

  if (moment <= 0.0d0 .or. moment >= 3.0d0) then
     ! Handle cases where the moment is infinite.
     if (present(isInfinite)) then
        isInfinite=.true.
        return
     else
        call Galacticus_Error_Report('Mass_Distribution_Hernquist_Density_Radial_Moment','requested radial density moment is infinite')
     end if
  else if (                                              &
       &   Values_Agree(moment,1.0d0,absTol=1.0d-3) .or. &
       &   Values_Agree(moment,2.0d0,absTol=1.0d-3)      &
       &  ) then
     ! Handle special cases.
     if (present(isInfinite)) isInfinite=.false.
     Mass_Distribution_Hernquist_Density_Radial_Moment=0.5d0
  else
     ! Compute the general case.
     if (present(isInfinite)) isInfinite=.false.
     Mass_Distribution_Hernquist_Density_Radial_Moment=0.5d0*Pi*(moment-1.0d0)*(moment-2.0d0)/sin(Pi*moment)
  end if
  Mass_Distribution_Hernquist_Density_Radial_Moment=Mass_Distribution_Hernquist_Density_Radial_Moment*self%densityNormalization&
       &*self%scaleLength**moment
  return
end function Mass_Distribution_Hernquist_Density_Radial_Moment

double precision function Mass_Distribution_Hernquist_Mass_Enc_By_Sphere(self,radius)
  !% Computes the mass enclosed within a sphere of given {\tt radius} for Hernquist mass distributions.
  use Numerical_Constants_Math
  implicit none
  class(massDistributionHernquist), intent(in   ), target :: self
  double precision                , intent(in   )         :: radius
  double precision                                        :: fractionalRadius
  
  fractionalRadius=radius/self%scaleLength
  Mass_Distribution_Hernquist_Mass_Enc_By_Sphere=2.0d0*Pi*self%densityNormalization*fractionalRadius**2/(1.0d0+fractionalRadius)**2
  return
end function Mass_Distribution_Hernquist_Mass_Enc_By_Sphere

double precision function Mass_Distribution_Hernquist_Potential(self,coordinates)
  !% Return the potential at the specified {\tt coordinates} in a Hernquist mass distribution.
  use Numerical_Constants_Physical
  use Coordinates
  implicit none
  class           (massDistributionHernquist), intent(in   ) :: self
  class           (coordinate               ), intent(in   ) :: coordinates
  type            (coordinateSpherical      )                :: position
  double precision                                           :: r

  ! Get position in spherical coordinate system.
  position=coordinates
  ! Compute the potential at this position.
  Mass_Distribution_Hernquist_Potential=-self%mass/(self%scaleLength+position%r())
  if (.not.self%isDimensionless()) Mass_Distribution_Hernquist_Potential=gravitationalConstantGalacticus&
       &*Mass_Distribution_Hernquist_Potential
  return
end function Mass_Distribution_Hernquist_Potential

double precision function Mass_Distribution_Hernquist_Half_Mass_Radius(self)
  !% Return the half-mass radius of a Hernquist mass distribution.
  implicit none
  class           (massDistributionHernquist), intent(in   ) :: self
  double precision                           , parameter     :: halfMassRadiusToScaleRadius=1.0d0/(sqrt(2.0d0)-1.0d0)

  Mass_Distribution_Hernquist_Half_Mass_Radius=halfMassRadiusToScaleRadius*self%scaleLength
  return
end function Mass_Distribution_Hernquist_Half_Mass_Radius
