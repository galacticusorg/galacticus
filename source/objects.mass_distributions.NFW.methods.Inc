!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu

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
