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

!% Defines methods for a Student t-distribution.

function distributionStudentTConstructor(degreesOfFreedom)
  !% Constructor for StudentT distribution class.
  implicit none
  type            (distributionStudentT)                :: distributionStudentTConstructor
  double precision                      , intent(in   ) :: degreesOfFreedom
  
  distributionStudentTConstructor%degreesOfFreedom=degreesOfFreedom
  return
end function distributionStudentTConstructor

double precision function distributionStudentTDensity(self,x)
  !% Return the density of a StudentT distribution.
  use FGSL
  implicit none
  class           (distributionStudentT), intent(in   ) :: self
  double precision                      , intent(in   ) :: x
  
  distributionStudentTDensity=FGSL_Ran_tDist_PDF(x,self%degreesOfFreedom)
  return
end function distributionStudentTDensity

double precision function distributionStudentTCumulative(self,x)
  !% Return the cumulative probability of a StudentT distribution.
  use FGSL
  implicit none
  class           (distributionStudentT), intent(in   ) :: self
  double precision                      , intent(in   ) :: x

  distributionStudentTCumulative=FGSL_CDF_tDist_P(x,self%degreesOfFreedom)
  return
end function distributionStudentTCumulative

double precision function distributionStudentTInverse(self,p)
  !% Return the cumulative probability of a StudentT distribution.
  use Galacticus_Error
  use FGSL
  implicit none
  class           (distributionStudentT), intent(in   ) :: self
  double precision                      , intent(in   ) :: p

  if (p < 0.0d0 .or. p > 1.0d0)                                      &
       & call Galacticus_Error_Report(                               &
       &                              'distributionStudentTInverse', &
       &                              'probability out of range'     &
       &                             )
  distributionStudentTInverse=FGSL_CDF_tDist_Pinv(p,self%degreesOfFreedom)
  return
end function distributionStudentTInverse
