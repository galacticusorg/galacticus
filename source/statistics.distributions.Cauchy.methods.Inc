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

!% Defines methods for a Cauchy distribution.

function distributionCauchyConstructor(median,scale)
  !% Constructor for Cauchy distribution class.
  implicit none
  type            (distributionCauchy)                :: distributionCauchyConstructor
  double precision                    , intent(in   ) :: median                       , scale
  
  distributionCauchyConstructor%median=median
  distributionCauchyConstructor%scale =scale
  return
end function distributionCauchyConstructor

double precision function distributionCauchyDensity(self,x)
  !% Return the density of a Cauchy distribution.
  use Numerical_Constants_Math
  implicit none
  class           (distributionCauchy), intent(in   ) :: self
  double precision                    , intent(in   ) :: x
  
  distributionCauchyDensity=1.0d0/Pi/self%scale/(1.0d0+((x-self%median)/self%scale)**2)
  return
end function distributionCauchyDensity

double precision function distributionCauchyCumulative(self,x)
  !% Return the cumulative probability of a Cauchy distribution.
  use Numerical_Constants_Math
  implicit none
  class           (distributionCauchy), intent(in   ) :: self
  double precision                    , intent(in   ) :: x

  distributionCauchyCumulative=0.5d0+atan((x-self%median)/self%scale)/Pi
  return
end function distributionCauchyCumulative

double precision function distributionCauchyInverse(self,p)
  !% Return the inverse of a Cauchy distribution.
  use Galacticus_Error
  use Numerical_Constants_Math
  implicit none
  class           (distributionCauchy), intent(in   ) :: self
  double precision                    , intent(in   ) :: p

  if (p < 0.0d0 .or. p > 1.0d0)                                    &
       & call Galacticus_Error_Report(                             &
       &                              'distributionCauchyInverse', &
       &                              'probability out of range'   &
       &                             )
  distributionCauchyInverse=self%median+self%scale*tan(Pi*(p-0.5d0))
  return
end function distributionCauchyInverse
