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

!% Defines methods for a uniform distribution.

function distributionUniformConstructor(limitLower,limitUpper)
  !% Constructor for uniform distribution class.
  implicit none
  type            (distributionUniform)                :: distributionUniformConstructor
  double precision                     , intent(in   ) :: limitLower                    , limitUpper
  
  distributionUniformConstructor%limitLower=limitLower
  distributionUniformConstructor%limitUpper=limitUpper
  return
end function distributionUniformConstructor

double precision function distributionUniformDensity(self,x)
  !% Return the density of a uniform distribution.
  implicit none
  class           (distributionUniform), intent(in   ) :: self
  double precision                     , intent(in   ) :: x
  
  if (x < self%limitLower .or. x > self%limitUpper) then
     distributionUniformDensity=0.0d0
  else
     distributionUniformDensity=1.0d0/(self%limitUpper-self%limitLower)
  end if
  return
end function distributionUniformDensity

double precision function distributionUniformCumulative(self,x)
  !% Return the cumulative probability of a uniform distribution.
  implicit none
  class           (distributionUniform), intent(in   ) :: self
  double precision                     , intent(in   ) :: x

  if      (x < self%limitLower) then
     distributionUniformCumulative=0.0d0
  else if (x > self%limitUpper) then
     distributionUniformCumulative=1.0d0
  else
     distributionUniformCumulative=(x-self%limitLower)/(self%limitUpper-self%limitLower)
  end if
  return
end function distributionUniformCumulative

double precision function distributionUniformInverse(self,p)
  !% Return the inverse of a uniform distribution.
  use Galacticus_Error
  implicit none
  class           (distributionUniform), intent(in   ) :: self
  double precision                     , intent(in   ) :: p

  if (p < 0.0d0 .or. p > 1.0d0)                                        &
       & call Galacticus_Error_Report(                                 &
       &                              'distributionUniformCumulative', &
       &                              'probability out of range'       &
       &                             )
  distributionUniformInverse=self%limitLower+p*(self%limitUpper-self%limitLower)
  return
end function distributionUniformInverse
