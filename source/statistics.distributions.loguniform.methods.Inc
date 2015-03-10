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

!% Defines methods for a log-uniform distribution.

function distributionLogUniformConstructor(limitLower,limitUpper)
  !% Constructor for uniform distribution class.
  implicit none
  type            (distributionLogUniform)                :: distributionLogUniformConstructor
  double precision                        , intent(in   ) :: limitLower                    , limitUpper
  
  distributionLogUniformConstructor%limitLower=limitLower
  distributionLogUniformConstructor%limitUpper=limitUpper
  return
end function distributionLogUniformConstructor

double precision function distributionLogUniformDensity(self,x)
  !% Return the density of a uniform distribution.
  implicit none
  class           (distributionLogUniform), intent(in   ) :: self
  double precision                        , intent(in   ) :: x
  
  if (x < self%limitLower .or. x > self%limitUpper) then
     distributionLogUniformDensity=0.0d0
  else
     distributionLogUniformDensity=1.0d0/x/(log(self%limitUpper)-log(self%limitLower))
  end if
  return
end function distributionLogUniformDensity

double precision function distributionLogUniformCumulative(self,x)
  !% Return the cumulative probability of a uniform distribution.
  implicit none
  class           (distributionLogUniform), intent(in   ) :: self
  double precision                        , intent(in   ) :: x

  if      (x < self%limitLower) then
     distributionLogUniformCumulative=0.0d0
  else if (x > self%limitUpper) then
     distributionLogUniformCumulative=1.0d0
  else
     distributionLogUniformCumulative=(log(x)-log(self%limitLower))/(log(self%limitUpper)-log(self%limitLower))
  end if
  return
end function distributionLogUniformCumulative

double precision function distributionLogUniformInverse(self,p)
  !% Return the inverse of a uniform distribution.
  use Galacticus_Error
  implicit none
  class           (distributionLogUniform), intent(in   ) :: self
  double precision                        , intent(in   ) :: p

  if (p < 0.0d0 .or. p > 1.0d0)                                           &
       & call Galacticus_Error_Report(                                    &
       &                              'distributionLogUniformCumulative', &
       &                              'probability out of range'          &
       &                             )
  distributionLogUniformInverse=exp(log(self%limitLower)+p*(log(self%limitUpper)-log(self%limitLower)))
  return
end function distributionLogUniformInverse
