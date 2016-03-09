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

!% Defines methods for a negative exponential distribution.

function distributionNegativeExponentialConstructor(rate)
  !% Constructor for negative exponential distribution class.
  implicit none
  type            (distributionNegativeExponential)                :: distributionNegativeExponentialConstructor
  double precision                                 , intent(in   ) :: rate
  
  distributionNegativeExponentialConstructor%rate=rate
  return
end function distributionNegativeExponentialConstructor

double precision function distributionNegativeExponentialDensity(self,x)
  !% Return the density of a negative exponential distribution.
  implicit none
  class           (distributionNegativeExponential), intent(in   ) :: self
  double precision                                 , intent(in   ) :: x
  
  if (x < 0.0d0) then
     distributionNegativeExponentialDensity=0.0d0
  else
     distributionNegativeExponentialDensity=self%rate*exp(-self%rate*x)
  end if
  return
end function distributionNegativeExponentialDensity

double precision function distributionNegativeExponentialCumulative(self,x)
  !% Return the cumulative probability of a negative exponential distribution.
  implicit none
  class           (distributionNegativeExponential), intent(in   ) :: self
  double precision                                 , intent(in   ) :: x

  if (x < 0.0d0) then
     distributionNegativeExponentialCumulative=0.0d0
  else
     distributionNegativeExponentialCumulative=1.0d0-exp(-self%rate*x)
  end if
  return
end function distributionNegativeExponentialCumulative

double precision function distributionNegativeExponentialInverse(self,p)
  !% Return the inverse of a negative exponential distribution.
  use Galacticus_Error
  implicit none
  class           (distributionNegativeExponential), intent(in   ) :: self
  double precision                                 , intent(in   ) :: p
  double precision                                 , parameter     :: pTiny=1.0d-6

  if (p < 0.0d0 .or. p > 1.0d0)                                                    &
       & call Galacticus_Error_Report(                                             &
       &                              'distributionNegativeExponentialCumulative', &
       &                              'probability out of range'                   &
       &                             )
  if (p > pTiny) then
     distributionNegativeExponentialInverse=-log(1.0d0-p)/self%rate
  else
     distributionNegativeExponentialInverse=+(            &
          &                                   +p**1/1.0d0 &
          &                                   +p**2/2.0d0 &
          &                                   +p**3/3.0d0 &
          &                                   +p**4/4.0d0 &
          &                                  )            &
          &                                 /self%rate
  end if

!! AJB HACK
if (distributionNegativeExponentialInverse <= 0.0d0) then
   write (0,*) "OUT OF RANGE NEG EXP ",p,distributionNegativeExponentialInverse
end if

  return
end function distributionNegativeExponentialInverse
