!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
!!    Andrew Benson <abenson@carnegiescience.edu>
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

!% Defines methods for a log-normal distribution with upper and lower limits.

function distributionLogNormalConstructor(mean,variance,limitLower,limitUpper)
  !% Constructor for log-normal distribution class.
  use Error_Functions
  implicit none
  type            (distributionLogNormal)                          :: distributionLogNormalConstructor
  double precision                       , intent(in   )           :: mean                            , variance
  double precision                       , intent(in   ), optional :: limitLower                      , limitUpper
  double precision                                                 :: meanNormal                      , varianceNormal

  varianceNormal=+log(               &
       &              +(             &
       &                +variance    &
       &                /mean    **2 &
       &               )             &
       &              +1.0d0         &
       &             )
  meanNormal    =+log(               &
       &              +mean          &
       &             )               &
       &         -    varianceNormal &
       &         /2.0d0
  if (present(limitLower)) then
     if (present(limitUpper)) then
        distributionLogNormalConstructor%distributionNormal=distributionNormal(meanNormal,varianceNormal,limitLower=log(limitLower),limitUpper=log(limitUpper))
     else
        distributionLogNormalConstructor%distributionNormal=distributionNormal(meanNormal,varianceNormal,limitLower=log(limitLower)                           )
     end if
  else
     if (present(limitUpper)) then
        distributionLogNormalConstructor%distributionNormal=distributionNormal(meanNormal,varianceNormal                           ,limitUpper=log(limitUpper))
     else
        distributionLogNormalConstructor%distributionNormal=distributionNormal(meanNormal,varianceNormal                                                      )
     end if
  end if
  return
end function distributionLogNormalConstructor

double precision function distributionLogNormalMinimum(self)
  !% Return the minimum possible value of a uniform distribution.
  implicit none
  class(distributionLogNormal), intent(in   ) :: self

  distributionLogNormalMinimum=exp(self%distributionNormal%minimum())
  return
end function distributionLogNormalMinimum

double precision function distributionLogNormalMaximum(self)
  !% Return the maximum possible value of a uniform distribution.
  implicit none
  class(distributionLogNormal), intent(in   ) :: self

  distributionLogNormalMaximum=exp(self%distributionNormal%maximum())
  return
end function distributionLogNormalMaximum

double precision function distributionLogNormalDensity(self,x)
  !% Return the density of a normal distribution.
  implicit none
  class           (distributionLogNormal), intent(in   ) :: self
  double precision                       , intent(in   ) :: x

  distributionLogNormalDensity=+self%distributionNormal%density(log(x)) &
       &                       /x
  return
end function distributionLogNormalDensity

double precision function distributionLogNormalCumulative(self,x)
  !% Return the cumulative probability of a normal distribution.
  implicit none
  class           (distributionLogNormal), intent(in   ) :: self
  double precision                       , intent(in   ) :: x

  distributionLogNormalCumulative=self%distributionNormal%cumulative(log(x))
  return
end function distributionLogNormalCumulative

double precision function distributionLogNormalInverse(self,p)
  !% Return the inverse of a normal distribution.
  implicit none
  class           (distributionLogNormal), intent(in   ) :: self
  double precision                       , intent(in   ) :: p

  distributionLogNormalInverse=exp(self%distributionNormal%inverse(p))
  return
end function distributionLogNormalInverse
