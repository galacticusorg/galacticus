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

!% Defines methods for a Gamma distribution with upper and lower limits.

function distributionGammaConstructor(shape,rate,limitLower,limitUpper)
  !% Constructor for Gamma distribution class.
  use Galacticus_Error
  implicit none
  type            (distributionGamma)                          :: distributionGammaConstructor
  double precision                   , intent(in   )           :: shape                        , rate
  double precision                   , intent(in   ), optional :: limitLower                   , limitUpper
  
  if (rate <= 0.0d0 .or. shape <= 0.0d0) call Galacticus_Error_Report('distributionGammaConstructor','rate>0 and shape>0 are required')
  distributionGammaConstructor%rate            =rate
  distributionGammaConstructor%shape           =shape
  distributionGammaConstructor%limitLowerExists=.false.
  distributionGammaConstructor%limitUpperExists=.false.
  distributionGammaConstructor%cdfAtLowerLimit =0.0d0
  distributionGammaConstructor%cdfAtUpperLimit =1.0d0
  if (present(limitLower)) then
     if (limitLower < 0.0d0) call Galacticus_Error_Report('distributionGammaConstructor','limitLower≥0 is required')
     distributionGammaConstructor%limitLower     =limitLower
     distributionGammaConstructor%cdfAtLowerLimit=distributionGammaCumulative(distributionGammaConstructor,limitLower)
  else
     distributionGammaConstructor%cdfAtLowerLimit=0.0d0
  end if
  if (present(limitUpper)) then
     if (limitUpper < 0.0d0) call Galacticus_Error_Report('distributionGammaConstructor','limitUpper≥0 is required')
     distributionGammaConstructor%limitUpper     =limitUpper
     distributionGammaConstructor%cdfAtUpperLimit=distributionGammaCumulative(distributionGammaConstructor,limitUpper)
  else
     distributionGammaConstructor%cdfAtUpperLimit=1.0d0
  end if
  distributionGammaConstructor%limitLowerExists=present(limitLower)
  distributionGammaConstructor%limitUpperExists=present(limitUpper)
  return
end function distributionGammaConstructor

double precision function distributionGammaDensity(self,x)
  !% Return the density of a Gamma distribution.
  use Gamma_Functions
  implicit none
  class           (distributionGamma), intent(in   ) :: self
  double precision                   , intent(in   ) :: x
  
  if     (                                                   &
       &                                x < 0.0d0            &
       &  .or.                                               &
       &   (self%limitLowerExists .and. x < self%limitLower) &
       &  .or.                                               &
       &   (self%limitUpperExists .and. x > self%limitUpper) &
       & ) then
     distributionGammaDensity=0.0d0
  else
     distributionGammaDensity=          &
          &  x**(self%shape-1.0d0)      &
          & *exp(-x*self%rate)          &
          & *self%rate**self%shape      &
          & /Gamma_Function(self%shape) &
          & /(                          &
          &   +self%cdfAtUpperLimit     &
          &   -self%cdfAtLowerLimit     &
          &  )
  end if
  return
end function distributionGammaDensity

double precision function distributionGammaCumulative(self,x)
  !% Return the cumulative probability of a Gamma distribution.
  use Gamma_Functions
  implicit none
  class           (distributionGamma), intent(in   ) :: self
  double precision                   , intent(in   ) :: x

  if      (                            x <= 0.0d0          ) then
     distributionGammaCumulative=0.0d0
  else if (self%limitLowerExists .and. x <  self%limitLower) then
     distributionGammaCumulative=0.0d0
  else if (self%limitUpperExists .and. x >  self%limitUpper) then
     distributionGammaCumulative=1.0d0
  else
     distributionGammaCumulative=                                           &
          & min(                                                            &
          &     1.0d0                                                     , &
          &     max(                                                        &
          &          0.0d0                                                , &
          &          (                                                      &
          &           +Gamma_Function_Incomplete_Complementary(             &
          &                                                     self%shape, &
          &                                                     self%rate   &
          &                                                    *x           &
          &                                                   )             &
          &           -self%cdfAtLowerLimit                                 &
          &          )                                                      &
          &         /(                                                      &
          &           +self%cdfAtUpperLimit                                 &
          &           -self%cdfAtLowerLimit                                 &
          &          )                                                      &
          &        )                                                        &
          &    )
  end if
  return
end function distributionGammaCumulative

double precision function distributionGammaInverse(self,p)
  !% Return the inverse of a Gamma distribution.
  use Galacticus_Error
  use Gamma_Functions
  implicit none
  class           (distributionGamma), intent(in   ) :: self
  double precision                   , intent(in   ) :: p

  if (p < 0.0d0 .or. p > 1.0d0)                                   &
       & call Galacticus_Error_Report(                            &
       &                              'distributionGammaInverse', &
       &                              'probability out of range'  &
       &                             )
  distributionGammaInverse=                                                        &
       & +Inverse_Gamma_Function_Incomplete_Complementary(                         &
       &                                                  self%shape             , &
       &                                                  +p                       &
       &                                                  *(                       &
       &                                                    +self%cdfAtUpperLimit  &
       &                                                    -self%cdfAtLowerLimit  &
       &                                                   )                       &
       &                                                  +  self%cdfAtLowerLimit  &
       &                                                 )                         &
       & /self%rate
  return
end function distributionGammaInverse
