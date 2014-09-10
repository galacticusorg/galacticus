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

!% Defines methods for a Voight distribution with upper and lower limits.

function distributionVoightConstructor(gamma,mu,sigma,limitLower,limitUpper)
  !% Constructor for Voight distribution class.
  use Error_Functions
  implicit none
  type            (distributionVoight)                          :: distributionVoightConstructor
  double precision                    , intent(in   )           :: gamma                        , mu        , &
       &                                                           sigma
  double precision                    , intent(in   ), optional :: limitLower                   , limitUpper
  double precision                                              :: cdfLower                     , cdfUpper

  distributionVoightConstructor%gamma=gamma
  distributionVoightConstructor%mu   =mu
  distributionVoightConstructor%sigma=sigma
  distributionVoightConstructor%limitLowerExists=.false.
  distributionVoightConstructor%limitUpperExists=.false.
  distributionVoightConstructor%cdfAtLowerLimit =0.0d0
  distributionVoightConstructor%cdfAtUpperLimit =1.0d0
  if (present(limitLower)) cdfLower=distributionVoightConstructor%cumulative(limitLower)
  if (present(limitUpper)) cdfUpper=distributionVoightConstructor%cumulative(limitUpper)
  distributionVoightConstructor%limitLowerExists=present(limitLower)
  distributionVoightConstructor%limitUpperExists=present(limitUpper)
  if (distributionVoightConstructor%limitLowerExists) then
     distributionVoightConstructor%limitLower     =limitLower
     distributionVoightConstructor%cdfAtLowerLimit=cdfLower
  end if
  if (distributionVoightConstructor%limitUpperExists) then
     distributionVoightConstructor%limitUpper     =limitUpper
     distributionVoightConstructor%cdfAtUpperLimit=cdfUpper
  end if
  return
end function distributionVoightConstructor

double precision function distributionVoightDensity(self,x)
  !% Return the density of a Voight distribution.
  use Numerical_Constants_Math
  use Error_Functions
  implicit none
  class           (distributionVoight), intent(in   ) :: self
  double precision                    , intent(in   ) :: x
  double precision                                    :: x0
  double complex                                      :: z   , w
  
  if     (                                                   &
       &   (self%limitLowerExists .and. x < self%limitLower) &
       &  .or.                                               &
       &   (self%limitUpperExists .and. x > self%limitUpper) &
       & ) then
     distributionVoightDensity=0.0d0
  else
     ! Compute the value of x relative to the mean of the Gaussian component.
     x0=x-self%mu
     ! Evaluate the Feddeeva function at w(z).
     z =dcmplx(x0,self%gamma)/sqrt(2.0d0)/self%sigma
     w =Faddeeva(z)
     ! Compute the density.
     distributionVoightDensity=real(w)/sqrt(2.0d0*Pi)/self%sigma/(self%cdfAtUpperLimit-self%cdfAtLowerLimit)
  end if
  return
end function distributionVoightDensity

double precision function distributionVoightCumulative(self,x)
  !% Return the cumulative probability of a Voight distribution.
  use Numerical_Constants_Math
  use Hypergeometric_Functions
  use Error_Functions
  implicit none
  class           (distributionVoight), intent(in   )               :: self
  double precision                    , intent(in   )               :: x
  double complex                      , parameter    , dimension(2) :: a=[dcmplx(1.0d0,0.0d0),dcmplx(1.0d0,0.0d0)]
  double complex                      , parameter    , dimension(2) :: b=[dcmplx(1.5d0,0.0d0),dcmplx(2.0d0,0.0d0)]
  double precision                                                  :: x0
  double complex                                                    :: z

  if      (self%limitLowerExists .and. x < self%limitLower) then
     distributionVoightCumulative=0.0d0
  else if (self%limitUpperExists .and. x > self%limitUpper) then
     distributionVoightCumulative=1.0d0
  else
     ! Compute the value of x relative to the mean of the Gaussian component.
     x0=x-self%mu
     ! Evaluate z.
     z =dcmplx(x0,self%gamma)/sqrt(2.0d0)/self%sigma
     ! Evaluate the cumulative distribution.
     distributionVoightCumulative=                 &
          & +(                                     &
          &   +real(                               &
          &         +0.5d0                         &
          &         +0.5d0                         &
          &         *Error_Function(z)             &
          &         +dcmplx(0.0d0,1.0d0)           &
          &         *z**2                          &
          &         /Pi                            &
          &         *Hypergeometric_pFq(a,b,-z**2) &
          &        )                               &
          &   -self%cdfAtLowerLimit                &
          &  )                                     &
          & /(                                     &
          &   +self%cdfAtUpperLimit                &
          &   -self%cdfAtLowerLimit                &
          &  )
  end if
  return
end function distributionVoightCumulative
