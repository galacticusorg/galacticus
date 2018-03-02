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

!% Defines methods for a distribution for cosmological density fields with normally distributed amplitude, but conditioned to have
!% not collapse on a larger scale \citep[][eqn.~9]{mo_analytic_1996}.

function distributionPeakBackgroundConstructor(varianceBackground,thresholdCollapse) result(self)
  !% Constructor for {\normalfont \ttfamily peakBackground} distribution class.
  use Error_Functions
  implicit none
  type            (distributionPeakBackground)                :: self
  double precision                            , intent(in   ) :: varianceBackground        , thresholdCollapse
  double precision                            , parameter     :: rangeTable          =7.0d0
  integer                                     , parameter     :: cdfCount            =1000
  integer                                                     :: i

  self%varianceBackground=varianceBackground
  self%thresholdCollapse =thresholdCollapse
  ! Compute normalization of the distribution.
  self%normalization     =+1.0d0                                                                  &
       &                  /Error_Function(thresholdCollapse/sqrt(2.0d0)/sqrt(varianceBackground))
  ! Tabulate the cumulative distribution function.
  call self%cdf%create(-rangeTable*sqrt(varianceBackground),thresholdCollapse,cdfCount)
  do i=1,cdfCount
     call self%cdf%populate(                                                                                            &
          &                 +self%normalization                                                                         &
          &                 *0.50d0                                                                                     &
          &                 *(                                                                                          &
          &                   +Error_Function((+self%cdf%x(i)                        )/sqrt(2.0d0*varianceBackground))  &
          &                   -Error_Function((+self%cdf%x(i)-2.0d0*thresholdCollapse)/sqrt(2.0d0*varianceBackground))  &
          &                  )                                                                                        , &
          &                 i                                                                                           &
          &                )
  end do
  call self%cdf%reverse(self%cdfInverse)
  return

contains

  double precision function cdfIntegrand(x)
    !% The integrand for the cumulative distribution function.
    implicit none
    double precision, intent(in   ) :: x

    cdfIntegrand=self%density(x)
    return
  end function cdfIntegrand
  
end function distributionPeakBackgroundConstructor

double precision function distributionPeakBackgroundMinimum(self)
  !% Return the minimum possible value of a uniform distribution.
  use Galacticus_Error
  implicit none
  class(distributionPeakBackground), intent(inout) :: self
  !GCC$ attributes unused :: self
  
  distributionPeakBackgroundMinimum=0.0d0
  call Galacticus_Error_Report('no minimum exists'//{introspection:location})
  return
end function distributionPeakBackgroundMinimum

double precision function distributionPeakBackgroundMaximum(self)
  !% Return the maximum possible value of a uniform distribution.
  use Galacticus_Error
  implicit none
  class(distributionPeakBackground), intent(inout) :: self
  !GCC$ attributes unused :: self
  
  call Galacticus_Error_Report('no maximum exists'//{introspection:location})
  return
end function distributionPeakBackgroundMaximum

double precision function distributionPeakBackgroundDensity(self,x)
  !% Return the density of a normal distribution.
  use Numerical_Constants_Math
  implicit none
  class           (distributionPeakBackground), intent(inout) :: self
  double precision                            , intent(in   ) :: x

  if (x < self%thresholdCollapse) then
     distributionPeakBackgroundDensity=+self%normalization                  &
          &                            /sqrt(                               &
          &                                  +2.0d0                         &
          &                                  *Pi                            &
          &                                 )                               &
          &                            *(                                   &
          &                              +exp(                              &
          &                                   -0.5d0                        &
          &                                   *x                        **2 &
          &                                   /  self%varianceBackground    &
          &                                  )                              &
          &                              -exp(                              &
          &                                   -0.5d0                        &
          &                                   *(                            &
          &                                     +x                          &
          &                                     -2.0d0                      &
          &                                     *self%thresholdCollapse     &
          &                                    )                        **2 &
          &                                   /  self%varianceBackground    &
          &                                  )                              &
          &                             )                                   &
          &                            /sqrt(                               &
          &                                  +   self%varianceBackground    &
          &                                 )
  else
     distributionPeakBackgroundDensity=+0.0d0
  end if
  return
end function distributionPeakBackgroundDensity

double precision function distributionPeakBackgroundCumulative(self,x)
  !% Return the cumulative probability of a normal distribution.
  implicit none
  class           (distributionPeakBackground), intent(inout) :: self
  double precision                            , intent(in   ) :: x

  if      (x < self%cdf%x(+1)) then
     distributionPeakBackgroundCumulative=0.0d0
  else if (x > self%cdf%x(-1)) then
     distributionPeakBackgroundCumulative=1.0d0
  else
     distributionPeakBackgroundCumulative=self%cdf%interpolate(x)
  end if
  return
end function distributionPeakBackgroundCumulative

double precision function distributionPeakBackgroundInverse(self,p)
  !% Return the inverse of a normal distribution.
  use Galacticus_Error
  implicit none
  class           (distributionPeakBackground), intent(inout) :: self
  double precision                            , intent(in   ) :: p

  if (p < 0.0d0 .or. p > 1.0d0)                                    &
       & call Galacticus_Error_Report(                             &
       &                              'probability out of range'// &
       &                              {introspection:location}     &
       &                             )
  distributionPeakBackgroundInverse=self%cdfInverse%interpolate(p)
  return
end function distributionPeakBackgroundInverse
