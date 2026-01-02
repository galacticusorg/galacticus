!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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

!!{
Contains a module which implements Gaussian distributions.
!!}

module Math_Distributions_Gaussian
  !!{
  Implements Gaussian distributions.
  !!}
  private
  public :: Gaussian_Distribution

contains

  double precision function Gaussian_Distribution(x,sigma)
    !!{
    Computes the Gaussian distribution with dispersion {\normalfont \ttfamily sigma} at argument {\normalfont \ttfamily x}.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision, intent(in   ) :: sigma, x

    ! Evaluate the Gaussian distribution.
    Gaussian_Distribution=exp(-0.5*(x/sigma)**2)/sqrt(2.0d0*Pi)/sigma
    return
  end function Gaussian_Distribution

end module Math_Distributions_Gaussian
