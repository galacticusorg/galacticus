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

!% Contains a module which implements exponential integrals.

module Exponential_Integrals
  !% Implements exponential integrals.
  use FGSL
  implicit none
  private
  public :: Sine_Integral, Cosine_Integral
  
contains
  
  double precision function Sine_Integral(x)
    !% Evaluate the $\hbox{Si}(x)\equiv\int_0^x \d t \sin(t)/t$ sine integral.
    implicit none
    double precision, intent(in) :: x

    Sine_Integral=FGSL_SF_Si(x)
    return
  end function Sine_Integral
  
  double precision function Cosine_Integral(x)
    !% Evaluate the $\hbox{Ci}(x)\equiv\int_0^x \d t \cos(t)/t$ cosine integral.
    implicit none
    double precision, intent(in) :: x

    Cosine_Integral=FGSL_SF_Ci(x)
    return
  end function Cosine_Integral
  
end module Exponential_Integrals
