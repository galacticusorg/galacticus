!! Copyright 2009, 2010, Andrew Benson <abenson@caltech.edu>
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






!% Contains a module which implements hypergeometric functions.

module Hypergeometric_Functions
  !% Implements hypergeometric functions.
  use FGSL
  private
  public :: Hypergeometric_2F1
  
contains
  
  double precision function Hypergeometric_2F1(a,b,x)
    !% Evaluate the $_2F_1(a_1,a_2;b_1;x)$ hypergeometric function.
    implicit none
    double precision, intent(in) :: a(2),b(1),x

    Hypergeometric_2F1=FGSL_SF_Hyperg_2F1(a(1),a(2),b(1),x)
    return
  end function Hypergeometric_2F1
  
end module Hypergeometric_Functions
