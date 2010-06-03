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






!% Contains a module which implements calculations of gamma functions.

module Gamma_Functions
  !% Implements calculations of gamma functions.
  private
  public :: Gamma_Function_Incomplete_Complementary

contains

  double precision function Gamma_Function_Incomplete_Complementary(exponent,argument)
    !% Computes the complementary incomplete gamma function.
    use FGSL
    implicit none
    double precision, intent(in) :: exponent,argument

    Gamma_Function_Incomplete_Complementary=FGSL_SF_Gamma_Inc_P(exponent,argument)
    return
  end function Gamma_Function_Incomplete_Complementary
  
end module Gamma_Functions
