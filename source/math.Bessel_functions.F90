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

!% Contains a module which implements calculations of Bessel functions.

module Bessel_Functions
  !% Implements calculations of Bessel functions.
  use FGSL
  implicit none
  private
  public :: Bessel_Function_K0, Bessel_Function_K1, Bessel_Function_I0, Bessel_Function_I1

contains

  double precision function Bessel_Function_K0(argument)
    !% Computes the $K_0$ Bessel function.
    implicit none
    double precision, intent(in   ) :: argument

    Bessel_Function_K0=FGSL_SF_Bessel_Kc0(argument)
    return
  end function Bessel_Function_K0

  double precision function Bessel_Function_K1(argument)
    !% Computes the $K_1$ Bessel function.
    implicit none
    double precision, intent(in   ) :: argument

    Bessel_Function_K1=FGSL_SF_Bessel_Kc1(argument)
    return
  end function Bessel_Function_K1

  double precision function Bessel_Function_I0(argument)
    !% Computes the $I_0$ Bessel function.
    implicit none
    double precision, intent(in   ) :: argument

    Bessel_Function_I0=FGSL_SF_Bessel_Ic0(argument)
    return
  end function Bessel_Function_I0

  double precision function Bessel_Function_I1(argument)
    !% Computes the $I_1$ Bessel function.
    implicit none
    double precision, intent(in   ) :: argument

    Bessel_Function_I1=FGSL_SF_Bessel_Ic1(argument)
    return
  end function Bessel_Function_I1

end module Bessel_Functions
