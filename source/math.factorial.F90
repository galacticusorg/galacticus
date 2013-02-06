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

!% Contains a module which implements calculations of factorials.

module Factorials
  !% Implements calculations of factorials
  implicit none
  private
  public :: Factorial, Logarithmic_Double_Factorial

contains

  double precision function Factorial(argument)
    !% Computes the factorial of {\tt argument}.
    use FGSL
    implicit none
    integer, intent(in) :: argument

    Factorial=FGSL_SF_Fact(argument)
    return
  end function Factorial

  double precision function Logarithmic_Double_Factorial(argument)
    !% Computes the natural logarithm of the double factorial, $k!!$.
    implicit none
    integer, intent(in) :: argument
    integer             :: i

    Logarithmic_Double_Factorial=0.0d0
    i=argument
    do while (i >= 2)
       Logarithmic_Double_Factorial=Logarithmic_Double_Factorial+dlog(dble(i))
       i=i-2
    end do
    return
  end function Logarithmic_Double_Factorial
  
end module Factorials
