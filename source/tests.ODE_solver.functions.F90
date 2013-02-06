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

!% Contains a module of ODEs for unit tests.

module Test_ODE_Solver_Functions
  !% Contains ODEs for unit tests.
  use FGSL
  use, intrinsic :: ISO_C_Binding
  implicit none
  private
  public :: ODE_Set_1, ODE_Set_2

contains
 
  function ODE_Set_1(x,y,dydx,parameterPointer) bind(c)
    !% A set of ODEs for unit tests.
    integer(c_int)                           :: ODE_Set_1
    real(c_double), value                    :: x
    real(c_double), dimension(1), intent(in) :: y
    real(c_double), dimension(1)             :: dydx
    type(c_ptr),    value                    :: parameterPointer

    dydx(1)=sin(x)
    ODE_Set_1=FGSL_Success
  end function ODE_Set_1
  
  function ODE_Set_2(x,y,dydx,parameterPointer) bind(c)
    !% A set of ODEs for unit tests.
    integer(c_int)                           :: ODE_Set_2
    real(c_double), value                    :: x
    real(c_double), dimension(2), intent(in) :: y
    real(c_double), dimension(2)             :: dydx
    type(c_ptr),    value                    :: parameterPointer

    dydx(1)=y(2)
    dydx(2)=-1.0d0*y(1)
    ODE_Set_2=FGSL_Success
  end function ODE_Set_2

end module Test_ODE_Solver_Functions
