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

!+    Contributions to this file made by:  Stéphane Mangeon, Andrew Benson.


!% Contains a module which implements a zero black hole binary recoil velocity.

module Black_Hole_Recoil_Velocities_Null
  !% Implements a zero black hole binary recoil velocity.
  implicit none
  private
  public :: Black_Hole_Binary_Recoil_Velocity_Null_Initialize

contains

  !# <blackHoleBinaryRecoilVelocityMethod>
  !#  <unitName>Black_Hole_Binary_Recoil_Velocity_Null_Initialize</unitName>
  !# </blackHoleBinaryRecoilVelocityMethod>
  subroutine Black_Hole_Binary_Recoil_Velocity_Null_Initialize(blackHoleBinaryRecoilVelocityMethod,Black_Hole_Binary_Recoil_Velocity_Get)
    !% Test if this method is to be used and set procedure pointer appropriately.
    use ISO_Varying_String
    implicit none
    type     (varying_string  ), intent(in   )          :: blackHoleBinaryRecoilVelocityMethod
    procedure(double precision), intent(inout), pointer :: Black_Hole_Binary_Recoil_Velocity_Get

    if (blackHoleBinaryRecoilVelocityMethod == 'null') Black_Hole_Binary_Recoil_Velocity_Get => Black_Hole_Binary_Recoil_Velocity_Null
    return
  end subroutine Black_Hole_Binary_Recoil_Velocity_Null_Initialize

  double precision function Black_Hole_Binary_Recoil_Velocity_Null(massBlackHole1,massBlackHole2,spinBlackHole1,spinBlackHole2)
    !% Returns a zero recoil velocity for black hole binary mergers.
    implicit none
    double precision, intent(in   ) :: massBlackHole1, massBlackHole2, spinBlackHole1, spinBlackHole2

    Black_Hole_Binary_Recoil_Velocity_Null=0.0d0
    return
  end function Black_Hole_Binary_Recoil_Velocity_Null

end module Black_Hole_Recoil_Velocities_Null
