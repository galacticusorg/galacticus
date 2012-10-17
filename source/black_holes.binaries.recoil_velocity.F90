!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements calculations of black hole binary recoil velocities.

module Black_Hole_Binary_Recoil_Velocities
  !% Implements calculations of black hole binary recoil velocities.
  use ISO_Varying_String
  implicit none
  private
  public :: Black_Hole_Binary_Recoil_Velocity

  ! Flag to indicate if this module has been initialized.  
  logical                                               :: blackHoleBinaryRecoilVelocityInitialized=.false.

  ! Name of mass movement method used.
  type(varying_string)                                  :: blackHoleBinaryRecoilVelocityMethod

  ! Pointer to the subroutine that returns descriptors for mass movement.
  procedure(Black_Hole_Binary_Recoil_Velocity), pointer :: Black_Hole_Binary_Recoil_Velocity_Get => null()
  
contains

  double precision function Black_Hole_Binary_Recoil_Velocity(massBlackHole1,massBlackHole2,spinBlackHole1,spinBlackHole2)
    !% Computes the recoil velocity of a black hole binary.
    use Galacticus_Error
    use Input_Parameters
    use Tree_Nodes
    !# <include directive="blackHoleBinaryRecoilVelocityMethod" type="moduleUse">
    include 'black_holes.binary.recoil_velocity.modules.inc'
    !# </include>
    implicit none
    double precision, intent(in) :: massBlackHole1,massBlackHole2,spinBlackHole1,spinBlackHole2

    if (.not.blackHoleBinaryRecoilVelocityInitialized) then
       !$omp critical(blackHoleBinaryRecoilVelocityInitialize)
       if (.not.blackHoleBinaryRecoilVelocityInitialized) then
          ! Get the binary black hole recoil velocity method parameter.
          !@ <inputParameter>
          !@   <name>blackHoleBinaryRecoilVelocityMethod</name>
          !@   <defaultValue>null</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be used for computing the recoil velocity of black hole binaries.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('blackHoleBinaryRecoilVelocityMethod',blackHoleBinaryRecoilVelocityMethod,defaultValue='null')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="blackHoleBinaryRecoilVelocityMethod" type="code" action="subroutine">
          !#  <subroutineArgs>blackHoleBinaryRecoilVelocityMethod,Black_Hole_Binary_Recoil_Velocity_Get</subroutineArgs>
          include 'black_holes.binaries.recoil_velocity.inc'
          !# </include>
          if (.not.associated(Black_Hole_Binary_Recoil_Velocity_Get)) call Galacticus_Error_Report('Black_Hole_Binary_Recoil_Velocity','method ' &
               &//char(blackHoleBinaryRecoilVelocityMethod)//' is unrecognized')
          ! Flag that the module is now initialized.
          blackHoleBinaryRecoilVelocityInitialized=.true.
       end if
       !$omp end critical(blackHoleBinaryRecoilVelocityInitialize)
    end if

    ! Call the routine to do the calculation.
    Black_Hole_Binary_Recoil_Velocity=Black_Hole_Binary_Recoil_Velocity_Get(massBlackHole1,massBlackHole2,spinBlackHole1,spinBlackHole2)

    return
  end function Black_Hole_Binary_Recoil_Velocity
  
end module Black_Hole_Binary_Recoil_Velocities
