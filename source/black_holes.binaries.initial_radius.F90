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

!% Contains a module which implements calculations of black hole binary initial separations.

module Black_Hole_Binary_Initial_Radii
  !% Implements calculations of black hole binary initial separations.
  use ISO_Varying_String
  implicit none
  private
  public :: Black_Hole_Binary_Initial_Radius

  ! Flag to indicate if this module has been initialized.
  logical                                              :: blackHoleBinaryInitialRadiiInitialized=.false.

  ! Name of mass movement method used.
  type     (varying_string                  )          :: blackHoleBinaryInitialRadiiMethod

  ! Pointer to the subroutine that returns descriptors for mass movement.
  procedure(Black_Hole_Binary_Initial_Radius), pointer :: Black_Hole_Binary_Initial_Radius_Get  =>null()

contains

  double precision function Black_Hole_Binary_Initial_Radius(thisNode,hostNode)
    !% Computes the initial radius of a newly formed black hole binary.
    use Galacticus_Error
    use Input_Parameters
    use Galacticus_Nodes
    !# <include directive="blackHoleBinaryInitialRadiiMethod" type="moduleUse">
    include 'black_holes.binary.initial_radius.modules.inc'
    !# </include>
    implicit none
    type(treeNode), intent(inout), pointer :: hostNode, thisNode

    if (.not.blackHoleBinaryInitialRadiiInitialized) then
       !$omp critical(blackHoleBinaryInitialRadiiInitialize)
       if (.not.blackHoleBinaryInitialRadiiInitialized) then
          ! Get the binary black hole initial radii method parameter.
          !@ <inputParameter>
          !@   <name>blackHoleBinaryInitialRadiiMethod</name>
          !@   <defaultValue>spheroidRadiusFraction</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be used for computing the initial separation of black hole binaries.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('blackHoleBinaryInitialRadiiMethod',blackHoleBinaryInitialRadiiMethod,defaultValue='spheroidRadiusFraction')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="blackHoleBinaryInitialRadiiMethod" type="functionCall" functionType="void">
          !#  <functionArgs>blackHoleBinaryInitialRadiiMethod,Black_Hole_Binary_Initial_Radius_Get</functionArgs>
          include 'black_holes.binaries.initial_radius.inc'
          !# </include>
          if (.not.associated(Black_Hole_Binary_Initial_Radius_Get)) call Galacticus_Error_Report('Black_Hole_Binary_Initial_Radius','method ' &
               &//char(blackHoleBinaryInitialRadiiMethod)//' is unrecognized')
          ! Flag that the module is now initialized.
          blackHoleBinaryInitialRadiiInitialized=.true.
       end if
       !$omp end critical(blackHoleBinaryInitialRadiiInitialize)
    end if

    ! Call the routine to do the calculation.
    Black_Hole_Binary_Initial_Radius=Black_Hole_Binary_Initial_Radius_Get(thisNode,hostNode)

    return
  end function Black_Hole_Binary_Initial_Radius

end module Black_Hole_Binary_Initial_Radii
