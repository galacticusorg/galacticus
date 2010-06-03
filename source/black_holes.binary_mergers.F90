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






!% Contains a module which implements calculations of black hole binary mergers.

module Black_Hole_Binary_Mergers
  !% Implements calculations of black hole binary mergers.
  use ISO_Varying_String
  private
  public :: Black_Hole_Binary_Merger

  ! Flag to indicate if this module has been initialized.  
  logical                                      :: blackHoleBinaryMergersInitialized=.false.

  ! Name of mass movement method used.
  type(varying_string)                         :: blackHoleBinaryMergersMethod

  ! Pointer to the subroutine that returns descriptors for mass movement.
  procedure(Black_Hole_Binary_Merger), pointer :: Black_Hole_Binary_Merger_Do => null()
  
contains

  subroutine Black_Hole_Binary_Merger(blackHoleMassA,blackHoleMassB,blackHoleSpinA,blackHoleSpinB,blackHoleMassFinal&
       &,blackHoleSpinFinal)
    !% Computes the effects of a black hole binary merger.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="blackHoleBinaryMergersMethod" type="moduleUse">
    include 'black_holes.binary_mergers.modules.inc'
    !# </include>
    implicit none
    double precision, intent(in)  :: blackHoleMassA,blackHoleMassB,blackHoleSpinA,blackHoleSpinB
    double precision, intent(out) :: blackHoleMassFinal,blackHoleSpinFinal

    !$omp critical(blackHoleBinaryMergersInitialize)
    if (.not.blackHoleBinaryMergersInitialized) then
       ! Do the binary black hole merger method parameter.
       !@ <inputParameter>
       !@   <name>blackHoleBinaryMergersMethod</name>
       !@   <defaultValue>Rezzolla2008</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of the method to be used for computing the effects of black hole binary mergers.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('blackHoleBinaryMergersMethod',blackHoleBinaryMergersMethod,defaultValue='Rezzolla2008')
       ! Include file that makes calls to all available method initialization routines.
       !# <include directive="blackHoleBinaryMergersMethod" type="code" action="subroutine">
       !#  <subroutineArgs>blackHoleBinaryMergersMethod,Black_Hole_Binary_Merger_Do</subroutineArgs>
       include 'black_holes.binary_mergers.inc'
       !# </include>
       if (.not.associated(Black_Hole_Binary_Merger_Do)) call Galacticus_Error_Report('Black_Hole_Binary_Merger','method ' &
            &//char(blackHoleBinaryMergersMethod)//' is unrecognized')
       ! Flag that the module is now initialized.
       blackHoleBinaryMergersInitialized=.true.
    end if
    !$omp end critical(blackHoleBinaryMergersInitialize)

    ! Call the routine to do the calculation.
    call Black_Hole_Binary_Merger_Do(blackHoleMassA,blackHoleMassB,blackHoleSpinA,blackHoleSpinB,blackHoleMassFinal&
       &,blackHoleSpinFinal)

    return
  end subroutine Black_Hole_Binary_Merger
  
end module Black_Hole_Binary_Mergers
