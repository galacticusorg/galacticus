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

!% Contains a module which implements calculations of star formation timescales for galactic disks.

module Star_Formation_Timescales_Disks
  !% Implements calculations of star formation timescales for galactic disks.
  use, intrinsic :: ISO_C_Binding
  use ISO_Varying_String
  use Tree_Nodes
  implicit none
  private
  public :: Star_Formation_Timescale_Disk
  
  ! Flag to indicate if this module has been initialized.  
  logical              :: starFormationTimescaleDisksInitialized=.false.

  ! Name of cooling rate available method used.
  type(varying_string) :: starFormationTimescaleDisksMethod

  ! Flag indicating if the work procedure is Fortran.
  logical              :: functionIsFortran

  ! Pointer to the function that actually does the calculation.
  procedure(Star_Formation_Timescale_Disk_FTemplate), pointer :: Star_Formation_Timescale_Disk_FGet => null()
  procedure(Star_Formation_Timescale_Disk_CTemplate), pointer :: Star_Formation_Timescale_Disk_CGet => null()
  abstract interface
     double precision function Star_Formation_Timescale_Disk_FTemplate(thisNode)
       import treeNode
       type(treeNode), intent(inout), pointer :: thisNode
     end function Star_Formation_Timescale_Disk_FTemplate
  end interface
  abstract interface
     function Star_Formation_Timescale_Disk_CTemplate(cNode   )
       import
       real(c_double)       :: Star_Formation_Timescale_Disk_CTemplate
       type(c_ptr),   value :: cNode
     end function Star_Formation_Timescale_Disk_CTemplate
  end interface

contains

  subroutine Star_Formation_Timescale_Disks_Initialize
    !% Initialize the disk star formation timecale module.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="starFormationTimescaleDisksMethod" type="moduleUse">
    include 'star_formation.timescales.disks.modules.inc'
    !# </include>
    implicit none
    type(c_funptr) :: cFunctionPointer=c_null_funptr
     !# <include directive="starFormationTimescaleDisksMethod" type="cBinding">
    include 'star_formation.timescales.disks.cBinding.inc'
    !# </include>

    ! Initialize if necessary.
    if (.not.starFormationTimescaleDisksInitialized) then
       !$omp critical(Star_Formation_Timescale_Disks_Initialization) 
       if (.not.starFormationTimescaleDisksInitialized) then
          ! Get the disk star formation timescale method parameter.
          !@ <inputParameter>
          !@   <name>starFormationTimescaleDisksMethod</name>
          !@   <defaultValue>KMT09</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be used for computing star formation timescales in disks.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@   <group>starFormation</group>
          !@ </inputParameter>
          call Get_Input_Parameter('starFormationTimescaleDisksMethod',starFormationTimescaleDisksMethod,defaultValue='KMT09')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="starFormationTimescaleDisksMethod" type="code" action="subroutine">
          !#  <subroutineArgs>
          !#   <fortran>starFormationTimescaleDisksMethod,Star_Formation_Timescale_Disk_FGet</fortran>
          !#   <c>char(starFormationTimescaleDisksMethod)//c_null_char,cFunctionPointer</c>
          !#  </subroutineArgs>
          include 'star_formation.timescales.disks.inc'
          !# </include>
          functionIsFortran=associated(Star_Formation_Timescale_Disk_FGet)
          if (.not.functionIsFortran) then
             ! Check if a C implementation was selected.
             if (c_associated(cFunctionPointer)) then
                ! One was, so transfer to the Fortran procedure pointer.
                call c_f_procpointer(cFunctionPointer,Star_Formation_Timescale_Disk_CGet)
             else
                call Galacticus_Error_Report('Star_Formation_Timescale_Disks'&
                     &,'method '//char(starFormationTimescaleDisksMethod)//' is unrecognized')
             end if
          end if
          starFormationTimescaleDisksInitialized=.true.
       end if
       !$omp end critical(Star_Formation_Timescale_Disks_Initialization) 
    end if
    return
  end subroutine Star_Formation_Timescale_Disks_Initialize

  double precision function Star_Formation_Timescale_Disk(thisNode)
    !% Returns the timescale (in Gyr) for star formation in the disk component of {\tt thisNode}.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    type(c_ptr)                            :: cNode

    ! Initialize the module.
    call Star_Formation_Timescale_Disks_Initialize

    ! Get the energy using the selected method.
    if (functionIsFortran) then
       Star_Formation_Timescale_Disk=Star_Formation_Timescale_Disk_FGet(thisNode)
    else
       cNode=c_loc(thisNode)
       Star_Formation_Timescale_Disk=Star_Formation_Timescale_Disk_CGet(cNode   )
    end if
    return
  end function Star_Formation_Timescale_Disk
  
end module Star_Formation_Timescales_Disks
