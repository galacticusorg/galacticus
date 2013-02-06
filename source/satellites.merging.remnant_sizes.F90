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

!% Contains a module which implements calculations of merger remnant sizes.

module Satellite_Merging_Remnant_Sizes
  !% Implements calculations of merger remnant sizes.
  use ISO_Varying_String
  implicit none
  private
  public :: Satellite_Merging_Remnant_Size

  ! Flag to indicate if this module has been initialized.  
  logical                                            :: satelliteMergingRemnantSizeInitialized=.false.

  ! Name of mass movement method used.
  type(varying_string)                               :: satelliteMergingRemnantSizeMethod

  ! Pointer to the subroutine that returns descriptors for mass movement.
  procedure(Satellite_Merging_Remnant_Size), pointer :: Satellite_Merging_Remnant_Size_Do => null()
  
contains

  !# <satelliteMergerTask>
  !#  <unitName>Satellite_Merging_Remnant_Size</unitName>
  !#  <after>Satellite_Merging_Mass_Movement_Store</after>
  !# </satelliteMergerTask>
  subroutine Satellite_Merging_Remnant_Size(thisNode)
    !% Computes the size of a merger remnant.
    use Galacticus_Nodes
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="satelliteMergingRemnantSizeMethod" type="moduleUse">
    include 'satellites.merging.remnant_sizes.modules.inc'
    !# </include>
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode
    
    if (.not.satelliteMergingRemnantSizeInitialized) then
       !$omp critical(satelliteMergingRemnantSizeInitialize)
       if (.not.satelliteMergingRemnantSizeInitialized) then
          ! Do the satellite merging remnant sizes method parameter.
          !@ <inputParameter>
          !@   <name>satelliteMergingRemnantSizeMethod</name>
          !@   <defaultValue>Covington2008</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be used for computing merger remnant sizes.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('satelliteMergingRemnantSizeMethod',satelliteMergingRemnantSizeMethod,defaultValue='Covington2008')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="satelliteMergingRemnantSizeMethod" type="functionCall" functionType="void">
          !#  <functionArgs>satelliteMergingRemnantSizeMethod,Satellite_Merging_Remnant_Size_Do</functionArgs>
          include 'satellites.merging.remnant_sizes.inc'
          !# </include>
          if (.not.associated(Satellite_Merging_Remnant_Size_Do)) call Galacticus_Error_Report('Satellite_Merging_Remnant_Size','method ' &
               &//char(satelliteMergingRemnantSizeMethod)//' is unrecognized')
          ! Flag that the module is now initialized.
          satelliteMergingRemnantSizeInitialized=.true.
       end if
       !$omp end critical(satelliteMergingRemnantSizeInitialize)
    end if

    ! Call the routine to do the calculation.
    call Satellite_Merging_Remnant_Size_Do(thisNode)

    return
  end subroutine Satellite_Merging_Remnant_Size
  
end module Satellite_Merging_Remnant_Sizes
