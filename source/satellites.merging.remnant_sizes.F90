!% Contains a module which implements calculations of merger remnant sizes.

module Satellite_Merging_Remnant_Sizes
  !% Implements calculations of merger remnant sizes.
  use ISO_Varying_String
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
    use Tree_Nodes
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="satelliteMergingRemnantSizeMethod" type="moduleUse">
    include 'satellites.merging.remnant_sizes.modules.inc'
    !# </include>
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode
    
    !$omp critical(satelliteMergingRemnantSizeInitialize)
    if (.not.satelliteMergingRemnantSizeInitialized) then
       ! Do the satellite merging remnant sizes method parameter.
       !@ <inputParameter>
       !@   <name>satelliteMergingRemnantSizeMethod</name>
       !@   <defaultValue>Cole2000</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of the method to be used for computing merger remnant sizes.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('satelliteMergingRemnantSizeMethod',satelliteMergingRemnantSizeMethod,defaultValue='Cole2000')
       ! Include file that makes calls to all available method initialization routines.
       !# <include directive="satelliteMergingRemnantSizeMethod" type="code" action="subroutine">
       !#  <subroutineArgs>satelliteMergingRemnantSizeMethod,Satellite_Merging_Remnant_Size_Do</subroutineArgs>
       include 'satellites.merging.remnant_sizes.inc'
       !# </include>
       if (.not.associated(Satellite_Merging_Remnant_Size_Do)) call Galacticus_Error_Report('Satellite_Merging_Remnant_Size','method ' &
            &//char(satelliteMergingRemnantSizeMethod)//' is unrecognized')
       ! Flag that the module is now initialized.
       satelliteMergingRemnantSizeInitialized=.true.
    end if
    !$omp end critical(satelliteMergingRemnantSizeInitialize)

    ! Call the routine to do the calculation.
    call Satellite_Merging_Remnant_Size_Do(thisNode)

    return
  end subroutine Satellite_Merging_Remnant_Size
  
end module Satellite_Merging_Remnant_Sizes
