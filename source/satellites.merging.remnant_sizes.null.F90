!% Contains a module which implements a null algorithm for merger remnant sizes.

module Satellite_Merging_Remnant_Sizes_Null
  !% Implements a null algorithm for merger remnant sizes.
  private
  public :: Satellite_Merging_Remnant_Sizes_Null_Initialize

contains

  !# <satelliteMergingRemnantSizeMethod>
  !#  <unitName>Satellite_Merging_Remnant_Sizes_Null_Initialize</unitName>
  !# </satelliteMergingRemnantSizeMethod>
  subroutine Satellite_Merging_Remnant_Sizes_Null_Initialize(satelliteMergingRemnantSizeMethod,Satellite_Merging_Remnant_Size_Do)
    !% Test if this method is to be used and set procedure pointer appropriately.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string),          intent(in)    :: satelliteMergingRemnantSizeMethod
    procedure(),          pointer, intent(inout) :: Satellite_Merging_Remnant_Size_Do
    
    if (satelliteMergingRemnantSizeMethod == 'null') Satellite_Merging_Remnant_Size_Do => Satellite_Merging_Remnant_Size_Null
    return
  end subroutine Satellite_Merging_Remnant_Sizes_Null_Initialize

  subroutine Satellite_Merging_Remnant_Size_Null(thisNode)
    !% A null implementation of merger remnant size. Does nothing.
    use Tree_Nodes
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    ! Do nothing.

    return
  end subroutine Satellite_Merging_Remnant_Size_Null

end module Satellite_Merging_Remnant_Sizes_Null
