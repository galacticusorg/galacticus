!% A null black hole method.

module Tree_Node_Methods_Black_Hole_Null
  !% A null black hole method.
  private
  public :: Tree_Node_Methods_Black_Hole_Null_Initialize
  
  ! Define procedure pointers.
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Black_Hole_Mass</methodName>
  !# </treeNodeMethodsPointer>
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Black_Hole_Spin</methodName>
  !# </treeNodeMethodsPointer>

contains

  !! Initialization routine.

  !# <treeNodeCreateInitialize>
  !#  <unitName>Tree_Node_Methods_Black_Hole_Null_Initialize</unitName>
  !#  <optionName>treeNodeMethodBlackHole</optionName>
  !# </treeNodeCreateInitialize>
  subroutine Tree_Node_Methods_Black_Hole_Null_Initialize(componentOption,componentTypeCount)
    !% Initializes the tree node black hole null methods module.
    use Tree_Node_Methods
    use ISO_Varying_String
    use Galacticus_Display
    implicit none
    type(varying_string), intent(in)    :: componentOption
    integer,              intent(inout) :: componentTypeCount
    type(varying_string)                :: message

    ! Check if this implementation is selected.
    if (componentOption == 'null') then

       ! Display message.
       message='Null black hole method selected'
       call Galacticus_Display_Message(message,verbosityInfo)

       ! Set up procedure pointers.
       Tree_Node_Black_Hole_Mass              => null()
       Tree_Node_Black_Hole_Mass_Set          => null()
       Tree_Node_Black_Hole_Mass_Rate_Adjust  => Tree_Node_Rate_Adjust_Dummy
       Tree_Node_Black_Hole_Mass_Rate_Compute => Tree_Node_Rate_Rate_Compute_Dummy
       Tree_Node_Black_Hole_Spin              => null()
       Tree_Node_Black_Hole_Spin_Set          => null()
       Tree_Node_Black_Hole_Spin_Rate_Adjust  => Tree_Node_Rate_Adjust_Dummy
       Tree_Node_Black_Hole_Spin_Rate_Compute => Tree_Node_Rate_Rate_Compute_Dummy

    end if
    return
  end subroutine Tree_Node_Methods_Black_Hole_Null_Initialize
  
end module Tree_Node_Methods_Black_Hole_Null
