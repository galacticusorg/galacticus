!% A null spheroid method.

module Tree_Node_Methods_Spheroid_Null
  !% A null spheroid method.
  private
  public :: Tree_Node_Methods_Spheroid_Null_Initialize
  
contains

  !! Initialization routine.

  !# <treeNodeCreateInitialize>
  !#  <unitName>Tree_Node_Methods_Spheroid_Null_Initialize</unitName>
  !#  <optionName>treeNodeMethodSpheroid</optionName>
  !# </treeNodeCreateInitialize>
  subroutine Tree_Node_Methods_Spheroid_Null_Initialize(componentOption,componentTypeCount)
    !% Initializes the tree node spheroid null methods module.
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
       message='Null spheroid method selected'
       call Galacticus_Display_Message(message,verbosityInfo)

    end if
    return
  end subroutine Tree_Node_Methods_Spheroid_Null_Initialize
  
end module Tree_Node_Methods_Spheroid_Null
