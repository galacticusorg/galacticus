!% A null dark matter profile method.

module Tree_Node_Methods_Dark_Matter_Profile_Null
  !% A null dark matter profile method.
  private
  public :: Tree_Node_Methods_Dark_Matter_Profile_Null_Initialize
  
  ! Define procedure pointers.
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Dark_Matter_Profile_Scale</methodName>
  !# </treeNodeMethodsPointer>
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Dark_Matter_Profile_Scale_Growth_Rate</methodName>
  !# </treeNodeMethodsPointer>

contains

  !! Initialization routine.

  !# <treeNodeCreateInitialize>
  !#  <unitName>Tree_Node_Methods_Dark_Matter_Profile_Null_Initialize</unitName>
  !#  <optionName default="null">treeNodeMethodDarkMatterProfile</optionName>
  !# </treeNodeCreateInitialize>
  subroutine Tree_Node_Methods_Dark_Matter_Profile_Null_Initialize(componentOption,componentTypeCount)
    !% Initializes the tree node dark matter profile null methods module.
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
       message='Null dark matter profile method selected'
       call Galacticus_Display_Message(message,verbosityInfo)

       ! Set up procedure pointers.
       Tree_Node_Dark_Matter_Profile_Scale                          => null()
       Tree_Node_Dark_Matter_Profile_Scale_Set                      => null()
       Tree_Node_Dark_Matter_Profile_Scale_Rate_Adjust              => Tree_Node_Rate_Adjust_Dummy
       Tree_Node_Dark_Matter_Profile_Scale_Rate_Compute             => Tree_Node_Rate_Rate_Compute_Dummy
       Tree_Node_Dark_Matter_Profile_Scale_Growth_Rate              => null()
       Tree_Node_Dark_Matter_Profile_Scale_Growth_Rate_Set          => null()
       Tree_Node_Dark_Matter_Profile_Scale_Growth_Rate_Rate_Adjust  => Tree_Node_Rate_Adjust_Dummy
       Tree_Node_Dark_Matter_Profile_Scale_Growth_Rate_Rate_Compute => Tree_Node_Rate_Rate_Compute_Dummy

    end if
    return
  end subroutine Tree_Node_Methods_Dark_Matter_Profile_Null_Initialize
  
end module Tree_Node_Methods_Dark_Matter_Profile_Null
