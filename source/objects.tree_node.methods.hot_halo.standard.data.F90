!% Contains a module that stores data that is shared between modules for the standard hot halo method.

module Tree_Node_Methods_Hot_Halo_Data
  !% Stores data that is shared between modules for the standard hot halo method.
  public

  ! The index used as a reference for this component.
  integer :: componentIndex=-1

  ! Flag to indicate if this method is selected.
  logical          :: methodSelected=.false.

end module Tree_Node_Methods_Hot_Halo_Data
