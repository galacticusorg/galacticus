!% Contains a module which implements a dark matter profile method that provides a scale radius.

module Tree_Node_Methods_Dark_Matter_Profile_Scales
  !% Implements a dark matter profile method that provides a scale radius.
  use Tree_Nodes
  use Tree_Node_Methods
  use Components
  private
  public :: Tree_Node_Methods_Profile_Scale_Initialize, Tree_Node_Methods_Profile_Scale_Initialize_Scale,&
       & Galacticus_Output_Tree_Profile_Scale, Galacticus_Output_Tree_Profile_Scale_Property_Count,&
       & Galacticus_Output_Tree_Profile_Scale_Names, Tree_Node_Methods_Profile_Scale_Dump, Tree_Node_Dark_Matter_Profile_Scale_Promote
  
  ! The index used as a reference for this component.
  integer            :: componentIndex=-1

  ! Property indices.
  integer, parameter :: propertyCount =1, dataCount=1, historyCount=0
  integer, parameter :: scaleIndex    =1
  integer, parameter :: scaleRateIndex=1

  ! Define procedure pointers.
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Dark_Matter_Profile_Scale</methodName>
  !# </treeNodeMethodsPointer>
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Dark_Matter_Profile_Scale_Growth_Rate</methodName>
  !# </treeNodeMethodsPointer>

  ! Flag to indicate if this method is selected.
  logical          :: methodSelected=.false.

  ! Parameters of the method.
  double precision :: darkMatterProfileMinimumConcentration

contains

  !! Initialization routine.

  !# <treeNodeCreateInitialize>
  !#  <unitName>Tree_Node_Methods_Profile_Scale_Initialize</unitName>
  !#  <optionName>treeNodeMethodDarkMatterProfile</optionName>
  !# </treeNodeCreateInitialize>
  subroutine Tree_Node_Methods_Profile_Scale_Initialize(componentOption,componentTypeCount)
    !% Initializes the tree node ``scale'' dark matter profile method module.
    use ISO_Varying_String
    use Input_Parameters
    use String_Handling
    use Galacticus_Display
    implicit none
    type(varying_string), intent(in)    :: componentOption
    integer,              intent(inout) :: componentTypeCount
    type(varying_string)                :: message

    ! Check if this implementation is selected.
    if (componentOption.eq.'scale') then
       ! Record that method is selected.
       methodSelected=.true.

       ! Increment the component count and store the value for later reference.
       componentTypeCount=componentTypeCount+1
       componentIndex    =componentTypeCount

       ! Display message.
       message='Scale dark matter profile method selected [component index '
       message=message//componentIndex//']'
       call Galacticus_Display_Message(message,verbosityInfo)

       ! Get parameters.
       !@ <inputParameter>
       !@   <name>darkMatterProfileMinimumConcentration</name>
       !@   <defaultValue>4</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The minimum concentration allowed for dark matter profiles.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('darkMatterProfileMinimumConcentration',darkMatterProfileMinimumConcentration,defaultValue=4.0d0)

       ! Set up procedure pointers.
       Tree_Node_Dark_Matter_Profile_Scale                          => Tree_Node_Dark_Matter_Profile_Scale_Scale
       Tree_Node_Dark_Matter_Profile_Scale_Set                      => null()
       Tree_Node_Dark_Matter_Profile_Scale_Rate_Adjust              => Tree_Node_Dark_Matter_Profile_Scale_Rate_Adjust_Scale
       Tree_Node_Dark_Matter_Profile_Scale_Rate_Compute             => Tree_Node_Dark_Matter_Profile_Scale_Rate_Compute_Scale
       Tree_Node_Dark_Matter_Profile_Scale_Growth_Rate              => Tree_Node_Dark_Matter_Profile_Scale_Growth_Rate_Scale
       Tree_Node_Dark_Matter_Profile_Scale_Growth_Rate_Set          => null()
       Tree_Node_Dark_Matter_Profile_Scale_Growth_Rate_Rate_Adjust  => null()
       Tree_Node_Dark_Matter_Profile_Scale_Growth_Rate_Rate_Compute => Tree_Node_Rate_Rate_Compute_Dummy
    end if
    return
  end subroutine Tree_Node_Methods_Profile_Scale_Initialize
  
  !! Methods.
  
  double precision function Tree_Node_Dark_Matter_Profile_Scale_Scale(thisNode)
    !% Return the node dark matter profile scale length.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    integer                                :: thisIndex
    
    ! Ensure the spin has been initialized.
    call Tree_Node_Methods_Profile_Scale_Initialize_Scale(thisNode)
    thisIndex=Tree_Node_Profile_Scale_Index(thisNode)
    Tree_Node_Dark_Matter_Profile_Scale_Scale=thisNode%components(thisIndex)%properties(scaleIndex,propertyValue)
    return
  end function Tree_Node_Dark_Matter_Profile_Scale_Scale

  subroutine Tree_Node_Dark_Matter_Profile_Scale_Rate_Adjust_Scale(thisNode,interrupt,interruptProcedure,rateAdjustment)
    !% Adjust the rate of scale radius growth.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    logical,                   intent(inout) :: interrupt
    procedure(),      pointer, intent(inout) :: interruptProcedure
    double precision,          intent(in)    :: rateAdjustment
    integer                                  :: thisIndex

    ! Apply the change in rate.
    thisIndex=Tree_Node_Profile_Scale_Index(thisNode)
    thisNode%components(thisIndex)%properties(scaleIndex,propertyDerivative)=thisNode%components(thisIndex)%properties(scaleIndex&
         &,propertyDerivative)+rateAdjustment
    return
  end subroutine Tree_Node_Dark_Matter_Profile_Scale_Rate_Adjust_Scale

  subroutine Tree_Node_Dark_Matter_Profile_Scale_Rate_Compute_Scale(thisNode,interrupt,interruptProcedure)
    !% Compute the rate of change of the scale radius.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    logical,                 intent(inout) :: interrupt
    procedure(),    pointer, intent(inout) :: interruptProcedure
    integer                                :: thisIndex
    double precision                       :: growthRate

    thisIndex=Tree_Node_Profile_Scale_Index(thisNode)
    growthRate=Tree_Node_Dark_Matter_Profile_Scale_Growth_Rate_Scale(thisNode)
    call Tree_Node_Dark_Matter_Profile_Scale_Rate_Adjust_Scale(thisNode,interrupt,interruptProcedure,growthRate)
    return
  end subroutine Tree_Node_Dark_Matter_Profile_Scale_Rate_Compute_Scale

  double precision function Tree_Node_Dark_Matter_Profile_Scale_Growth_Rate_Scale(thisNode)
    !% Return the node dark matter profile scale length rate of growth (assumed to be zero).
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    integer                                  :: thisIndex
    double precision                         :: deltaTime

    ! Get component index for this node.
    thisIndex=Tree_Node_Profile_Scale_Index(thisNode)
    ! Return the pre-computed growth rate.
    Tree_Node_Dark_Matter_Profile_Scale_Growth_Rate_Scale=thisNode%components(thisIndex)%data(scaleRateIndex)
    return
  end function Tree_Node_Dark_Matter_Profile_Scale_Growth_Rate_Scale

  !! Tasks.

  !# <mergerTreeInitializeTask>
  !#  <unitName>Tree_Node_Methods_Profile_Scale_Initialize_Scale</unitName>
  !# </mergerTreeInitializeTask>
  subroutine Tree_Node_Methods_Profile_Scale_Initialize_Scale(thisNode)
    !% Initialize the scale radius of {\tt thisNode}.
    use Dark_Matter_Profiles_Concentrations
    use Dark_Matter_Halo_Scales
    implicit none
    type(treeNode),  pointer, intent(inout) :: thisNode
    type(treeNode),  pointer                :: relatedNode
    integer                                 :: thisIndex,relatedIndex
    double precision                        :: concentration,deltaTime

    if (methodSelected) then
       if (.not.thisNode%componentExists(componentIndex)) then
          thisIndex=Tree_Node_Profile_Scale_Index(thisNode)
          ! Set the scale radius of the halo.
          concentration=max(Dark_Matter_Profile_Concentration(thisNode),darkMatterProfileMinimumConcentration)
          thisNode%components(thisIndex)%properties(scaleIndex,propertyValue)=Dark_Matter_Halo_Virial_Radius(thisNode)/concentration
          ! Check if this node is the primary progenitor.
          if (thisNode%isPrimaryProgenitor()) then
             ! It is, so compute the scale radius growth rate.
             ! First ensure that parent node has scale radius set.
             call Tree_Node_Methods_Profile_Scale_Initialize_Scale(thisNode%parentNode)
             ! Now compute the growth rate.
             deltaTime=Tree_Node_Time(thisNode%parentNode)-Tree_Node_Time(thisNode)
             if (deltaTime > 0.0d0) then
                thisNode%components(thisIndex)%data(scaleRateIndex)=(Tree_Node_Dark_Matter_Profile_Scale_Scale(thisNode%parentNode)&
                     &-Tree_Node_Dark_Matter_Profile_Scale_Scale(thisNode))/deltaTime
             else
                thisNode%components(thisIndex)%data(scaleRateIndex)=0.0d0
             end if
          else
             ! It is not, so set scale radius growth rate to zero.
             thisNode%components(thisIndex)%data(scaleRateIndex)=0.0d0
          end if
       end if
    end if
    return
  end subroutine Tree_Node_Methods_Profile_Scale_Initialize_Scale
  
  !# <nodePromotionTask>
  !#  <unitName>Tree_Node_Dark_Matter_Profile_Scale_Promote</unitName>
  !# </nodePromotionTask>
  subroutine Tree_Node_Dark_Matter_Profile_Scale_Promote(thisNode)
    !% Ensure that {\tt thisNode} is ready for promotion to its parent. In this case, we simply update the growth rate of {\tt thisNode}
    !% to be that of its parent.
    use Galacticus_Error
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    type(treeNode), pointer                :: parentNode
    integer                                :: thisIndex,parentIndex

    if (methodSelected) then
       parentNode => thisNode%parentNode
       if (Tree_Node_Time(thisNode) /= Tree_Node_Time(parentNode)) call Galacticus_Error_Report('Tree_Node_Dark_Matter_Profile_Scale_Promote','thisNode&
            & has not been evolved to its parent')
       ! Get component indices.
       thisIndex  =Tree_Node_Profile_Scale_Index(thisNode)
       parentIndex=Tree_Node_Profile_Scale_Index(parentNode)
       ! Adjust the scale radius to that of the parent node.
       thisNode%components(thisIndex)%properties(scaleIndex,propertyValue)=parentNode%components(parentIndex)%properties(scaleIndex,propertyValue)
       ! Adjust the growth rate to that of the parent node.
       thisNode%components(thisIndex)%data      (scaleRateIndex          )=parentNode%components(parentIndex)%data      (scaleRateIndex          )
    end if
    return
  end subroutine Tree_Node_Dark_Matter_Profile_Scale_Promote
  
  !! Auxiliary functions.

  integer function Tree_Node_Profile_Scale_Index(thisNode)
    !% Ensure the profile component exists and return its position in the components array.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    
    call thisNode%createComponent(componentIndex,propertyCount,dataCount,historyCount)
    Tree_Node_Profile_Scale_Index=thisNode%componentIndex(componentIndex)
    return
  end function Tree_Node_Profile_Scale_Index

  !! Output functions.
  
  !# <mergerTreeOutputNames>
  !#  <unitName>Galacticus_Output_Tree_Profile_Scale_Names</unitName>
  !#  <sortName>Galacticus_Output_Tree_Profile_Scale</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Tree_Profile_Scale_Names(integerProperty,integerPropertyNames,integerPropertyComments,doubleProperty&
       &,doublePropertyNames,doublePropertyComments,time)
    !% Set names of scale properties to be written to the \glc\ output file.
    implicit none
    double precision, intent(in)                  :: time
    integer,          intent(inout)               :: integerProperty,doubleProperty
    character(len=*), intent(inout), dimension(:) :: integerPropertyNames,integerPropertyComments,doublePropertyNames &
         &,doublePropertyComments
    
    if (methodSelected) then
       doubleProperty=doubleProperty+1
       doublePropertyNames   (doubleProperty)='darkMatterScaleRadius'
       doublePropertyComments(doubleProperty)='Scale radius of the dark matter profile [Mpc].'
    end if
    return
  end subroutine Galacticus_Output_Tree_Profile_scale_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Galacticus_Output_Tree_Profile_Scale_Property_Count</unitName>
  !#  <sortName>Galacticus_Output_Tree_Profile_Scale</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Galacticus_Output_Tree_Profile_Scale_Property_Count(integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of spin properties to be written to the the \glc\ output file.
    implicit none
    double precision, intent(in)    :: time
    integer,          intent(inout) :: integerPropertyCount,doublePropertyCount

    if (methodSelected) doublePropertyCount=doublePropertyCount+1
    return
  end subroutine Galacticus_Output_Tree_Profile_Scale_Property_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Galacticus_Output_Tree_Profile_Scale</unitName>
  !#  <sortName>Galacticus_Output_Tree_Profile_Scale</sortName>
  !# </mergerTreeOutputTask>
  subroutine Galacticus_Output_Tree_Profile_Scale(thisNode,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time)
    !% Store spin properties in the \glc\ output file buffers.
    use Tree_Nodes
    use Tree_Node_Methods
    implicit none
    double precision, intent(in)             :: time
    type(treeNode),   intent(inout), pointer :: thisNode
    integer,          intent(inout)          :: integerProperty,integerBufferCount,doubleProperty,doubleBufferCount
    integer,          intent(inout)          :: integerBuffer(:,:)
    double precision, intent(inout)          :: doubleBuffer(:,:)

    if (methodSelected) then
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=Tree_Node_Dark_Matter_Profile_Scale(thisNode)
    end if
    return
  end subroutine Galacticus_Output_Tree_Profile_Scale

  !# <nodeDumpTask>
  !#  <unitName>Tree_Node_Methods_Profile_Scale_Dump</unitName>
  !# </nodeDumpTask>
  subroutine Tree_Node_Methods_Profile_Scale_Dump(thisNode)
    !% Dump all properties of {\tt thisNode} to screen.
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    if (methodSelected) then
       if (thisNode%componentExists(componentIndex)) then
          write (0,'(1x,a)'           ) 'dark matter profile properties:'
          write (0,'(2x,a50,1x,e12.6)') 'scale:',Tree_Node_Dark_Matter_Profile_Scale(thisNode)
       else
          write (0,'(1x,a)'           ) 'no dark matter profile component'
       end if
    end if
    return
  end subroutine Tree_Node_Methods_Profile_Scale_Dump

end module Tree_Node_Methods_Dark_Matter_Profile_Scales
