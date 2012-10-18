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

!% Contains a module with the standard implementation of basic tree node methods.

module Tree_Node_Methods_Basic_Non_Evolving
  !% The standard implementation of basic tree node methods.
  use Tree_Nodes
  use Components
  implicit none
  private
  public :: Tree_Node_Methods_Basic_Initialize_Non_Evolving, Tree_Node_Basic_Promote_Non_Evolving,&
       & Galacticus_Output_Tree_Basic_Names_Non_Evolving, Galacticus_Output_Tree_Basic_Property_Count_Non_Evolving, Galacticus_Output_Tree_Basic_Non_Evolving,&
       & Tree_Node_Methods_Basic_Dump_Non_Evolving, Basic_Non_Evolving_Scale_Set
  
  ! The index used as a reference for this component.
  integer :: componentIndex=-1

  ! Property indices.
  integer, parameter :: propertyCount=1, dataCount=2, historyCount=0
  integer, parameter :: timeIndex     =1
  integer, parameter :: isolationIndex=1
  integer, parameter :: massIndex     =2

  ! Define procedure pointers.
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Mass</methodName>
  !# </treeNodeMethodsPointer>
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Time</methodName>
  !# </treeNodeMethodsPointer>
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Time_Last_Isolated</methodName>
  !# </treeNodeMethodsPointer>

  ! Flag to indicate if this method is selected.
  logical :: methodSelected=.false.

contains


  !# <treeNodeCreateInitialize>
  !#  <unitName>Tree_Node_Methods_Basic_Initialize_Non_Evolving</unitName>
  !#  <optionName>treeNodeMethodBasic</optionName>
  !# </treeNodeCreateInitialize>
  subroutine Tree_Node_Methods_Basic_Initialize_Non_Evolving(componentOption,componentTypeCount)
    !% Initializes the tree node basic methods module.
    use ISO_Varying_String
    use String_Handling
    use Galacticus_Display
    implicit none
    type(varying_string), intent(in)    :: componentOption
    integer,              intent(inout) :: componentTypeCount
    type(varying_string)                :: message

    ! Check if this implementation is selected.
    if (componentOption == 'nonEvolving') then
       ! Record that method is selected.
       methodSelected=.true.
       
       ! Increment the component count and store the value for later reference.
       componentTypeCount=componentTypeCount+1
       componentIndex=componentTypeCount

       ! Display message.
       message='Non-evolving basic method selected [component index '
       message=message//componentIndex//']'
       call Galacticus_Display_Message(message,verbosityInfo)

       ! Set up procedure pointers.
       Tree_Node_Mass              => Tree_Node_Mass_Basic
       Tree_Node_Mass_Set          => Tree_Node_Mass_Set_Basic
       Tree_Node_Mass_Rate_Adjust  => null()
       Tree_Node_Mass_Rate_Compute => Tree_Node_Rate_Rate_Compute_Dummy

       Tree_Node_Time              => Tree_Node_Time_Basic
       Tree_Node_Time_Set          => Tree_Node_Time_Set_Basic
       Tree_Node_Time_Rate_Adjust  => Tree_Node_Time_Rate_Adjust_Basic
       Tree_Node_Time_Rate_Compute => Tree_Node_Time_Rate_Compute_Basic

       Tree_Node_Time_Last_Isolated              => Tree_Node_Time_Last_Isolated_Basic
       Tree_Node_Time_Last_Isolated_Set          => null()
       Tree_Node_Time_Last_Isolated_Rate_Adjust  => null()
       Tree_Node_Time_Last_Isolated_Rate_Compute => Tree_Node_Rate_Rate_Compute_Dummy
      end if
    return
  end subroutine Tree_Node_Methods_Basic_Initialize_Non_Evolving

  double precision function Tree_Node_Mass_Basic(thisNode,instance)
    !% Return the node mass.
    implicit none
    integer, intent(in), optional :: instance
    type(treeNode), pointer, intent(inout) :: thisNode
    integer                                :: thisIndex

    thisIndex=Tree_Node_Basic_Index(thisNode)
    Tree_Node_Mass_Basic=thisNode%components(thisIndex)%instance(1)%data(massIndex)
    return
  end function Tree_Node_Mass_Basic

  subroutine Tree_Node_Mass_Set_Basic(thisNode,mass,instance)
    !% Set the node mass.
    implicit none
    integer, intent(in), optional :: instance
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: mass
    integer                                  :: thisIndex

    thisIndex=Tree_Node_Basic_Index(thisNode)
    thisNode%components(thisIndex)%instance(1)%data(massIndex)=mass
    return
  end subroutine Tree_Node_Mass_Set_Basic

  double precision function Tree_Node_Time_Basic(thisNode,instance)
    !% Return the node time.
    implicit none
    integer, intent(in), optional :: instance
    type(treeNode), pointer, intent(inout) :: thisNode
    integer                                :: thisIndex

    thisIndex=Tree_Node_Basic_Index(thisNode)
    Tree_Node_Time_Basic=thisNode%components(thisIndex)%instance(1)%properties(timeIndex,propertyValue)
    return
  end function Tree_Node_Time_Basic

  subroutine Tree_Node_Time_Set_Basic(thisNode,time,instance)
    !% Set the node time.
    implicit none
    integer, intent(in), optional :: instance
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: time
    integer                                  :: thisIndex

    thisIndex=Tree_Node_Basic_Index(thisNode)
    thisNode%components(thisIndex)%instance(1)%properties(timeIndex,propertyValue)=time
    return
  end subroutine Tree_Node_Time_Set_Basic

  subroutine Tree_Node_Time_Rate_Adjust_Basic(thisNode,interrupt,interruptProcedure,rateAdjustment,instance)
    !% Return the node time rate of change.
    implicit none
    integer, intent(in), optional :: instance
    type(treeNode),   pointer, intent(inout) :: thisNode
    logical,                   intent(inout) :: interrupt
    procedure(), pointer, intent(inout) :: interruptProcedure
    double precision,          intent(in)    :: rateAdjustment
    integer                                  :: thisIndex
    
    thisIndex=Tree_Node_Basic_Index(thisNode)
    thisNode%components(thisIndex)%instance(1)%properties(timeIndex,propertyDerivative)=thisNode%components(thisIndex)%instance(1)%properties(timeIndex&
         &,propertyDerivative)+rateAdjustment
    return
  end subroutine Tree_Node_Time_Rate_Adjust_Basic

  subroutine Tree_Node_Time_Rate_Compute_Basic(thisNode,interrupt,interruptProcedure)
    !% Compute the node time rate of change.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    logical,                 intent(inout) :: interrupt
    procedure(), pointer, intent(inout) :: interruptProcedure
    
    call Tree_Node_Time_Rate_Adjust_Basic(thisNode,interrupt,interruptProcedure,1.0d0) ! By definition.
    return
  end subroutine Tree_Node_Time_Rate_Compute_Basic

  double precision function Tree_Node_Time_Last_Isolated_Basic(thisNode,instance)
    !% Returns the time at which {\tt thisNode} was last an isolated halo.
    implicit none
    integer, intent(in), optional :: instance
    type(treeNode), pointer, intent(inout) :: thisNode
    integer                                :: thisIndex

    if (thisNode%isSatellite()) then
       thisIndex=Tree_Node_Basic_Index(thisNode)
       Tree_Node_Time_Last_Isolated_Basic=thisNode%components(thisIndex)%instance(1)%data(isolationIndex)
    else
       Tree_Node_Time_Last_Isolated_Basic=Tree_Node_Time_Basic(thisNode)
    end if
    return
  end function Tree_Node_Time_Last_Isolated_Basic

  !# <scaleSetTask>
  !#  <unitName>Basic_Non_Evolving_Scale_Set</unitName>
  !# </scaleSetTask>
  subroutine Basic_Non_Evolving_Scale_Set(thisNode)
    !% Set scales for properties of {\tt thisNode}.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision, parameter              :: timeScale=1.0d-3
    integer                                  :: thisIndex
 
    ! Determine if method is active and a basic component exists.
    if (methodSelected.and.thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Basic_Index(thisNode)

       ! Set scale for time.
       thisNode%components(thisIndex)%instance(1)%properties(timeIndex,propertyScale)=timeScale

    end if
    return
  end subroutine Basic_Non_Evolving_Scale_Set

  !# <nodePromotionTask>
  !#  <unitName>Tree_Node_Basic_Promote_Non_Evolving</unitName>
  !# </nodePromotionTask>
  subroutine Tree_Node_Basic_Promote_Non_Evolving(thisNode)
    !% Ensure that {\tt thisNode} is ready for promotion to its parent. In this case, we simply update the mass of {\tt thisNode}
    !% to be that of its parent.
    use Galacticus_Error
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    type(treeNode), pointer                :: parentNode

    if (methodSelected) then
       parentNode => thisNode%parentNode
       if (Tree_Node_Time(thisNode) /= Tree_Node_Time(parentNode)) call Galacticus_Error_Report('Tree_Node_Basic_Promote','thisNode&
            & has not been evolved to its parent')
       ! Adjust the mass to that of the parent node.
       call Tree_Node_Mass_Set_Basic(thisNode,Tree_Node_Mass_Basic(parentNode))
    end if
    return
  end subroutine Tree_Node_Basic_Promote_Non_Evolving

  integer function Tree_Node_Basic_Index(thisNode)
    !% Ensure the basic component exists and return its position in the components array.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode

    call thisNode%createComponent(componentIndex,propertyCount,dataCount,historyCount)
    Tree_Node_Basic_Index=thisNode%componentIndex(componentIndex)
    return
  end function Tree_Node_Basic_Index

  !# <mergerTreeOutputNames>
  !#  <unitName>Galacticus_Output_Tree_Basic_Names_Non_Evolving</unitName>
  !#  <sortName>Galacticus_Output_Tree_Basic_Non_Evolving</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Tree_Basic_Names_Non_Evolving(integerProperty,integerPropertyNames,integerPropertyComments&
       &,integerPropertyUnitsSI,doubleProperty,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time)
    !% Set names of basic properties to be written to the \glc\ output file.
    use Numerical_Constants_Astronomical
    implicit none
    double precision, intent(in)                  :: time
    integer,          intent(inout)               :: integerProperty,doubleProperty
    character(len=*), intent(inout), dimension(:) :: integerPropertyNames,integerPropertyComments,doublePropertyNames &
         &,doublePropertyComments
    double precision, intent(inout), dimension(:) :: integerPropertyUnitsSI,doublePropertyUnitsSI

    if (methodSelected) then
       doubleProperty=doubleProperty+1
       !@ <outputProperty>
       !@   <name>nodeMass</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>Total mass of the node, assuming univeral baryon fraction.</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@ </outputProperty>
       doublePropertyNames   (doubleProperty)='nodeMass'
       doublePropertyComments(doubleProperty)='Total mass of the node, assuming univeral baryon fraction.'
       doublePropertyUnitsSI (doubleProperty)=massSolar
       doubleProperty=doubleProperty+1
       !@ <outputProperty>
       !@   <name>nodeTimeLastIsolated</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>Time at which node was last an isolated halo.</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@ </outputProperty>
       doublePropertyNames   (doubleProperty)='nodeTimeLastIsolated'
       doublePropertyComments(doubleProperty)='Time at which node was last an isolated halo.'
       doublePropertyUnitsSI (doubleProperty)=gigaYear
    end if
    return
  end subroutine Galacticus_Output_Tree_Basic_Names_Non_Evolving

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Galacticus_Output_Tree_Basic_Property_Count_Non_Evolving</unitName>
  !#  <sortName>Galacticus_Output_Tree_Basic_Non_Evolving</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Galacticus_Output_Tree_Basic_Property_Count_Non_Evolving(integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of basic properties to be written to the the \glc\ output file.
    implicit none
    double precision, intent(in)    :: time
    integer,          intent(inout) :: integerPropertyCount,doublePropertyCount

    if (methodSelected) doublePropertyCount=doublePropertyCount+2
    return
  end subroutine Galacticus_Output_Tree_Basic_Property_Count_Non_Evolving

  !# <mergerTreeOutputTask>
  !#  <unitName>Galacticus_Output_Tree_Basic_Non_Evolving</unitName>
  !#  <sortName>Galacticus_Output_Tree_Basic_Non_Evolving</sortName>
  !# </mergerTreeOutputTask>
  subroutine Galacticus_Output_Tree_Basic_Non_Evolving(thisNode,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time)
    !% Store basic properties in the \glc\ output file buffers.
    use Tree_Nodes
    use Kind_Numbers
    implicit none
    double precision,        intent(in)             :: time
    type(treeNode),          intent(inout), pointer :: thisNode
    integer,                 intent(inout)          :: integerProperty,integerBufferCount,doubleProperty,doubleBufferCount
    integer(kind=kind_int8), intent(inout)          :: integerBuffer(:,:)
    double precision,        intent(inout)          :: doubleBuffer(:,:)

    if (methodSelected) then
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=Tree_Node_Mass(thisNode)
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=Tree_Node_Time_Last_Isolated_Basic(thisNode)
    end if
    return
  end subroutine Galacticus_Output_Tree_Basic_Non_Evolving

  !# <nodeDumpTask>
  !#  <unitName>Tree_Node_Methods_Basic_Dump_Non_Evolving</unitName>
  !# </nodeDumpTask>
  subroutine Tree_Node_Methods_Basic_Dump_Non_Evolving(thisNode)
    !% Dump all properties of {\tt thisNode} to screen.
    implicit none
    type(treeNode),   intent(inout), pointer     :: thisNode

    if (methodSelected) then
       if (thisNode%componentExists(componentIndex)) then
          write (0,'(1x,a)'           ) 'non-evolving node component -> properties:'
          write (0,'(2x,a50,1x,e12.6)') 'node mass:'              ,Tree_Node_Mass(thisNode)
          write (0,'(2x,a50,1x,e12.6)') 'node time:'              ,Tree_Node_Time(thisNode)
          write (0,'(2x,a50,1x,e12.6)') 'node time last isolated:',Tree_Node_Time_Last_Isolated_Basic(thisNode)
       else
          write (0,'(1x,a)'           ) 'basic node component -> nonexistant'
       end if
    end if
    return
  end subroutine Tree_Node_Methods_Basic_Dump_Non_Evolving

end module Tree_Node_Methods_Basic_Non_Evolving
