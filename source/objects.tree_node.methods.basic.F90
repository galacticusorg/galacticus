!! Copyright 2009, 2010, Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


!% Contains a module of basic tree node methods.

module Tree_Node_Methods_Basic
  !% Implement basic tree node methods.
  use Tree_Nodes
  use Tree_Node_Methods
  use Components
  private
  public :: Tree_Node_Methods_Basic_Initialize, Tree_Node_Basic_Promote, Halo_Mass_Accretion_Rate, Tree_Node_Mass_Stop_Accretion,&
       & Galacticus_Output_Tree_Basic_Names, Galacticus_Output_Tree_Basic_Property_Count, Galacticus_Output_Tree_Basic,&
       & Tree_Node_Methods_Basic_Dump
  
  ! The index used as a reference for this component.
  integer :: componentIndex=-1

  ! Property indices.
  integer, parameter :: propertyCount=2, dataCount=2, historyCount=0
  integer, parameter :: massIndex     =1
  integer, parameter :: timeIndex     =2
  integer, parameter :: rateIndex     =1
  integer, parameter :: isolationIndex=2

  ! Define procedure pointers.
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Mass</methodName>
  !# </treeNodeMethodsPointer>
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Time</methodName>
  !# </treeNodeMethodsPointer>
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Mass_Accretion_Rate</methodName>
  !# </treeNodeMethodsPointer>
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Time_Last_Isolated</methodName>
  !# </treeNodeMethodsPointer>

  ! Flag to indicate if this method is selected.
  logical :: methodSelected=.false.

contains


  !# <treeNodeCreateInitialize>
  !#  <unitName>Tree_Node_Methods_Basic_Initialize</unitName>
  !#  <optionName default="standard">treeNodeMethodBasic</optionName>
  !# </treeNodeCreateInitialize>
  subroutine Tree_Node_Methods_Basic_Initialize(componentOption,componentTypeCount)
    !% Initializes the tree node basic methods module.
    use ISO_Varying_String
    implicit none
    type(varying_string), intent(in)    :: componentOption
    integer,              intent(inout) :: componentTypeCount

    ! Check if this implementation is selected.
    if (componentOption.eq.'standard') then
       ! Record that method is selected.
       methodSelected=.true.
       
       ! Increment the component count and store the value for later reference.
       componentTypeCount=componentTypeCount+1
       componentIndex=componentTypeCount

       ! Set up procedure pointers.
       Tree_Node_Mass              => Tree_Node_Mass_Basic
       Tree_Node_Mass_Set          => Tree_Node_Mass_Set_Basic
       Tree_Node_Mass_Rate_Adjust  => Tree_Node_Mass_Rate_Adjust_Basic
       Tree_Node_Mass_Rate_Compute => Tree_Node_Mass_Rate_Compute_Basic
 
       Tree_Node_Mass_Accretion_Rate              => Tree_Node_Mass_Accretion_Rate_Basic
       Tree_Node_Mass_Accretion_Rate_Set          => null()
       Tree_Node_Mass_Accretion_Rate_Rate_Adjust  => null()
       Tree_Node_Mass_Accretion_Rate_Rate_Compute => Tree_Node_Rate_Rate_Compute_Dummy

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
  end subroutine Tree_Node_Methods_Basic_Initialize


  double precision function Tree_Node_Mass_Basic(thisNode)
    !% Return the node mass.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    integer                                :: thisIndex

    thisIndex=Tree_Node_Basic_Index(thisNode)
    Tree_Node_Mass_Basic=thisNode%components(thisIndex)%properties(massIndex,propertyValue)
    return
  end function Tree_Node_Mass_Basic

  subroutine Tree_Node_Mass_Set_Basic(thisNode,mass)
    !% Set the node mass.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: mass
    integer                                  :: thisIndex

    thisIndex=Tree_Node_Basic_Index(thisNode)
    thisNode%components(thisIndex)%properties(massIndex,propertyValue)=mass
    return
  end subroutine Tree_Node_Mass_Set_Basic

  subroutine Tree_Node_Mass_Rate_Adjust_Basic(thisNode,interrupt,interruptProcedure,rateAdjustment)
    !% Return the node mass rate of change.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    logical,                   intent(inout) :: interrupt
    procedure(),      pointer, intent(inout) :: interruptProcedure
    double precision,          intent(in)    :: rateAdjustment
    integer                                  :: thisIndex

    ! Use the stored accretion rate value.
    thisIndex=Tree_Node_Basic_Index(thisNode)
    thisNode%components(thisIndex)%properties(massIndex,propertyDerivative)=thisNode%components(thisIndex)%properties(massIndex&
         &,propertyDerivative)+rateAdjustment
    return
  end subroutine Tree_Node_Mass_Rate_Adjust_Basic

  subroutine Tree_Node_Mass_Rate_Compute_Basic(thisNode,interrupt,interruptProcedure)
    !% Compute the node mass rate of change.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    logical,                 intent(inout) :: interrupt
    procedure(),    pointer, intent(inout) :: interruptProcedure
    integer                                :: thisIndex

    thisIndex=Tree_Node_Basic_Index(thisNode)
    call Tree_Node_Mass_Rate_Adjust_Basic(thisNode,interrupt,interruptProcedure,thisNode%components(thisIndex)%data(rateIndex))
    return
  end subroutine Tree_Node_Mass_Rate_Compute_Basic

  double precision function Tree_Node_Mass_Accretion_Rate_Basic(thisNode)
    !% Returns the mass accretion rate for {\tt thisNode}.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    integer                                :: thisIndex

    thisIndex=Tree_Node_Basic_Index(thisNode)
    Tree_Node_Mass_Accretion_Rate_Basic=thisNode%components(thisIndex)%data(rateIndex)
    return
  end function Tree_Node_Mass_Accretion_Rate_Basic

  double precision function Tree_Node_Time_Basic(thisNode)
    !% Return the node time.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    integer                                :: thisIndex

    thisIndex=Tree_Node_Basic_Index(thisNode)
    Tree_Node_Time_Basic=thisNode%components(thisIndex)%properties(timeIndex,propertyValue)
    return
  end function Tree_Node_Time_Basic

  subroutine Tree_Node_Time_Set_Basic(thisNode,time)
    !% Set the node time.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: time
    integer                                  :: thisIndex

    thisIndex=Tree_Node_Basic_Index(thisNode)
    thisNode%components(thisIndex)%properties(timeIndex,propertyValue)=time
    return
  end subroutine Tree_Node_Time_Set_Basic

  subroutine Tree_Node_Time_Rate_Adjust_Basic(thisNode,interrupt,interruptProcedure,rateAdjustment)
    !% Return the node time rate of change.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    logical,                   intent(inout) :: interrupt
    procedure(),      pointer, intent(inout) :: interruptProcedure
    double precision,          intent(in)    :: rateAdjustment
    integer                                  :: thisIndex
    
    thisIndex=Tree_Node_Basic_Index(thisNode)
    thisNode%components(thisIndex)%properties(timeIndex,propertyDerivative)=thisNode%components(thisIndex)%properties(timeIndex&
         &,propertyDerivative)+rateAdjustment
    return
  end subroutine Tree_Node_Time_Rate_Adjust_Basic

  subroutine Tree_Node_Time_Rate_Compute_Basic(thisNode,interrupt,interruptProcedure)
    !% Compute the node time rate of change.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    logical,                 intent(inout) :: interrupt
    procedure(),    pointer, intent(inout) :: interruptProcedure
    
    call Tree_Node_Time_Rate_Adjust_Basic(thisNode,interrupt,interruptProcedure,1.0d0) ! By definition.
    return
  end subroutine Tree_Node_Time_Rate_Compute_Basic

  double precision function Tree_Node_Time_Last_Isolated_Basic(thisNode)
    !% Returns the time at which {\tt thisNode} was last an isolated halo.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    integer                                :: thisIndex

    if (thisNode%isSatellite()) then
       thisIndex=Tree_Node_Basic_Index(thisNode)
       Tree_Node_Time_Last_Isolated_Basic=thisNode%components(thisIndex)%data(isolationIndex)
    else
       Tree_Node_Time_Last_Isolated_Basic=Tree_Node_Time_Basic(thisNode)
    end if
    return
  end function Tree_Node_Time_Last_Isolated_Basic


  !# <mergerTreeInitializeTask>
  !#  <unitName>Halo_Mass_Accretion_Rate</unitName>
  !# </mergerTreeInitializeTask>
  subroutine Halo_Mass_Accretion_Rate(thisNode)
    !% Set the mass accretion rate for {\tt thisNode}.
    implicit none
    type(treeNode),  pointer, intent(inout) :: thisNode
    integer                                 :: thisIndex
    double precision                        :: deltaTime

    if (methodSelected) then
       ! Get component index for this node.
       thisIndex=Tree_Node_Basic_Index(thisNode)
       ! Compute the mass accretion rate for this node if it's the primary progenitor.
       if (thisNode%isPrimaryProgenitor()) then
          deltaTime=Tree_Node_Time(thisNode%parentNode)-Tree_Node_Time(thisNode)
          if (deltaTime > 0.0d0) thisNode%components(thisIndex)%data(rateIndex)=Unresolved_Mass(thisNode%parentNode)/deltaTime
       end if
    end if
    return
  end subroutine Halo_Mass_Accretion_Rate

  !# <nodeMergerTask>
  !#  <unitName>Tree_Node_Mass_Stop_Accretion</unitName>
  !# </nodeMergerTask>
  subroutine Tree_Node_Mass_Stop_Accretion(thisNode)
    !% Switch off accretion of new mass onto this node once it becomes a satellite.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    integer                                :: thisIndex

    if (methodSelected) then
       thisIndex=Tree_Node_Basic_Index(thisNode)
       ! Shut down mass accretion onto the halo now that it is a satellite.
       thisNode%components(thisIndex)%data(rateIndex)=0.0d0
       ! Record the time at which the node became a satellite - used for computing halo scales etc.
       thisNode%components(thisIndex)%data(isolationIndex)=Tree_Node_Time(thisNode)
    end if
    return
  end subroutine Tree_Node_Mass_Stop_Accretion

  !# <nodePromotionTask>
  !#  <unitName>Tree_Node_Basic_Promote</unitName>
  !# </nodePromotionTask>
  subroutine Tree_Node_Basic_Promote(thisNode)
    !% Ensure that {\tt thisNode} is ready for promotion to its parent. In this case, we simply update the mass of {\tt thisNode}
    !% to be that of its parent.
    use Galacticus_Error
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    type(treeNode), pointer                :: parentNode
    integer                                :: thisIndex,parentIndex

    if (methodSelected) then
       parentNode => thisNode%parentNode
       if (Tree_Node_Time(thisNode) /= Tree_Node_Time(parentNode)) call Galacticus_Error_Report('Tree_Node_Basic_Promote','thisNode&
            & has not been evolved to its parent')
       ! Adjust the mass to that of the parent node.
       call Tree_Node_Mass_Set_Basic(thisNode,Tree_Node_Mass_Basic(parentNode))
       ! Adjust the accretion rate to that of the parent node.
       thisIndex  =Tree_Node_Basic_Index(thisNode)
       parentIndex=Tree_Node_Basic_Index(parentNode)
       thisNode%components(thisIndex)%data(rateIndex)=parentNode%components(parentIndex)%data(rateIndex)
    end if
    return
  end subroutine Tree_Node_Basic_Promote
  

  integer function Tree_Node_Basic_Index(thisNode)
    !% Ensure the basic component exists and return its position in the components array.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode

    call thisNode%createComponent(componentIndex,propertyCount,dataCount,historyCount)
    Tree_Node_Basic_Index=thisNode%componentIndex(componentIndex)
    return
  end function Tree_Node_Basic_Index

  double precision function Unresolved_Mass(thisNode)
    !% Return the unresolved mass for {\tt thisNode}.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    type(treeNode), pointer                :: childNode

    Unresolved_Mass=Tree_Node_Mass(thisNode)
    childNode => thisNode%childNode
    do while (associated(childNode))
       Unresolved_Mass=Unresolved_Mass-Tree_Node_Mass(childNode)
       childNode => childNode%siblingNode
    end do
    return
  end function Unresolved_Mass

  
  !# <mergerTreeOutputNames>
  !#  <unitName>Galacticus_Output_Tree_Basic_Names</unitName>
  !#  <sortName>Galacticus_Output_Tree_Basic</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Tree_Basic_Names(integerProperty,integerPropertyNames,integerPropertyComments,doubleProperty&
       &,doublePropertyNames,doublePropertyComments,time)
    !% Set names of basic properties to be written to the \glc\ output file.
    implicit none
    double precision, intent(in)                  :: time
    integer,          intent(inout)               :: integerProperty,doubleProperty
    character(len=*), intent(inout), dimension(:) :: integerPropertyNames,integerPropertyComments,doublePropertyNames &
         &,doublePropertyComments
    
    if (methodSelected) then
       doubleProperty=doubleProperty+1
       doublePropertyNames   (doubleProperty)='nodeMass'
       doublePropertyComments(doubleProperty)='Total mass of the node, assuming univeral baryon fraction.'
       doubleProperty=doubleProperty+1
       doublePropertyNames   (doubleProperty)='nodeTimeLastIsolated'
       doublePropertyComments(doubleProperty)='Time at which node was last an isolated halo.'
    end if
    return
  end subroutine Galacticus_Output_Tree_Basic_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Galacticus_Output_Tree_Basic_Property_Count</unitName>
  !#  <sortName>Galacticus_Output_Tree_Basic</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Galacticus_Output_Tree_Basic_Property_Count(integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of basic properties to be written to the the \glc\ output file.
    implicit none
    double precision, intent(in)    :: time
    integer,          intent(inout) :: integerPropertyCount,doublePropertyCount

    if (methodSelected) doublePropertyCount=doublePropertyCount+2
    return
  end subroutine Galacticus_Output_Tree_Basic_Property_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Galacticus_Output_Tree_Basic</unitName>
  !#  <sortName>Galacticus_Output_Tree_Basic</sortName>
  !# </mergerTreeOutputTask>
  subroutine Galacticus_Output_Tree_Basic(thisNode,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time)
    !% Store basic properties in the \glc\ output file buffers.
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
       doubleBuffer(doubleBufferCount,doubleProperty)=Tree_Node_Mass(thisNode)
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=Tree_Node_Time_Last_Isolated_Basic(thisNode)
    end if
    return
  end subroutine Galacticus_Output_Tree_Basic

  !# <nodeDumpTask>
  !#  <unitName>Tree_Node_Methods_Basic_Dump</unitName>
  !# </nodeDumpTask>
  subroutine Tree_Node_Methods_Basic_Dump(thisNode)
    !% Dump all properties of {\tt thisNode} to screen.
    implicit none
    type(treeNode),   intent(inout), pointer     :: thisNode

    if (methodSelected) then
       if (thisNode%componentExists(componentIndex)) then
          write (0,'(1x,a)'           ) 'basic node component -> properties:'
          write (0,'(2x,a50,1x,e12.6)') 'node mass:'              ,Tree_Node_Mass(thisNode)
          write (0,'(2x,a50,1x,e12.6)') 'node time:'              ,Tree_Node_Time(thisNode)
          write (0,'(2x,a50,1x,e12.6)') 'node time last isolated:',Tree_Node_Time_Last_Isolated_Basic(thisNode)
       else
          write (0,'(1x,a)'           ) 'basic node component -> nonexistant'
       end if
    end if
    return
  end subroutine Tree_Node_Methods_Basic_Dump

end module Tree_Node_Methods_Basic
