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


!% Contains a module which implements a dark matter profile method that provides a scale radius.

module Tree_Node_Methods_Dark_Matter_Profile_Scales
  !% Implements a dark matter profile method that provides a scale radius.
  use Tree_Nodes
  use Components
  private
  public :: Tree_Node_Methods_Profile_Scale_Initialize, Tree_Node_Methods_Profile_Scale_Initialize_Scale,&
       & Galacticus_Output_Tree_Profile_Scale, Galacticus_Output_Tree_Profile_Scale_Property_Count,&
       & Galacticus_Output_Tree_Profile_Scale_Names, Tree_Node_Methods_Profile_Scale_Dump,&
       & Tree_Node_Dark_Matter_Profile_Scale_Promote, Tree_Node_Methods_Profile_Scale_Merger_Tree_Output, Profile_Scale_Scale_Set
  
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

  ! Flag indicating whether scale radius data should be output when full merger trees are output.
  logical          :: mergerTreeStructureOutputDarkMatterScaleRadius

contains


  !# <treeNodeCreateInitialize>
  !#  <unitName>Tree_Node_Methods_Profile_Scale_Initialize</unitName>
  !#  <optionName default="scale">treeNodeMethodDarkMatterProfile</optionName>
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
       !@ <inputParameter>
       !@   <name>mergerTreeStructureOutputDarkMatterScaleRadius</name>
       !@   <defaultValue>false</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Determines whether or not dark matter halo scale radius is included in outputs of merger trees.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeStructureOutputDarkMatterScaleRadius',mergerTreeStructureOutputDarkMatterScaleRadius,defaultValue=.false.)

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
    procedure(), pointer, intent(inout) :: interruptProcedure
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
    procedure(), pointer, intent(inout) :: interruptProcedure
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

    ! Get component index for this node.
    thisIndex=Tree_Node_Profile_Scale_Index(thisNode)
    ! Return the pre-computed growth rate.
    Tree_Node_Dark_Matter_Profile_Scale_Growth_Rate_Scale=thisNode%components(thisIndex)%data(scaleRateIndex)
    return
  end function Tree_Node_Dark_Matter_Profile_Scale_Growth_Rate_Scale


  !# <mergerTreeInitializeTask>
  !#  <unitName>Tree_Node_Methods_Profile_Scale_Initialize_Scale</unitName>
  !# </mergerTreeInitializeTask>
  subroutine Tree_Node_Methods_Profile_Scale_Initialize_Scale(thisNode)
    !% Initialize the scale radius of {\tt thisNode}.
    use Dark_Matter_Profiles_Concentrations
    use Dark_Matter_Halo_Scales
    implicit none
    type(treeNode),  pointer, intent(inout) :: thisNode
    integer                                 :: thisIndex
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
  

  integer function Tree_Node_Profile_Scale_Index(thisNode)
    !% Ensure the profile component exists and return its position in the components array.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    
    call thisNode%createComponent(componentIndex,propertyCount,dataCount,historyCount)
    Tree_Node_Profile_Scale_Index=thisNode%componentIndex(componentIndex)
    return
  end function Tree_Node_Profile_Scale_Index

  !# <scaleSetTask>
  !#  <unitName>Profile_Scale_Scale_Set</unitName>
  !# </scaleSetTask>
  subroutine Profile_Scale_Scale_Set(thisNode)
    !% Set scales for properties of {\tt thisNode}.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    integer                                  :: thisIndex
 
    ! Determine if method is active and a black hole component exists.
    if (methodSelected.and.thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Profile_Scale_Index(thisNode)

       ! Set scale for the scale radius.
       thisNode%components(thisIndex)%properties(scaleIndex,propertyScale)=thisNode%components(thisIndex)%properties(scaleIndex,propertyValue)

    end if
    return
  end subroutine Profile_Scale_Scale_Set
  
  !# <mergerTreeOutputNames>
  !#  <unitName>Galacticus_Output_Tree_Profile_Scale_Names</unitName>
  !#  <sortName>Galacticus_Output_Tree_Profile_Scale</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Tree_Profile_Scale_Names(integerProperty,integerPropertyNames,integerPropertyComments&
       &,integerPropertyUnitsSI,doubleProperty,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time)
    !% Set names of scale properties to be written to the \glc\ output file.
    use Numerical_Constants_Astronomical
    implicit none
    double precision, intent(in)                  :: time
    integer,          intent(inout)               :: integerProperty,doubleProperty
    character(len=*), intent(inout), dimension(:) :: integerPropertyNames,integerPropertyComments,doublePropertyNames &
         &,doublePropertyComments
    double precision, intent(inout), dimension(:) :: integerPropertyUnitsSI,doublePropertyUnitsSI

    if (methodSelected) then
       doubleProperty=doubleProperty+1
       doublePropertyNames   (doubleProperty)='darkMatterScaleRadius'
       doublePropertyComments(doubleProperty)='Scale radius of the dark matter profile [Mpc].'
       doublePropertyUnitsSI (doubleProperty)=megaParsec
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
    use Kind_Numbers
    implicit none
    double precision,        intent(in)             :: time
    type(treeNode),          intent(inout), pointer :: thisNode
    integer,                 intent(inout)          :: integerProperty,integerBufferCount,doubleProperty,doubleBufferCount
    integer(kind=kind_int8), intent(inout)          :: integerBuffer(:,:)
    double precision,        intent(inout)          :: doubleBuffer(:,:)

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

  !# <mergerTreeStructureOutputTask>
  !#  <unitName>Tree_Node_Methods_Profile_Scale_Merger_Tree_Output</unitName>
  !# </mergerTreeStructureOutputTask>
  subroutine Tree_Node_Methods_Profile_Scale_Merger_Tree_Output(baseNode,nodeProperty,treeGroup)
    !% Write the scale radius property to a full merger tree output.
    use IO_HDF5
    use Tree_Nodes
    use Numerical_Constants_Astronomical
    implicit none
    type(treeNode),   intent(in),    pointer      :: baseNode
    double precision, intent(inout), dimension(:) :: nodeProperty
    type(hdf5Object), intent(inout)               :: treeGroup
    type(treeNode),                  pointer      :: thisNode
    integer                                       :: nodeCount
    type(hdf5Object)                              :: nodeDataset

    ! Check if scale radius is to be included in merger tree outputs.
    if (methodSelected.and.mergerTreeStructureOutputDarkMatterScaleRadius) then
       
       ! Extract node scale radius and output to file.
       nodeCount=0
       thisNode => baseNode
       do while (associated(thisNode))
          nodeCount=nodeCount+1
          nodeProperty(nodeCount)=Tree_Node_Dark_Matter_Profile_Scale_Scale(thisNode)
          ! <gfortan 4.6> explicitly specify the target as thisNode since we can't use the "_Same_Node" tree walking procedures.
          call thisNode%walkTree(thisNode)
       end do
       call treeGroup%writeDataset(nodeProperty,'darkMatterScaleRadius','Scale radius of the dark matter profile [Mpc].',datasetReturned=nodeDataset)
       call nodeDataset%writeAttribute(megaParsec,"unitsInSI")
       call nodeDataset%close()
       
    end if

    return
  end subroutine Tree_Node_Methods_Profile_Scale_Merger_Tree_Output

end module Tree_Node_Methods_Dark_Matter_Profile_Scales
