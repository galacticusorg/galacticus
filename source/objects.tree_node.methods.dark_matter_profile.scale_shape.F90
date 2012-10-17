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

!% Contains a module which implements a dark matter profile method that provides a scale radius.

module Tree_Node_Methods_Dark_Matter_Profile_Scale_Shapes
  !% Implements a dark matter profile method that provides a scale radius and a shape parameter.
  use Tree_Nodes
  use Components
  implicit none
  private
  public :: Tree_Node_Methods_Profile_ScaleShape_Initialize, Tree_Node_Methods_Profile_ScaleShape_Initialize_Rates,&
       & Galacticus_Output_Tree_Profile_ScaleShape, Galacticus_Output_Tree_Profile_ScaleShape_Property_Count,&
       & Galacticus_Output_Tree_Profile_ScaleShape_Names, Tree_Node_Methods_Profile_ScaleShape_Dump,&
       & Tree_Node_Dark_Matter_Profile_ScaleShape_Promote, Tree_Node_Methods_Profile_ScaleShape_Merger_Tree_Output,&
       & Profile_ScaleShape_Scale_Set, Tree_Node_Methods_Profile_Scale_Shape_Plausibility
  
  ! The index used as a reference for this component.
  integer            :: componentIndex=-1

  ! Property indices.
  integer, parameter :: propertyCount =2, dataCount=2, historyCount=0
  integer, parameter :: scaleIndex    =1
  integer, parameter :: scaleRateIndex=1
  integer, parameter :: shapeIndex    =2
  integer, parameter :: shapeRateIndex=2

  ! Define procedure pointers.
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Dark_Matter_Profile_Scale</methodName>
  !# </treeNodeMethodsPointer>
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Dark_Matter_Profile_Scale_Growth_Rate</methodName>
  !# </treeNodeMethodsPointer>
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Dark_Matter_Profile_Shape</methodName>
  !# </treeNodeMethodsPointer>
  !# <treeNodeMethodsPointer>
  !#  <methodName>Tree_Node_Dark_Matter_Profile_Shape_Growth_Rate</methodName>
  !# </treeNodeMethodsPointer>

  ! Flag to indicate if this method is selected.
  logical          :: methodSelected=.false.

  ! Parameters of the method.
  double precision :: darkMatterProfileMinimumConcentration,darkMatterProfileMaximumConcentration

  ! Flag indicating whether scale radius and shape data should be output when full merger trees are output.
  logical          :: mergerTreeStructureOutputDarkMatterProfileProperties

contains

  !# <treeNodeCreateInitialize>
  !#  <unitName>Tree_Node_Methods_Profile_ScaleShape_Initialize</unitName>
  !#  <optionName>treeNodeMethodDarkMatterProfile</optionName>
  !# </treeNodeCreateInitialize>
  subroutine Tree_Node_Methods_Profile_ScaleShape_Initialize(componentOption,componentTypeCount)
    !% Initializes the tree node ``scale + shape'' dark matter profile method module.
    use ISO_Varying_String
    use Input_Parameters
    use String_Handling
    use Galacticus_Display
    implicit none
    type(varying_string), intent(in)    :: componentOption
    integer,              intent(inout) :: componentTypeCount
    type(varying_string)                :: message

    ! Check if this implementation is selected.
    if (componentOption == 'scale+shape') then
       ! Record that method is selected.
       methodSelected=.true.

       ! Increment the component count and store the value for later reference.
       componentTypeCount=componentTypeCount+1
       componentIndex    =componentTypeCount

       ! Display message.
       message='Scale+shape dark matter profile method selected [component index '
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
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('darkMatterProfileMinimumConcentration',darkMatterProfileMinimumConcentration,defaultValue=4.0d0)
       !@ <inputParameter>
       !@   <name>darkMatterProfileMaximumConcentration</name>
       !@   <defaultValue>100</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The maximum concentration allowed for dark matter profiles.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('darkMatterProfileMaximumConcentration',darkMatterProfileMaximumConcentration,defaultValue=100.0d0)
       !@ <inputParameter>
       !@   <name>mergerTreeStructureOutputDarkMatterProfileProperties</name>
       !@   <defaultValue>false</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Determines whether or not dark matter halo scale radius is included in outputs of merger trees.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@   <group>output</group>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeStructureOutputDarkMatterProfileProperties',mergerTreeStructureOutputDarkMatterProfileProperties,defaultValue=.false.)

       ! Set up procedure pointers.
       Tree_Node_Dark_Matter_Profile_Scale                          => Tree_Node_Dark_Matter_Profile_Scale_ScaleShape
       Tree_Node_Dark_Matter_Profile_Scale_Set                      => null()
       Tree_Node_Dark_Matter_Profile_Scale_Rate_Adjust              => Tree_Node_Dark_Matter_Profile_Scale_Rate_Adjust_ScaleShape
       Tree_Node_Dark_Matter_Profile_Scale_Rate_Compute             => Tree_Node_Dark_Matter_Profile_Scale_Rate_Compute_ScaleShape
       Tree_Node_Dark_Matter_Profile_Scale_Growth_Rate              => Tree_Node_Dark_Matter_Profile_Scale_Growth_Rate_ScaleShape
       Tree_Node_Dark_Matter_Profile_Scale_Growth_Rate_Set          => null()
       Tree_Node_Dark_Matter_Profile_Scale_Growth_Rate_Rate_Adjust  => null()
       Tree_Node_Dark_Matter_Profile_Scale_Growth_Rate_Rate_Compute => Tree_Node_Rate_Rate_Compute_Dummy
       Tree_Node_Dark_Matter_Profile_Shape                          => Tree_Node_Dark_Matter_Profile_Shape_ScaleShape
       Tree_Node_Dark_Matter_Profile_Shape_Set                      => null()
       Tree_Node_Dark_Matter_Profile_Shape_Rate_Adjust              => Tree_Node_Dark_Matter_Profile_Shape_Rate_Adjust_ScaleShape
       Tree_Node_Dark_Matter_Profile_Shape_Rate_Compute             => Tree_Node_Dark_Matter_Profile_Shape_Rate_Compute_ScaleShape
       Tree_Node_Dark_Matter_Profile_Shape_Growth_Rate              => Tree_Node_Dark_Matter_Profile_Shape_Growth_Rate_ScaleShape
       Tree_Node_Dark_Matter_Profile_Shape_Growth_Rate_Set          => null()
       Tree_Node_Dark_Matter_Profile_Shape_Growth_Rate_Rate_Adjust  => null()
       Tree_Node_Dark_Matter_Profile_Shape_Growth_Rate_Rate_Compute => Tree_Node_Rate_Rate_Compute_Dummy
    end if
    return
  end subroutine Tree_Node_Methods_Profile_ScaleShape_Initialize 
 
  double precision function Tree_Node_Dark_Matter_Profile_Scale_ScaleShape(thisNode,instance)
    !% Return the node dark matter profile scale length.
    use Dark_Matter_Halo_Scales
    implicit none
    integer,        optional, intent(in   ) :: instance
    type(treeNode), pointer,  intent(inout) :: thisNode
    integer                                 :: thisIndex
    double precision                        :: scaleLengthMaximum,scaleLengthMinimum

    ! Ensure the spin has been initialized.
    call Tree_Node_Methods_Profile_ScaleShape_Initialize_ScaleShape(thisNode)
    thisIndex=Tree_Node_Profile_ScaleShape_Index(thisNode)
    scaleLengthMaximum=Dark_Matter_Halo_Virial_Radius(thisNode)/darkMatterProfileMinimumConcentration
    scaleLengthMinimum=Dark_Matter_Halo_Virial_Radius(thisNode)/darkMatterProfileMaximumConcentration
    Tree_Node_Dark_Matter_Profile_Scale_ScaleShape=                                                 &
         & min(                                                                                     &
         &     scaleLengthMaximum                                                                 , &
         &     max(                                                                                 &
         &         scaleLengthMinimum                                                             , &
         &         thisNode%components(thisIndex)%instance(1)%properties(scaleIndex,propertyValue)  &
         &        )                                                                                 &
         &    )
    return
  end function Tree_Node_Dark_Matter_Profile_Scale_ScaleShape

  subroutine Tree_Node_Dark_Matter_Profile_Scale_Rate_Adjust_ScaleShape(thisNode,interrupt,interruptProcedure,rateAdjustment,instance)
    !% Adjust the rate of scale radius growth.
    implicit none
    integer, intent(in), optional :: instance
    type(treeNode),   pointer, intent(inout) :: thisNode
    logical,                   intent(inout) :: interrupt
    procedure(),      pointer, intent(inout) :: interruptProcedure
    double precision,          intent(in)    :: rateAdjustment
    integer                                  :: thisIndex

    ! Apply the change in rate.
    thisIndex=Tree_Node_Profile_ScaleShape_Index(thisNode)
    thisNode%components(thisIndex)%instance(1)%properties(scaleIndex,propertyDerivative)=thisNode%components(thisIndex)%instance(1)%properties(scaleIndex&
         &,propertyDerivative)+rateAdjustment
    return
  end subroutine Tree_Node_Dark_Matter_Profile_Scale_Rate_Adjust_ScaleShape

  subroutine Tree_Node_Dark_Matter_Profile_Scale_Rate_Compute_ScaleShape(thisNode,interrupt,interruptProcedure)
    !% Compute the rate of change of the scale radius.
    use Dark_Matter_Halo_Scales  
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    logical,                 intent(inout) :: interrupt
    procedure(),    pointer, intent(inout) :: interruptProcedure
    integer                                :: thisIndex
    double precision                       :: growthRate,concentration

    thisIndex=Tree_Node_Profile_ScaleShape_Index(thisNode)
    ! Find the concentration of this halo.
    concentration=Dark_Matter_Halo_Virial_Radius(thisNode)/Tree_Node_Dark_Matter_Profile_Scale(thisNode)
    ! Find the growth rate and limit to ensure minimum and maximum concentrations are not exceeded.
    if (concentration <= darkMatterProfileMinimumConcentration) growthRate=min(growthRate,Dark_Matter_Halo_Virial_Radius_Growth_Rate(thisNode)/darkMatterProfileMinimumConcentration)
    if (concentration >= darkMatterProfileMaximumConcentration) growthRate=max(growthRate,Dark_Matter_Halo_Virial_Radius_Growth_Rate(thisNode)/darkMatterProfileMaximumConcentration)
    ! Adjust the growth rate.
    growthRate=Tree_Node_Dark_Matter_Profile_Scale_Growth_Rate_ScaleShape(thisNode)
    call Tree_Node_Dark_Matter_Profile_Scale_Rate_Adjust_ScaleShape(thisNode,interrupt,interruptProcedure,growthRate)
    return
  end subroutine Tree_Node_Dark_Matter_Profile_Scale_Rate_Compute_ScaleShape

  double precision function Tree_Node_Dark_Matter_Profile_Scale_Growth_Rate_ScaleShape(thisNode,instance)
    !% Return the node dark matter profile scale length rate of growth (assumed to be zero).
    implicit none
    integer, intent(in), optional :: instance
    type(treeNode),   pointer, intent(inout) :: thisNode
    integer                                  :: thisIndex

    ! Get component index for this node.
    thisIndex=Tree_Node_Profile_ScaleShape_Index(thisNode)
    ! Return the pre-computed growth rate.
    Tree_Node_Dark_Matter_Profile_Scale_Growth_Rate_ScaleShape=thisNode%components(thisIndex)%instance(1)%data(scaleRateIndex)
    return
  end function Tree_Node_Dark_Matter_Profile_Scale_Growth_Rate_ScaleShape

  double precision function Tree_Node_Dark_Matter_Profile_Shape_ScaleShape(thisNode,instance)
    !% Return the node dark matter profile shape parameter.
    implicit none
    integer, intent(in), optional :: instance
    type(treeNode), pointer, intent(inout) :: thisNode
    integer                                :: thisIndex
    
    ! Ensure the spin has been initialized.
    call Tree_Node_Methods_Profile_ScaleShape_Initialize_ScaleShape(thisNode)
    thisIndex=Tree_Node_Profile_ScaleShape_Index(thisNode)
    Tree_Node_Dark_Matter_Profile_Shape_ScaleShape=thisNode%components(thisIndex)%instance(1)%properties(shapeIndex,propertyValue)
    return
  end function Tree_Node_Dark_Matter_Profile_Shape_ScaleShape

  subroutine Tree_Node_Dark_Matter_Profile_Shape_Rate_Adjust_ScaleShape(thisNode,interrupt,interruptProcedure,rateAdjustment,instance)
    !% Adjust the rate of shape parameter growth.
    implicit none
    integer, intent(in), optional :: instance
    type(treeNode),   pointer, intent(inout) :: thisNode
    logical,                   intent(inout) :: interrupt
    procedure(),      pointer, intent(inout) :: interruptProcedure
    double precision,          intent(in)    :: rateAdjustment
    integer                                  :: thisIndex

    ! Apply the change in rate.
    thisIndex=Tree_Node_Profile_ScaleShape_Index(thisNode)
    thisNode%components(thisIndex)%instance(1)%properties(shapeIndex,propertyDerivative)=thisNode%components(thisIndex)%instance(1)%properties(shapeIndex&
         &,propertyDerivative)+rateAdjustment
    return
  end subroutine Tree_Node_Dark_Matter_Profile_Shape_Rate_Adjust_ScaleShape

  subroutine Tree_Node_Dark_Matter_Profile_Shape_Rate_Compute_ScaleShape(thisNode,interrupt,interruptProcedure)
    !% Compute the rate of change of the shape parameter.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    logical,                 intent(inout) :: interrupt
    procedure(),    pointer, intent(inout) :: interruptProcedure
    integer                                :: thisIndex
    double precision                       :: growthRate

    thisIndex=Tree_Node_Profile_ScaleShape_Index(thisNode)
    growthRate=Tree_Node_Dark_Matter_Profile_Scale_Growth_Rate_ScaleShape(thisNode)
    call Tree_Node_Dark_Matter_Profile_Shape_Rate_Adjust_ScaleShape(thisNode,interrupt,interruptProcedure,growthRate)
    return
  end subroutine Tree_Node_Dark_Matter_Profile_Shape_Rate_Compute_ScaleShape

  double precision function Tree_Node_Dark_Matter_Profile_Shape_Growth_Rate_ScaleShape(thisNode,instance)
    !% Return the node dark matter profile shape parameter rate of growth.
    implicit none
    integer,          optional, intent(in   ) :: instance
    type(treeNode),   pointer,  intent(inout) :: thisNode
    integer                                   :: thisIndex

    ! Get component index for this node.
    thisIndex=Tree_Node_Profile_ScaleShape_Index(thisNode)
    ! Return the pre-computed growth rate.
    Tree_Node_Dark_Matter_Profile_Shape_Growth_Rate_ScaleShape=thisNode%components(thisIndex)%instance(1)%data(shapeRateIndex)
    return
  end function Tree_Node_Dark_Matter_Profile_Shape_Growth_Rate_ScaleShape

  subroutine Tree_Node_Methods_Profile_ScaleShape_Initialize_ScaleShape(thisNode)
    !% Initialize the scale radius and shape parameter of {\tt thisNode}.
    use Dark_Matter_Profiles_Concentrations
    use Dark_Matter_Halo_Scales
     use Dark_Matter_Profiles_Shapes
    implicit none
    type(treeNode),  pointer, intent(inout) :: thisNode
    integer                                 :: thisIndex
    double precision                        :: concentration,deltaTime

    if (methodSelected) then
       if (.not.thisNode%componentExists(componentIndex)) then
          thisIndex=Tree_Node_Profile_ScaleShape_Index(thisNode)
          ! Set the scale radius parameter of the halo.
          concentration=max(Dark_Matter_Profile_Concentration(thisNode),darkMatterProfileMinimumConcentration)
          thisNode%components(thisIndex)%instance(1)%properties(scaleIndex,propertyValue)=Dark_Matter_Halo_Virial_Radius(thisNode)/concentration
          ! Set the shape parameter of the halo.
          thisNode%components(thisIndex)%instance(1)%properties(shapeIndex,propertyValue)=Dark_Matter_Profile_Shape     (thisNode)
       end if
    end if
    return
  end subroutine Tree_Node_Methods_Profile_ScaleShape_Initialize_ScaleShape

  !# <radiusSolverPlausibility>
  !#  <unitName>Tree_Node_Methods_Profile_Scale_Shape_Plausibility</unitName>
  !# </radiusSolverPlausibility>
  subroutine Tree_Node_Methods_Profile_Scale_Shape_Plausibility(thisNode,galaxyIsPhysicallyPlausible)
    !% Determines whether the dark matter profile is physically plausible for radius solving tasks.
    use Dark_Matter_Halo_Scales
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    logical,                 intent(inout) :: galaxyIsPhysicallyPlausible

    ! Return immediately if our method is not selected.
    if (.not.methodSelected) return
    ! Determine the plausibility of the current disk.
    if (Tree_Node_Dark_Matter_Profile_Scale(thisNode) <= 0.0d0) galaxyIsPhysicallyPlausible=.false.
    return
  end subroutine Tree_Node_Methods_Profile_Scale_Shape_Plausibility

  !# <mergerTreeInitializeTask>
  !#  <unitName>Tree_Node_Methods_Profile_ScaleShape_Initialize_Rates</unitName>
  !# </mergerTreeInitializeTask>
  subroutine Tree_Node_Methods_Profile_ScaleShape_Initialize_Rates(thisNode)
    !% Initialize the scale radius and shape parameter of {\tt thisNode}.
    use Dark_Matter_Profiles_Concentrations
    use Dark_Matter_Halo_Scales
    implicit none
    type(treeNode),  pointer, intent(inout) :: thisNode
    integer                                 :: thisIndex
    double precision                        :: deltaTime

    if (methodSelected) then
       ! Ensure that current node has its scale set.
       call Tree_Node_Methods_Profile_ScaleShape_Initialize_ScaleShape(thisNode)
       thisIndex=Tree_Node_Profile_ScaleShape_Index(thisNode)
       ! Check if this node is the primary progenitor.
       if (thisNode%isPrimaryProgenitor()) then
          ! It is, so compute the scale radius and shape parameter growth rates.
          ! First ensure that parent node has scale radius and shape parameter set.
          call Tree_Node_Methods_Profile_ScaleShape_Initialize_ScaleShape(thisNode%parentNode)
          ! Now compute the growth rate.
          deltaTime=Tree_Node_Time(thisNode%parentNode)-Tree_Node_Time(thisNode)
          if (deltaTime > 0.0d0) then
             thisNode%components(thisIndex)%instance(1)%data(scaleRateIndex)=(Tree_Node_Dark_Matter_Profile_Scale_ScaleShape(thisNode%parentNode)&
                  &-Tree_Node_Dark_Matter_Profile_Scale_ScaleShape(thisNode))/deltaTime
             thisNode%components(thisIndex)%instance(1)%data(shapeRateIndex)=(Tree_Node_Dark_Matter_Profile_Shape_ScaleShape(thisNode%parentNode)&
                  &-Tree_Node_Dark_Matter_Profile_Shape_ScaleShape(thisNode))/deltaTime
          else
             thisNode%components(thisIndex)%instance(1)%data(scaleRateIndex)=0.0d0
             thisNode%components(thisIndex)%instance(1)%data(shapeRateIndex)=0.0d0
          end if
       else
          ! It is not, so set scale radius and shape parameter growth rates to zero.
          thisNode%components(thisIndex)%instance(1)%data(scaleRateIndex)=0.0d0
          thisNode%components(thisIndex)%instance(1)%data(shapeRateIndex)=0.0d0
       end if
    end if
    return
  end subroutine Tree_Node_Methods_Profile_ScaleShape_Initialize_Rates
     
  !# <nodePromotionTask>
  !#  <unitName>Tree_Node_Dark_Matter_Profile_ScaleShape_Promote</unitName>
  !# </nodePromotionTask>
  subroutine Tree_Node_Dark_Matter_Profile_ScaleShape_Promote(thisNode)
    !% Ensure that {\tt thisNode} is ready for promotion to its parent. In this case, we simply update the growth rate of {\tt thisNode}
    !% to be that of its parent.
    use Galacticus_Error
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    type(treeNode), pointer                :: parentNode
    integer                                :: thisIndex,parentIndex

    if (methodSelected) then
       parentNode => thisNode%parentNode
       if (Tree_Node_Time(thisNode) /= Tree_Node_Time(parentNode)) call Galacticus_Error_Report('Tree_Node_Dark_Matter_Profile_ScaleShape_Promote','thisNode&
            & has not been evolved to its parent')
       ! Get component indices.
       thisIndex  =Tree_Node_Profile_ScaleShape_Index(thisNode)
       parentIndex=Tree_Node_Profile_ScaleShape_Index(parentNode)
       ! Adjust the scale radius to that of the parent node.
       thisNode%components(thisIndex)%instance(1)%properties(scaleIndex,propertyValue)=parentNode%components(parentIndex)%instance(1)%properties(scaleIndex,propertyValue)
       ! Adjust the growth rate to that of the parent node.
       thisNode%components(thisIndex)%instance(1)%data      (scaleRateIndex          )=parentNode%components(parentIndex)%instance(1)%data      (scaleRateIndex          )
       ! Adjust the shape parameter to that of the parent node.
       thisNode%components(thisIndex)%instance(1)%properties(shapeIndex,propertyValue)=parentNode%components(parentIndex)%instance(1)%properties(shapeIndex,propertyValue)
       ! Adjust the growth rate to that of the parent node.
       thisNode%components(thisIndex)%instance(1)%data      (shapeRateIndex          )=parentNode%components(parentIndex)%instance(1)%data      (shapeRateIndex          )
    end if
    return
  end subroutine Tree_Node_Dark_Matter_Profile_ScaleShape_Promote
  
  integer function Tree_Node_Profile_ScaleShape_Index(thisNode)
    !% Ensure the profile component exists and return its position in the components array.
    implicit none
    type(treeNode), pointer, intent(inout) :: thisNode
    
    call thisNode%createComponent(componentIndex,propertyCount,dataCount,historyCount)
    Tree_Node_Profile_ScaleShape_Index=thisNode%componentIndex(componentIndex)
    return
  end function Tree_Node_Profile_ScaleShape_Index

  !# <scaleSetTask>
  !#  <unitName>Profile_ScaleShape_Scale_Set</unitName>
  !# </scaleSetTask>
  subroutine Profile_ScaleShape_Scale_Set(thisNode)
    !% Set scales for properties of {\tt thisNode}.
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    integer                                  :: thisIndex
 
    ! Determine if method is active and a profile component exists.
    if (methodSelected.and.thisNode%componentExists(componentIndex)) then
       thisIndex=Tree_Node_Profile_ScaleShape_Index(thisNode)

       ! Set scale for the scale radius.
       thisNode%components(thisIndex)%instance(1)%properties(scaleIndex,propertyScale)=thisNode%components(thisIndex)%instance(1)%properties(scaleIndex,propertyValue)
       ! Set scale for the shape parameter.
       thisNode%components(thisIndex)%instance(1)%properties(shapeIndex,propertyScale)=thisNode%components(thisIndex)%instance(1)%properties(shapeIndex,propertyValue)

    end if
    return
  end subroutine Profile_ScaleShape_Scale_Set
  
  !# <mergerTreeOutputNames>
  !#  <unitName>Galacticus_Output_Tree_Profile_ScaleShape_Names</unitName>
  !#  <sortName>Galacticus_Output_Tree_Profile_ScaleShape</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Tree_Profile_ScaleShape_Names(integerProperty,integerPropertyNames,integerPropertyComments&
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
       !@ <outputPropertyGroup>
       !@   <name>darkMatterProfile</name>
       !@   <description>Dark matter profile properities</description>
       !@   <outputType>nodeData</outputType>
       !@ </outputPropertyGroup>
       doubleProperty=doubleProperty+1
       !@ <outputProperty>
       !@   <name>darkMatterScaleRadius</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>Scale radius of the dark matter profile [Mpc].</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@   <group>darkMatterProfile</group>
       !@ </outputProperty>
       doublePropertyNames   (doubleProperty)='darkMatterScaleRadius'
       doublePropertyComments(doubleProperty)='Scale radius of the dark matter profile [Mpc].'
       doublePropertyUnitsSI (doubleProperty)=megaParsec
       doubleProperty=doubleProperty+1
       !@ <outputProperty>
       !@   <name>darkMatterShapeParameter</name>
       !@   <datatype>real</datatype>
       !@   <cardinality>0..1</cardinality>
       !@   <description>Shape parameter of the dark matter profile.</description>
       !@   <label>???</label>
       !@   <outputType>nodeData</outputType>
       !@   <group>darkMatterProfile</group>
       !@ </outputProperty>
       doublePropertyNames   (doubleProperty)='darkMatterShapeParameter'
       doublePropertyComments(doubleProperty)='Shape parameter of the dark matter profile.'
       doublePropertyUnitsSI (doubleProperty)=0.0d0
    end if
    return
  end subroutine Galacticus_Output_Tree_Profile_ScaleShape_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Galacticus_Output_Tree_Profile_ScaleShape_Property_Count</unitName>
  !#  <sortName>Galacticus_Output_Tree_Profile_ScaleShape</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Galacticus_Output_Tree_Profile_ScaleShape_Property_Count(integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of spin properties to be written to the the \glc\ output file.
    implicit none
    double precision, intent(in)    :: time
    integer,          intent(inout) :: integerPropertyCount,doublePropertyCount

    if (methodSelected) doublePropertyCount=doublePropertyCount+2
    return
  end subroutine Galacticus_Output_Tree_Profile_ScaleShape_Property_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Galacticus_Output_Tree_Profile_ScaleShape</unitName>
  !#  <sortName>Galacticus_Output_Tree_Profile_ScaleShape</sortName>
  !# </mergerTreeOutputTask>
  subroutine Galacticus_Output_Tree_Profile_ScaleShape(thisNode,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
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
       doubleProperty=doubleProperty+1
       doubleBuffer(doubleBufferCount,doubleProperty)=Tree_Node_Dark_Matter_Profile_Shape(thisNode)
    end if
    return
  end subroutine Galacticus_Output_Tree_Profile_ScaleShape
  
  !# <nodeDumpTask>
  !#  <unitName>Tree_Node_Methods_Profile_ScaleShape_Dump</unitName>
  !# </nodeDumpTask>
  subroutine Tree_Node_Methods_Profile_ScaleShape_Dump(thisNode)
    !% Dump all properties of {\tt thisNode} to screen.
    implicit none
    type(treeNode), intent(inout), pointer :: thisNode

    if (methodSelected) then
       if (thisNode%componentExists(componentIndex)) then
          write (0,'(1x,a)'           ) 'dark matter profile properties:'
          write (0,'(2x,a50,1x,e12.6)') 'scale:',Tree_Node_Dark_Matter_Profile_Scale(thisNode)
          write (0,'(2x,a50,1x,e12.6)') 'shape:',Tree_Node_Dark_Matter_Profile_Shape(thisNode)
       else
          write (0,'(1x,a)'           ) 'no dark matter profile component'
       end if
    end if
    return
  end subroutine Tree_Node_Methods_Profile_ScaleShape_Dump

  !# <mergerTreeStructureOutputTask>
  !#  <unitName>Tree_Node_Methods_Profile_ScaleShape_Merger_Tree_Output</unitName>
  !# </mergerTreeStructureOutputTask>
  subroutine Tree_Node_Methods_Profile_ScaleShape_Merger_Tree_Output(baseNode,nodeProperty,treeGroup)
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
    
    ! Check if scale radius and shape parameter are to be included in merger tree outputs.
    if (methodSelected.and.mergerTreeStructureOutputDarkMatterProfileProperties) then
       
       ! Extract node scale radius and output to file.
       nodeCount=0
       thisNode => baseNode
       do while (associated(thisNode))
          nodeCount=nodeCount+1
          nodeProperty(nodeCount)=Tree_Node_Dark_Matter_Profile_Scale_ScaleShape(thisNode)
          call thisNode%walkTree(thisNode)
       end do
       call treeGroup%writeDataset(nodeProperty,'darkMatterScaleRadius','Scale radius of the dark matter profile [Mpc].',datasetReturned=nodeDataset)
       call nodeDataset%writeAttribute(megaParsec,"unitsInSI")
       
       ! Extract node shape parameter and output to file.
       nodeCount=0
       thisNode => baseNode
       do while (associated(thisNode))
          nodeCount=nodeCount+1
          nodeProperty(nodeCount)=Tree_Node_Dark_Matter_Profile_Shape_ScaleShape(thisNode)
          call thisNode%walkTree(thisNode)
       end do
       call treeGroup%writeDataset(nodeProperty,'darkMatterShapeParameter','Shape parameter of the dark matter profile.')
       
    end if

    return
  end subroutine Tree_Node_Methods_Profile_ScaleShape_Merger_Tree_Output

end module Tree_Node_Methods_Dark_Matter_Profile_Scale_Shapes
