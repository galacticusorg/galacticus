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

!% Contains a module which implements a dark matter profile method that provides a scale radius and a shape parameter.

module Node_Component_Dark_Matter_Profile_Scale_Shape
  !% Implements a dark matter profile method that provides a scale radius and a shape parameter.
  use Galacticus_Nodes
  implicit none
  private
  public :: Node_Component_Dark_Matter_Profile_Scale_Shape_Rate_Compute, Node_Component_Dark_Matter_Profile_Scale_Shape_Tree_Initialize, &
       &    Node_Component_Dark_Matter_Profile_Scale_Shape_Promote     , Node_Component_Dark_Matter_Profile_Scale_Shape_Scale_Set      , &
       &    Node_Component_Dark_Matter_Profile_Scale_Shape_Tree_Output

  !# <component>
  !#  <class>darkMatterProfile</class>
  !#  <name>scaleShape</name>
  !#  <isDefault>no</isDefault>
  !#  <extends>
  !#   <class>darkMatterProfile</class>
  !#   <name>scale</name>
  !#  </extends>
  !#  <methods>
  !#   <method>
  !#     <name>shape</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="true" />
  !#     <output unitsInSI="0.0d0" comment="Shape parameter of the dark matter profile."/>
  !#   </method>
  !#   <method>
  !#     <name>shapeGrowthRate</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#   </method>
  !#  </methods>
  !# </component>

  ! Flag indicating whether scale radius and shape data should be output when full merger trees are output.
  logical :: mergerTreeStructureOutputDarkMatterProfileShape

  ! Record of whether the module has been initialized.
  logical :: moduleInitialized=.false.

contains

  subroutine Node_Component_Dark_Matter_Profile_Scale_Shape_Initialize()
    !% Initializes the ``scale'' implementation of the dark matter halo profile component.
    use Input_Parameters
    implicit none

    ! Check if this implementation is selected.
    !$omp critical (Node_Component_Dark_Matter_Profile_Scale_Shape_Initialize)
    if (.not.moduleInitialized) then
       !@ <inputParameter>
       !@   <name>mergerTreeStructureOutputDarkMatterProfileShape</name>
       !@   <defaultValue>false</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Determines whether or not dark matter halo shape parameter is included in outputs of merger trees.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@   <group>output</group>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeStructureOutputDarkMatterProfileShape',mergerTreeStructureOutputDarkMatterProfileShape,defaultValue=.false.)
       ! Record that the module is now initialize.
       moduleInitialized=.true.
    end if
    !$omp end critical (Node_Component_Dark_Matter_Profile_Scale_Shape_Initialize)
    return
  end subroutine Node_Component_Dark_Matter_Profile_Scale_Shape_Initialize

  !# <rateComputeTask>
  !#  <unitName>Node_Component_Dark_Matter_Profile_Scale_Shape_Rate_Compute</unitName>
  !# </rateComputeTask>
  subroutine Node_Component_Dark_Matter_Profile_Scale_Shape_Rate_Compute(thisNode,interrupt,interruptProcedure)
    !% Compute the rate of change of the scale radius.
    use Dark_Matter_Halo_Scales  
    implicit none
    type (treeNode                      ), pointer, intent(inout) :: thisNode
    logical,                                        intent(inout) :: interrupt
    procedure(),                           pointer, intent(inout) :: interruptProcedure
    class(nodeComponentDarkMatterProfile), pointer                :: thisDarkMatterProfileComponent

    ! Get the dark matter profile component.
    thisDarkMatterProfileComponent => thisNode%darkMatterProfile()
    ! Ensure that it is of the scale+shape class.
    select type (thisDarkMatterProfileComponent)
       class is (nodeComponentDarkMatterProfileScaleShape)
       call thisDarkMatterProfileComponent%shapeRate(thisDarkMatterProfileComponent%shapeGrowthRate())
    end select
    return
  end subroutine Node_Component_Dark_Matter_Profile_Scale_Shape_Rate_Compute

  subroutine Node_Component_Dark_Matter_Profile_Scale_Shape_Initialize_Shape(thisNode)
    !% Initialize the shape parameter of {\tt thisNode}.
    use Dark_Matter_Profiles_Shapes
    implicit none
    type (treeNode                      ), pointer, intent(inout) :: thisNode
    class(nodeComponentDarkMatterProfile), pointer                :: thisDarkMatterProfileComponent

    ! Ensure that the module is initialized.
    call Node_Component_Dark_Matter_Profile_Scale_Shape_Initialize()
    ! Get the dark matter profile component.
    thisDarkMatterProfileComponent => thisNode%darkMatterProfile(autoCreate=.true.)
    ! Set the shape parameter of the halo.
    call thisDarkMatterProfileComponent%shapeSet(Dark_Matter_Profile_Shape(thisNode))
    return
  end subroutine Node_Component_Dark_Matter_Profile_Scale_Shape_Initialize_Shape

  !# <mergerTreeInitializeTask>
  !#  <unitName>Node_Component_Dark_Matter_Profile_Scale_Shape_Tree_Initialize</unitName>
  !#  <sortName>darkMatterProfile</sortName>
  !# </mergerTreeInitializeTask>
  subroutine Node_Component_Dark_Matter_Profile_Scale_Shape_Tree_Initialize(thisNode)
    !% Initialize the scale radius of {\tt thisNode}.
    implicit none
    type (treeNode                      ), pointer, intent(inout) :: thisNode
    class(nodeComponentDarkMatterProfile), pointer                :: thisDarkMatterProfileComponent,parentDarkMatterProfileComponent
    class(nodeComponentBasic            ), pointer                :: thisBasicComponent,parentBasicComponent
    double precision                                              :: deltaTime

    ! Get the dark matter profile component.
    thisDarkMatterProfileComponent => thisNode%darkMatterProfile()
    if (defaultDarkMatterProfileComponent%scaleShapeIsActive()) then
       ! Ensure that current node has its shape set.
       call Node_Component_Dark_Matter_Profile_Scale_Shape_Initialize_Shape(thisNode)
       ! Check if this node is the primary progenitor.
       if (thisNode%isPrimaryProgenitor()) then
          ! It is, so compute the shape parameter growth rate.
          ! First ensure that parent node has scale radius set.
          call Node_Component_Dark_Matter_Profile_Scale_Shape_Initialize_Shape(thisNode%parent)
          ! Now compute the growth rate.
          thisBasicComponent   => thisNode       %basic()
          parentBasicComponent => thisNode%parent%basic()
          deltaTime=parentBasicComponent%time()-thisBasicComponent%time()
          if (deltaTime > 0.0d0) then
             parentDarkMatterProfileComponent => thisNode%parent%darkMatterProfile()
             call thisDarkMatterProfileComponent%shapeGrowthRateSet(                                           &
                  &                                                 (                                          &
                  &                                                   parentDarkMatterProfileComponent%shape() &
                  &                                                  -thisDarkMatterProfileComponent  %shape() &
                  &                                                 )                                          &
                  &                                                 /deltaTime                                 &
                  &                                                )
          else
             call thisDarkMatterProfileComponent%shapeGrowthRateSet(0.0d0)
          end if
       else
          ! It is not, so set shape parameter growth rates to zero.
          call    thisDarkMatterProfileComponent%shapeGrowthRateSet(0.0d0)
       end if
    end if
    return
  end subroutine Node_Component_Dark_Matter_Profile_Scale_Shape_Tree_Initialize

  !# <nodePromotionTask>
  !#  <unitName>Node_Component_Dark_Matter_Profile_Scale_Shape_Promote</unitName>
  !# </nodePromotionTask>
  subroutine Node_Component_Dark_Matter_Profile_Scale_Shape_Promote(thisNode)
    !% Ensure that {\tt thisNode} is ready for promotion to its parent. In this case, we simply update the growth rate of {\tt thisNode}
    !% to be that of its parent.
    use Galacticus_Error
    implicit none
    type (treeNode                      ), pointer, intent(inout) :: thisNode
    class(nodeComponentDarkMatterProfile), pointer                :: thisDarkMatterProfileComponent,parentDarkMatterProfileComponent
    class(nodeComponentBasic            ), pointer                :: thisBasicComponent,parentBasicComponent

    ! Get the dark matter profile component.
    thisDarkMatterProfileComponent => thisNode%darkMatterProfile()
    ! Ensure it is of the scale+shape class.
    select type (thisDarkMatterProfileComponent)
       class is (nodeComponentDarkMatterProfileScaleShape)
       parentDarkMatterProfileComponent => thisNode%parent%darkMatterProfile()
       thisBasicComponent               => thisNode       %basic            ()
       parentBasicComponent             => thisNode%parent%basic            ()
       if (thisBasicComponent%time() /= parentBasicComponent%time()) call Galacticus_Error_Report('Node_Component_Dark_Matter_Profile_Scale_Shape_Promote','thisNode&
            & has not been evolved to its parent')
       ! Adjust the shape parameter to that of the parent node.
       call thisDarkMatterProfileComponent%shapeSet          (parentDarkMatterProfileComponent%shape          ())
       ! Adjust the growth rate to that of the parent node.
       call thisDarkMatterProfileComponent%shapeGrowthRateSet(parentDarkMatterProfileComponent%shapeGrowthRate())
    end select
    return
  end subroutine Node_Component_Dark_Matter_Profile_Scale_Shape_Promote

  !# <scaleSetTask>
  !#  <unitName>Node_Component_Dark_Matter_Profile_Scale_Shape_Scale_Set</unitName>
  !# </scaleSetTask>
  subroutine Node_Component_Dark_Matter_Profile_Scale_Shape_Scale_Set(thisNode)
    !% Set scales for properties of {\tt thisNode}.
    implicit none
    type (treeNode                      ), pointer, intent(inout) :: thisNode
    class(nodeComponentDarkMatterProfile), pointer                :: thisDarkMatterProfileComponent

    ! Get the dark matter profile component.
    thisDarkMatterProfileComponent => thisNode%darkMatterProfile()
    ! Ensure it is of the scale+shape class.
    select type (thisDarkMatterProfileComponent)
       class is (nodeComponentDarkMatterProfileScale)        
          ! Set scale for the scale radius.
       call thisDarkMatterProfileComponent%shapeScale(thisDarkMatterProfileComponent%shape())
    end select
    return
  end subroutine Node_Component_Dark_Matter_Profile_Scale_Shape_Scale_Set

  !# <mergerTreeStructureOutputTask>
  !#  <unitName>Node_Component_Dark_Matter_Profile_Scale_Shape_Tree_Output</unitName>
  !# </mergerTreeStructureOutputTask>
  subroutine Node_Component_Dark_Matter_Profile_Scale_Shape_Tree_Output(baseNode,nodeProperty,treeGroup)
    !% Write the scale radius property to a full merger tree output.
    use IO_HDF5
    use Numerical_Constants_Astronomical
    implicit none
    type (treeNode                      ), intent(in   ), pointer      :: baseNode
    double precision                     , intent(inout), dimension(:) :: nodeProperty
    type (hdf5Object                    ), intent(inout)               :: treeGroup
    type (treeNode                      ),                pointer      :: thisNode
    integer                                                            :: nodeCount
    class(nodeComponentDarkMatterProfile),                pointer      :: baseDarkMatterProfileComponent,thisDarkMatterProfileComponent
    type (hdf5Object                    )                              :: nodeDataset

    ! Check if scale radius is to be included in merger tree outputs.
    if (mergerTreeStructureOutputDarkMatterProfileShape) then
       ! Get the dark matter profile component.
       baseDarkMatterProfileComponent => baseNode%darkMatterProfile()
       ! Ensure it is of the scale+shape class.
       select type (baseDarkMatterProfileComponent)
          class is (nodeComponentDarkMatterProfileScaleShape)                  
             ! Extract node shape parameter and output to file.
          nodeCount=0
          thisNode => baseNode
          do while (associated(thisNode))
             thisDarkMatterProfileComponent => thisNode%darkMatterProfile()
             nodeCount=nodeCount+1
             nodeProperty(nodeCount)=thisDarkMatterProfileComponent%shape()
             call thisNode%walkTree(thisNode)
          end do
          call treeGroup%writeDataset(nodeProperty,'darkMatterShapeParameter','Shape parameter of the dark matter profile.')
       end select
    end if

    return
  end subroutine Node_Component_Dark_Matter_Profile_Scale_Shape_Tree_Output

end module Node_Component_Dark_Matter_Profile_Scale_Shape
