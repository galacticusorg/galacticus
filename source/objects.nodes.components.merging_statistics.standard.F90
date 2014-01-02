!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements the standard merging statistics component.

module Node_Component_Merging_Statistics_Standard
  !% Implements the standard merging statistics component.
  use Galacticus_Nodes
  implicit none
  private
  public :: Node_Component_Merging_Statistics_Standard_Merger_Tree_Init, Node_Component_Merging_Statistics_Standard_Node_Merger, &
       &    Node_Component_Merging_Statistics_Standard_Node_Promotion  , Node_Component_Merging_Statistics_Standard_Satellite_Merger

  !# <component>
  !#  <class>mergingStatistics</class>
  !#  <name>standard</name>
  !#  <isDefault>no</isDefault>
  !#  <properties>
  !#   <property>
  !#     <name>galaxyMajorMergerTime</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#     <classDefault>-1.0d0</classDefault>
  !#     <output unitsInSI="gigaYear" comment="Time of the last major merger."/>
  !#   </property>
  !#   <property>
  !#     <name>nodeMajorMergerTime</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#     <classDefault>-1.0d0</classDefault>
  !#     <output unitsInSI="gigaYear" comment="Time of the last node major merger."/>
  !#   </property>
  !#   <property>
  !#     <name>nodeFormationTime</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#     <classDefault>-1.0d0</classDefault>
  !#     <output unitsInSI="gigaYear" comment="Formation time of the node."/>
  !#   </property>
  !#   <property>
  !#     <name>nodeHierarchyLevel</name>
  !#     <type>integer</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#     <output unitsInSI="0.0d0" comment="Initial level of the node in the tree hierarchy."/>
  !#   </property>
  !#  </properties>
  !# </component>

  ! Parameters controlling the statistics gathered.
  double precision :: nodeFormationMassFraction        , nodeMajorMergerFraction

  ! Record of whether this module has been initialized.
  logical          :: moduleInitialized        =.false.

contains

  subroutine Node_Component_Merging_Statistics_Standard_Initialize()
    !% Initializes the standard merging statistics component.
    use Input_Parameters
    implicit none

    ! Test whether module is already initialize.
    !$omp critical (Node_Component_Merging_Statistics_Standard_Initialize)
    if (.not.moduleInitialized) then
       !@ <inputParameter>
       !@   <name>nodeMajorMergerFraction</name>
       !@   <defaultValue>0.25</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The mass ratio ($M_2/M_1$ where $M_2 &lt; M_1$) of merging halos above which the merger should be considered to be ``major''.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('nodeMajorMergerFraction',nodeMajorMergerFraction,defaultValue=0.25d0)
       !@ <inputParameter>
       !@   <name>nodeFormationMassFraction</name>
       !@   <defaultValue>0.5</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The mass fraction in the main branch progenitor used to define the formation time of each halo.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('nodeFormationMassFraction',nodeFormationMassFraction,defaultValue=0.5d0)
       ! Record that the module is now initialized.
       moduleInitialized=.true.
    end if
    !$omp end critical (Node_Component_Merging_Statistics_Standard_Initialize)
    return
  end subroutine Node_Component_Merging_Statistics_Standard_Initialize

  !# <mergerTreeInitializeTask>
  !#  <unitName>Node_Component_Merging_Statistics_Standard_Merger_Tree_Init</unitName>
  !# </mergerTreeInitializeTask>
  subroutine Node_Component_Merging_Statistics_Standard_Merger_Tree_Init(thisNode)
    !% Initialize the merging statistics component by creating components in nodes and computing formation times.
    implicit none
    type   (treeNode                      ), intent(inout), pointer :: thisNode
    type   (treeNode                      )               , pointer :: descendentNode
    class  (nodeComponentMergingStatistics)               , pointer :: thisMergingStatisticsComponent,descendentMergingStatistics
    integer                                                         :: hierarchyLevel

    ! Return immediately if this class is not active.
    if (.not.defaultMergingStatisticsComponent%standardIsActive()) return

    ! Ensure that the module is initialized.
    call Node_Component_Merging_Statistics_Standard_Initialize()

    ! Create a merger statistics component and initialize it.
    thisMergingStatisticsComponent => thisNode%mergingStatistics()
    select type (thisMergingStatisticsComponent)
    class is (nodeComponentMergingStatisticsStandard)
       ! Node has been initialized already.
       return
    end select
    thisMergingStatisticsComponent => thisNode%mergingStatistics(autoCreate=.true.)
    select type (thisMergingStatisticsComponent)
    class is (nodeComponentMergingStatisticsStandard)
       hierarchyLevel=0
       descendentNode => thisNode
       do while (associated(descendentNode))
          if (associated(descendentNode%parent).and..not.descendentNode%isPrimaryProgenitor()) hierarchyLevel=hierarchyLevel+1
          if (.not.associated(thisNode,descendentNode)) then
             descendentMergingStatistics => descendentNode%mergingStatistics()
             select type (descendentMergingStatistics)
             class is (nodeComponentMergingStatisticsStandard)
                hierarchyLevel=hierarchyLevel+descendentMergingStatistics%nodeHierarchyLevel()
                exit
             end select
          end if
          descendentNode => descendentNode%parent
       end do
       descendentNode => thisNode
       do while (associated(descendentNode))
          call Node_Component_Merging_Statistics_Standard_Merger_Tree_Init_Set(descendentNode,hierarchyLevel)
          if (descendentNode%isPrimaryProgenitor()) then
             descendentNode => descendentNode%parent
          else
             descendentNode => null()
          end if
       end do
    end select
    return
  end subroutine Node_Component_Merging_Statistics_Standard_Merger_Tree_Init

  subroutine Node_Component_Merging_Statistics_Standard_Merger_Tree_Init_Set(thisNode,hierarchyLevel)
    !% Set the initial properties of the standard merging statistics component in a node.
    use Dark_Matter_Halo_Formation_Times
    implicit none
    type   (treeNode                      ), intent(inout), pointer :: thisNode
    integer                                , intent(in   )          :: hierarchyLevel
    class  (nodeComponentMergingStatistics),                pointer :: thisMergingStatistics

    thisMergingStatistics => thisNode%mergingStatistics(autoCreate=.true.)
    call thisMergingStatistics%   nodeHierarchyLevelSet(hierarchyLevel)
    call thisMergingStatistics%galaxyMajorMergerTimeSet(        -1.0d0)
    call thisMergingStatistics%  nodeMajorMergerTimeSet(        -1.0d0)
    call thisMergingStatistics%    nodeFormationTimeSet(Dark_Matter_Halo_Formation_Time(thisNode,nodeFormationMassFraction))
    return
  end subroutine Node_Component_Merging_Statistics_Standard_Merger_Tree_Init_Set

  !# <nodeMergerTask>
  !#  <unitName>Node_Component_Merging_Statistics_Standard_Node_Merger</unitName>
  !# </nodeMergerTask>
  subroutine Node_Component_Merging_Statistics_Standard_Node_Merger(thisNode)
    !% Record any major merger of {\tt thisNode}.
    implicit none
    type (treeNode                      ), intent(inout), pointer :: thisNode
    class(nodeComponentMergingStatistics)               , pointer :: parentMergingStatisticsComponent
    class(nodeComponentBasic            )               , pointer :: parentBasicComponent            , thisBasicComponent

    ! Return immediately if this class is not active.
    if (.not.defaultMergingStatisticsComponent%standardIsActive()) return

    thisBasicComponent               => thisNode       %basic            ()
    parentBasicComponent             => thisNode%parent%basic            ()
    parentMergingStatisticsComponent => thisNode%parent%mergingStatistics()
    ! Record the merger time if this is a major merger.
    if (thisBasicComponent%mass() >= nodeMajorMergerFraction*parentBasicComponent%mass()) &
            &  call parentMergingStatisticsComponent%nodeMajorMergerTimeSet(thisBasicComponent%time())
    return
  end subroutine Node_Component_Merging_Statistics_Standard_Node_Merger

  !# <nodePromotionTask>
  !#  <unitName>Node_Component_Merging_Statistics_Standard_Node_Promotion</unitName>
  !# </nodePromotionTask>
  subroutine Node_Component_Merging_Statistics_Standard_Node_Promotion(thisNode)
    !% Ensure that {\tt thisNode} is ready for promotion to its parent. In this case, we simply update the node merger time.
    implicit none
    type (treeNode                      ), intent(inout), pointer :: thisNode
    class(nodeComponentMergingStatistics)               , pointer :: parentMergingStatisticsComponent, thisMergingStatisticsComponent

    ! Return immediately if this class is not active.
    if (.not.defaultMergingStatisticsComponent%standardIsActive()) return

    ! Get the merging statistics components.
    parentMergingStatisticsComponent => thisNode%parent%mergingStatistics()
    thisMergingStatisticsComponent   => thisNode       %mergingStatistics()
    if (parentMergingStatisticsComponent%nodeMajorMergerTime() > thisMergingStatisticsComponent%nodeMajorMergerTime()) &
         & call thisMergingStatisticsComponent%nodeMajorMergerTimeSet(parentMergingStatisticsComponent%nodeMajorMergerTime())
    call thisMergingStatisticsComponent%nodeFormationTimeSet(parentMergingStatisticsComponent%nodeFormationTime())
    return
  end subroutine Node_Component_Merging_Statistics_Standard_Node_Promotion

  !# <satelliteMergerTask>
  !#  <unitName>Node_Component_Merging_Statistics_Standard_Satellite_Merger</unitName>
  !#  <after>Satellite_Merging_Mass_Movement_Store</after>
  !# </satelliteMergerTask>
  subroutine Node_Component_Merging_Statistics_Standard_Satellite_Merger(thisNode)
    !% Record properties of a merging event for {\tt thisNode}.
    use Satellite_Merging_Mass_Movements_Descriptors
    implicit none
    type (treeNode                      ), intent(inout), pointer :: thisNode
    type (treeNode                      )               , pointer :: hostNode
    class(nodeComponentBasic            )               , pointer :: hostBasicComponent
    class(nodeComponentMergingStatistics)               , pointer :: hostMergingStatisticsComponent

    ! Return immediately if this class is not active.
    if (.not.defaultMergingStatisticsComponent%standardIsActive()) return

    ! Record the time of this merger if it is a major merger.
    if (thisMergerIsMajor) then
       ! Find the node to merge with.
       hostNode                       => thisNode%mergesWith       ()
       hostBasicComponent             => hostNode%basic            ()
       hostMergingStatisticsComponent => hostNode%mergingStatistics()
       ! Record the merger time.
       call hostMergingStatisticsComponent%galaxyMajorMergerTimeSet(hostBasicComponent%time())
    end if
    return
  end subroutine Node_Component_Merging_Statistics_Standard_Satellite_Merger

end module Node_Component_Merging_Statistics_Standard
