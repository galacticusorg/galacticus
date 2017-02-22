!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
!!    Andrew Benson <abenson@carnegiescience.edu>
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
  public :: Node_Component_Merging_Statistics_Standard_Merger_Tree_Init, Node_Component_Merging_Statistics_Standard_Node_Merger      , &
       &    Node_Component_Merging_Statistics_Standard_Node_Promotion  , Node_Component_Merging_Statistics_Standard_Satellite_Merging, &
       &    Node_Component_Merging_Statistics_Standard_Reset_Hierarchy , Node_Component_Merging_Statistics_Standard_Initialize

  !# <component>
  !#  <class>mergingStatistics</class>
  !#  <name>standard</name>
  !#  <isDefault>false</isDefault>
  !#  <properties>
  !#   <property>
  !#     <name>galaxyMajorMergerTime</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#     <classDefault>-1.0d0</classDefault>
  !#     <output unitsInSI="gigaYear" comment="Time of the last major merger."/>
  !#   </property>
  !#   <property>
  !#     <name>nodeMajorMergerTime</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#     <classDefault>-1.0d0</classDefault>
  !#     <output unitsInSI="gigaYear" comment="Time of the last node major merger."/>
  !#   </property>
  !#   <property>
  !#     <name>nodeFormationTime</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#     <classDefault>-1.0d0</classDefault>
  !#     <output unitsInSI="gigaYear" comment="Formation time of the node."/>
  !#   </property>
  !#   <property>
  !#     <name>nodeHierarchyLevel</name>
  !#     <type>integer</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" isDeferred="get" />
  !#     <output unitsInSI="0.0d0" comment="Initial level of the node in the tree hierarchy."/>
  !#     <classDefault>-1.0d0</classDefault>
  !#   </property>
  !#   <property>
  !#     <name>nodeHierarchyLevelMaximum</name>
  !#     <type>integer</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" isDeferred="get" />
  !#     <output unitsInSI="0.0d0" comment="Maximum level of the node in the tree hierarchy."/>
  !#     <classDefault>-1.0d0</classDefault>
  !#   </property>
  !#   <property>
  !#     <name>massWhenFirstIsolated</name>
  !#     <type>double</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#   </property>
  !#  </properties>
  !# </component>

  ! Parameters controlling the statistics gathered.
  double precision                                         :: nodeFormationMassFraction        , nodeMajorMergerFraction, &
       &                                                      hierarchyLevelResetFactor

  ! Record of whether this module has been initialized.
  logical                                                  :: moduleInitialized        =.false.

  ! Queriable merging statistics object.
  type            (nodeComponentMergingStatisticsStandard) :: mergingStatistics

contains

  !# <nodeComponentInitializationTask>
  !#  <unitName>Node_Component_Merging_Statistics_Standard_Initialize</unitName>
  !# </nodeComponentInitializationTask>
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
       !@   <type>double</type>
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
       !@   <type>double</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('nodeFormationMassFraction',nodeFormationMassFraction,defaultValue=0.5d0)
       !@ <inputParameter>
       !@   <name>hierarchyLevelResetFactor</name>
       !@   <defaultValue>$10^{30}$</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The factor by which a node's mass must increase before the previous maximum hierarchy level is forgotten.
       !@   </description>
       !@   <type>double</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('hierarchyLevelResetFactor',hierarchyLevelResetFactor,defaultValue=huge(1.0d0))
       ! Bind the hierarchy level get functions.
       call mergingStatistics%nodeHierarchyLevelFunction       (Node_Component_Merging_statistics_Standard_Hierarchy_Level)
       call mergingStatistics%nodeHierarchyLevelMaximumFunction(Node_Component_Merging_statistics_Standard_HLM            )
       ! Record that the module is now initialized.
       moduleInitialized=.true.
    end if
    !$omp end critical (Node_Component_Merging_Statistics_Standard_Initialize)
    return
  end subroutine Node_Component_Merging_Statistics_Standard_Initialize

  integer function Node_Component_Merging_Statistics_Standard_Hierarchy_Level(self)
    !% Return the hierarchy level of a node, computing it if necessary.
    implicit none
    class  (nodeComponentMergingStatisticsStandard), intent(inout) :: self
    type   (treeNode                              ), pointer       :: hostNode
    integer                                                        :: nodeHierarchyLevel

    if (self%nodeHierarchyLevelValue() < 0.0d0) then
       nodeHierarchyLevel =  0
       hostNode       => self%hostNode
       do while (hostNode%isSatellite())
          nodeHierarchyLevel =  nodeHierarchyLevel       +1
          hostNode           => hostNode          %parent
       end do
       call self%nodeHierarchyLevelSet(nodeHierarchyLevel)
    end if
    Node_Component_Merging_statistics_Standard_Hierarchy_Level=self%nodeHierarchyLevelValue()
    return
  end function Node_Component_Merging_Statistics_Standard_Hierarchy_Level

  integer function Node_Component_Merging_Statistics_Standard_HLM(self)
    !% Return the maximum hierarchy level of a node, computing it if necessary.
    implicit none
    class  (nodeComponentMergingStatisticsStandard), intent(inout) :: self
    class  (nodeComponentMergingStatistics        ), pointer       :: mergingStatisticsProgenitor
    class  (nodeEvent                             ), pointer       :: event
    type   (treeNode                              ), pointer       :: nodeProgenitor
    class  (nodeComponentBasic                    ), pointer       :: basicSelf                  , basicProgenitor
    integer                                                        :: nodeHierarchyLevel         , nodeHierarchyLevelProgenitor

    if (self%nodeHierarchyLevelMaximumValue() < 0.0d0) then
       ! Get current hierarchy level.
       nodeHierarchyLevel=self%nodeHierarchyLevel()
       ! Check for a subhalo promotion event to this node.
       nodeHierarchyLevelProgenitor =  -1
       event                        => self%hostNode%event
       do while (associated(event))
          ! Select only subhalo promotion events.
          select type (event)
          class is (nodeEventSubhaloPromotion)           
             ! If the event task is not associated, this indicates that our node is the node being promoted to (not the node
             ! being promoted). Select these cases, as in these instances we want to propagate the maximum hierarchy level of
             ! the node being promoted.
             if (.not.associated(event%task)) then
                ! We have a subhalo promotion. Set the progenitor maximum hierarchy level to the greater of the progenitor's
                ! maximum hierarchy level and 1 (since we know it is a subhalo prior to promotion to our node).
                mergingStatisticsProgenitor  => event   %node%mergingStatistics(autoCreate=.true.)
                nodeHierarchyLevelProgenitor =  max(mergingStatisticsProgenitor%nodeHierarchyLevelMaximum(),1)
             end if
          end select
          event => event%next
       end do
       ! If no subhalo promotion event found, check for a direct, primary progenitor.
       if (nodeHierarchyLevelProgenitor < 0) then
          if (associated(self%hostNode%firstChild)) then
             ! Found a direct, primary progenitor. Our maximum hierarchy level is just equal to that of our progenitor.
             mergingStatisticsProgenitor  => self                       %hostNode%firstChild%mergingStatistics        (autoCreate=.true.)
             nodeHierarchyLevelProgenitor =  mergingStatisticsProgenitor                    %nodeHierarchyLevelMaximum(                 )
          else
             ! No progenitor found, set maximum hierarchy level to zero.
             nodeHierarchyLevelProgenitor =  0
          end if
       end if
       ! Check for cases where a node has grown sufficiently in mass that its maximum hierarchy level should be reset.
       nodeProgenitor => self%hostNode
       do while (associated(nodeProgenitor%firstChild))
          nodeProgenitor => nodeProgenitor%firstChild
       end do
       basicProgenitor => nodeProgenitor         %basic()
       basicSelf       => self          %hostNode%basic()
       if (basicSelf%mass() > hierarchyLevelResetFactor*basicProgenitor%mass()) nodeHierarchyLevelProgenitor=nodeHierarchyLevel
       ! Assign the computed maximum hierarchy level.
       call self%nodeHierarchyLevelMaximumSet(max(nodeHierarchyLevel,nodeHierarchyLevelProgenitor))       
    end if
    Node_Component_Merging_statistics_Standard_HLM=self%nodeHierarchyLevelMaximumValue()
    return
  end function Node_Component_Merging_Statistics_Standard_HLM

  !# <mergerTreeInitializeTask>
  !#  <unitName>Node_Component_Merging_Statistics_Standard_Merger_Tree_Init</unitName>
  !# </mergerTreeInitializeTask>
  subroutine Node_Component_Merging_Statistics_Standard_Merger_Tree_Init(node)
    !% Initialize the merging statistics component by creating components in nodes and computing formation times.
    use Dark_Matter_Halo_Formation_Times
    implicit none
    type   (treeNode                      ), intent(inout), pointer :: node
    type   (treeNode                      )               , pointer :: nodeHost
    class  (nodeComponentMergingStatistics)               , pointer :: mergingStatistics
    class  (nodeComponentBasic            )               , pointer :: basic
    integer                                                         :: nodeHierarchyLevel

    ! Return immediately if this class is not active.
    if (.not.defaultMergingStatisticsComponent%standardIsActive()) return

    ! Ensure that the module is initialized.
    call Node_Component_Merging_Statistics_Standard_Initialize()
    ! Find the initial hierarchy level.
    nodeHierarchyLevel =  0
    nodeHost       => node
    do while (nodeHost%isSatellite())
       nodeHierarchyLevel =  nodeHierarchyLevel       +1
       nodeHost           => nodeHost          %parent
    end do    
    ! Create a merger statistics component and initialize it. 
    mergingStatistics => node%mergingStatistics(autoCreate=.true.)
    basic             => node%basic            (                 )
    call mergingStatistics%       nodeHierarchyLevelSet(nodeHierarchyLevel)
    call mergingStatistics%nodeHierarchyLevelMaximumSet(nodeHierarchyLevel)
    call mergingStatistics%    massWhenFirstIsolatedSet(      basic%mass())
    call mergingStatistics%    galaxyMajorMergerTimeSet(            -1.0d0)
    call mergingStatistics%      nodeMajorMergerTimeSet(            -1.0d0)
    call mergingStatistics%        nodeFormationTimeSet(Dark_Matter_Halo_Formation_Time(node,nodeFormationMassFraction))
    return
  end subroutine Node_Component_Merging_Statistics_Standard_Merger_Tree_Init

  !# <nodeMergerTask>
  !#  <unitName>Node_Component_Merging_Statistics_Standard_Node_Merger</unitName>
  !# </nodeMergerTask>
  subroutine Node_Component_Merging_Statistics_Standard_Node_Merger(node)
    !% Record any major merger of {\normalfont \ttfamily node}.
    implicit none
    type (treeNode                      ), intent(inout), pointer :: node
    class(nodeComponentMergingStatistics)               , pointer :: mergingStatisticsParent, mergingStatistics
    class(nodeComponentBasic            )               , pointer :: parentBasicComponent            , basic
    
    ! Return immediately if this class is not active.
    if (.not.defaultMergingStatisticsComponent%standardIsActive()) return

    basic                            => node       %basic            ()
    parentBasicComponent             => node%parent%basic            ()
    mergingStatistics       => node       %mergingStatistics()
    mergingStatisticsParent => node%parent%mergingStatistics()
    ! Record the merger time if this is a major merger.
    if (basic%mass() >= nodeMajorMergerFraction*parentBasicComponent%mass()) &
         &  call mergingStatisticsParent%nodeMajorMergerTimeSet(basic%time())
    ! Increment the hierarchy level of the merging node.
    call Node_Component_Merging_Statistics_Standard_Reset_Hierarchy(node)
    call mergingStatistics%nodeHierarchyLevelSet       (    mergingStatistics%nodeHierarchyLevel       ()+1)
    call mergingStatistics%nodeHierarchyLevelMaximumSet(                                                                   &
         &                                                           max(                                                               &
         &                                                               mergingStatistics%nodeHierarchyLevel       ()   , &
         &                                                               mergingStatistics%nodeHierarchyLevelMaximum()     &
         &                                                              )                                                               &
         &                                                          )
    return
  end subroutine Node_Component_Merging_Statistics_Standard_Node_Merger

  !# <nodePromotionTask>
  !#  <unitName>Node_Component_Merging_Statistics_Standard_Node_Promotion</unitName>
  !# </nodePromotionTask>
  subroutine Node_Component_Merging_Statistics_Standard_Node_Promotion(node)
    !% Ensure that {\normalfont \ttfamily node} is ready for promotion to its parent. In this case, we simply update the node merger time.
    implicit none
    type (treeNode                      ), intent(inout), pointer :: node
    class(nodeComponentMergingStatistics)               , pointer :: mergingStatisticsParent, mergingStatistics
    class(nodeComponentBasic            )               , pointer :: basicParent
    
    ! Return immediately if this class is not active.
    if (.not.defaultMergingStatisticsComponent%standardIsActive()) return

    ! Get the merging statistics components.
    mergingStatisticsParent => node%parent%mergingStatistics()
    mergingStatistics       => node       %mergingStatistics()
    basicParent             => node%parent%basic            ()
    call Node_Component_Merging_Statistics_Standard_Reset_Hierarchy(node)
    if (mergingStatisticsParent%nodeMajorMergerTime() > mergingStatistics%nodeMajorMergerTime()) &
         & call mergingStatistics%nodeMajorMergerTimeSet(mergingStatisticsParent%nodeMajorMergerTime())
    call mergingStatistics%nodeHierarchyLevelSet(mergingStatisticsParent%nodeHierarchyLevel())
    call mergingStatistics%nodeFormationTimeSet (mergingStatisticsParent%nodeFormationTime ())
    return
  end subroutine Node_Component_Merging_Statistics_Standard_Node_Promotion

  !# <postEvolveTask>
  !# <unitName>Node_Component_Merging_Statistics_Standard_Reset_Hierarchy</unitName>
  !# </postEvolveTask>
  subroutine Node_Component_Merging_Statistics_Standard_Reset_Hierarchy(node)
    !% Reset the maximum node hierarchy level if the node has grown sufficiently in mass.
    implicit none
    type            (treeNode                      ), intent(inout), pointer :: node
    class           (nodeComponentBasic            )               , pointer :: basic
    class           (nodeComponentMergingStatistics)               , pointer :: mergingStatistics

    ! Return immediately if this class is not active.
    if (.not.defaultMergingStatisticsComponent%standardIsActive()) return
    ! Test if previous mass at promotion has been exceeded by a sufficient factor.
    basic    => node%basic            ()
    mergingStatistics => node%mergingStatistics()
    if (basic%mass() > hierarchyLevelResetFactor*mergingStatistics%massWhenFirstIsolated()) &
         & call mergingStatistics%nodeHierarchyLevelMaximumSet(mergingStatistics%nodeHierarchyLevel())
    return
  end subroutine Node_Component_Merging_Statistics_Standard_Reset_Hierarchy
    
  !# <satelliteMergerTask>
  !#  <unitName>Node_Component_Merging_Statistics_Standard_Satellite_Merging</unitName>
  !#  <after>Satellite_Merging_Mass_Movement_Store</after>
  !# </satelliteMergerTask>
  subroutine Node_Component_Merging_Statistics_Standard_Satellite_Merging(node)
    !% Record properties of a merging event for {\normalfont \ttfamily node}.
    use Satellite_Merging_Mass_Movements_Descriptors
    implicit none
    type (treeNode                      ), intent(inout), pointer :: node
    type (treeNode                      )               , pointer :: nodeHost
    class(nodeComponentBasic            )               , pointer :: basicHost
    class(nodeComponentMergingStatistics)               , pointer :: mergingStatisticsHost

    ! Return immediately if this class is not active.
    if (.not.defaultMergingStatisticsComponent%standardIsActive()) return

    ! Record the time of this merger if it is a major merger.
    if (thisMergerIsMajor) then
       ! Find the node to merge with.
       nodeHost              => node    %mergesWith       ()
       basicHost             => nodeHost%basic            ()
       mergingStatisticsHost => nodeHost%mergingStatistics()
       ! Record the merger time.
       call mergingStatisticsHost%galaxyMajorMergerTimeSet(basicHost%time())
    end if
    return
  end subroutine Node_Component_Merging_Statistics_Standard_Satellite_Merging

end module Node_Component_Merging_Statistics_Standard
