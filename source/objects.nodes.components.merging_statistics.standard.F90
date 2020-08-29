!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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
  use :: Dark_Matter_Halo_Mass_Accretion_Histories, only : darkMatterHaloMassAccretionHistoryClass
  use :: Galacticus_Nodes                         , only : nodeComponentMergingStatisticsStandard
  use :: Satellite_Merging_Mass_Movements         , only : mergerMassMovementsClass
  implicit none
  private
  public :: Node_Component_Merging_Statistics_Standard_Merger_Tree_Init     , Node_Component_Merging_Statistics_Standard_Node_Merger        , &
       &    Node_Component_Merging_Statistics_Standard_Reset_Hierarchy      , Node_Component_Merging_Statistics_Standard_Initialize         , &
       &    Node_Component_Merging_Statistics_Standard_Thread_Initialize    , Node_Component_Merging_Statistics_Standard_Thread_Uninitialize, &
       &    Node_Component_Merging_Statistics_Standard_Host_Change

  !# <component>
  !#  <class>mergingStatistics</class>
  !#  <name>standard</name>
  !#  <isDefault>false</isDefault>
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
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" isDeferred="get" />
  !#     <output unitsInSI="0.0d0" comment="The current level of the node in the tree hierarchy (0 for a (non-sub-)halo; 1 for a sub-halo; 2 for a sub-sub-halo; etc.)."/>
  !#     <classDefault>-1.0d0</classDefault>
  !#   </property>
  !#   <property>
  !#     <name>nodeHierarchyLevelMaximum</name>
  !#     <type>integer</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" isDeferred="get" />
  !#     <output unitsInSI="0.0d0" comment="Maximum level that the node ever reached in the tree hierarchy (0 for a (non-sub-)halo; 1 for a sub-halo; 2 for a sub-sub-halo; etc.)."/>
  !#     <classDefault>-1.0d0</classDefault>
  !#   </property>
  !#   <property>
  !#     <name>nodeHierarchyLevelDepth</name>
  !#     <type>integer</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#     <output unitsInSI="0.0d0" comment="Maximum level of the node in the tree hierarchy that could possibly ever be reached (0 for a (non-sub-)halo; 1 for a sub-halo; 2 for a sub-sub-halo; etc.)."/>
  !#     <classDefault>-1.0d0</classDefault>
  !#   </property>
  !#   <property>
  !#     <name>massWhenFirstIsolated</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#   </property>
  !#  </properties>
  !# </component>

  ! Classes used.
  class          (darkMatterHaloMassAccretionHistoryClass), pointer :: darkMatterHaloMassAccretionHistory_
  class          (mergerMassMovementsClass               ), pointer :: mergerMassMovements_
  !$omp threadprivate(darkMatterHaloMassAccretionHistory_,mergerMassMovements_)

  ! Parameters controlling the statistics gathered.
  double precision                                                  :: nodeFormationMassFraction          , nodeMajorMergerFraction, &
       &                                                               hierarchyLevelResetFactor

  ! Queriable merging statistics object.
  type            (nodeComponentMergingStatisticsStandard)          :: mergingStatistics

contains

  !# <nodeComponentInitializationTask>
  !#  <unitName>Node_Component_Merging_Statistics_Standard_Initialize</unitName>
  !# </nodeComponentInitializationTask>
  subroutine Node_Component_Merging_Statistics_Standard_Initialize(parameters_)
    !% Initializes the standard merging statistics component.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters_

    !# <inputParameter>
    !#   <name>nodeMajorMergerFraction</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>0.25d0</defaultValue>
    !#   <description>The mass ratio ($M_2/M_1$ where $M_2 &lt; M_1$) of merging halos above which the merger should be considered to be ``major''.</description>
    !#   <source>parameters_</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>nodeFormationMassFraction</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>0.5d0</defaultValue>
    !#   <description>The mass fraction in the main branch progenitor used to define the formation time of each halo.</description>
    !#   <source>parameters_</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>hierarchyLevelResetFactor</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>1.0d100</defaultValue>
    !#   <description>The factor by which a node's mass must increase before the previous maximum hierarchy level is forgotten.</description>
    !#   <source>parameters_</source>
    !#   <type>real</type>
    !# </inputParameter>
    ! Bind the hierarchy level get functions.
    call mergingStatistics%nodeHierarchyLevelFunction       (Node_Component_Merging_statistics_Standard_Hierarchy_Level)
    call mergingStatistics%nodeHierarchyLevelMaximumFunction(Node_Component_Merging_statistics_Standard_HLM            )
    return
  end subroutine Node_Component_Merging_Statistics_Standard_Initialize

  !# <nodeComponentThreadInitializationTask>
  !#  <unitName>Node_Component_Merging_Statistics_Standard_Thread_Initialize</unitName>
  !# </nodeComponentThreadInitializationTask>
  subroutine Node_Component_Merging_Statistics_Standard_Thread_Initialize(parameters_)
    !% Initializes the tree node standard merging statistics module.
    use :: Events_Hooks    , only : nodePromotionEvent               , subhaloPromotionEvent, openMPThreadBindingAtLevel, dependencyExact, &
         &                          dependencyDirectionBefore        , satelliteMergerEvent , postEvolveEvent           , dependencyRegEx, &
         &                          dependencyDirectionAfter
    use :: Galacticus_Nodes, only : defaultMergingStatisticsComponent
    use :: Input_Parameters, only : inputParameter                   , inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters_
    type(dependencyExact), dimension(1)  :: dependenciesSubhaloPromotion
    type(dependencyRegEx), dimension(1)  :: dependenciesSatelliteMerger

    if (defaultMergingStatisticsComponent%standardIsActive()) then
       !# <objectBuilder class="darkMatterHaloMassAccretionHistory" name="darkMatterHaloMassAccretionHistory_" source="parameters_"/>
       !# <objectBuilder class="mergerMassMovements"                name="mergerMassMovements_"                source="parameters_"/>
       dependenciesSubhaloPromotion(1)=dependencyExact(dependencyDirectionBefore,'mergerTreeNodeEvolver')
       dependenciesSatelliteMerger (1)=dependencyRegEx(dependencyDirectionAfter ,'^remnantStructure:'   )
       call nodePromotionEvent   %attach(defaultMergingStatisticsComponent,nodePromotion       ,openMPThreadBindingAtLevel,label='nodeComponentMergingStatisticsStandard'                                          )
       call subhaloPromotionEvent%attach(defaultMergingStatisticsComponent,nodeSubhaloPromotion,openMPThreadBindingAtLevel,label='nodeComponentMergingStatisticsStandard',dependencies=dependenciesSubhaloPromotion)
       call satelliteMergerEvent %attach(defaultMergingStatisticsComponent,satelliteMerger     ,openMPThreadBindingAtLevel,label='nodeComponentMergingStatisticsStandard',dependencies=dependenciesSatelliteMerger )
       call postEvolveEvent      %attach(defaultMergingStatisticsComponent,postEvolve          ,openMPThreadBindingAtLevel,label='nodeComponentMergingStatisticsStandard'                                          )
    end if
    return
  end subroutine Node_Component_Merging_Statistics_Standard_Thread_Initialize

  !# <nodeComponentThreadUninitializationTask>
  !#  <unitName>Node_Component_Merging_Statistics_Standard_Thread_Uninitialize</unitName>
  !# </nodeComponentThreadUninitializationTask>
  subroutine Node_Component_Merging_Statistics_Standard_Thread_Uninitialize()
    !% Uninitializes the tree node standard merging statistics module.
    use :: Events_Hooks    , only : nodePromotionEvent               , subhaloPromotionEvent, satelliteMergerEvent, postEvolveEvent
    use :: Galacticus_Nodes, only : defaultMergingStatisticsComponent
    implicit none

    if (defaultMergingStatisticsComponent%standardIsActive()) then
       !# <objectDestructor name="darkMatterHaloMassAccretionHistory_"/>
       !# <objectDestructor name="mergerMassMovements_"               />
       call nodePromotionEvent   %detach(defaultMergingStatisticsComponent,nodePromotion       )
       call subhaloPromotionEvent%detach(defaultMergingStatisticsComponent,nodeSubhaloPromotion)
       call satelliteMergerEvent %detach(defaultMergingStatisticsComponent,satelliteMerger     )
       call postEvolveEvent      %detach(defaultMergingStatisticsComponent,postEvolve          )
    end if
    return
  end subroutine Node_Component_Merging_Statistics_Standard_Thread_Uninitialize

  integer function Node_Component_Merging_Statistics_Standard_Hierarchy_Level(self)
    !% Return the hierarchy level of a node, computing it if necessary.
    use :: Galacticus_Nodes, only : nodeComponentMergingStatisticsStandard, treeNode
    implicit none
    class  (nodeComponentMergingStatisticsStandard), intent(inout) :: self
    type   (treeNode                              ), pointer       :: hostNode
    integer                                                        :: nodeHierarchyLevel

    if (self%nodeHierarchyLevelValue() < 0) then
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
    use :: Galacticus_Nodes, only : nodeComponentBasic       , nodeComponentMergingStatistics, nodeComponentMergingStatisticsStandard, nodeEvent, &
          &                         nodeEventSubhaloPromotion, treeNode
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
    use :: Dark_Matter_Halo_Formation_Times, only : Dark_Matter_Halo_Formation_Time
    use :: Galacticus_Nodes                , only : defaultMergingStatisticsComponent, nodeComponentBasic, nodeComponentMergingStatistics, treeNode
    use :: Merger_Tree_Walkers             , only : mergerTreeWalkerAllNodes
    implicit none
    type   (treeNode                      ), intent(inout), pointer :: node
    type   (treeNode                      )               , pointer :: nodeHost          , nodeWork
    class  (nodeComponentMergingStatistics)               , pointer :: mergingStatistics
    class  (nodeComponentBasic            )               , pointer :: basic
    type   (mergerTreeWalkerAllNodes      )                         :: treeWalker
    integer                                                         :: nodeHierarchyLevel, nodeHierarchyLevelDepth

    ! Return immediately if this class is not active.
    if (.not.defaultMergingStatisticsComponent%standardIsActive()) return
    ! Perform our own depth-first walk of the tree. This is necessary as we must set merging statistic properties based on the
    ! initial state of the tree, while for nodes that exist at later times this function may only be called on that node after the
    ! tree is already evolved.
    treeWalker=mergerTreeWalkerAllNodes(node%hostTree,spanForest=.false.)
    do while (treeWalker%next(nodeWork))
       ! Create a merger statistics component and initialize it.
       mergingStatistics => nodeWork%mergingStatistics(autoCreate=.true.)
       ! If merging statistics are already set in this node, then they have been set in all nodes in this tree and we can simply
       ! return.
       if (mergingStatistics%nodeFormationTime() >= 0.0d0) return
       nodeHierarchyLevel =  0
       nodeHost           => nodeWork
       do while (nodeHost%isSatellite())
          nodeHierarchyLevel =  nodeHierarchyLevel       +1
          nodeHost           => nodeHost          %parent
       end do
       ! Find the maximum possible hierarchy level.
       nodeHierarchyLevelDepth =  0
       nodeHost                => nodeWork
       do while (associated(nodeHost))
          if (.not.nodeHost%isPrimaryProgenitor() .and. associated(nodeHost%parent)) &
               & nodeHierarchyLevelDepth =  nodeHierarchyLevelDepth       +1
          nodeHost                       => nodeHost               %parent
       end do
       ! Set the initial merging statistics.
       basic => nodeWork%basic()
       call mergingStatistics%       nodeHierarchyLevelSet(nodeHierarchyLevel     )
       call mergingStatistics%nodeHierarchyLevelMaximumSet(nodeHierarchyLevel     )
       call mergingStatistics%  nodeHierarchyLevelDepthSet(nodeHierarchyLevelDepth)
       call mergingStatistics%    massWhenFirstIsolatedSet(           basic%mass())
       call mergingStatistics%    galaxyMajorMergerTimeSet(                 -1.0d0)
       call mergingStatistics%      nodeMajorMergerTimeSet(                 -1.0d0)
       call mergingStatistics%        nodeFormationTimeSet(Dark_Matter_Halo_Formation_Time(nodeWork,nodeFormationMassFraction,darkMatterHaloMassAccretionHistory_))
    end do
    return
  end subroutine Node_Component_Merging_Statistics_Standard_Merger_Tree_Init

  !# <nodeMergerTask>
  !#  <unitName>Node_Component_Merging_Statistics_Standard_Node_Merger</unitName>
  !# </nodeMergerTask>
  subroutine Node_Component_Merging_Statistics_Standard_Node_Merger(node)
    !% Record any major merger of {\normalfont \ttfamily node}.
    use :: Galacticus_Nodes, only : defaultMergingStatisticsComponent, nodeComponentBasic, nodeComponentMergingStatistics, treeNode
    implicit none
    type (treeNode                      ), intent(inout), pointer :: node
    class(nodeComponentMergingStatistics)               , pointer :: mergingStatisticsParent
    class(nodeComponentBasic            )               , pointer :: parentBasicComponent   , basic

    ! Return immediately if this class is not active.
    if (.not.defaultMergingStatisticsComponent%standardIsActive()) return
    ! Record the merger time if this is a major merger.
    basic                   => node       %basic            ()
    parentBasicComponent    => node%parent%basic            ()
    mergingStatisticsParent => node%parent%mergingStatistics()
    if (basic%mass() >= nodeMajorMergerFraction*parentBasicComponent%mass()) &
         &  call mergingStatisticsParent%nodeMajorMergerTimeSet(basic%time())
    ! Increment the hierarchy level of the merging node.
    call Node_Component_Merging_Statistics_Standard_Reset_Hierarchy(node)
    call hierarchyLevelIncrement(node)
    return
    
  contains

    recursive subroutine hierarchyLevelIncrement(node)
      !% Increment the hierarchy level of the given node, and then call our self on any satellite nodes.
      implicit none
      type (treeNode                      ), intent(inout), target  :: node
      type (treeNode                      ),                pointer :: satelliteNode
      class(nodeComponentMergingStatistics)               , pointer :: mergingStatistics

      ! Increment the hierarchy level.
      mergingStatistics => node%mergingStatistics()
      call mergingStatistics%nodeHierarchyLevelSet       (    mergingStatistics%nodeHierarchyLevel       ()+1)
      call mergingStatistics%nodeHierarchyLevelMaximumSet(                                                      &
           &                                              max(                                                  &
           &                                                  mergingStatistics%nodeHierarchyLevel       ()   , &
           &                                                  mergingStatistics%nodeHierarchyLevelMaximum()     &
           &                                                 )                                                  &
           &                                             )
      ! Increment the hierarchy level of any satellites.
      satelliteNode => node%firstSatellite
      do while (associated(satelliteNode))
         call hierarchyLevelIncrement(satelliteNode)
         satelliteNode => satelliteNode%sibling
      end do
      return
    end subroutine hierarchyLevelIncrement
    
  end subroutine Node_Component_Merging_Statistics_Standard_Node_Merger

  !# <satelliteHostChangeTask>
  !#  <unitName>Node_Component_Merging_Statistics_Standard_Host_Change</unitName>
  !# </satelliteHostChangeTask>
  recursive subroutine Node_Component_Merging_Statistics_Standard_Host_Change(node)
    !% Handle cases where a satellite switches host node.
    use :: Galacticus_Nodes, only : defaultMergingStatisticsComponent, nodeComponentMergingStatistics, treeNode
    implicit none
    type   (treeNode                      ), intent(inout), target  :: node
    type   (treeNode                      ),                pointer :: satelliteNode     , nodeHost
    class  (nodeComponentMergingStatistics)               , pointer :: mergingStatistics
    integer                                                         :: nodeHierarchyLevel

    ! Return immediately if this class is not active.
    if (.not.defaultMergingStatisticsComponent%standardIsActive()) return
    ! Recompute the hierarchy level of the satellite node.
    nodeHierarchyLevel =  0
    nodeHost           => node
    mergingStatistics  => node%mergingStatistics()
    do while (nodeHost%isSatellite())
       nodeHierarchyLevel =  nodeHierarchyLevel       +1
       nodeHost           => nodeHost          %parent
    end do
    call mergingStatistics%nodeHierarchyLevelSet(nodeHierarchyLevel)
    ! Call this function on any satellites of the satellite node.
    satelliteNode => node%firstSatellite
    do while (associated(satelliteNode))
       call Node_Component_Merging_Statistics_Standard_Host_Change(satelliteNode)
       satelliteNode => satelliteNode%sibling
    end do
    return
  end subroutine Node_Component_Merging_Statistics_Standard_Host_Change

  subroutine nodePromotion(self,node)
    !% Ensure that {\normalfont \ttfamily node} is ready for promotion to its parent. In this case, we simply update the node merger time.
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentMergingStatistics, treeNode
    implicit none
    class(*                             ), intent(inout)          :: self
    type (treeNode                      ), intent(inout), target  :: node
    class(nodeComponentMergingStatistics)               , pointer :: mergingStatisticsParent, mergingStatistics
    class(nodeComponentBasic            )               , pointer :: basicParent
    !$GLC attributes unused :: self
    
    mergingStatistics       => node       %mergingStatistics()
    mergingStatisticsParent => node%parent%mergingStatistics()
    basicParent             => node%parent%basic            ()
    call Node_Component_Merging_Statistics_Standard_Reset_Hierarchy(node)
    if (mergingStatisticsParent%nodeMajorMergerTime() > mergingStatistics%nodeMajorMergerTime()) &
         & call mergingStatistics%nodeMajorMergerTimeSet(mergingStatisticsParent%nodeMajorMergerTime())
    call        mergingStatistics% nodeHierarchyLevelSet(mergingStatisticsParent%nodeHierarchyLevel ())
    call        mergingStatistics%  nodeFormationTimeSet(mergingStatisticsParent%nodeFormationTime  ())
    return
  end subroutine nodePromotion

  subroutine nodeSubhaloPromotion(self,node,nodePromotion)
    !% Reset the mass-when-first-isolated property of the merging statistics component in the event of the subhalo promotion.
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentMergingStatistics, treeNode
    implicit none
    class(*                             ), intent(inout)          :: self
    type (treeNode                      ), intent(inout), pointer :: node             , nodePromotion
    class(nodeComponentMergingStatistics)               , pointer :: mergingStatistics
    class(nodeComponentBasic            )               , pointer :: basicParent
    !$GLC attributes unused :: self

    mergingStatistics => node         %mergingStatistics()
    basicParent       => nodePromotion%basic            ()
    call mergingStatistics%massWhenFirstIsolatedSet(basicParent%mass())
    return
  end subroutine nodeSubhaloPromotion

  subroutine postEvolve(self,node)
    !% Handle post differential evoloution.
    use :: Galacticus_Nodes, only : treeNode
    implicit none
    class(*       ), intent(inout) :: self
    type (treeNode), intent(inout) :: node
    !$GLC attributes unused :: self

    call Node_Component_Merging_Statistics_Standard_Reset_Hierarchy(node)
    return
  end subroutine postEvolve

  subroutine Node_Component_Merging_Statistics_Standard_Reset_Hierarchy(node)
    !% Reset the maximum node hierarchy level if the node has grown sufficiently in mass.
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentMergingStatistics, treeNode
    implicit none
    type (treeNode                      ), intent(inout) :: node
    class(nodeComponentBasic            ), pointer       :: basic
    class(nodeComponentMergingStatistics), pointer       :: mergingStatistics

    ! Test if previous mass at promotion has been exceeded by a sufficient factor.
    basic             => node%basic            ()
    mergingStatistics => node%mergingStatistics()
    if (basic%mass() > hierarchyLevelResetFactor*mergingStatistics%massWhenFirstIsolated()) &
         & call mergingStatistics%nodeHierarchyLevelMaximumSet(mergingStatistics%nodeHierarchyLevel())
    return
  end subroutine Node_Component_Merging_Statistics_Standard_Reset_Hierarchy

  subroutine satelliteMerger(self,node)
    !% Record properties of a merging event for {\normalfont \ttfamily node}.
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentMergingStatistics, treeNode
    implicit none
    class  (*                             ), intent(inout) :: self
    type   (treeNode                      ), intent(inout) :: node
    type   (treeNode                      ), pointer       :: nodeHost
    class  (nodeComponentBasic            ), pointer       :: basicHost
    class  (nodeComponentMergingStatistics), pointer       :: mergingStatisticsHost
    integer                                                :: destinationGasSatellite, destinationGasHost       , &
         &                                                    destinationStarsHost   , destinationStarsSatellite
    logical                                                :: mergerIsMajor
    !$GLC attributes unused :: self

    ! Record the time of this merger if it is a major merger.
    call mergerMassMovements_%get(node,destinationGasSatellite,destinationStarsSatellite,destinationGasHost,destinationStarsHost,mergerIsMajor)
    if (mergerIsMajor) then
       ! Find the node to merge with.
       nodeHost              => node    %mergesWith       ()
       basicHost             => nodeHost%basic            ()
       mergingStatisticsHost => nodeHost%mergingStatistics()
       ! Record the merger time.
       call mergingStatisticsHost%galaxyMajorMergerTimeSet(basicHost%time())
    end if
    return
  end subroutine satelliteMerger

end module Node_Component_Merging_Statistics_Standard
