!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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

  !!{
  Implements a node operator class that computes quantities related to a node's position within the halo/subhalo hierarchy.
  !!}

  !![
  <nodeOperator name="nodeOperatorHierarchy">
   <description>A node operator class that computes quantities related to a node's position within the halo/subhalo hierarchy.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorHierarchy
     !!{
     A node operator class that computes quantities related to a node's position within the halo/subhalo hierarchy.
     !!}
     private
     integer          :: nodeHierarchyLevelID       , nodeHierarchyLevelDepthID, &
          &              nodeHierarchyLevelMaximumID, massWhenFirstIsolatedID
     double precision :: factorMassReset
   contains
     !![
     <methods>
       <method method="increment" description="Increment the hierarchy level of a node."    />
       <method method="reset"     description="Reset the maximum hierarchy level of a node."/>
     </methods>
     !!]
     final     ::                              hierarchyDestructor
     procedure :: nodeTreeInitialize        => hierarchyNodeTreeInitialize
     procedure :: nodesMerge                => hierarchyNodesMerge
     procedure :: nodePromote               => hierarchyNodePromote
     procedure :: differentialEvolutionPost => hierarchyDifferentialEvolutionPost
     procedure :: increment                 => hierarchyIncrement
     procedure :: reset                     => hierarchyReset
     procedure :: autoHook                  => hierarchyAutoHook
  end type nodeOperatorHierarchy
  
  interface nodeOperatorHierarchy
     !!{
     Constructors for the \refClass{nodeOperatorHierarchy} node operator class.
     !!}
     module procedure hierarchyConstructorParameters
     module procedure hierarchyConstructorInternal
  end interface nodeOperatorHierarchy
  
contains

  function hierarchyConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorHierarchy} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nodeOperatorHierarchy)                :: self
    type            (inputParameters      ), intent(inout) :: parameters
    double precision                                       :: factorMassReset
    
    !![
    <inputParameter>
      <name>factorMassReset</name>
      <defaultValue>1.0d100</defaultValue>
      <description>The factor by which a node's mass must increase before the previous maximum hierarchy level is forgotten.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=nodeOperatorHierarchy(factorMassReset)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function hierarchyConstructorParameters

  function hierarchyConstructorInternal(factorMassReset) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorHierarchy} node operator class.
    !!}
    use :: Galacticus_Nodes, only : defaultBasicComponent
    implicit none
    type            (nodeOperatorHierarchy)                :: self
    double precision                       , intent(in   ) :: factorMassReset
    !![
    <constructorAssign variables="factorMassReset"/>
    !!]
    
    !![
    <addMetaProperty component="basic" name="nodeHierarchyLevel"        type="integer" id="self%nodeHierarchyLevelID"        isCreator="yes"                 />
    <addMetaProperty component="basic" name="nodeHierarchyLevelDepth"   type="integer" id="self%nodeHierarchyLevelDepthID"   isCreator="yes"                 />
    <addMetaProperty component="basic" name="nodeHierarchyLevelMaximum" type="integer" id="self%nodeHierarchyLevelMaximumID" isCreator="yes"                 />
    <addMetaProperty component="basic" name="massWhenFirstIsolated"                    id="self%massWhenFirstIsolatedID"     isCreator="yes" isEvolvable="no"/>
    !!]
    return
  end function hierarchyConstructorInternal

  subroutine hierarchyAutoHook(self)
    !!{
    Attach to various event hooks.
    !!}
    use :: Events_Hooks, only : satelliteHostChangeEvent , subhaloPromotionEvent, openMPThreadBindingAtLevel, dependencyExact, &
         &                      dependencyDirectionBefore
    implicit none
    class(nodeOperatorHierarchy), intent(inout) :: self
    type (dependencyExact      ), dimension(1)  :: dependenciesSubhaloPromotion
    
    dependenciesSubhaloPromotion(1)=dependencyExact(dependencyDirectionBefore,'mergerTreeNodeEvolver')
    call satelliteHostChangeEvent%attach(self,satelliteHostChange ,openMPThreadBindingAtLevel,label='hierarchy'                                          )
    call    subhaloPromotionEvent%attach(self,nodeSubhaloPromotion,openMPThreadBindingAtLevel,label='hierarchy',dependencies=dependenciesSubhaloPromotion)
    return
  end subroutine hierarchyAutoHook

  subroutine hierarchyDestructor(self)
    !!{
    Destructor for the \refClass{nodeOperatorHierarchy} node operator class.
    !!}
    use :: Events_Hooks, only : satelliteHostChangeEvent, subhaloPromotionEvent
    implicit none
    type(nodeOperatorHierarchy), intent(inout) :: self

    if (satelliteHostChangeEvent%isAttached(self,satelliteHostChange )) call satelliteHostChangeEvent%detach(self,satelliteHostChange )
    if (   subhaloPromotionEvent%isAttached(self,nodeSubhaloPromotion)) call    subhaloPromotionEvent%detach(self,nodeSubhaloPromotion)
    return
  end subroutine hierarchyDestructor

  subroutine hierarchyNodeTreeInitialize(self,node)
    !!{
    Initialize hierarchy level data.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class   (nodeOperatorHierarchy   ), intent(inout), target  :: self
    type    (treeNode                ), intent(inout), target  :: node
    type    (treeNode                )               , pointer :: nodeHost               , nodeParent
    class   (nodeComponentBasic      )               , pointer :: basicParent            , basic
    integer                                                    :: nodeHierarchyLevelDepth, nodeHierarchyLevel
    
    ! Find the maximum possible hierarchy level.
    if (.not.associated(node%firstChild)) then
       ! No primary progenitor, so we must compute node hierarchy level depth.
       nodeHierarchyLevelDepth =  0
       nodeHost                => node
       do while (associated(nodeHost))
          if (.not.nodeHost%isPrimaryProgenitor() .and. associated(nodeHost%parent)) then
             nodeHierarchyLevelDepth =  nodeHierarchyLevelDepth                 +1
             basicParent             => nodeHost%parent%basic(autoCreate=.true.)
             if (basicParent%integerRank0MetaPropertyGet(self%nodeHierarchyLevelDepthID) > -1) then
                nodeHierarchyLevelDepth= nodeHierarchyLevelDepth+basicParent%integerRank0MetaPropertyGet(self%nodeHierarchyLevelDepthID)
                exit
             end if
          end if
          nodeHost => nodeHost%parent
       end do
       ! Set this depth in all nodes along the branch.
       nodeParent => node
       do while (associated(nodeParent))
          basicParent => nodeParent%basic(autoCreate=.true.)
          call basicParent%integerRank0MetaPropertySet(self%nodeHierarchyLevelDepthID,nodeHierarchyLevelDepth)
          if (nodeParent%isPrimaryProgenitor()) then
             nodeParent => nodeParent%parent
          else
             nodeParent => null()
          end if
       end do
    end if
    ! Find the current hierarchy level.
    nodeHierarchyLevel =  0
    nodeHost           => node
    do while (nodeHost%isSatellite())
       nodeHierarchyLevel =  nodeHierarchyLevel       +1
       nodeHost           => nodeHost          %parent
    end do
    basic => node%basic()
    call basic%integerRank0MetaPropertySet(self%nodeHierarchyLevelID       ,      nodeHierarchyLevel  )
    call basic%integerRank0MetaPropertySet(self%nodeHierarchyLevelMaximumID,      nodeHierarchyLevel  )
    call basic%  floatRank0MetaPropertySet(self%massWhenFirstIsolatedID    ,basic%mass              ())
    return
  end subroutine hierarchyNodeTreeInitialize

  subroutine hierarchyNodesMerge(self,node)
    !!{
    Update hierarchy levels in response to a node merger.
    !!}
    implicit none
    class(nodeOperatorHierarchy), intent(inout) :: self
    type (treeNode             ), intent(inout) :: node

    call self%reset    (node)
    call self%increment(node)
    return
  end subroutine hierarchyNodesMerge
 
  recursive subroutine hierarchyIncrement(self,node)
    !!{
    Increment the hierarchy level of the given node, and then call our self on any satellite nodes.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(nodeOperatorHierarchy), intent(inout)          :: self
    type (treeNode             ), intent(inout), target  :: node
    type (treeNode             )               , pointer :: nodeSatellite
    class(nodeComponentBasic   )               , pointer :: basic
    
    basic => node%basic()
    call basic%integerRank0MetaPropertySet(self%nodeHierarchyLevelID       ,    basic%integerRank0MetaPropertyGet(self%nodeHierarchyLevelID       )+1)
    call basic%integerRank0MetaPropertySet(self%nodeHierarchyLevelMaximumID,                                                                            &
         &                                                                  max(                                                                        &
         &                                                                      basic%integerRank0MetaPropertyGet(self%nodeHierarchyLevelID       )   , &
         &                                                                      basic%integerRank0MetaPropertyGet(self%nodeHierarchyLevelMaximumID)     &
         &                                                                     )                                                                        &
         &                                )
    ! Increment the hierarchy level of any satellites.
    nodeSatellite => node%firstSatellite
    do while (associated(nodeSatellite))
       call self%increment(nodeSatellite)
       nodeSatellite => nodeSatellite%sibling
    end do
    return
  end subroutine hierarchyIncrement

  recursive subroutine satelliteHostChange(self,node)
    !!{
    Handle cases where a satellite switches host node.
    !!}
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class  (*                 ), intent(inout)          :: self
    type   (treeNode          ), intent(inout), target  :: node
    type   (treeNode          ),                pointer :: nodeSatellite     , nodeHost
    class  (nodeComponentBasic)               , pointer :: basic
    integer                                             :: nodeHierarchyLevel

    select type (self)
    class is (nodeOperatorHierarchy)
       ! Recompute the hierarchy level of the satellite node.
       nodeHierarchyLevel =  0
       nodeHost           => node
       basic  => node%basic()
       do while (nodeHost%isSatellite())
          nodeHierarchyLevel =  nodeHierarchyLevel       +1
          nodeHost           => nodeHost          %parent
       end do
       call basic%integerRank0MetaPropertySet(self%nodeHierarchyLevelID,nodeHierarchyLevel)
       ! Call this function on any satellites of the satellite node.
       nodeSatellite => node%firstSatellite
       do while (associated(nodeSatellite))
          call satelliteHostChange(self,nodeSatellite)
          nodeSatellite => nodeSatellite%sibling
       end do
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine satelliteHostChange

  subroutine nodeSubhaloPromotion(self,node,nodePromotion)
    !!{
    Reset the mass-when-first-isolated property of the merging statistics component in the event of the subhalo promotion.
    !!}
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class(*                 ), intent(inout)          :: self
    type (treeNode          ), intent(inout), pointer :: node , nodePromotion
    class(nodeComponentBasic), pointer                :: basic, basicParent

    select type (self)
    class is (nodeOperatorHierarchy)
       basic       => node         %basic()
       basicParent => nodePromotion%basic()
       call basic%floatRank0MetaPropertySet(self%massWhenFirstIsolatedID,basicParent%mass())
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine nodeSubhaloPromotion

  subroutine hierarchyNodePromote(self,node)
    !!{
    Ensure that {\normalfont \ttfamily node} is ready for promotion to its parent.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(nodeOperatorHierarchy), intent(inout) :: self
    type (treeNode             ), intent(inout) :: node
    class(nodeComponentBasic   ), pointer       :: basic, basicParent
    
    basic       => node       %basic()
    basicParent => node%parent%basic()
    call self %reset                      (node                                                                                        )
    call basic%integerRank0MetaPropertySet(self%nodeHierarchyLevelID,basicParent%integerRank0MetaPropertyGet(self%nodeHierarchyLevelID))
    return
  end subroutine hierarchyNodePromote

  subroutine hierarchyDifferentialEvolutionPost(self,node)
    !!{
    Handle post differential evolution.
    !!}
    implicit none
    class(nodeOperatorHierarchy), intent(inout) :: self
    type (treeNode             ), intent(inout) :: node

    call self%reset(node)
    return
  end subroutine hierarchyDifferentialEvolutionPost

  subroutine hierarchyReset(self,node)
    !!{
    Reset the maximum node hierarchy level if the node has grown sufficiently in mass.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(nodeOperatorHierarchy), intent(inout) :: self
    type (treeNode             ), intent(inout) :: node
    class(nodeComponentBasic   ), pointer       :: basic

    ! Test if previous mass has been exceeded by a sufficient factor to reset the maximum hierarchy level.
    basic => node%basic()
    if (basic%mass() > self%factorMassReset*basic%floatRank0MetaPropertyGet(self%massWhenFirstIsolatedID)) &
         & call basic%integerRank0MetaPropertySet(self%nodeHierarchyLevelMaximumID,basic%integerRank0MetaPropertyGet(self%nodeHierarchyLevelID))
    return
  end subroutine hierarchyReset
