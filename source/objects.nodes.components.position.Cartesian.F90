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
Contains a module which implements a position component in Cartesian coordinates.
!!}

module Node_Component_Position_Cartesian
  !!{
  Implements a position component in Cartesian coordinates.
  !!}
  implicit none
  private
  public :: Node_Component_Position_Cartesian_Thread_Initialize, Node_Component_Position_Cartesian_Thread_Uninitialize, &
       &    Node_Component_Position_Cartesian_Inter_Tree_Insert 
  !![
  <component>
   <class>position</class>
   <name>cartesian</name>
   <isDefault>false</isDefault>
   <properties>
    <property>
      <name>position</name>
      <type>double</type>
      <rank>1</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" />
      <output labels="[X,Y,Z]" unitsInSI="megaParsec" comment="Position of the node (in physical coordinates)."/>
      <classDefault>[0.0d0,0.0d0,0.0d0]</classDefault>
    </property>
    <property>
      <name>velocity</name>
      <type>double</type>
      <rank>1</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="true" />
      <output labels="[X,Y,Z]" unitsInSI="kilo" comment="Velocity of the node (in physical coordinates)."/>
      <classDefault>[0.0d0,0.0d0,0.0d0]</classDefault>
    </property>
    <property>
      <name>positionHistory</name>
      <type>history</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="false" />
    </property>
   </properties>
  </component>
  !!]

  ! A threadprivate object used to track to which thread events are attached.
  integer :: thread
  !$omp threadprivate(thread)

contains

  !![
  <nodeComponentThreadInitializationTask>
   <unitName>Node_Component_Position_Cartesian_Thread_Initialize</unitName>
  </nodeComponentThreadInitializationTask>
  !!]
  subroutine Node_Component_Position_Cartesian_Thread_Initialize(parameters_)
    !!{
    Initializes the tree node scale dark matter profile module.
    !!}
    use :: Events_Hooks    , only : nodePromotionEvent      , openMPThreadBindingAtLevel
    use :: Galacticus_Nodes, only : defaultPositionComponent
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters_
    !$GLC attributes unused :: parameters_

    if (defaultPositionComponent%cartesianIsActive()) &
         call nodePromotionEvent%attach(thread,nodePromotion,openMPThreadBindingAtLevel,label='nodeComponentPositionCartesian')
    return
  end subroutine Node_Component_Position_Cartesian_Thread_Initialize

  !![
  <nodeComponentThreadUninitializationTask>
   <unitName>Node_Component_Position_Cartesian_Thread_Uninitialize</unitName>
  </nodeComponentThreadUninitializationTask>
  !!]
  subroutine Node_Component_Position_Cartesian_Thread_Uninitialize()
    !!{
    Uninitializes the tree node scale dark matter profile module.
    !!}
    use :: Events_Hooks    , only : nodePromotionEvent
    use :: Galacticus_Nodes, only : defaultPositionComponent
    implicit none

    if (defaultPositionComponent%cartesianIsActive() .and. nodePromotionEvent%isAttached(thread,nodePromotion)) &
         & call nodePromotionEvent%detach(thread,nodePromotion)
    return
  end subroutine Node_Component_Position_Cartesian_Thread_Uninitialize

  subroutine nodePromotion(self,node)
    !!{
    Ensure that {\normalfont \ttfamily node} is ready for promotion to its parent. In this case, update the position of {\normalfont \ttfamily
    node} to that of the parent.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentPosition, nodeComponentPositionCartesian, treeNode
    implicit none
    class(*                    ), intent(inout)          :: self
    type (treeNode             ), intent(inout), target  :: node
    class(nodeComponentPosition)               , pointer :: positionParent, position
    !$GLC attributes unused :: self
    
    position       => node       %position()
    positionParent => node%parent%position()
    select type (positionParent)
    class is (nodeComponentPositionCartesian)
       call position%       positionSet(positionParent%position       ())
       call position%       velocitySet(positionParent%velocity       ())
       call position%positionHistorySet(positionParent%positionHistory())
    end select
    return
  end subroutine nodePromotion

  !![
  <interTreePositionInsert>
   <unitName>Node_Component_Position_Cartesian_Inter_Tree_Insert</unitName>
  </interTreePositionInsert>
  !!]
  subroutine Node_Component_Position_Cartesian_Inter_Tree_Insert(node,replaceNode)
    !!{
    A satellite node is being moved between trees, and being added as a new satellite. Its (future-)histories will have been
    assigned to the {\normalfont \ttfamily replaceNode} so must be transferred.
    !!}
    use :: Galacticus_Nodes, only : defaultPositionComponent, nodeComponentBasic, nodeComponentPosition, treeNode
    use :: Histories       , only : history
    implicit none
    type (treeNode             ), intent(inout), pointer :: node               , replaceNode
    class(nodeComponentPosition)               , pointer :: position           , replacePosition
    class(nodeComponentBasic   )               , pointer :: basic
    type (history              )                         :: historyPosition    , replaceHistoryPosition, &
         &                                                  moveHistoryPosition

    ! Return immediately if the cartesian position implementation is not active.
    if (.not.defaultPositionComponent%cartesianIsActive()) return
    ! Get the basic component of the pulled node.
    basic                  =>           node%basic          ()
    ! Get the position components to both nodes.
    position               =>           node%position       ()
    replacePosition        =>    replaceNode%position       ()
     ! Transfer subhalo mass history.
    historyPosition        =        position%positionHistory()
    replaceHistoryPosition = replacePosition%positionHistory()
    ! Cut off history in node being replaced subsequent to current time.
    call replaceHistoryPosition%trimForward       (basic%time(),   moveHistoryPosition)
    ! Append removed history to pulled node.
    call historyPosition       %append            (                moveHistoryPosition)
    ! Set the histories.
    call        position       %positionHistorySet(                    historyPosition)
    call replacePosition       %positionHistorySet(             replaceHistoryPosition)
    return
  end subroutine Node_Component_Position_Cartesian_Inter_Tree_Insert

end module Node_Component_Position_Cartesian
