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

!% Contains a module which implements a component class that tracks the maximum host mass seen by each halo.

module Node_Component_Host_History_Standard
  !% Implements a component class that tracks the maximum host mass seen by each halo.
  implicit none
  private
  public :: Node_Component_Host_History_Standard_Merger_Tree_Init   , Node_Component_Host_History_Standard_Thread_Initialize, &
       &    Node_Component_Host_History_Standard_Thread_Uninitialize, Node_Component_Host_History_Standard_Update_History

  !# <component>
  !#  <class>hostHistory</class>
  !#  <name>standard</name>
  !#  <isDefault>false</isDefault>
  !#  <properties>
  !#   <property>
  !#     <name>hostMassMaximum</name>
  !#     <type>real</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#     <classDefault>-1.0d0</classDefault>
  !#     <output unitsInSI="massSolar" comment="The maximum mass of halo in which this node has ever been hosted. (Or $-1$ if this node has never been a subhalo.)"/>
  !#   </property>
  !#  </properties>
  !# </component>

contains

  !# <mergerTreeInitializeTask>
  !#  <unitName>Node_Component_Host_History_Standard_Merger_Tree_Init</unitName>
  !# </mergerTreeInitializeTask>
  subroutine Node_Component_Host_History_Standard_Merger_Tree_Init(thisNode)
    !% Initialize the standard host history component by creating components in nodes and assigning
    !% host mass for satellites.
    use :: Galacticus_Nodes, only : defaultHostHistoryComponent, nodeComponentBasic, nodeComponentHostHistory, nodeComponentHostHistoryStandard, &
          &                         treeNode
    implicit none
    type (treeNode                ), pointer, intent(inout) :: thisNode
    class(nodeComponentHostHistory), pointer                :: thisHostHistory
    class(nodeComponentBasic      ), pointer                :: hostBasic

    ! Return immediately if this class is not active.
    if (.not.defaultHostHistoryComponent%standardIsActive()) return

    ! Create a host history component and initialize it.
    thisHostHistory => thisNode%hostHistory(autoCreate=.true.)
    select type (thisHostHistory)
    class is (nodeComponentHostHistoryStandard)
       if (thisNode%isSatellite()) then
          hostBasic => thisNode%parent%basic()
          call thisHostHistory%hostMassMaximumSet(hostBasic%mass())
       end if
    end select
    return
  end subroutine Node_Component_Host_History_Standard_Merger_Tree_Init
  
  !# <nodeComponentThreadInitializationTask>
  !#  <unitName>Node_Component_Host_History_Standard_Thread_Initialize</unitName>
  !# </nodeComponentThreadInitializationTask>
  subroutine Node_Component_Host_History_Standard_Thread_Initialize(parameters_)
    !% Initializes the tree node scale dark matter profile module.
    use :: Events_Hooks    , only : nodePromotionEvent         , postEvolveEvent, openMPThreadBindingAtLevel
    use :: Galacticus_Nodes, only : defaultHostHistoryComponent
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters_
    !$GLC attributes unused :: parameters_
    
    if (defaultHostHistoryComponent%standardIsActive()) then
       call nodePromotionEvent%attach(defaultHostHistoryComponent,nodePromotion,openMPThreadBindingAtLevel,label='nodeComponentHostHistoryStandard')
       call postEvolveEvent   %attach(defaultHostHistoryComponent,postEvolve   ,openMPThreadBindingAtLevel,label='nodeComponentHostHistoryStandard')
    end if
    return
  end subroutine Node_Component_Host_History_Standard_Thread_Initialize

  !# <nodeComponentThreadUninitializationTask>
  !#  <unitName>Node_Component_Host_History_Standard_Thread_Uninitialize</unitName>
  !# </nodeComponentThreadUninitializationTask>
  subroutine Node_Component_Host_History_Standard_Thread_Uninitialize()
    !% Uninitializes the tree node scale dark matter profile module.
    use :: Events_Hooks    , only : nodePromotionEvent         , postEvolveEvent
    use :: Galacticus_Nodes, only : defaultHostHistoryComponent
    implicit none

    if (defaultHostHistoryComponent%standardIsActive()) then
       call nodePromotionEvent%detach(defaultHostHistoryComponent,nodePromotion)
       call postEvolveEvent   %detach(defaultHostHistoryComponent,postEvolve   )
    end if
    return
  end subroutine Node_Component_Host_History_Standard_Thread_Uninitialize

  subroutine nodePromotion(self,node)
    !% Call the history update function in the event of a node promotion.
    use :: Galacticus_Nodes, only : treeNode
    implicit none
    class(*       ), intent(inout)         :: self
    type (treeNode), intent(inout), target :: node
    !$GLC attributes unused :: self
    
    call Node_Component_Host_History_Standard_Update_History(node)
    return
  end subroutine nodePromotion
  
  subroutine postEvolve(self,node)
    !% Wrapper function for the history update function.
    use :: Galacticus_Nodes, only : treeNode
    implicit none
    class(*       ), intent(inout) :: self
    type (treeNode), intent(inout) :: node
    !$GLC attributes unused :: self
    
    call Node_Component_Host_History_Standard_Update_History(node)
    return
  end subroutine postEvolve
  
  !# <nodeMergerTask>
  !#  <unitName>Node_Component_Host_History_Standard_Update_History</unitName>
  !# </nodeMergerTask>
  subroutine Node_Component_Host_History_Standard_Update_History(thisNode)
    !% Record any major merger of {\normalfont \ttfamily thisNode}.
    use :: Galacticus_Nodes, only : defaultHostHistoryComponent, nodeComponentBasic, nodeComponentHostHistory, treeNode
    implicit none
    type (treeNode                ), intent(inout) :: thisNode
    class(nodeComponentHostHistory), pointer       :: thisHostHistory
    class(nodeComponentBasic      ), pointer       :: hostBasic

    ! Return immediately if this class is not active.
    if (.not.defaultHostHistoryComponent%standardIsActive()) return
    ! Return immediately if thisNode is not a satellite.
    if (.not.thisNode                   %isSatellite     ()) return
    ! Set the maximum host mass to the larger of the current host mass and the previous maximum host mass.
    thisHostHistory => thisNode       %hostHistory()
    hostBasic       => thisNode%parent%basic      ()
    call thisHostHistory%hostMassMaximumSet(                                       &
         &                                  max(                                   &
         &                                      hostBasic      %mass           (), &
         &                                      thisHostHistory%hostMassMaximum()  &
         &                                     )                                   &
         &                                 )
    return
  end subroutine Node_Component_Host_History_Standard_Update_History

end module Node_Component_Host_History_Standard
