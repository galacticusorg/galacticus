!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
!!    Andrew Benson <abenson@obs.carnegiescience.edu>
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
  use Galacticus_Nodes
  implicit none
  private
  public :: Node_Component_Host_History_Standard_Merger_Tree_Init, Node_Component_Host_History_Standard_Update_History

  !# <component>
  !#  <class>hostHistory</class>
  !#  <name>standard</name>
  !#  <isDefault>no</isDefault>
  !#  <properties>
  !#   <property>
  !#     <name>hostMassMaximum</name>
  !#     <type>double</type>
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

  !# <nodeMergerTask>
  !#  <unitName>Node_Component_Host_History_Standard_Update_History</unitName>
  !# </nodeMergerTask>
  !# <nodePromotionTask>
  !#  <unitName>Node_Component_Host_History_Standard_Update_History</unitName>
  !# </nodePromotionTask>
  !# <postEvolveTask>
  !# <unitName>Node_Component_Host_History_Standard_Update_History</unitName>
  !# </postEvolveTask>
  subroutine Node_Component_Host_History_Standard_Update_History(thisNode)
    !% Record any major merger of {\tt thisNode}.
    implicit none
    type (treeNode                ), pointer, intent(inout) :: thisNode
    class(nodeComponentHostHistory), pointer                :: thisHostHistory
    class(nodeComponentBasic      ), pointer                :: hostBasic
    
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
