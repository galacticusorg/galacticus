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

!% Contains a module which implements the standard indices component.

module Node_Component_Indices_Standard
  !% Implements the standard indices component.
  use Galacticus_Nodes
  implicit none
  private
  public :: Node_Component_Indices_Standard_Merger_Tree_Init

  !# <component>
  !#  <class>indices</class>
  !#  <name>standard</name>
  !#  <isDefault>no</isDefault>
  !#  <properties>
  !#   <property>
  !#     <name>branchTip</name>
  !#     <type>longInteger</type>
  !#     <rank>0</rank>
  !#     <attributes isSettable="true" isGettable="true" isEvolvable="false" />
  !#     <output unitsInSI="0.0d0" comment="Index of the node at the tip of the branch in which this galaxy formed."/>
  !#   </property>
  !#  </properties>
  !# </component>

contains

  !# <mergerTreeInitializeTask>
  !#  <unitName>Node_Component_Indices_Standard_Merger_Tree_Init</unitName>
  !# </mergerTreeInitializeTask>
  subroutine Node_Component_Indices_Standard_Merger_Tree_Init(thisNode)
    !% Initialize the indices component by creating components in nodes and storing indices.
    implicit none
    type (treeNode            ), intent(inout), pointer :: thisNode
    class(nodeComponentIndices)               , pointer :: thisIndices

    ! Return immediately if this class is not active.
    if (.not.defaultIndicesComponent%standardIsActive()) return

    ! Return immediately if this node is not a branch tip.
    if (associated(thisNode%firstChild)                ) return

    ! Create an indices component and initialize it.
    thisIndices => thisNode%indices(autoCreate=.true.)
    select type (thisIndices)
    class is (nodeComponentIndicesStandard)
       call thisIndices%branchTipSet(thisNode%index())
    end select
    return
  end subroutine Node_Component_Indices_Standard_Merger_Tree_Init

end module Node_Component_Indices_Standard
