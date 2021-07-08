!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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
Contains a module which implements the standard indices component.
!!}

module Node_Component_Indices_Standard
  !!{
  Implements the standard indices component.
  !!}
  implicit none
  private
  public :: Node_Component_Indices_Standard_Merger_Tree_Init

  !![
  <component>
   <class>indices</class>
   <name>standard</name>
   <isDefault>false</isDefault>
   <properties>
    <property>
      <name>branchTip</name>
      <type>longInteger</type>
      <rank>0</rank>
      <attributes isSettable="true" isGettable="true" isEvolvable="false" />
      <output unitsInSI="0.0d0" comment="Index of the node at the tip of the branch in which this galaxy formed."/>
    </property>
   </properties>
  </component>
  !!]

contains

  !![
  <mergerTreeInitializeTask>
   <unitName>Node_Component_Indices_Standard_Merger_Tree_Init</unitName>
  </mergerTreeInitializeTask>
  !!]
  subroutine Node_Component_Indices_Standard_Merger_Tree_Init(node)
    !!{
    Initialize the indices component by creating components in nodes and storing indices.
    !!}
    use :: Galacticus_Nodes, only : defaultIndicesComponent, nodeComponentIndices, nodeComponentIndicesStandard, treeNode
    implicit none
    type (treeNode            ), intent(inout), pointer :: node
    class(nodeComponentIndices)               , pointer :: indices

    ! Return immediately if this class is not active.
    if (.not.defaultIndicesComponent%standardIsActive()) return

    ! Return immediately if this node is not a branch tip.
    if (associated(node%firstChild)                    ) return

    ! Create an indices component and initialize it.
    indices => node%indices(autoCreate=.true.)
    select type (indices)
    class is (nodeComponentIndicesStandard)
       call indices%branchTipSet(node%index())
    end select
    return
  end subroutine Node_Component_Indices_Standard_Merger_Tree_Init

end module Node_Component_Indices_Standard
