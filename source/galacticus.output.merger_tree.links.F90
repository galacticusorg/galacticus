!! Copyright 2009, 2010, Andrew Benson <abenson@caltech.edu>
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






!% Contains a module which handles outputting of tree link data to the \glc\ output file.

module Galacticus_Output_Trees_Links
  !% Handles outputting of tree link data to the \glc\ output file.
  private
  public :: Galacticus_Output_Tree_Links, Galacticus_Output_Tree_Links_Property_Count, Galacticus_Output_Tree_Links_Names

  ! Number of link properties.
  integer, parameter   :: linkPropertyCount=6

contains

  !# <mergerTreeOutputNames>
  !#  <unitName>Galacticus_Output_Tree_Links_Names</unitName>
  !#  <sortName>Galacticus_Output_Tree_Links</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Tree_Links_Names(integerProperty,integerPropertyNames,integerPropertyComments,doubleProperty&
       &,doublePropertyNames,doublePropertyComments,time)
    !% Set the names of link properties to be written to the \glc\ output file.
    implicit none
    double precision, intent(in)                  :: time
    integer,          intent(inout)               :: integerProperty,doubleProperty
    character(len=*), intent(inout), dimension(:) :: integerPropertyNames,integerPropertyComments,doublePropertyNames &
         &,doublePropertyComments
    
    integerProperty=integerProperty+1
    integerPropertyNames   (integerProperty)='nodeIndex'
    integerPropertyComments(integerProperty)='Tree-unique ID for this node.'
    integerProperty=integerProperty+1
    integerPropertyNames   (integerProperty)='parentNode'
    integerPropertyComments(integerProperty)='ID of parent node.'
    integerProperty=integerProperty+1
    integerPropertyNames   (integerProperty)='childNode'
    integerPropertyComments(integerProperty)='ID of primary child node.'
    integerProperty=integerProperty+1
    integerPropertyNames   (integerProperty)='siblingNode'
    integerPropertyComments(integerProperty)='ID of sibling node.'
    integerProperty=integerProperty+1
    integerPropertyNames   (integerProperty)='satelliteNode'
    integerPropertyComments(integerProperty)='ID of first satellite node.'
    integerProperty=integerProperty+1
    integerPropertyNames   (integerProperty)='nodeIsIsolated'
    integerPropertyComments(integerProperty)='Is the node isolated (0|1)?'
    return
  end subroutine Galacticus_Output_Tree_Links_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Galacticus_Output_Tree_Links_Property_Count</unitName>
  !#  <sortName>Galacticus_Output_Tree_Links</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Galacticus_Output_Tree_Links_Property_Count(integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of link properties to be written to the \glc\ output file.
    implicit none
    double precision, intent(in)    :: time
    integer,          intent(inout) :: integerPropertyCount,doublePropertyCount

    integerPropertyCount=integerPropertyCount+linkPropertyCount
    return
  end subroutine Galacticus_Output_Tree_Links_Property_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Galacticus_Output_Tree_Links</unitName>
  !#  <sortName>Galacticus_Output_Tree_Links</sortName>
  !# </mergerTreeOutputTask>
  subroutine Galacticus_Output_Tree_Links(thisNode,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time)
    !% Store link properties in the \glc\ output file buffers.
    use Tree_Nodes
    use Tree_Node_Methods
    implicit none
    double precision, intent(in)             :: time
    type(treeNode),   intent(inout), pointer :: thisNode
    integer,          intent(inout)          :: integerProperty,integerBufferCount,doubleProperty,doubleBufferCount
    integer,          intent(inout)          :: integerBuffer(:,:)
    double precision, intent(inout)          :: doubleBuffer(:,:)

    integerProperty=integerProperty+1
    integerBuffer(integerBufferCount,integerProperty)=thisNode%index()
    integerProperty=integerProperty+1
    integerBuffer(integerBufferCount,integerProperty)=thisNode%parentNode%index()
    integerProperty=integerProperty+1
    integerBuffer(integerBufferCount,integerProperty)=thisNode%childNode%index()
    integerProperty=integerProperty+1
    integerBuffer(integerBufferCount,integerProperty)=thisNode%siblingNode%index()
    integerProperty=integerProperty+1
    integerBuffer(integerBufferCount,integerProperty)=thisNode%satelliteNode%index()
    integerProperty=integerProperty+1
    select case (thisNode%isSatellite())
    case (.true.)
       integerBuffer(integerBufferCount,integerProperty)=0
    case (.false.)
       integerBuffer(integerBufferCount,integerProperty)=1
    end select
    return
  end subroutine Galacticus_Output_Tree_Links

end module Galacticus_Output_Trees_Links
