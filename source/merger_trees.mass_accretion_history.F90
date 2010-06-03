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






!% Contains a module which outputs mass accretion histories of merger trees.

module Merger_Tree_Mass_Accretion_History
  !% Outputs mass accretion histories of merger trees.
  private
  public :: Merger_Tree_Mass_Accretion_History_Output
  
  ! Flag indicating if module is initialized.
  logical :: accretionHistoryModuleInitialized=.false.
  
  ! Flag indicating if output is required.
  logical :: massAccretionHistoryOutput

  ! HDF5 group index.
  integer :: accretionGroupID
  
contains

  !# <mergerTreePreEvolveTask>
  !#   <unitName>Merger_Tree_Mass_Accretion_History_Output</unitName>
  !# </mergerTreePreEvolveTask>
  subroutine Merger_Tree_Mass_Accretion_History_Output(thisTree)
    !% Output the mass accretion history of {\tt thisTree}.
    use Merger_Trees
    use Tree_Nodes
    use Tree_Node_Methods
    use Input_Parameters
    use Memory_Management
    use Galacticus_HDF5_Groups
    use HDF5
    use ISO_Varying_String
    use String_Handling
    implicit none
    type(mergerTree), intent(in)                :: thisTree
    type(treeNode),   pointer                   :: thisNode
    integer,          allocatable, dimension(:) :: accretionHistoryNodeIndex
    double precision, allocatable, dimension(:) :: accretionHistoryNodeMass,accretionHistoryNodeTime
    integer                                     :: accretionHistoryCount,accretionDataID,treeGroupID
    type(varying_string)                        :: groupName,groupComment

    ! Check if module is initialized.
    if (.not.accretionHistoryModuleInitialized) then
       ! Get parameter specifying if output is required.
       !@ <inputParameter>
       !@   <name>massAccretionHistoryOutput</name>
       !@   <defaultValue>false</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Specifies whether or not to output mass accretion histories for the main branches of merger trees.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('massAccretionHistoryOutput',massAccretionHistoryOutput,defaultValue=.false.)
       ! Create an output group if necessary.
       if (massAccretionHistoryOutput) then
          groupName   ='massAccretionHistories'
          groupComment='Mass accretion histories of main branches in merger trees.'
          accretionGroupID=Galacticus_Output_Make_Group(groupName,groupComment)
       end if
       ! Flag that module is initialized.
       accretionHistoryModuleInitialized=.true.
    end if

    ! Output the mass accretion history.
    if (massAccretionHistoryOutput) then
       ! Count up number of entries expected for accretion history.
       accretionHistoryCount=0
       thisNode => thisTree%baseNode
       do while (associated(thisNode))
          accretionHistoryCount=accretionHistoryCount+1
          thisNode => thisNode%childNode
       end do
       ! Allocate storage space.
       call Alloc_Array(accretionHistoryNodeIndex,accretionHistoryCount,'accretionHistoryNodeIndex')
       call Alloc_Array(accretionHistoryNodeTime ,accretionHistoryCount,'accretionHistoryNodeTime' )
       call Alloc_Array(accretionHistoryNodeMass ,accretionHistoryCount,'accretionHistoryNodeMass' )
       ! Extract accretion history.
       accretionHistoryCount=0
       thisNode => thisTree%baseNode
       do while (associated(thisNode))
          accretionHistoryCount=accretionHistoryCount+1
          accretionHistoryNodeIndex(accretionHistoryCount)=               thisNode%index()
          accretionHistoryNodeTime (accretionHistoryCount)=Tree_Node_Time(thisNode)
          accretionHistoryNodeMass (accretionHistoryCount)=Tree_Node_Mass(thisNode)
          thisNode => thisNode%childNode
       end do
       ! Output to HDF5 file.
       groupName   ='mergerTree'
       groupName   =groupName//thisTree%index
       groupComment='Mass accretion history for main branch of merger tree.'
       treeGroupID =Galacticus_Output_Make_Group(groupName,groupComment,accretionGroupID)
       accretionDataID=0
       call Galacticus_Output_Dataset(treeGroupID,accretionDataID,'nodeIndex','Index of the node.'         ,accretionHistoryNodeIndex)
       accretionDataID=0
       call Galacticus_Output_Dataset(treeGroupID,accretionDataID,'nodeTime' ,'Time at node [Gyr].'        ,accretionHistoryNodeTime )
       accretionDataID=0
       call Galacticus_Output_Dataset(treeGroupID,accretionDataID,'nodeMass' ,'Mass of the node [M_Solar].',accretionHistoryNodeMass )
       ! Deallocate storage space.
       call Dealloc_Array(accretionHistoryNodeIndex)
       call Dealloc_Array(accretionHistoryNodeTime )
       call Dealloc_Array(accretionHistoryNodeMass )
    end if

    return
  end subroutine Merger_Tree_Mass_Accretion_History_Output
  
end module Merger_Tree_Mass_Accretion_History
