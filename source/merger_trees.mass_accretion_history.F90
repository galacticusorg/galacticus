!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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
  use IO_HDF5
  implicit none
  private
  public :: Merger_Tree_Mass_Accretion_History_Output, Merger_Tree_Mass_Accretion_History_Close

  ! Flag indicating if module is initialized.
  logical          :: accretionHistoryModuleInitialized=.false.

  ! Flag indicating if output is required.
  logical          :: massAccretionHistoryOutput

  ! Accretion group object.
  type(hdf5Object) :: accretionGroup

contains

  !# <mergerTreePreEvolveTask>
  !#   <unitName>Merger_Tree_Mass_Accretion_History_Output</unitName>
  !# </mergerTreePreEvolveTask>
  subroutine Merger_Tree_Mass_Accretion_History_Output(thisTree)
    !% Output the mass accretion history of {\tt thisTree}.
    use Merger_Trees
    use Galacticus_Nodes
    use Input_Parameters
    use Memory_Management
    use Galacticus_HDF5
    use ISO_Varying_String
    use String_Handling
    use Numerical_Constants_Astronomical
    implicit none
    type(mergerTree),          intent(in) , target       :: thisTree
    type(treeNode),            pointer                   :: thisNode
    integer(kind=kind_int8),   allocatable, dimension(:) :: accretionHistoryNodeIndex
    double precision,          allocatable, dimension(:) :: accretionHistoryNodeMass,accretionHistoryNodeTime
    class(nodeComponentBasic), pointer                   :: thisBasicComponent
    type(mergerTree),          pointer                   :: currentTree
    integer(kind=kind_int8)                              :: accretionHistoryCount
    type(varying_string)                                 :: groupName
    type(hdf5Object)                                     :: treeGroup,accretionDataset

    ! Check if module is initialized.
    if (.not.accretionHistoryModuleInitialized) then
       !$omp critical(accretionHistoryModuleInitialize)
       if (.not.accretionHistoryModuleInitialized) then
          ! Get parameter specifying if output is required.
          !@ <inputParameter>
          !@   <name>massAccretionHistoryOutput</name>
          !@   <defaultValue>false</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Specifies whether or not to output mass accretion histories for the main branches of merger trees.
          !@   </description>
          !@   <type>boolean</type>
          !@   <cardinality>1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('massAccretionHistoryOutput',massAccretionHistoryOutput,defaultValue=.false.)
          ! Create an output group if necessary.
          !$omp critical (HDF5_Access)
          if (massAccretionHistoryOutput) accretionGroup=galacticusOutputFile%openGroup('massAccretionHistories','Mass&
               & accretion histories of main branches in merger trees.')
          !$omp end critical (HDF5_Access)
          ! Flag that module is initialized.
          accretionHistoryModuleInitialized=.true.
       end if
       !$omp end critical(accretionHistoryModuleInitialize)
    end if

    ! Output the mass accretion history.
    if (massAccretionHistoryOutput) then
       ! Iterate over trees.
       currentTree => thisTree
       do while (associated(currentTree))
          ! Count up number of entries expected for accretion history.
          accretionHistoryCount=0
          thisNode => currentTree%baseNode
          do while (associated(thisNode))
             accretionHistoryCount=accretionHistoryCount+1
             thisNode => thisNode%firstChild
          end do
          ! Allocate storage space.
          call Alloc_Array(accretionHistoryNodeIndex,[int(accretionHistoryCount)])
          call Alloc_Array(accretionHistoryNodeTime ,[int(accretionHistoryCount)])
          call Alloc_Array(accretionHistoryNodeMass ,[int(accretionHistoryCount)])
          ! Extract accretion history.
          accretionHistoryCount=0
          thisNode => currentTree%baseNode
          do while (associated(thisNode))
             thisBasicComponent => thisNode%basic()
             accretionHistoryCount=accretionHistoryCount+1
             accretionHistoryNodeIndex(accretionHistoryCount)=               thisNode%index()
             accretionHistoryNodeTime (accretionHistoryCount)=thisBasicComponent%time()
             accretionHistoryNodeMass (accretionHistoryCount)=thisBasicComponent%mass()
             thisNode => thisNode%firstChild
          end do
          ! Output to HDF5 file.
          groupName   ='mergerTree'
          groupName   =groupName//currentTree%index
          !$omp critical (HDF5_Access)
          treeGroup=accretionGroup%openGroup(char(groupName),'Mass accretion history for main branch of merger tree.')
          call treeGroup%writeDataset(accretionHistoryNodeIndex,'nodeIndex','Index of the node.'         )
          call treeGroup%writeDataset(accretionHistoryNodeTime ,'nodeTime' ,'Time at node [Gyr].'        ,datasetReturned=accretionDataset)
          call accretionDataset%writeAttribute(gigaYear ,"unitsInSI")
          call accretionDataset%close()
          call treeGroup%writeDataset(accretionHistoryNodeMass ,'nodeMass' ,'Mass of the node [M⊙].',datasetReturned=accretionDataset)
          call accretionDataset%writeAttribute(massSolar,"unitsInSI")
          call accretionDataset%close()
          call treeGroup       %close()
          !$omp end critical (HDF5_Access)
          ! Deallocate storage space.
          call Dealloc_Array(accretionHistoryNodeIndex)
          call Dealloc_Array(accretionHistoryNodeTime )
          call Dealloc_Array(accretionHistoryNodeMass )
          ! Move to the next tree.
          currentTree => currentTree%nextTree
       end do
    end if

    return
  end subroutine Merger_Tree_Mass_Accretion_History_Output

  !# <hdfPreCloseTask>
  !#   <unitName>Merger_Tree_Mass_Accretion_History_Close</unitName>
  !# </hdfPreCloseTask>
  subroutine Merger_Tree_Mass_Accretion_History_Close()
    !% Close the mass accretion history group before closing the HDF5 file.
    implicit none

    if (massAccretionHistoryOutput) call accretionGroup%close()
    return
  end subroutine Merger_Tree_Mass_Accretion_History_Close

end module Merger_Tree_Mass_Accretion_History
