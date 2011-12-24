!! Copyright 2009, 2010, 2011 Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


!% Contains a module which outputs mass accretion histories of merger trees.

module Merger_Tree_Mass_Accretion_History
  !% Outputs mass accretion histories of merger trees.
  implicit none
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
    use Input_Parameters
    use Memory_Management
    use IO_HDF5
    use Galacticus_HDF5
    use ISO_Varying_String
    use String_Handling
    use Numerical_Constants_Astronomical
    use Kind_Numbers
    implicit none
    type(mergerTree),        intent(in)                :: thisTree
    type(treeNode),          pointer                   :: thisNode
    integer(kind=kind_int8), allocatable, dimension(:) :: accretionHistoryNodeIndex
    double precision,        allocatable, dimension(:) :: accretionHistoryNodeMass,accretionHistoryNodeTime
    type(hdf5Object),        save                      :: accretionGroup
    integer(kind=kind_int8)                            :: accretionHistoryCount
    type(varying_string)                               :: groupName
    type(hdf5Object)                                   :: treeGroup,accretionDataset

    ! Check if module is initialized.
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
       if (massAccretionHistoryOutput) accretionGroup=IO_HDF5_Open_Group(galacticusOutputFile,'massAccretionHistories','Mass&
            & accretion histories of main branches in merger trees.')
       !$omp end critical (HDF5_Access)
       ! Flag that module is initialized.
       accretionHistoryModuleInitialized=.true.
    end if
    !$omp end critical(accretionHistoryModuleInitialize)

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
       call Alloc_Array(accretionHistoryNodeIndex,[int(accretionHistoryCount)])
       call Alloc_Array(accretionHistoryNodeTime ,[int(accretionHistoryCount)])
       call Alloc_Array(accretionHistoryNodeMass ,[int(accretionHistoryCount)])
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
       !$omp critical (HDF5_Access)
       treeGroup=IO_HDF5_Open_Group(accretionGroup,char(groupName),'Mass accretion history for main branch of merger tree.')
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
    end if

    return
  end subroutine Merger_Tree_Mass_Accretion_History_Output
  
end module Merger_Tree_Mass_Accretion_History
