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






!% Contains a module which implements reading of merger trees from an HDF5 file.

module Merger_Tree_Read
  !% Implements reading of merger trees from an HDF5 file.
  use Merger_Trees
  use ISO_Varying_String
  use IO_HDF5
  use HDF5
  private
  public :: Merger_Tree_Read_Initialize

  ! The name and identifier of the file from which to read merger trees.
  integer(kind=HID_T)  :: mergerTreeFileID
  type(varying_string) :: mergerTreeReadFileName

  ! Identifier for the merger trees group.
  integer(kind=HID_T)  :: mergerTreesGroupID

  ! Index of the next merger tree to read.
  integer              :: nextTreeToRead=1

contains

  !# <mergerTreeConstructMethod>
  !#  <unitName>Merger_Tree_Read_Initialize</unitName>
  !# </mergerTreeConstructMethod>
  subroutine Merger_Tree_Read_Initialize(mergerTreeConstructMethod,Merger_Tree_Construct)
    !% Initializes the merger tree reading module.
    use Input_Parameters
    use Galacticus_Error
    implicit none
    type(varying_string),          intent(in)    :: mergerTreeConstructMethod
    procedure(),          pointer, intent(inout) :: Merger_Tree_Construct
    integer                                      :: errorCode

    ! Check if our method is to be used.
    if (mergerTreeConstructMethod == 'read') then
       ! Assign pointer to our merger tree construction subroutine.
       Merger_Tree_Construct => Merger_Tree_Read_Do
       ! Read parameters for halo mass sampling.
       !@ <inputParameter>
       !@   <name>mergerTreeReadFileName</name>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The minimum mass of merger tree base halos to consider when building merger trees, in units of $M_\odot$.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('mergerTreeReadFileName',mergerTreeReadFileName)

       ! Read basic data from the merger tree file.
       ! Ensure HDF5 system is initialized.
       call IO_HDF5_Initialize
       ! Open the file.
       call h5fopen_f(char(mergerTreeReadFileName),H5F_ACC_RDONLY_F,mergerTreeFileID,errorCode)
       if (errorCode < 0) call Galacticus_Error_Report('Merger_Tree_Read_Initialize','failed to open input file')
       ! Open the merger trees group.
       call h5gopen_f(mergerTreeFileID,"mergerTrees",mergerTreesGroupID,errorCode)
       if (errorCode < 0) call Galacticus_Error_Report('Merger_Tree_Read_Initialize','failed to open mergerTrees group')

    end if
    return
  end subroutine Merger_Tree_Read_Initialize

  subroutine Merger_Tree_Read_Do(thisTree)
    !% Read a merger tree from file.
    use Tree_Nodes
    use Tree_Node_Methods
    use Cosmology_Functions
    use Galacticus_Error
    use String_Handling
    use H5Lt
    implicit none
    type(mergerTree),      intent(inout)             :: thisTree
    integer,               allocatable, dimension(:) :: nodeIndex,childNode,parentNode,siblingNode
    double precision,      allocatable, dimension(:) :: nodeMass,nodeRedshift,nodeTime
    type(treeNodeList),    allocatable, dimension(:) :: thisNodeList ! memoryManagementIgnore (force memory management system to ignore this)
    integer(kind=HID_T)                              :: mergerTreeGroupID
    integer(kind=HSIZE_T)                            :: datasetDimensions(1)
    integer(kind=SIZE_T)                             :: datasetTypeSize
    integer                                          :: errorCode,datasetTypeClass,nodeCount,iNode,jNode,treeIndex(1)
    double precision                                 :: volumeWeight(1)
    type(varying_string)                             :: datasetName,message
    logical                                          :: foundParent,foundChild,foundSibling

    !$omp critical(mergerTreeReadTree)
    ! Open the merger tree group.
    datasetName="mergerTree"
    datasetName=datasetName//nextTreeToRead
    call h5gopen_f(mergerTreesGroupID,char(datasetName),mergerTreeGroupID,errorCode)
    if (errorCode < 0)  then
       ! Could not find this group - assume we've read all of the trees.
       ! Close the merger trees group.
       call h5gclose_f(mergerTreesGroupID,errorCode)
       if (errorCode < 0) call Galacticus_Error_Report('Merger_Tree_Read_Do','failed to close mergerTrees group')
       ! Close the file.
       call h5fclose_f(mergerTreeFileID,errorCode)
       if (errorCode < 0) call Galacticus_Error_Report('Merger_Tree_Read_Do','failed to close input file')
       ! Uninitialize the HDF5 system.
       call IO_HDF5_Uninitialize
       return
    else
       ! Increment the tree to read index.
       nextTreeToRead=nextTreeToRead+1
    end if

    ! Get the dimensions of the dataset.
    call h5ltget_dataset_info_f(mergerTreeGroupID,"nodeIndex",datasetDimensions,datasetTypeClass,datasetTypeSize,errorCode)
    if (errorCode < 0) call Galacticus_Error_Report('Merger_Tree_Read_Do','could not get info on nodeIndex dataset')
    nodeCount=datasetDimensions(1)
    
    ! Allocate temporary arrays to hold the merger tree data.
    allocate(nodeIndex   (nodeCount))
    allocate(childNode   (nodeCount))
    allocate(parentNode  (nodeCount))
    allocate(siblingNode (nodeCount))
    allocate(nodeMass    (nodeCount))
    allocate(nodeRedshift(nodeCount))
    allocate(nodeTime    (nodeCount))
    allocate(thisNodeList(nodeCount))

    ! Read data from the file.
    ! nodeIndex
    call h5ltget_dataset_info_f(mergerTreeGroupID,"nodeIndex",datasetDimensions,datasetTypeClass,datasetTypeSize,errorCode)
    if (errorCode < 0) call Galacticus_Error_Report('Merger_Tree_Read_Do','could not get info on nodeIndex dataset')
    if (datasetDimensions(1) /= nodeCount) call Galacticus_Error_Report('Merger_Tree_Read_Do','nodeIndex dataset has incorrect dimension')
    call h5ltread_dataset_int_f(mergerTreeGroupID,"nodeIndex",nodeIndex,datasetDimensions,errorCode)
    if (errorCode < 0) call Galacticus_Error_Report('Merger_Tree_Read_Do','failed to read nodeIndex dataset')
    ! childNode
    call h5ltget_dataset_info_f(mergerTreeGroupID,"childNode",datasetDimensions,datasetTypeClass,datasetTypeSize,errorCode)
    if (errorCode < 0) call Galacticus_Error_Report('Merger_Tree_Read_Do','could not get info on childNode dataset')
    if (datasetDimensions(1) /= nodeCount) call Galacticus_Error_Report('Merger_Tree_Read_Do','childNode dataset has incorrect dimension')
    call h5ltread_dataset_int_f(mergerTreeGroupID,"childNode",childNode,datasetDimensions,errorCode)
    if (errorCode < 0) call Galacticus_Error_Report('Merger_Tree_Read_Do','failed to read childNode dataset')
    ! parentNode
    call h5ltget_dataset_info_f(mergerTreeGroupID,"parentNode",datasetDimensions,datasetTypeClass,datasetTypeSize,errorCode)
    if (errorCode < 0) call Galacticus_Error_Report('Merger_Tree_Read_Do','could not get info on parentNode dataset')
    if (datasetDimensions(1) /= nodeCount) call Galacticus_Error_Report('Merger_Tree_Read_Do','parentNode dataset has incorrect dimension')
    call h5ltread_dataset_int_f(mergerTreeGroupID,"parentNode",parentNode,datasetDimensions,errorCode)
    if (errorCode < 0) call Galacticus_Error_Report('Merger_Tree_Read_Do','failed to read parentNode dataset')
    ! siblingNode
    call h5ltget_dataset_info_f(mergerTreeGroupID,"siblingNode",datasetDimensions,datasetTypeClass,datasetTypeSize,errorCode)
    if (errorCode < 0) call Galacticus_Error_Report('Merger_Tree_Read_Do','could not get info on siblingNode dataset')
    if (datasetDimensions(1) /= nodeCount) call Galacticus_Error_Report('Merger_Tree_Read_Do','siblingNode dataset has incorrect dimension')
    call h5ltread_dataset_int_f(mergerTreeGroupID,"siblingNode",siblingNode,datasetDimensions,errorCode)
    if (errorCode < 0) call Galacticus_Error_Report('Merger_Tree_Read_Do','failed to read siblingNode dataset')
    ! nodeMass
    call h5ltget_dataset_info_f(mergerTreeGroupID,"nodeMass",datasetDimensions,datasetTypeClass,datasetTypeSize,errorCode)
    if (errorCode < 0) call Galacticus_Error_Report('Merger_Tree_Read_Do','could not get info on nodeMass dataset')
    if (datasetDimensions(1) /= nodeCount) call Galacticus_Error_Report('Merger_Tree_Read_Do','nodeMass dataset has incorrect dimension')
    call h5ltread_dataset_double_f(mergerTreeGroupID,"nodeMass",nodeMass,datasetDimensions,errorCode)
    if (errorCode < 0) call Galacticus_Error_Report('Merger_Tree_Read_Do','failed to read nodeMass dataset')
    ! nodeRedshift
    call h5ltget_dataset_info_f(mergerTreeGroupID,"nodeRedshift",datasetDimensions,datasetTypeClass,datasetTypeSize,errorCode)
    if (errorCode < 0) call Galacticus_Error_Report('Merger_Tree_Read_Do','could not get info on nodeRedshift dataset')
    if (datasetDimensions(1) /= nodeCount) call Galacticus_Error_Report('Merger_Tree_Read_Do','nodeRedshift dataset has incorrect dimension')
    call h5ltread_dataset_double_f(mergerTreeGroupID,"nodeRedshift",nodeRedshift,datasetDimensions,errorCode)
    if (errorCode < 0) call Galacticus_Error_Report('Merger_Tree_Read_Do','failed to read nodeRedshift dataset')
    ! treeIndex
    datasetDimensions(1)=1
    call h5ltread_dataset_int_f(mergerTreeGroupID,"treeIndex",treeIndex,datasetDimensions,errorCode)
    if (errorCode < 0) call Galacticus_Error_Report('Merger_Tree_Read_Do','failed to read treeIndex dataset')
    thisTree%index=treeIndex(1)
    ! volumeWeight
    datasetDimensions(1)=1
    call h5ltread_dataset_double_f(mergerTreeGroupID,"volumeWeight",volumeWeight,datasetDimensions,errorCode)
    if (errorCode < 0) call Galacticus_Error_Report('Merger_Tree_Read_Do','failed to read volumeWeight dataset')
    thisTree%volumeWeight=volumeWeight(1)

    ! Allocate tree node objects in our list.
    do iNode=1,nodeCount
       ! Create the node.
       call thisTree%createNode(thisNodeList(iNode)%node)
       ! Translate parent/child/sibling indices into array positions.
       foundParent =(parentNode (iNode) < 0)
       foundSibling=(siblingNode(iNode) < 0)
       foundChild  =(childNode  (iNode) < 0)
       do jNode=1,nodeCount
          if (nodeIndex(jNode) == parentNode (iNode)) then
             parentNode (iNode)=jNode
             foundParent       =.true.
          end if
          if (nodeIndex(jNode) == siblingNode(iNode)) then
             siblingNode(iNode)=jNode
             foundSibling      =.true.
          end if
          if (nodeIndex(jNode) == childNode  (iNode)) then
             childNode  (iNode)=jNode
             foundChild        =.true.
          end if
       end do
       if (.not.foundParent)  then
          message='failed to find parent node: '
          message=message//parentNode(iNode)//' at '//iNode
          call Galacticus_Error_Report('Merger_Tree_Read_Do',message)
       end if
       if (.not.foundSibling) then
          message='failed to find sibling node: '
          message=message//siblingNode(iNode)//' at '//iNode
          call Galacticus_Error_Report('Merger_Tree_Read_Do',message)
       end if
       if (.not.foundChild)   then
          message='failed to find child node: '
          message=message//childNode(iNode)//' at '//iNode
          call Galacticus_Error_Report('Merger_Tree_Read_Do',message)
       end if
    end do

    ! Set pointers and properties of the nodes.
    do iNode=1,nodeCount
       call thisNodeList(iNode)%node%indexSet(nodeIndex(iNode))
       if (parentNode(iNode) > 0) then
          thisNodeList(iNode)%node%parentNode  => thisNodeList(parentNode(iNode))%node
       else
          thisNodeList(iNode)%node%parentNode  => null()
          thisTree%baseNode => thisNodeList(iNode)%node
       end if
       if (siblingNode(iNode) > 0) then
          thisNodeList(iNode)%node%siblingNode => thisNodeList(siblingNode(iNode))%node
       else
          thisNodeList(iNode)%node%siblingNode => null()
       end if
       if (childNode(iNode) > 0) then
          thisNodeList(iNode)%node%childNode   => thisNodeList(childNode(iNode))%node
       else
          thisNodeList(iNode)%node%childNode   => null()
       end if
       call Tree_Node_Mass_Set(thisNodeList(iNode)%node,                                             nodeMass    (iNode)  )
       call Tree_Node_Time_Set(thisNodeList(iNode)%node,Cosmology_Age(Expansion_Factor_from_Redshift(nodeRedshift(iNode))))       
    end do

    ! Deallocate the temporary arrays.
    deallocate(nodeIndex,childNode,parentNode,siblingNode,nodeMass,nodeRedshift,nodeTime,thisNodeList)

    ! Close the merger tree group.
    call h5gclose_f(mergerTreeGroupID,errorCode)
    if (errorCode < 0) call Galacticus_Error_Report('Merger_Tree_Read_Do','failed to close mergerTree group')
    !$omp end critical(mergerTreeReadTree)

    return
  end subroutine Merger_Tree_Read_Do

end module Merger_Tree_Read
