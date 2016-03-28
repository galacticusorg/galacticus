!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
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

  !% Contains a module which implements a merger tree operator which outputs mass accretion
  !% histories.

  use IO_HDF5
  
  !# <mergerTreeOperator name="mergerTreeOperatorMassAccretionHistory">
  !#  <description>
  !#   A merger tree operator which outputs mass accretion histories.
  !# </description>
  !# </mergerTreeOperator>
  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorMassAccretionHistory
     !% A merger tree operator class which outputs mass accretion histories.
     private
     type(hdf5Object) :: outputGroup
   contains
     final     ::             massAccretionHistoryDestructor
     procedure :: operate  => massAccretionHistoryOperate
     procedure :: finalize => massAccretionHistoryFinalize
  end type mergerTreeOperatorMassAccretionHistory
  
  interface mergerTreeOperatorMassAccretionHistory
     !% Constructors for the mass accretion history merger tree operator class.
     module procedure massAccretionHistoryConstructorParameters
     module procedure massAccretionHistoryConstructorInternal
  end interface mergerTreeOperatorMassAccretionHistory

contains

  function massAccretionHistoryConstructorParameters(parameters)
    !% Constructor for the mass accretion history merger tree operator class which takes a
    !% parameter set as input.
    implicit none
    type(mergerTreeOperatorMassAccretionHistory)                :: massAccretionHistoryConstructorParameters
    type(inputParameters                       ), intent(in   ) :: parameters
    type(varying_string                        )                :: outputGroupName
    !# <inputParameterList label="allowedParameterNames" />
   
    call parameters%checkParameters(allowedParameterNames)
    !# <inputParameter>
    !#   <name>outputGroupName</name>
    !#   <source>parameters</source>
    !#   <defaultValue>var_str('massAccretionHistories')</defaultValue>
    !#   <description>The name of the \gls{hdf5} group to output mass accretion histories to.</description>
    !#   <type>string</type>
    !#   <cardinality>1</cardinality>
    !#   <group>output</group>
    !# </inputParameter>
    massAccretionHistoryConstructorParameters=massAccretionHistoryConstructorInternal(char(outputGroupName))
    return
  end function massAccretionHistoryConstructorParameters

  function massAccretionHistoryConstructorInternal(outputGroupName)
    !% Internal constructor for the mass accretion history merger tree operator class.
    use Galacticus_HDF5
    implicit none
    type     (mergerTreeOperatorMassAccretionHistory)                :: massAccretionHistoryConstructorInternal
    character(len=*                                 ), intent(in   ) :: outputGroupName
    
    ! Create the output group.
    !$omp critical (HDF5_Access)
    massAccretionHistoryConstructorInternal%outputGroup=galacticusOutputFile%openGroup(outputGroupName,'Mass accretion histories of main branches in merger trees.')
    !$omp end critical (HDF5_Access)
    return
  end function massAccretionHistoryConstructorInternal

  elemental subroutine massAccretionHistoryDestructor(self)
    !% Destructor for the mass accretion history merger tree operator function class.
    implicit none
    type(mergerTreeOperatorMassAccretionHistory), intent(inout) :: self

    ! Nothing to do.
    return
  end subroutine massAccretionHistoryDestructor

  subroutine massAccretionHistoryOperate(self,tree)
    !% Output the mass accretion history for a merger tree.
    use, intrinsic :: ISO_C_Binding
    use Galacticus_Nodes
    use Input_Parameters
    use Memory_Management
    use ISO_Varying_String
    use String_Handling
    use Numerical_Constants_Astronomical
    use Galacticus_Error
    implicit none
    class           (mergerTreeOperatorMassAccretionHistory), intent(inout)               :: self
    type            (mergerTree                            ), intent(inout), target       :: tree
    type            (treeNode                              )               , pointer      :: node
    integer         (kind=kind_int8                        ), allocatable  , dimension(:) :: accretionHistoryNodeIndex
    double precision                                        , allocatable  , dimension(:) :: accretionHistoryNodeMass , accretionHistoryNodeTime
    class           (nodeComponentBasic                    )               , pointer      :: basic
    type            (mergerTree                            )               , pointer      :: treeCurrent
    integer         (c_size_t                              )                              :: accretionHistoryCount
    type            (varying_string                        )                              :: groupName
    type            (hdf5Object                            )                              :: accretionDataset         , treeGroup

    ! Iterate over trees.
    treeCurrent => tree
    do while (associated(treeCurrent))
       ! Count up number of entries expected for accretion history.
       accretionHistoryCount =  0
       node                  => treeCurrent%baseNode
       do while (associated(node))
          accretionHistoryCount =  accretionHistoryCount           +1
          node                  => node                 %firstChild
       end do
       ! Allocate storage space.
       call Alloc_Array(accretionHistoryNodeIndex,[int(accretionHistoryCount)])
       call Alloc_Array(accretionHistoryNodeTime ,[int(accretionHistoryCount)])
       call Alloc_Array(accretionHistoryNodeMass ,[int(accretionHistoryCount)])
       ! Extract accretion history.
       accretionHistoryCount =  0
       node                  => treeCurrent%baseNode
       do while (associated(node))
          accretionHistoryCount                            =  accretionHistoryCount+1
          basic                                            => node%basic      ()
          accretionHistoryNodeIndex(accretionHistoryCount) =  node %index     ()
          accretionHistoryNodeTime (accretionHistoryCount) =  basic%time      ()
          accretionHistoryNodeMass (accretionHistoryCount) =  basic%mass      ()
          node                                             => node %firstChild
       end do
       ! Output to HDF5 file.
       groupName='mergerTree'
       groupName=groupName//treeCurrent%index
       !$omp critical (HDF5_Access)
       if (self%outputGroup%hasGroup(char(groupName))) call Galacticus_Error_Report('Merger_Tree_Mass_Accretion_History_Output','duplicate tree index detected - mass accretion history can not be output'//char(10)//'  HELP: This can happen if reading merger trees which contain multiple root nodes from file. To avoid this problem, force tree indices to be reset to the index of the root node by adding the following to your input parameter file:'//char(10)//'  <mergerTreeReadTreeIndexToRootNodeIndex value="true" />>')
       treeGroup=self%outputGroup%openGroup(char(groupName)                      ,'Mass accretion history for main branch of merger tree.'                                 )
       call treeGroup       %writeDataset  (accretionHistoryNodeIndex,'nodeIndex','Index of the node.'                                                                     )
       call treeGroup       %writeDataset  (accretionHistoryNodeTime ,'nodeTime' ,'Time at node [Gyr].'                                   ,datasetReturned=accretionDataset)
       call accretionDataset%writeAttribute(gigaYear                 ,'unitsInSI'                                                                                          )
       call accretionDataset%close         (                                                                                                                               )
       call treeGroup       %writeDataset  (accretionHistoryNodeMass ,'nodeMass' ,'Mass of the node [M⊙].'                                ,datasetReturned=accretionDataset)
       call accretionDataset%writeAttribute(massSolar,"unitsInSI")
       call accretionDataset%close         (                                                                                                                               )
       call treeGroup       %close         (                                                                                                                               )
       !$omp end critical (HDF5_Access)
       ! Deallocate storage space.
       call Dealloc_Array(accretionHistoryNodeIndex)
       call Dealloc_Array(accretionHistoryNodeTime )
       call Dealloc_Array(accretionHistoryNodeMass )
       ! Move to the next tree.
       treeCurrent => treeCurrent%nextTree
    end do    
    return
  end subroutine massAccretionHistoryOperate
  
  subroutine massAccretionHistoryFinalize(self)
    !% Close the mass accretion history group before closing the HDF5 file.
    implicit none
    class(mergerTreeOperatorMassAccretionHistory), intent(inout) :: self

    !$omp critical (HDF5_Access)
    if (self%outputGroup%isOpen()) call self%outputGroup%close()
    !$omp end critical (HDF5_Access)
    return
  end subroutine massAccretionHistoryFinalize
