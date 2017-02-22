!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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

!% Contains a module which dumps information on merger tree structure useful for rendering 3D views of merger trees.

module Merger_Trees_Render
  !% Implements dumping of information on merger tree structure useful for rendering 3D views of merger trees.
  use Kind_Numbers
  implicit none
  private
  public :: Merger_Trees_Render_Dump

  ! Counters for output file names.
  integer(kind=kind_int8) :: treeIndexPrevious=-1
  integer                 :: outputCounter

contains

  subroutine Merger_Trees_Render_Dump(tree)
    !% Dumps information on merger tree structure useful for rendering 3D views of merger trees.
    use Galacticus_Nodes
    use Dark_Matter_Halo_Scales
    use Cosmology_Functions
    use File_Utilities
    use IO_HDF5
    use Numerical_Constants_Astronomical
    use Memory_Management
    implicit none
    type            (mergerTree        ), intent(inout)                          :: tree
    type            (treeNode          )                               , pointer :: node
    class           (nodeComponentBasic)                               , pointer :: basic
    integer         (kind=kind_int8    ), allocatable  , dimension(:  )          :: childIndex     , nodeIndex   , parentIndex
    double precision                    , allocatable  , dimension(:  )          :: expansionFactor, radiusVirial, time
    double precision                    , allocatable  , dimension(:,:)          :: position
    integer                                                                      :: iNode          , nodesInTree
    character       (len=39            )                                         :: fileName
    type            (hdf5Object        )                                         :: fileObject     , treeDataset

    ! Reset output incremental counter if this tree is not the same as the previous one.
    if (tree%index /= treeIndexPrevious) then
       treeIndexPrevious=tree%index
       outputCounter    =-1
    end if

    ! Increment the output counter.
    outputCounter=outputCounter+1

    ! Construct a file name for the output.
    write (fileName,'(a7,i16.16,a1,i10.10,a5)') "render:",tree%index,":",outputCounter,".hdf5"

    ! Count the number of nodes in the tree.
    nodesInTree=0
    node => tree%baseNode
    do while (associated(node))
       nodesInTree=nodesInTree+1
       call node%walkTreeWithSatellites(node)
    end do

    ! Allocate arrays for temporary storage.
    call allocateArray(nodeIndex      ,[  nodesInTree])
    call allocateArray(parentIndex    ,[  nodesInTree])
    call allocateArray(childIndex     ,[  nodesInTree])
    call allocateArray(time           ,[  nodesInTree])
    call allocateArray(expansionFactor,[  nodesInTree])
    call allocateArray(radiusVirial   ,[  nodesInTree])
    call allocateArray(position       ,[3,nodesInTree])

    ! Populate arrays with data.
    node => tree%baseNode
    iNode=0
    do while (associated(node))
       iNode=iNode+1
       basic              => node           %basic()
       nodeIndex      (iNode) =  node           %index()
       parentIndex    (iNode) =  node%parent    %index()
       childIndex     (iNode) =  node%firstChild%index()
       time           (iNode) =                                 basic%time()
       expansionFactor(iNode) =  cosmologyFunctionsDefault%expansionFactor              (basic%time())
       radiusVirial   (iNode) =  darkMatterHaloScale_%virialRadius(node        )
       call Tree_Node_Position(node,position(:,iNode))
       call node%walkTreeWithSatellites(node)
    end do

    ! Open an HDF5 file.
    call fileObject%openFile(fileName,overWrite=.true.,objectsOverwritable=.true.)

    ! Write the datasets.
    call fileObject%writeDataset(nodeIndex      ,"nodeIndex"      ,"Node index []"                                  )
    call fileObject%writeDataset(parentIndex    ,"parentIndex"    ,"Parent index []"                                )
    call fileObject%writeDataset(childIndex     ,"childIndex"     ,"Child index []"                                 )
    call fileObject%writeDataset(expansionFactor,"expansionFactor","Expansion factor []"                            )

    call fileObject%writeDataset(time           ,"time"           ,"Time [Gyr]"         ,datasetReturned=treeDataset)
    call treeDataset%writeAttribute(gigaYear  ,"unitsInSI")
    call treeDataset%close()

    call fileObject%writeDataset(radiusVirial   ,"radiusVirial"   ,"Virial radius [Mpc]",datasetReturned=treeDataset)
    call treeDataset%writeAttribute(megaParsec,"unitsInSI")
    call treeDataset%close()

    call fileObject%writeDataset(position       ,"position"       ,"Position [Mpc]"     ,datasetReturned=treeDataset)
    call treeDataset%writeAttribute(megaParsec,"unitsInSI")
    call treeDataset%close()

    ! Close the output file.
    call fileObject%close()

    ! Deallocate temporary arrays.
    call deallocateArray(nodeIndex      )
    call deallocateArray(parentIndex    )
    call deallocateArray(childIndex     )
    call deallocateArray(time           )
    call deallocateArray(expansionFactor)
    call deallocateArray(radiusVirial   )
    call deallocateArray(position       )

    return
  end subroutine Merger_Trees_Render_Dump

end module Merger_Trees_Render
