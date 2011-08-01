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

  subroutine Merger_Trees_Render_Dump(thisTree)
    !% Dumps information on merger tree structure useful for rendering 3D views of merger trees.
    use Merger_Trees
    use Tree_Nodes
    use Dark_Matter_Halo_Scales
    use Cosmology_Functions
    use File_Utilities
    use IO_HDF5
    use Numerical_Constants_Astronomical
    use Memory_Management
    implicit none
    type(mergerTree), intent(inout)               :: thisTree
    type(treeNode),   pointer                     :: thisNode
    integer,          dimension(:  ), allocatable :: nodeIndex,parentIndex,childIndex
    double precision, dimension(:  ), allocatable :: time,expansionFactor,radiusVirial
    double precision, dimension(:,:), allocatable :: position
    integer                                       :: nodesInTree,iNode
    character(len=39)                             :: fileName
    type(hdf5Object)                              :: fileObject,treeDataset

    ! Reset output incremental counter if this tree is not the same as the previous one.
    if (thisTree%index /= treeIndexPrevious) then
       treeIndexPrevious=thisTree%index
       outputCounter    =-1
    end if

    ! Increment the output counter.
    outputCounter=outputCounter+1

    ! Construct a file name for the output.
    write (fileName,'(a7,i16.16,a1,i10.10,a5)') "render:",thisTree%index,":",outputCounter,".hdf5"

    ! Count the number of nodes in the tree.
    nodesInTree=0
    thisNode => thisTree%baseNode
    do while (associated(thisNode))
       nodesInTree=nodesInTree+1
       call thisNode%walkTreeWithSatellites(thisNode)
    end do

    ! Allocate arrays for temporary storage.
    call Alloc_Array(nodeIndex      ,[  nodesInTree])
    call Alloc_Array(parentIndex    ,[  nodesInTree])
    call Alloc_Array(childIndex     ,[  nodesInTree])
    call Alloc_Array(time           ,[  nodesInTree])
    call Alloc_Array(expansionFactor,[  nodesInTree])
    call Alloc_Array(radiusVirial   ,[  nodesInTree])
    call Alloc_Array(position       ,[3,nodesInTree])

    ! Populate arrays with data.
    thisNode => thisTree%baseNode
    iNode=0
    do while (associated(thisNode))
       iNode=iNode+1
       nodeIndex      (iNode)=thisNode           %index()
       parentIndex    (iNode)=thisNode%parentNode%index()
       childIndex     (iNode)=thisNode%childNode %index()
       time           (iNode)=                               Tree_Node_Time(thisNode)
       expansionFactor(iNode)=Expansion_Factor              (Tree_Node_Time(thisNode))
       radiusVirial   (iNode)=Dark_Matter_Halo_Virial_Radius(               thisNode )
       call Tree_Node_Position(thisNode,position(:,iNode))
       call thisNode%walkTreeWithSatellites(thisNode)
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
    call Dealloc_Array(nodeIndex      )
    call Dealloc_Array(parentIndex    )
    call Dealloc_Array(childIndex     )
    call Dealloc_Array(time           )
    call Dealloc_Array(expansionFactor)
    call Dealloc_Array(radiusVirial   )
    call Dealloc_Array(position       )

    return
  end subroutine Merger_Trees_Render_Dump

end module Merger_Trees_Render
