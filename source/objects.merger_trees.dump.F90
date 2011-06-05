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


!% Contains a module which implements dumping of the structure of a merger tree to a file for plotting with \href{http://www.graphviz.org/}{\sc dot}.

module Merger_Trees_Dump
  !% Implements dumping of the structure of a merger tree to a file for plotting with \href{http://www.graphviz.org/}{\sc dot}.
  use Kind_Numbers
  private
  public :: Merger_Tree_Dump

  ! Internal indices used for labelling outputs.
  integer(kind=kind_int8) :: outputCount
  integer(kind=kind_int8) :: treeIndexPrevious=-1
  !$omp threadprivate(outputCount,treeIndexPrevious)

contains

  subroutine Merger_Tree_Dump(treeIndex,baseNode,highlightNodes)
    !% Dumps the tree structure to a file in a format suitable for processing with \href{http://www.graphviz.org/}{\sc dot}. Nodes
    !% are shown as circles if isolated or rectangles if satellites. Isolated nodes are connected (by white lines) to their
    !% descendent halo, while satellites are connected (by red lines) to their host halo. Optionally, a list of node indices to
    !% highlight can be specified---any such nodes will be filled in green.
    use File_Utilities
    use Tree_Nodes
    use Input_Parameters
    implicit none
    integer(kind=kind_int8), intent(in)                         :: treeIndex
    type(treeNode),          intent(in), pointer                :: baseNode
    integer(kind=kind_int8), intent(in), dimension(:), optional :: highlightNodes
    type(treeNode),                      pointer                :: thisNode
    integer                                                     :: fileUnit
    character(len=  20)                                         :: color,style,treeIndexFormatted,outputCountFormatted
    character(len=1024)                                         :: fileName

    ! If the tree index differs from the previous one, then reset the output count.
    if (treeIndex /= treeIndexPrevious) then
       treeIndexPrevious=treeIndex
       outputCount=0
    end if

    ! Increment the count of the number of outputs.
    outputCount=outputCount+1
    
    ! Open an output file and write the GraphViz opening.
    write (treeIndexFormatted,'(i8)') treeIndex
    write (outputCountFormatted,'(i8)') outputCount
    write (fileName,'(a,a,a,a,a)') 'mergerTreeDump:',trim(adjustl(treeIndexFormatted)),':',trim(adjustl(outputCountFormatted)),'.gv'
    fileUnit=File_Units_Get()
    open(fileUnit,file=fileName,status='unknown',form='formatted')
    write (fileUnit,*) 'digraph Tree {'
    write (fileUnit,*) 'bgcolor=black;'
    write (fileUnit,*) 'size="8,11";'

    ! Loop over all nodes.
    thisNode => baseNode
    do while (associated(thisNode))
       ! Write each node, setting the node shape to a box for subhalos and a circle for halos. Node label consists of the node
       ! index plus the redshift, separated by a colon.
       ! Determine node color.
       if (present(highlightNodes)) then
          if (any(highlightNodes == thisNode%index())) then
             color='green'
             style='filled'
          else
             color='white'
             style='filled'
          end if
       else
          color='white'
          style='filled'
       end if
       if (thisNode%isSatellite()) then
          write (fileUnit,'(a,i16.16,a,i16.16,a,f5.2,a,a,a,a,a)') '"',thisNode%index(),'" [shape=box   , label="',thisNode%index(),':',Tree_Node_Time(thisNode),'", color=',trim(color),', style=',trim(style),'];'
          ! Make a link to the hosting node.
          write (fileUnit,'(a,i16.16,a,i16.16,a)') '"',thisNode%index(),'" -> "',thisNode%parentNode%index(),'" [color=red];'
       else
          write (fileUnit,'(a,i16.16,a,i16.16,a,f5.2,a,a,a,a,a)') '"',thisNode%index(),'" [shape=circle, label="',thisNode%index(),':',Tree_Node_Time(thisNode),'", color=',trim(color),', style=',trim(style),'];'
          ! Make a link to the descendent node using a white line.
          if (associated(thisNode%parentNode)) write (fileUnit,'(a,i16.16,a,i16.16,a)') '"',thisNode%index(),'" -> "',thisNode%parentNode%index(),'" [color=white];'
       endif

       ! Walk the tree, including satellite nodes.
       call thisNode%walkTreeWithSatellites(thisNode)

    end do

    ! Close the file.
    write (fileUnit,*) '}'
    close(fileUnit)

    return
  end subroutine Merger_Tree_Dump
  
end module Merger_Trees_Dump
