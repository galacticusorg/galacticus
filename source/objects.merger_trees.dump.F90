!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements dumping of the structure of a merger tree to a file for plotting with \href{http://www.graphviz.org/}{\sc dot}.

module Merger_Trees_Dump
  !% Implements dumping of the structure of a merger tree to a file for plotting with \href{http://www.graphviz.org/}{\sc dot}.
  use Kind_Numbers
  implicit none
  private
  public :: Merger_Tree_Dump

  ! Internal indices used for labelling outputs.
  integer(kind=kind_int8) :: outputCount
  integer(kind=kind_int8) :: treeIndexPrevious=-1
  !$omp threadprivate(outputCount,treeIndexPrevious)

contains

  subroutine Merger_Tree_Dump(treeIndex,baseNode,highlightNodes,backgroundColor,nodeColor,edgeColor,highlightColor,nodeStyle&
       &,highlightStyle ,edgeStyle ,labelNodes,scaleNodesByLogMass,edgeLengthsToTimes,path)
    !% Dumps the tree structure to a file in a format suitable for processing with \href{http://www.graphviz.org/}{\sc dot}. Nodes
    !% are shown as circles if isolated or rectangles if satellites. Isolated nodes are connected to their descendent halo, while
    !% satellites are connected (by red lines) to their host halo. Optionally, a list of node indices to highlight can be
    !% specified.
    use Galacticus_Nodes
    use ISO_Varying_String
    implicit none
    integer(kind=kind_int8),   intent(in)                         :: treeIndex
    type(treeNode),            intent(in), pointer                :: baseNode
    integer(kind=kind_int8),   intent(in), dimension(:), optional :: highlightNodes
    character(len=*),          intent(in),               optional :: backgroundColor,nodeColor,highlightColor,edgeColor,nodeStyle&
         &,highlightStyle,edgeStyle
    logical,                   intent(in),               optional :: labelNodes,scaleNodesByLogMass,edgeLengthsToTimes
    type(varying_string),      intent(in),               optional :: path
    type(treeNode),                        pointer                :: thisNode
    class(nodeComponentBasic),             pointer                :: thisBasicComponent,parentBasicComponent
    logical                                                       :: labelNodesActual,scaleNodesByLogMassActual,edgeLengthsToTimesActual
    integer                                                       :: fileUnit
    double precision                                              :: timeDifference,nodeMassMinimum,nodeMassMaximum,nodeMass,timeMinimum,timeMaximum
    character(len=  20)                                           :: color,style,treeIndexFormatted,outputCountFormatted&
         &,backgroundColorActual,nodeColorActual,highlightColorActual,edgeColorActual,nodeStyleActual,highlightStyleActual,edgeStyleActual
    character(len=1024)                                           :: fileName
    type(varying_string)                                          :: fullFileName

    ! If the tree index differs from the previous one, then reset the output count.
    if (treeIndex /= treeIndexPrevious) then
       treeIndexPrevious=treeIndex
       outputCount=0
    end if

    ! Increment the count of the number of outputs.
    outputCount=outputCount+1
    
    ! Get optional arguments or set defaults.
    if (present(backgroundColor    )) then
       backgroundColorActual    =backgroundColor
    else
       backgroundColorActual    ='white'
    end if
    if (present(nodeColor          )) then
       nodeColorActual          =nodeColor
    else
       nodeColorActual          ='black'
    end if
    if (present(highlightColor     )) then
       highlightColorActual     =highlightColor
    else
       highlightColorActual     ='green'
    end if
    if (present(edgeColor          )) then
       edgeColorActual          =edgeColor
    else
       edgeColorActual          ='black'
    end if
    if (present(nodeStyle          )) then
       nodeStyleActual          =nodeStyle
    else
       nodeStyleActual          ='filled'
    end if
    if (present(highlightStyle     )) then
       highlightStyleActual     =highlightStyle
    else
       highlightStyleActual     ='filled'
    end if
    if (present(edgeStyle          )) then
       edgeStyleActual          =edgeStyle
    else
       edgeStyleActual          ='solid'
    end if
    if (present(labelNodes         )) then
       labelNodesActual         =labelNodes
    else
       labelNodesActual         =.true.
    end if
    if (present(scaleNodesByLogMass)) then
       scaleNodesByLogMassActual=scaleNodesByLogMass
    else
       scaleNodesByLogMassActual=.false.
    end if
    if (present(edgeLengthsToTimes )) then
       edgeLengthsToTimesActual =edgeLengthsToTimes
    else
       edgeLengthsToTimesActual =.false.
    end if

    ! If sizes are to be scaled by mass, find the range of masses in the tree. If edges are set to time intervals, find minimum
    ! and maximum times in tree.
    if (scaleNodesByLogMassActual.or.edgeLengthsToTimesActual) then
       thisNode           => baseNode
       thisBasicComponent => thisNode%basic()
       nodeMassMinimum=thisBasicComponent%mass()
       nodeMassMaximum=thisBasicComponent%mass()
       timeMinimum    =thisBasicComponent%time()
       timeMaximum    =thisBasicComponent%time()
       do while (associated(thisNode))
          thisBasicComponent => thisNode%basic()
          if (thisBasicComponent%mass() < nodeMassMinimum) nodeMassMinimum=thisBasicComponent%mass()
          if (thisBasicComponent%mass() > nodeMassMaximum) nodeMassMaximum=thisBasicComponent%mass()
          if (thisBasicComponent%time() < timeMinimum    ) timeMinimum    =thisBasicComponent%time()
          if (thisBasicComponent%time() > timeMaximum    ) timeMaximum    =thisBasicComponent%time()
          call thisNode%walkTreeWithSatellites(thisNode)
       end do
       nodeMassMinimum=log(nodeMassMinimum)
       nodeMassMaximum=log(nodeMassMaximum)
       timeMinimum    =log(timeMinimum)
       timeMaximum    =log(timeMaximum)
    end if

    ! Open an output file and write the GraphViz opening.
    write (treeIndexFormatted  ,'(i16)') treeIndex
    write (outputCountFormatted,'(i08)') outputCount
    write (fileName,'(a,a,a,a,a)') 'mergerTreeDump:',trim(adjustl(treeIndexFormatted)),':',trim(adjustl(outputCountFormatted)),'.gv'
    if (present(path)) then
       fullFileName=path//"/"//trim(fileName)
    else
       fullFileName=trim(fileName)
    end if
    open(newunit=fileUnit,file=char(fullFileName),status='unknown',form='formatted')
    write (fileUnit,'(a)'    ) 'digraph Tree {'
    write (fileUnit,'(a,a,a)') 'bgcolor=',trim(backgroundColorActual),';'
    write (fileUnit,'(a)'    ) 'size="8,11";'

    ! Loop over all nodes.
    thisNode => baseNode
    do while (associated(thisNode))
       ! Get the basic component.
       thisBasicComponent => thisNode%basic()
       ! Write each node, setting the node shape to a box for subhalos and a circle for halos. Node label consists of the node
       ! index plus the redshift, separated by a colon.
       ! Determine node color.
       if (present(highlightNodes)) then
          if (any(highlightNodes == thisNode%index())) then
             color=highlightColorActual
             style=highlightStyleActual
          else
             color=nodeColorActual
             style=nodeStyleActual
          end if
       else
          color=nodeColorActual
          style=nodeStyleActual
       end if
       if (thisNode%isSatellite()) then
          write (fileUnit,'(a,i16.16,a,a,a,a,a)') '"',thisNode%index(),'" [shape=box   , color=',trim(color),', style=',trim(style),'];'
          ! Make a link to the hosting node.
          write (fileUnit,'(a,i16.16,a,i16.16,a,a,a)') '"',thisNode%index(),'" -> "',thisNode%parent%index(),'" [color=red, style=',trim(edgeStyleActual),'];'
       else
          write (fileUnit,'(a,i16.16,a,a,a,a,a,f10.6,a)') '"',thisNode%index(),'" [shape=circle, color=',trim(color),', style=',trim(style),'];'
          ! Make a link to the descendent node.
          if (associated(thisNode%parent)) then
             write (fileUnit,'(a,i16.16,a,i16.16,a,a,a,a,$)') '"',thisNode%index(),'" -> "',thisNode%parent%index(),'" [color="',trim(edgeColorActual),'", style=',trim(edgeStyleActual)
             if (edgeLengthsToTimesActual) then
                parentBasicComponent => thisNode%parent%basic()
                timeDifference=1000.0d0*log10(parentBasicComponent%time()/thisBasicComponent%time())/(timeMaximum-timeMinimum)
                write (fileUnit,'(a,f10.6,$)') ', minlen=',timeDifference
             end if
             write (fileUnit,'(a)') '];'
          end if
       end if

       ! Set size of node if requested.
       if (scaleNodesByLogMassActual.and.nodeMassMaximum > nodeMassMinimum) then
          nodeMass=10.0d0*(log(thisBasicComponent%mass())-nodeMassMinimum)/(nodeMassMaximum-nodeMassMinimum)+1.0d0
          write (fileUnit,'(a,i16.16,a,f10.6,a)') '"',thisNode%index(),'" [width=',nodeMass,'];'
       end if

       ! Remove label is requested.
       if (.not.labelNodesActual) then
          write (fileUnit,'(a,i16.16,a)') '"',thisNode%index(),'" [label=""];'
       end if

       ! Walk the tree, including satellite nodes.
       call thisNode%walkTreeWithSatellites(thisNode)

    end do

    ! Close the file.
    write (fileUnit,*) '}'
    close(fileUnit)

    return
  end subroutine Merger_Tree_Dump
  
end module Merger_Trees_Dump
