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

!% Contains a module which implements dumping of the structure of a merger tree to a file for plotting with \href{http://www.graphviz.org/}{\normalfont \scshape dot}.

module Merger_Trees_Dump
  !% Implements dumping of the structure of a merger tree to a file for plotting with \href{http://www.graphviz.org/}{\normalfont \scshape dot}.
  use, intrinsic :: ISO_C_Binding
  use Kind_Numbers
  implicit none
  private
  public :: Merger_Tree_Dump

  ! Internal indices used for labelling outputs.
  integer(kind=c_size_t ) :: outputCount
  integer(kind=kind_int8) :: treeIndexPrevious=-1
  !$omp threadprivate(outputCount,treeIndexPrevious)
contains

  subroutine Merger_Tree_Dump(treeIndex,baseNode,highlightNodes,backgroundColor,nodeColor,edgeColor,highlightColor,nodeStyle&
       &,highlightStyle ,edgeStyle ,labelNodes,labelUnique,scaleNodesByLogMass,edgeLengthsToTimes,timeRange,path)
    !% Dumps the tree structure to a file in a format suitable for processing with \href{http://www.graphviz.org/}{\normalfont \scshape dot}. Nodes
    !% are shown as circles if isolated or rectangles if satellites. Isolated nodes are connected to their descendent halo, while
    !% satellites are connected (by red lines) to their host halo. Optionally, a list of node indices to highlight can be
    !% specified.
    use Galacticus_Nodes
    use ISO_Varying_String
    implicit none
    integer         (kind=kind_int8    )              , intent(in   )                    :: treeIndex
    type            (treeNode          )              , intent(in   )          , pointer :: baseNode
    integer         (kind=kind_int8    ), dimension(:), intent(in   ), optional          :: highlightNodes
    character       (len=*             )              , intent(in   ), optional          :: backgroundColor          , edgeColor           , &
         &                                                                                  edgeStyle                , highlightColor      , &
         &                                                                                  highlightStyle           , nodeColor           , &
         &                                                                                  nodeStyle
    logical                                           , intent(in   ), optional          :: edgeLengthsToTimes       , labelNodes          , &
         &                                                                                  scaleNodesByLogMass      , labelUnique
    double precision                    , dimension(2), intent(in   ), optional          :: timeRange
    type            (varying_string    )              , intent(in   ), optional          :: path
    type            (treeNode          )                                       , pointer :: thisNode
    class           (nodeComponentBasic)                                       , pointer :: parentBasicComponent     , thisBasicComponent
    logical                                                                              :: edgeLengthsToTimesActual , labelNodesActual    , &
         &                                                                                  scaleNodesByLogMassActual, labelUniqueActual
    integer                                                                              :: fileUnit
    double precision                                                                     :: nodeMass                 , nodeMassMaximum     , &
         &                                                                                  nodeMassMinimum          , timeDifference      , &
         &                                                                                  timeMaximum              , timeMinimum
    character       (len=  16          )                                                 :: label
    character       (len=  20          )                                                 :: backgroundColorActual    , color               , &
         &                                                                                  edgeColorActual          , edgeStyleActual     , &
         &                                                                                  highlightColorActual     , highlightStyleActual, &
         &                                                                                  nodeColorActual          , nodeStyleActual     , &
         &                                                                                  outputCountFormatted     , style               , &
         &                                                                                  treeIndexFormatted
    character       (len=1024          )                                                 :: fileName
    type            (varying_string    )                                                 :: fullFileName

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
    if (present(labelUnique        )) then
       labelUniqueActual        =labelUnique
    else
       labelUniqueActual        =.false.
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
          thisNode => thisNode%walkTreeWithSatellites()
       end do
       nodeMassMinimum=log(nodeMassMinimum)
       nodeMassMaximum=log(nodeMassMaximum)
       timeMinimum    =log(timeMinimum)
       timeMaximum    =log(timeMaximum)
    else
       nodeMassMinimum=0.0d0
       nodeMassMaximum=0.0d0
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
       ! Skip the node if it lies outside of the specified time range.
       if     (                                             &
            &   .not.present(timeRange)                     &
            &  .or.                                         &
            &   (                                           &
            &     thisBasicComponent%time() >= timeRange(1) &
            &    .and.                                      &
            &     thisBasicComponent%time() <= timeRange(2) &
            &   )                                           &
            & ) then
          ! Write each node, setting the node shape to a box for subhalos and a circle for halos. Node label consists of the node
          ! index plus the redshift, separated by a colon.
          ! Determine node color.
          if (present(highlightNodes)) then
             if     (                                                                           &
                  &   (     labelUniqueActual .and. any(highlightNodes == thisNode%uniqueID())) &
                  &  .or.                                                                       &
                  &   (.not.labelUniqueActual .and. any(highlightNodes == thisNode%index   ())) &
                  & ) then
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
          ! Create node labels.
          if (labelUniqueActual) then
             write (      label,'(i16.16)') thisNode       %uniqueID()
          else
             write (      label,'(i16.16)') thisNode       %index   ()
          end if
          ! Create node.
          if (thisNode%isSatellite()) then
             write (fileUnit,'(a,i16.16,a,a,a,a,a)') '"',thisNode%uniqueID(),'" [shape=box   , color=',trim(color),', style=',trim(style),'];'
             ! Make a link to the hosting node.
             write (fileUnit,'(a,i16.16,a,i16.16,a,a,a)') '"',thisNode%uniqueID(),'" -> "',thisNode%parent%uniqueID(),'" [color=red, style=',trim(edgeStyleActual),'];'
          else
             write (fileUnit,'(a,i16.16,a,a,a,a,a)') '"',thisNode%uniqueID(),'" [shape=circle, color=',trim(color),', style=',trim(style),'];'
             ! Make a link to the descendent node.
             if (associated(thisNode%parent)) then
                write (fileUnit,'(a,i16.16,a,i16.16,a,a,a,a,$)') '"',thisNode%uniqueID(),'" -> "',thisNode%parent%uniqueID(),'" [color="',trim(edgeColorActual),'", style=',trim(edgeStyleActual)
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
             write (fileUnit,'(a,i16.16,a,f10.6,a)') '"',thisNode%uniqueID(),'" [width=',nodeMass,'];'
          end if

          ! Remove label if requested.
          if (.not.labelNodesActual) then
             write (fileUnit,'(a,i16.16,a)'    ) '"',thisNode%uniqueID(),'" [label=""];'
          else
             write (fileUnit,'(a,i16.16,a,a,a)') '"',thisNode%uniqueID(),'" [label="',label,'"];'
          end if

       end if

       ! Walk the tree, including satellite nodes.
       thisNode => thisNode%walkTreeWithSatellites()

    end do

    ! Close the file.
    write (fileUnit,*) '}'
    close(fileUnit)

    return
  end subroutine Merger_Tree_Dump

end module Merger_Trees_Dump
