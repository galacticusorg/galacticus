!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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

!!{
Contains a module which implements dumping of the structure of a merger tree to a file for plotting with \href{http://www.graphviz.org/}{\normalfont \scshape dot}.
!!}

module Merger_Trees_Dump
  !!{
  Implements dumping of the structure of a merger tree to a file for plotting with \href{http://www.graphviz.org/}{\normalfont \scshape dot}.
  !!}
  use, intrinsic :: ISO_C_Binding, only : c_size_t
  use            :: Kind_Numbers , only : kind_int8
  implicit none
  private
  public :: Merger_Tree_Dump

  ! Internal indices used for labeling outputs.
  integer(kind=c_size_t ) :: outputCount
  integer(kind=kind_int8) :: treeIndexPrevious=-1
  !$omp threadprivate(outputCount,treeIndexPrevious)
contains

  subroutine Merger_Tree_Dump(tree,highlightNodes,backgroundColor,nodeColor,edgeColor,highlightColor,nodeStyle&
       &,highlightStyle ,edgeStyle ,labelNodes,labelUnique,scaleNodesByLogMass,edgeLengthsToTimes,useNodeLabels,timeRange,path)
    !!{
    Dumps the tree structure to a file in a format suitable for processing with \href{http://www.graphviz.org/}{\normalfont \scshape dot}. Nodes
    are shown as circles if isolated or rectangles if satellites. Isolated nodes are connected to their descendant halo, while
    satellites are connected (by red lines) to their host halo. Optionally, a list of node indices to highlight can be
    specified.
    !!}
    use :: Galacticus_Nodes   , only : mergerTree              , nodeComponentBasic, treeNode
    use :: ISO_Varying_String , only : assignment(=)           , char              , operator(//)  , trim, &
          &                            varying_string
    use :: Merger_Tree_Walkers, only : mergerTreeWalkerAllNodes
    use :: Nodes_Labels       , only : nodeLabelNames          , nodeLabelIsPresent, nodeLabelCount
    implicit none
    type            (mergerTree              )              , intent(in   )                    :: tree
    integer         (kind=kind_int8          ), dimension(:), intent(in   ), optional          :: highlightNodes
    character       (len=*                   )              , intent(in   ), optional          :: backgroundColor          , edgeColor           , &
         &                                                                                        edgeStyle                , highlightColor      , &
         &                                                                                        highlightStyle           , nodeColor           , &
         &                                                                                        nodeStyle
    logical                                                 , intent(in   ), optional          :: edgeLengthsToTimes       , labelNodes          , &
         &                                                                                        scaleNodesByLogMass      , labelUnique         , &
         &                                                                                        useNodeLabels
    double precision                          , dimension(2), intent(in   ), optional          :: timeRange
    type            (varying_string          )              , intent(in   ), optional          :: path
    type            (varying_string          ), dimension(:), allocatable                      :: labelNames
    type            (treeNode                )                                       , pointer :: node
    class           (nodeComponentBasic      )                                       , pointer :: basicParent              , basic
    type            (mergerTreeWalkerAllNodes)                                                 :: treeWalker
    logical                                                                                    :: edgeLengthsToTimesActual , labelNodesActual    , &
         &                                                                                        scaleNodesByLogMassActual, labelUniqueActual   , &
         &                                                                                        useNodeLabelsActual      , first
    integer                                                                                    :: fileUnit                 , i
    double precision                                                                           :: nodeMass                 , nodeMassMaximum     , &
         &                                                                                        nodeMassMinimum          , timeDifference      , &
         &                                                                                        timeMaximum              , timeMinimum
    character       (len=  16                )                                                 :: label
    character       (len=  20                )                                                 :: backgroundColorActual    , color               , &
         &                                                                                        edgeColorActual          , edgeStyleActual     , &
         &                                                                                        highlightColorActual     , highlightStyleActual, &
         &                                                                                        nodeColorActual          , nodeStyleActual     , &
         &                                                                                        outputCountFormatted     , style               , &
         &                                                                                        treeIndexFormatted
    character       (len=1024                )                                                 :: fileName
    type            (varying_string          )                                                 :: fullFileName             , label_

    ! If the tree index differs from the previous one, then reset the output count.
    if (tree%index /= treeIndexPrevious) then
       treeIndexPrevious=tree%index
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
    if (present(useNodeLabels      )) then
       useNodeLabelsActual      =useNodeLabels
    else
      useNodeLabelsActual       =.false.
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
    ! Get node label names.
    if (useNodeLabelsActual) call nodeLabelNames(labelNames)
    ! If sizes are to be scaled by mass, find the range of masses in the tree. If edges are set to time intervals, find minimum
    ! and maximum times in tree.
    if (scaleNodesByLogMassActual.or.edgeLengthsToTimesActual) then
       basic           =>                          tree%nodeBase%basic()
       nodeMassMinimum =                           basic        %mass ()
       nodeMassMaximum =                           basic        %mass ()
       timeMinimum     =                           basic        %time ()
       timeMaximum     =                           basic        %time ()
       treeWalker      =  mergerTreeWalkerAllNodes(tree                 )
       do while (treeWalker%next(node))
          basic => node%basic()
          if (basic%mass() < nodeMassMinimum) nodeMassMinimum=basic%mass()
          if (basic%mass() > nodeMassMaximum) nodeMassMaximum=basic%mass()
          if (basic%time() < timeMinimum    ) timeMinimum    =basic%time()
          if (basic%time() > timeMaximum    ) timeMaximum    =basic%time()
       end do
       nodeMassMinimum=log(nodeMassMinimum)
       nodeMassMaximum=log(nodeMassMaximum)
       timeMinimum    =log(timeMinimum    )
       timeMaximum    =log(timeMaximum    )
    else
       nodeMassMinimum=0.0d0
       nodeMassMaximum=0.0d0
       timeMinimum    =0.0d0
       timeMaximum    =0.0d0
    end if
    ! Open an output file and write the GraphViz opening.
    write (treeIndexFormatted  ,'(i16)') tree%index
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
    treeWalker=mergerTreeWalkerAllNodes(tree)
    do while (treeWalker%next(node))
       ! Get the basic component.
       basic => node%basic()
       ! Skip the node if it lies outside of the specified time range.
       if     (                                             &
            &   .not.present(timeRange)                     &
            &  .or.                                         &
            &   (                                           &
            &     basic%time() >= timeRange(1) &
            &    .and.                                      &
            &     basic%time() <= timeRange(2) &
            &   )                                           &
            & ) then
          ! Write each node, setting the node shape to a box for subhalos and a circle for halos. Node label consists of the node
          ! index plus the redshift, separated by a colon.
          ! Determine node color.
          if (present(highlightNodes)) then
             if     (                                                                           &
                  &   (     labelUniqueActual .and. any(highlightNodes == node%uniqueID())) &
                  &  .or.                                                                       &
                  &   (.not.labelUniqueActual .and. any(highlightNodes == node%index   ())) &
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
          if      (useNodeLabelsActual) then
             label_=""
             first=.true.
             do i=1,nodeLabelCount()
                if (nodeLabelIsPresent(i,node)) then
                   if (.not.first) label_=label_//"; "
                   label_=label_//labelNames(i)
                   first=.false.
                end if
             end do
          else if (labelUniqueActual  ) then
             write (      label,'(i16.16)') node       %uniqueID()
             label_=trim(label)
          else
             write (      label,'(i16.16)') node       %index   ()
             label_=trim(label)
          end if
          ! Create node.
          if (node%isSatellite()) then
             write (fileUnit,'(a,i16.16,a,a,a,a,a)') '"',node%uniqueID(),'" [shape=box   , color=',trim(color),', style=',trim(style),'];'
             ! Make a link to the hosting node.
             write (fileUnit,'(a,i16.16,a,i16.16,a,a,a)') '"',node%uniqueID(),'" -> "',node%parent%uniqueID(),'" [color=red, style=',trim(edgeStyleActual),'];'
          else
             write (fileUnit,'(a,i16.16,a,a,a,a,a)') '"',node%uniqueID(),'" [shape=circle, color=',trim(color),', style=',trim(style),'];'
             ! Make a link to the descendant node.
             if (associated(node%parent)) then
                write (fileUnit,'(a,i16.16,a,i16.16,a,a,a,a,$)') '"',node%uniqueID(),'" -> "',node%parent%uniqueID(),'" [color="',trim(edgeColorActual),'", style=',trim(edgeStyleActual)
                if (edgeLengthsToTimesActual) then
                   basicParent => node%parent%basic()
                   timeDifference=1000.0d0*log10(basicParent%time()/basic%time())/(timeMaximum-timeMinimum)
                   write (fileUnit,'(a,f10.6,$)') ', minlen=',timeDifference
                end if
                write (fileUnit,'(a)') '];'
             end if
          end if
          ! Set size of node if requested.
          if (scaleNodesByLogMassActual.and.nodeMassMaximum > nodeMassMinimum) then
             nodeMass=10.0d0*(log(basic%mass())-nodeMassMinimum)/(nodeMassMaximum-nodeMassMinimum)+1.0d0
             write (fileUnit,'(a,i16.16,a,f10.6,a)') '"',node%uniqueID(),'" [width=',nodeMass,'];'
          end if
          ! Remove label if requested.
          if (.not.labelNodesActual) then
             write (fileUnit,'(a,i16.16,a)'    ) '"',node%uniqueID(),'" [label=""];'
          else
             write (fileUnit,'(a,i16.16,a,a,a)') '"',node%uniqueID(),'" [label="',char(label_),'"];'
          end if
       end if
    end do
    ! Close the file.
    write (fileUnit,*) '}'
    close(fileUnit)
    return
  end subroutine Merger_Tree_Dump

end module Merger_Trees_Dump
