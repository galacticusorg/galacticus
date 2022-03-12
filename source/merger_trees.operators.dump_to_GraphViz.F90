!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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
Contains a module which implements a dump to \gls{graphviz} operator on merger trees.
!!}

  !![
  <mergerTreeOperator name="mergerTreeOperatorDumpToGraphViz">
   <description>
      A merger tree operator class which dumps the full structure of each merger tree to a file using the \gls{graphviz}
      format. All trees with root node basic mass between {\normalfont \ttfamily [massMinimum]} and {\normalfont \ttfamily
      [massMaximum]} will be dumped to a file named ``{\normalfont \ttfamily mergerTreeDump:\textless
      treeIndex\textgreater:1.gv}'' in the directory specified by {\normalfont \ttfamily [path]}. If {\normalfont \ttfamily
      [scaleNodesByLogMass]}$=${\normalfont \ttfamily true} then the size of each \gls{graphviz} node is scaled in proportion to
      the logarithm of the halo mass. If {\normalfont \ttfamily [edgeLengthsToTimes]}$=${\normalfont \ttfamily true} then the
      lengths of edges in the \gls{graphviz} graph are scaled in proportion to the time difference between the connected nodes.
   </description>
  </mergerTreeOperator>
  !!]
  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorDumpToGraphViz
     !!{
     A dump to \gls{graphviz} merger tree operator class.
     !!}
     private
     type            (varying_string) :: path
     double precision                 :: massMinimum        , massMaximum
     logical                          :: scaleNodesByLogMass, edgeLengthsToTimes
   contains
     procedure :: operatePreEvolution => dumpToGraphVizOperatePreEvolution
  end type mergerTreeOperatorDumpToGraphViz

  interface mergerTreeOperatorDumpToGraphViz
     !!{
     Constructors for the dump to \gls{graphviz} merger tree operator class.
     !!}
     module procedure dumpToGraphVizConstructorParameters
     module procedure dumpToGraphVizConstructorInternal
  end interface mergerTreeOperatorDumpToGraphViz

contains

  function dumpToGraphVizConstructorParameters(parameters)
    !!{
    Constructor for the dump-to-\gls{graphviz} merger tree operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(mergerTreeOperatorDumpToGraphViz)                :: dumpToGraphVizConstructorParameters
    type(inputParameters                 ), intent(inout) :: parameters

    !![
    <inputParameter>
      <name>path</name>
      <defaultValue>var_str('.')</defaultValue>
      <source>parameters</source>
      <variable>dumpToGraphVizConstructorParameters%path</variable>
      <description>Specifies the directory to which merger tree structure should be dumped.</description>
    </inputParameter>
    <inputParameter>
      <name>massMinimum</name>
      <defaultValue>0.0d0</defaultValue>
      <source>parameters</source>
      <variable>dumpToGraphVizConstructorParameters%massMinimum</variable>
      <description>Specifies the minimum root mass for which merger tree structure should be dumped.</description>
    </inputParameter>
    <inputParameter>
      <name>massMaximum</name>
      <defaultValue>huge(0.0d0)</defaultValue>
      <source>parameters</source>
      <variable>dumpToGraphVizConstructorParameters%massMaximum</variable>
      <description>Specifies the minimum root mass for which merger tree structure should be dumped.</description>
    </inputParameter>
    <inputParameter>
      <name>scaleNodesByLogMass</name>
      <defaultValue>.true.</defaultValue>
      <source>parameters</source>
      <variable>dumpToGraphVizConstructorParameters%scaleNodesByLogMass</variable>
      <description>Specifies whether or not node sizes should be scaled by the logarithm of their mass.</description>
    </inputParameter>
    <inputParameter>
      <name>edgeLengthsToTimes</name>
      <defaultValue>.true.</defaultValue>
      <source>parameters</source>
      <variable>dumpToGraphVizConstructorParameters%edgeLengthsToTimes</variable>
      <description>Specifies whether or not the lengths of edges in the graph should be scaled to time differences between nodes.</description>
    </inputParameter>
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function dumpToGraphVizConstructorParameters

  function dumpToGraphVizConstructorInternal(path,massMinimum,massMaximum)
    !!{
    Internal constructor for the dump-to-\gls{graphviz} merger tree operator class.
    !!}
    implicit none
    type            (mergerTreeOperatorDumpToGraphViz)                :: dumpToGraphVizConstructorInternal
    type            (varying_string                  ), intent(in   ) :: path
    double precision                                  , intent(in   ) :: massMinimum                      , massMaximum

    dumpToGraphVizConstructorInternal%path       =path
    dumpToGraphVizConstructorInternal%massMinimum=massMinimum
    dumpToGraphVizConstructorInternal%massMaximum=massMaximum
    return
  end function dumpToGraphVizConstructorInternal

  subroutine dumpToGraphVizOperatePreEvolution(self,tree)
    !!{
    Output the structure of {\normalfont \ttfamily tree}.
    !!}
    use :: Galacticus_Nodes , only : mergerTree      , nodeComponentBasic
    use :: Merger_Trees_Dump, only : Merger_Tree_Dump
    implicit none
    class(mergerTreeOperatorDumpToGraphViz), intent(inout), target :: self
    type (mergerTree                      ), intent(inout), target :: tree
    type (mergerTree                      ), pointer               :: treeCurrent
    class(nodeComponentBasic              ), pointer               :: basicBase

    ! Iterate over trees.
    treeCurrent => tree
    do while (associated(treeCurrent))
       ! Dump the tree.
       basicBase => treeCurrent%nodeBase%basic()
       if     (                                                                   &
            &   basicBase%mass() >= self%massMinimum                              &
            &  .and.                                                              &
            &   basicBase%mass() <  self%massMaximum                              &
            & )                                                                   &
            & call Merger_Tree_Dump(                                              &
            &                       treeCurrent                                 , &
            &                       scaleNodesByLogMass=self%scaleNodesByLogMass, &
            &                       edgeLengthsToTimes =self%edgeLengthsToTimes , &
            &                       nodeStyle          ='solid'                 , &
            &                       path               =self%path                 &
            &                      )
       ! Move to the next tree.
       treeCurrent => treeCurrent%nextTree
    end do
    return
  end subroutine dumpToGraphVizOperatePreEvolution
