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
Implements a dump to \gls{graphviz} operator on merger trees.
!!}

  use :: Cosmology_Functions, only : cosmologyFunctionsClass
  
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
     class           (cosmologyFunctionsClass), pointer :: cosmologyFunctions_ => null()
     type            (varying_string         )          :: path
     double precision                                   :: massMinimum                  , massMaximum       , &
          &                                                redshiftMinimum              , redshiftMaximum   , &
          &                                                timeMinimum                  , timeMaximum
     logical                                            :: scaleNodesByLogMass          , edgeLengthsToTimes, &
          &                                                useNodeLabels
   contains
     final     ::                        dumpToGraphVizDestructor
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

  function dumpToGraphVizConstructorParameters(parameters) result(self)
    !!{
    Constructor for the dump-to-\gls{graphviz} merger tree operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (mergerTreeOperatorDumpToGraphViz)                :: self
    type            (inputParameters                 ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass         ), pointer       :: cosmologyFunctions_
    type            (varying_string                  )                :: path
    double precision                                                  :: massMinimum        , massMaximum       , &
         &                                                               redshiftMinimum    , redshiftMaximum
    logical                                                           :: scaleNodesByLogMass, edgeLengthsToTimes, &
         &                                                               useNodeLabels
    
    !![
    <inputParameter>
      <name>path</name>
      <defaultValue>var_str('.')</defaultValue>
      <source>parameters</source>
      <description>Specifies the directory to which merger tree structure should be dumped.</description>
    </inputParameter>
    <inputParameter>
      <name>massMinimum</name>
      <defaultValue>0.0d0</defaultValue>
      <source>parameters</source>
      <description>Specifies the minimum root mass for which merger tree structure should be dumped.</description>
    </inputParameter>
    <inputParameter>
      <name>massMaximum</name>
      <defaultValue>huge(0.0d0)</defaultValue>
      <source>parameters</source>
      <description>Specifies the minimum root mass for which merger tree structure should be dumped.</description>
    </inputParameter>
    <inputParameter>
      <name>redshiftMinimum</name>
      <defaultValue>0.0d0</defaultValue>
      <source>parameters</source>
      <description>Specifies the minimum redshift for which merger tree structure should be dumped.</description>
    </inputParameter>
    <inputParameter>
      <name>redshiftMaximum</name>
      <defaultValue>100.0d0</defaultValue>
      <source>parameters</source>
      <description>Specifies the maximum redshift for which merger tree structure should be dumped.</description>
    </inputParameter>
    <inputParameter>
      <name>scaleNodesByLogMass</name>
      <defaultValue>.true.</defaultValue>
      <source>parameters</source>
      <description>Specifies whether or not node sizes should be scaled by the logarithm of their mass.</description>
    </inputParameter>
    <inputParameter>
      <name>edgeLengthsToTimes</name>
      <defaultValue>.true.</defaultValue>
      <source>parameters</source>
      <description>Specifies whether or not the lengths of edges in the graph should be scaled to time differences between nodes.</description>
    </inputParameter>
    <inputParameter>
      <name>useNodeLabels</name>
      <defaultValue>.false.</defaultValue>
      <source>parameters</source>
      <description>If true, label nodes in the graph with any node labels. Otherwise, label using node IDs.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    !!]
    self=mergerTreeOperatorDumpToGraphViz(                                                                                                  &
         &                                path                                                                                            , &
         &                                massMinimum                                                                                     , &
         &                                massMaximum                                                                                     , &
         &                                cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshiftMaximum)), &
         &                                cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshiftMinimum)), &
         &                                scaleNodesByLogMass                                                                             , &
         &                                edgeLengthsToTimes                                                                              , &
         &                                useNodeLabels                                                                                   , &
         &                                cosmologyFunctions_                                                                               &
         &                               )
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    !!]
    return
  end function dumpToGraphVizConstructorParameters

  function dumpToGraphVizConstructorInternal(path,massMinimum,massMaximum,timeMinimum,timeMaximum,scaleNodesByLogMass,edgeLengthsToTimes,useNodeLabels,cosmologyFunctions_) result(self)
    !!{
    Internal constructor for the dump-to-\gls{graphviz} merger tree operator class.
    !!}
    implicit none
    type            (mergerTreeOperatorDumpToGraphViz)                        :: self
    type            (varying_string                  ), intent(in   )         :: path
    double precision                                  , intent(in   )         :: massMinimum        , massMaximum       , &
         &                                                                       timeMinimum        , timeMaximum
    logical                                           , intent(in   )         :: scaleNodesByLogMass, edgeLengthsToTimes, &
         &                                                                       useNodeLabels
    class           (cosmologyFunctionsClass         ), intent(in   ), target :: cosmologyFunctions_
    !![
    <constructorAssign variables="path, massMinimum, massMaximum, timeMinimum, timeMaximum, scaleNodesByLogMass, edgeLengthsToTimes, useNodeLabels, *cosmologyFunctions_"/>
    !!]

    self%redshiftMinimum=self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(timeMaximum))
    self%redshiftMaximum=self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(timeMinimum))
    return
  end function dumpToGraphVizConstructorInternal

  subroutine dumpToGraphVizDestructor(self)
    !!{
    Destructor for the dump-to-\gls{graphviz} merger tree operator class.
    !!}
    implicit none
    type(mergerTreeOperatorDumpToGraphViz), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"/>
    !!]
    return
  end subroutine dumpToGraphVizDestructor

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
       if     (                                                                    &
            &   basicBase%mass() >= self%massMinimum                               &
            &  .and.                                                               &
            &   basicBase%mass() <  self%massMaximum                               &
            & )                                                                    &
            & call Merger_Tree_Dump(                                               &
            &                       treeCurrent                                  , &
            &                       scaleNodesByLogMass= self%scaleNodesByLogMass, &
            &                       edgeLengthsToTimes = self%edgeLengthsToTimes , &
            &                       useNodeLabels      = self%useNodeLabels      , &
            &                       timeRange          =[                          &
            &                                            self%timeMinimum        , &
            &                                            self%timeMaximum          &
            &                                           ]                        , &
            &                       nodeStyle          ='solid'                  , &
            &                       path               =self%path                  &
            &                      )
       ! Move to the next tree.
       treeCurrent => treeCurrent%nextTree
    end do
    return
  end subroutine dumpToGraphVizOperatePreEvolution
