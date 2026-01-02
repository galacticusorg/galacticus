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
  Implements a merger tree outputter class which performs analyzes on the trees.
  !!}

  use :: Output_Analyses, only : outputAnalysis, outputAnalysisClass

  !![
  <mergerTreeOutputter name="mergerTreeOutputterAnalyzer">
   <description>A merger tree outputter class which performs analyzes on the trees.</description>
  </mergerTreeOutputter>
  !!]
  type, extends(mergerTreeOutputterClass) :: mergerTreeOutputterAnalyzer
     !!{
     Implementation of a merger tree outputter class which performs analyzes on the trees.
     !!}
     private
     class(outputAnalysisClass), pointer :: outputAnalysis_ => null()
   contains
     final     ::               analyzerDestructor
     procedure :: outputTree => analyzerOutputTree
     procedure :: outputNode => analyzerOutputNode
     procedure :: finalize   => analyzerFinalize
     procedure :: reduce     => analyzerReduce
  end type mergerTreeOutputterAnalyzer

  interface mergerTreeOutputterAnalyzer
     !!{
     Constructors for the \refClass{mergerTreeOutputterAnalyzer} merger tree outputter.
     !!}
     module procedure analyzerConstructorParameters
     module procedure analyzerConstructorInternal
  end interface mergerTreeOutputterAnalyzer

contains

  function analyzerConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerTreeOutputterAnalyzer} merger tree outputter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (mergerTreeOutputterAnalyzer)                :: self
    type   (inputParameters            ), intent(inout) :: parameters
    class  (outputAnalysisClass        ), pointer       :: outputAnalysis_

    !![
    <objectBuilder class="outputAnalysis" name="outputAnalysis_" source="parameters"/>
    !!]
    self=mergerTreeOutputterAnalyzer(outputAnalysis_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="outputAnalysis_"/>
    !!]
    return
  end function analyzerConstructorParameters

  function analyzerConstructorInternal(outputAnalysis_) result(self)
    !!{
    Internal constructor for the \refClass{mergerTreeOutputterAnalyzer} merger tree outputter class.
    !!}
    implicit none
    type (mergerTreeOutputterAnalyzer)                        :: self
    class(outputAnalysisClass        ), intent(in   ), target :: outputAnalysis_
    !![
    <constructorAssign variables="*outputAnalysis_"/>
    !!]

    return
  end function analyzerConstructorInternal

  subroutine analyzerDestructor(self)
    !!{
    Destructor  for the {\normalfont \ttfamily analyzer} merger tree outputter class.
    !!}
    implicit none
    type(mergerTreeOutputterAnalyzer), intent(inout) :: self

    !![
    <objectDestructor name="self%outputAnalysis_"/>
    !!]
    return
  end subroutine analyzerDestructor

  subroutine analyzerOutputTree(self,tree,indexOutput,time)
    !!{
    Write properties of nodes in {\normalfont \ttfamily tree} to the \glc\ output file.
    !!}
    use :: Calculations_Resets, only : Calculations_Reset
    use :: Galacticus_Nodes   , only : mergerTree              , nodeComponentBasic, treeNode
    use :: Merger_Tree_Walkers, only : mergerTreeWalkerAllNodes
    implicit none
    class           (mergerTreeOutputterAnalyzer), intent(inout)          :: self
    type            (mergerTree                 ), intent(inout), target  :: tree
    integer         (c_size_t                   ), intent(in   )          :: indexOutput
    double precision                             , intent(in   )          :: time
    type            (treeNode                   )               , pointer :: node
    class           (nodeComponentBasic         )               , pointer :: basic
    type            (mergerTree                 )               , pointer :: treeCurrent
    type            (mergerTreeWalkerAllNodes   )                         :: treeWalker

    ! Iterate over trees.
    treeCurrent => tree
    do while (associated(treeCurrent))
       ! Iterate over nodes.
       treeWalker=mergerTreeWalkerAllNodes(treeCurrent,spanForest=.false.)
       ! Inform analyzers that a new tree is being output.
       call self%outputAnalysis_%newTree(treeCurrent,indexOutput)
       do while (treeWalker%next(node))
          ! Reset calculations (necessary in case the last node to be evolved is the first one we output, in which case
          ! calculations would not be automatically reset because the node unique ID will not have changed).
          call Calculations_Reset(node)
          ! Get the basic component.
          basic => node%basic()
          if (basic%time() == time) call self%outputAnalysis_%analyze(node,indexOutput)
       end do
       ! Skip to the next tree.
       treeCurrent => treeCurrent%nextTree
    end do
    return
  end subroutine analyzerOutputTree

  subroutine analyzerOutputNode(self,node,indexOutput)
    !!{
    Perform no output.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class  (mergerTreeOutputterAnalyzer), intent(inout) :: self
    type   (treeNode                   ), intent(inout) :: node
    integer(c_size_t                   ), intent(in   ) :: indexOutput
    !$GLC attributes unused :: self, node, indexOutput

    call Error_Report('output of single nodes is not supported'//{introspection:location})
    return
  end subroutine analyzerOutputNode

  subroutine analyzerReduce(self,reduced)
    !!{
    Reduce over the outputter.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(mergerTreeOutputterAnalyzer), intent(inout) :: self
    class(mergerTreeOutputterClass   ), intent(inout) :: reduced

    select type (reduced)
    type is (mergerTreeOutputterAnalyzer)
       call self%outputAnalysis_%reduce(reduced%outputAnalysis_)
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine analyzerReduce

  subroutine analyzerFinalize(self)
    !!{
    Finalize merger tree output by finalizing analyses.
    !!}
    use :: Output_HDF5, only : outputFileIsOpen
    implicit none
    class  (mergerTreeOutputterAnalyzer), intent(inout) :: self

    if (outputFileIsOpen) call self%outputAnalysis_%finalize()
    return
  end subroutine analyzerFinalize
