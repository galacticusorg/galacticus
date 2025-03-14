!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
Implements a merger tree filter which passes if any node in the tree passes the given galactic filter.
!!}

  use :: Galactic_Filters, only : galacticFilterClass
  
  !![
  <mergerTreeFilter name="mergerTreeFilterAnyNode">
   <description>A merger tree filter which passes if any node in the tree passes the given galactic filter.</description>
  </mergerTreeFilter>
  !!]
  type, extends(mergerTreeFilterClass) :: mergerTreeFilterAnyNode
     !!{
     A merger tree filter class which passes if any node in the tree passes the given galactic filter.
     !!}
     private
     class  (galacticFilterClass), pointer :: galacticFilter_ => null()
     type   (varying_string     )          :: label                    , labelDescription
     integer                               :: labelID
     logical                               :: labelBranch
   contains
     final     ::           anyNodeDestructor
     procedure :: passes => anyNodePasses
  end type mergerTreeFilterAnyNode

  interface mergerTreeFilterAnyNode
     !!{
     Constructors for the ``anyNode'' merger tree filter class.
     !!}
     module procedure anyNodeConstructorParameters
     module procedure anyNodeConstructorInternal
  end interface mergerTreeFilterAnyNode

contains
  
  function anyNodeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``anyNode'' merger tree filter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (mergerTreeFilterAnyNode)                :: self
    type   (inputParameters        ), intent(inout) :: parameters
    class  (galacticFilterClass    ), pointer       :: galacticFilter_
    type   (varying_string         )                :: label          , labelDescription
    logical                                         :: labelBranch

    !![
    <inputParameter>
      <name>label</name>
      <source>parameters</source>
      <defaultValue>var_str(' ')</defaultValue>
      <description>A label to apply to nodes that pass the filter.</description>
    </inputParameter>
    <inputParameter>
      <name>labelBranch</name>
      <source>parameters</source>
      <defaultValue>.false.</defaultValue>
      <description>If true apply the label to the entire branch of any node that passes the filter.</description>
    </inputParameter>
    !!]
    if (label == '') label=' '
    if (trim(label) /= '') then
       !![
       <inputParameter>
         <name>labelDescription</name>
         <source>parameters</source>
         <description>A description of the label.</description>
       </inputParameter>
       !!]
    end if
    !![
    <objectBuilder class="galacticFilter" name="galacticFilter_" source="parameters"/>
    !!]
    self=mergerTreeFilterAnyNode(label,labelBranch,labelDescription,galacticFilter_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="galacticFilter_"/>
    !!]
    return
  end function anyNodeConstructorParameters

  function anyNodeConstructorInternal(label,labelBranch,labelDescription,galacticFilter_) result(self)
    !!{
    Internal constructor for the ``anyNode'' merger tree filter class.
    !!}
    use :: Nodes_Labels, only : nodeLabelRegister
    implicit none
    type   (mergerTreeFilterAnyNode)                        :: self
    type   (varying_string         ), intent(in   )         :: label          , labelDescription
    class  (galacticFilterClass    ), intent(in   ), target :: galacticFilter_
    logical                         , intent(in   )         :: labelBranch
    !![
    <constructorAssign variables="label, labelBranch, labelDescription, *galacticFilter_"/>
    !!]

    if (trim(label) /= '') then
       self%labelID=nodeLabelRegister(char(label),char(labelDescription))
    else
       self%labelID=-1
    end if
    return
  end function anyNodeConstructorInternal

  subroutine anyNodeDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily anyNode} merger tree filter class.
    !!}
    implicit none
    type(mergerTreeFilterAnyNode), intent(inout) :: self

    !![
    <objectDestructor name="self%galacticFilter_"/>
    !!]
    return
  end subroutine anyNodeDestructor

  logical function anyNodePasses(self,tree)
    !!{
    Implement a merger tree filter which passes if any node in the tree passes the given merger tree filter.
    !!}
    use :: Galacticus_Nodes   , only : treeNode
    use :: Merger_Tree_Walkers, only : mergerTreeWalkerIsolatedNodes
    use :: Nodes_Labels       , only : nodeLabelSet
    implicit none
    class(mergerTreeFilterAnyNode      ), intent(inout) :: self
    type (mergerTree                   ), intent(in   ) :: tree
    type (treeNode                     ), pointer       :: node      , nodeBranch
    type (mergerTreeWalkerIsolatedNodes)                :: treeWalker
    
    anyNodePasses=.false.
    treeWalker   =mergerTreeWalkerIsolatedNodes(tree,spanForest=.true.)
    do while (treeWalker%next(node))
       if (self%galacticFilter_%passes(node)) then
          anyNodePasses=.true.
          if (self%labelID > 0) then
             if (self%labelBranch) then
                ! Apply the label to the entire branch to which this node belongs.
                nodeBranch => node
                do while (nodeBranch%isPrimaryProgenitor())
                   nodeBranch => nodeBranch%parent
                end do
                do while (associated(nodeBranch))
                   call nodeLabelSet(self%labelID,nodeBranch)
                   nodeBranch => nodeBranch%firstChild
                end do
             else
                ! Apply the label only to the node itself.
                call nodeLabelSet(self%labelID,node)
             end if
          else
             exit
          end if
       end if
    end do
    return
  end function anyNodePasses
