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
  Implements a merger tree operator which makes mass growth along
  branch monotonically increasing.
  !!}

  !![
  <mergerTreeOperator name="mergerTreeOperatorMonotonizeMassGrowth">
   <description>
    A merger tree operator class which enforces monotonic growth a halo mass along each branch of each merger tree. It does
    this by searching the tree for nodes which are less massive than the sum of the masses of their immediate progenitors, and
    increasing the mass of such nodes to equal the sum of the masses of their immediate progenitors.
   </description>
  </mergerTreeOperator>
  !!]
  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorMonotonizeMassGrowth
     !!{
     A merger tree operator class makes mass growth along branch monotonically increasing.
     !!}
     private
   contains
     procedure :: operatePreEvolution => monotonizeMassGrowthOperatePreEvolution
  end type mergerTreeOperatorMonotonizeMassGrowth

  interface mergerTreeOperatorMonotonizeMassGrowth
     !!{
     Constructors for the mass growth monotonizing merger tree operator class.
     !!}
     module procedure monotonizeMassGrowthConstructorParameters
  end interface mergerTreeOperatorMonotonizeMassGrowth

contains

  function monotonizeMassGrowthConstructorParameters(parameters) result(self)
    !!{
    Constructor for the mass growth monotonizing merger tree operator class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(mergerTreeOperatorMonotonizeMassGrowth)                :: self
    type(inputParameters                       ), intent(inout) :: parameters

    self=mergerTreeOperatorMonotonizeMassGrowth()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function monotonizeMassGrowthConstructorParameters

  subroutine monotonizeMassGrowthOperatePreEvolution(self,tree)
    !!{
    Perform a mass growth monotonizing operation on a merger tree.
    !!}
    use :: Galacticus_Nodes   , only : mergerTree                   , nodeComponentBasic, treeNode
    use :: Merger_Tree_Walkers, only : mergerTreeWalkerIsolatedNodes
    implicit none
    class           (mergerTreeOperatorMonotonizeMassGrowth), intent(inout), target :: self
    type            (mergerTree                            ), intent(inout), target :: tree
    type            (treeNode                              ), pointer               :: nodeProgenitor , node
    class           (nodeComponentBasic                    ), pointer               :: basicProgenitor, basic
    type            (mergerTree                            ), pointer               :: treeCurrent
    type            (mergerTreeWalkerIsolatedNodes         )                        :: treeWalker
    logical                                                                         :: didModifyTree
    double precision                                                                :: massProgenitor
    !$GLC attributes unused :: self

    ! Iterate over trees.
    treeCurrent => tree
    do while (associated(treeCurrent))
       didModifyTree=.true.
       do while (didModifyTree)
          didModifyTree=.false.
          ! Walk the tree.
          treeWalker=mergerTreeWalkerIsolatedNodes(treeCurrent)
          do while (treeWalker%next(node))
             ! Find nodes that have children.
             if (associated(node%firstChild)) then
                ! Find the mass of all progenitor nodes.
                massProgenitor =  0.0d0
                nodeProgenitor => node%firstChild
                do while (associated(nodeProgenitor))
                   basicProgenitor =>  nodeProgenitor %basic  ()
                   massProgenitor  =  +massProgenitor            &
                        &             +basicProgenitor%mass   ()
                   nodeProgenitor  =>  nodeProgenitor %sibling
                end do
                ! Find nodes which are less massive than the sum of their progenitors.
                basic => node%basic()
                if (basic%mass() < massProgenitor) then
                   call basic%massSet(massProgenitor)
                   didModifyTree=.true.
                end if
             end if
          end do
       end do
       ! Move to the next tree.
       treeCurrent => treeCurrent%nextTree
    end do
    return
  end subroutine monotonizeMassGrowthOperatePreEvolution

