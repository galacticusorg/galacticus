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
Implements a deforestation operator on merger trees (i.e. removes all but the most massive tree in a forest).
!!}

  !![
  <mergerTreeOperator name="mergerTreeOperatorDeforest">
   <description>Provides a deforestation operator for merger trees. Given a forest, this operator will destroy all but the first tree in the forest.</description>
  </mergerTreeOperator>
  !!]
  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorDeforest
     !!{
     A deforestation merger tree operator class.
     !!}
     private
   contains
     procedure :: operatePreEvolution => deforestOperatePreEvolution
  end type mergerTreeOperatorDeforest

  interface mergerTreeOperatorDeforest
     !!{
     Constructors for the deforestation merger tree operator class.
     !!}
     module procedure deforestConstructorParameters
  end interface mergerTreeOperatorDeforest

contains

  function deforestConstructorParameters(parameters) result(self)
    !!{
    Constructor for the deforestation merger tree operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(mergerTreeOperatorDeforest)                :: self
    type(inputParameters           ), intent(inout) :: parameters

    self=mergerTreeOperatorDeforest()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function deforestConstructorParameters

  subroutine deforestOperatePreEvolution(self,tree)
    !!{
    Perform a deforestation operation on a merger tree.
    !!}
    use :: Galacticus_Nodes, only : mergerTree, nodeComponentBasic, treeNode
    implicit none
    class           (mergerTreeOperatorDeforest), intent(inout), target :: self
    type            (mergerTree                ), intent(inout), target :: tree
    type            (treeNode                  ), pointer               :: nodeBase       , nodeNext
    type            (mergerTree                ), pointer               :: currentTree
    class           (nodeComponentBasic        ), pointer               :: basic
    double precision                                                    :: massRootMaximum
    integer                                                             :: treeIndex      , massRootMaximumIndex
    !$GLC attributes unused :: self

    ! Iterate over trees to find the most massive.
    currentTree          => tree
    massRootMaximum      =  0.0d0
    treeIndex            =  0
    massRootMaximumIndex =  -1
    do while (associated(currentTree))
       treeIndex=treeIndex+1
       basic => currentTree%nodeBase%basic()
       if (basic%mass() > massRootMaximum) then
          massRootMaximum     =basic    %mass()
          massRootMaximumIndex=treeIndex
       end if
       ! Move to the next tree.
       currentTree => currentTree%nextTree
    end do
     ! Iterate over trees.
    currentTree => tree
    treeIndex   =  0
    do while (associated(currentTree))
       treeIndex=treeIndex+1
       if (treeIndex /= massRootMaximumIndex) then
          ! Get root node of the tree.
          nodeBase => currentTree%nodeBase%firstChild
          do while (associated(nodeBase))
             nodeNext => nodeBase%sibling
             call nodeBase%destroyBranch()
             nodeBase => nodeNext
          end do
       end if
       ! Move to the next tree.
       currentTree => currentTree%nextTree
    end do
    return
  end subroutine deforestOperatePreEvolution

