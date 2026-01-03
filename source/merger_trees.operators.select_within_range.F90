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
Implements a select-within-range operator on the base nodes of merger trees.
!!}

  !![
  <mergerTreeOperator name="mergerTreeOperatorSelectWithinRange">
   <description>Provides a select-within-range operator on the base nodes of merger trees.</description>
  </mergerTreeOperator>
  !!]
  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorSelectWithinRange
     !!{
     A select-within-range merger tree operator class.
     !!}
     private
     double precision :: baseMassMinimum, baseMassMaximum
   contains
     procedure :: operatePreEvolution => selectWithinRangeOperatePreEvolution
  end type mergerTreeOperatorSelectWithinRange

  interface mergerTreeOperatorSelectWithinRange
     !!{
     Constructors for the select-within-range merger tree operator class.
     !!}
     module procedure selectWithinRangeConstructorParameters
     module procedure selectWithinRangeConstructorInternal
  end interface mergerTreeOperatorSelectWithinRange

contains

  function selectWithinRangeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the select-within-range merger tree operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (mergerTreeOperatorSelectWithinRange)                :: self
    type            (inputParameters                    ), intent(inout) :: parameters
    double precision                                                     :: baseMassMinimum, baseMassMaximum
    
    !![
    <inputParameter>
      <name>baseMassMinimum</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>Base node mass below which trees should be ignored.</description>
    </inputParameter>
    <inputParameter>
      <name>baseMassMaximum</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>Base node mass above which trees should be ignored.</description>
    </inputParameter>
    !!]
    self=mergerTreeOperatorSelectWithinRange(baseMassMinimum,baseMassMaximum)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function selectWithinRangeConstructorParameters

  function selectWithinRangeConstructorInternal(baseMassMinimum,baseMassMaximum) result(self)
    !!{
    Internal constructor for the select-within-range merger tree operator class.
    !!}
    implicit none
    type            (mergerTreeOperatorSelectWithinRange)                :: self
    double precision                                     , intent(in   ) :: baseMassMinimum, baseMassMaximum
    !![
    <constructorAssign variables="baseMassMinimum, baseMassMaximum"/>
    !!]
 
    return
  end function selectWithinRangeConstructorInternal

  subroutine selectWithinRangeOperatePreEvolution(self,tree)
    !!{
    Perform a select-within-range operation on a merger tree.
    !!}
    use :: Galacticus_Nodes              , only : mergerTree                    , nodeComponentBasic, treeNode
    use :: Merger_Trees_Pruning_Utilities, only : Merger_Tree_Prune_Clean_Branch
    implicit none
    class  (mergerTreeOperatorSelectWithinRange), intent(inout), target :: self
    type   (mergerTree                         ), intent(inout), target :: tree
    type   (treeNode                           ), pointer               :: nodeBase     , nodeNext
    class  (nodeComponentBasic                 ), pointer               :: basicNodeBase
    type   (mergerTree                         ), pointer               :: currentTree

    ! Iterate over trees.
    currentTree => tree
    do while (associated(currentTree))
       ! Get root node of the tree.
       nodeBase      => currentTree%nodeBase
       basicNodeBase => nodeBase   %basic   ()
       if     (                                             &
            &   basicNodeBase%mass() < self%baseMassMinimum &
            &  .or.                                         &
            &   basicNodeBase%mass() > self%baseMassMaximum &
            & ) then
          ! Tree is outside range. Destroy all but the base node. (Leaving just the base node
          ! makes the tree inert - i.e. it can not do anything.)
          nodeBase => nodeBase%firstChild
          do while (associated(nodeBase))
             nodeNext => nodeBase%sibling
             call Merger_Tree_Prune_Clean_Branch(nodeBase)
             call nodeBase%destroyBranch()
             deallocate(nodeBase)
             nodeBase => nodeNext
          end do
       end if
       ! Move to the next tree.
       currentTree => currentTree%nextTree
    end do
    return
  end subroutine selectWithinRangeOperatePreEvolution
