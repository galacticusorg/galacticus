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
  Implements the standard class for initializing merger trees.
  !!}

  use :: Nodes_Operators, only : nodeOperatorClass
  
  !![
  <mergerTreeInitializor name="mergerTreeInitializorStandard">
   <description>The standard merger tree initializer.</description>
  </mergerTreeInitializor>
  !!]
  type, extends(mergerTreeInitializorClass) :: mergerTreeInitializorStandard
     !!{
     Implementation of the standard merger tree initializer.
     !!}
     private
     class(nodeOperatorClass), pointer :: nodeOperator_ => null()
   contains
     final     ::               standardDestructor
     procedure :: initialize => standardInitialize
  end type mergerTreeInitializorStandard

  interface mergerTreeInitializorStandard
     !!{
     Constructors for the \refClass{mergerTreeInitializorStandard} merger tree initializer.
     !!}
     module procedure standardConstructorParameters
     module procedure standardConstructorInternal
  end interface mergerTreeInitializorStandard

contains

  function standardConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerTreeInitializorStandard} merger tree initializer class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (mergerTreeInitializorStandard)                :: self
    type (inputParameters              ), intent(inout) :: parameters
    class(nodeOperatorClass            ), pointer       :: nodeOperator_


    !![
    <objectBuilder class="nodeOperator" name="nodeOperator_" source="parameters"/>
    !!]
    self=mergerTreeInitializorStandard(nodeOperator_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="nodeOperator_"/>
    !!]
    return
  end function standardConstructorParameters

   function standardConstructorInternal(nodeOperator_) result(self)
     !!{
     Internal constructor for the \refClass{mergerTreeInitializorStandard} merger tree initializer class.
     !!}
     implicit none
     type (mergerTreeInitializorStandard)                        :: self
     class(nodeOperatorClass            ), intent(in   ), target :: nodeOperator_
     !![
     <constructorAssign variables="*nodeOperator_"/>
     !!]

     return
   end function standardConstructorInternal

  subroutine standardDestructor(self)
    !!{
    Destructor for the \refClass{mergerTreeInitializorStandard} merger tree initializer class.
    !!}
    implicit none
    type(mergerTreeInitializorStandard), intent(inout) :: self

    !![
    <objectDestructor name="self%nodeOperator_"/>
    !!]
    return
  end subroutine standardDestructor

  subroutine standardInitialize(self,tree,timeEnd)
    !!{
    Walk through all nodes of a tree and call any routines that requested to perform initialization tasks.
    !!}
    use :: Galacticus_Nodes   , only : nodeComponentBasic, treeNode
    use :: Merger_Tree_Walkers, only : mergerTreeWalkerAllNodes
    !![
    <include directive="mergerTreeInitializeTask" type="moduleUse">
    !!]
    include 'merger_trees.initialize.tasks.modules.inc'
    !![
    </include>
    !!]
    implicit none
    class           (mergerTreeInitializorStandard), intent(inout) :: self
    type            (mergerTree                   ), intent(inout) :: tree
    double precision                               , intent(in   ) :: timeEnd
    type            (treeNode                     ), pointer       :: node
    class           (nodeComponentBasic           ), pointer       :: basic
    type            (mergerTreeWalkerAllNodes     )                :: treeWalker

    if (tree%initializedUntil >= timeEnd) return
    treeWalker=mergerTreeWalkerAllNodes(tree,spanForest=.false.)
    do while (treeWalker%next(node))
       ! Initialize only nodes that exist before the end time.
       basic => node%basic()
       if (basic%time() > tree%initializedUntil .and. basic%time() <= timeEnd) then
          ! Initialize using node operators first as we need, for example, spin to be initialized before some initialization
          ! within node components.
          call self%nodeOperator_%nodeInitialize(node)
          ! Call subroutines to perform any necessary initialization of this node.
          !![
          <include directive="mergerTreeInitializeTask" type="functionCall" functionType="void">
           <functionArgs>node</functionArgs>
          !!]
	  include 'merger_trees.initialize.tasks.inc'
          !![
          </include>
          !!]
       end if
    end do
    tree%initializedUntil=timeEnd
    return
  end subroutine standardInitialize
