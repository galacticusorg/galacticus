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
  Implements a node operator class that cleans any remaining stubs of the merger tree left in place by subsampling.  
  !!}
  
  !![
  <nodeOperator name="nodeOperatorCleanSubsampleStubs">
   <description>
    A node operator class that cleans any remaining stubs of the merger tree left in place by subsampling.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorCleanSubsampleStubs
     !!{
     A node operator class that cleans any remaining stubs of the merger tree left in place by subsampling.
     !!}
     private
     double precision :: factorMassGrowthConsolidate
   contains
     procedure :: nodeTreeInitialize => cleanSubsampleStubsNodeTreeInitialize
  end type nodeOperatorCleanSubsampleStubs
  
  interface nodeOperatorCleanSubsampleStubs
     !!{
     Constructors for the \refClass{nodeOperatorCleanSubsampleStubs} node operator class.
     !!}
     module procedure cleanSubsampleStubsConstructorParameters
     module procedure cleanSubsampleStubsConstructorInternal
  end interface nodeOperatorCleanSubsampleStubs
  
contains
  
  function cleanSubsampleStubsConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorCleanSubsampleStubs} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nodeOperatorCleanSubsampleStubs)                :: self
    type            (inputParameters                ), intent(inout) :: parameters
    double precision                                                 :: factorMassGrowthConsolidate

    !![
    <inputParameter>
      <name>factorMassGrowthConsolidate</name>
      <source>parameters</source>
      <description>The maximum factor by which the mass is allowed to grow between child and parent when consolidating nodes. A non-positive value prevents consolidation.</description>
      <defaultValue>0.0d0</defaultValue>
    </inputParameter>
    !!]
    self=nodeOperatorCleanSubsampleStubs(factorMassGrowthConsolidate)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function cleanSubsampleStubsConstructorParameters

  function cleanSubsampleStubsConstructorInternal(factorMassGrowthConsolidate) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorCleanSubsampleStubs} node operator class.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nodeOperatorCleanSubsampleStubs)                :: self
    double precision                                 , intent(in   ) :: factorMassGrowthConsolidate
    !![
    <constructorAssign variables="factorMassGrowthConsolidate"/>
    !!]

    return
  end function cleanSubsampleStubsConstructorInternal

  subroutine cleanSubsampleStubsNodeTreeInitialize(self,node)
    !!{
    Clean the tree of any stub branches left over from subsampling the merger tree.
    !!}
    use :: Merger_Tree_Walkers, only : mergerTreeWalkerIsolatedNodes
    use :: Galacticus_Nodes   , only : nodeComponentBasic
    implicit none
    class  (nodeOperatorCleanSubsampleStubs), intent(inout), target  :: self
    type   (treeNode                       ), intent(inout), target  :: node
    type   (treeNode                       )               , pointer :: nodeWork      , nodeChild  , &
         &                                                              nodeGrandchild, nodeParent , &
         &                                                              nodeNext
    class  (nodeComponentBasic             )               , pointer :: basic         , basicParent
    type   (mergerTreeWalkerIsolatedNodes  )                         :: treeWalker
    integer(c_size_t                       )                         :: countNodes
    logical                                                          :: nodesRemain

    ! Only operator once the root node is reached so that we can be sure that none of the stubs are needed.
    if (associated(node%parent)) return
    treeWalker=mergerTreeWalkerIsolatedNodes(node%hostTree,spanForest=.false.)
    nodesRemain  =treeWalker%next(nodeWork)
    do while (nodesRemain)
       nodesRemain=treeWalker%next(nodeNext)
       if (nodeWork%subsamplingWeight() < 0.0d0) then
          ! Decouple the node from the tree.
          nodeParent => nodeWork  %parent
          nodeChild  => nodeParent%firstChild
          do while (.not.associated(nodeChild%sibling,nodeWork))
             nodeChild => nodeChild%sibling
          end do
          nodeChild%sibling => nodeWork%sibling
          ! Destroy and deallocate the node.
          call nodeWork%destroy()
          deallocate(nodeWork)
          ! Determine if we can consolidate any nodes down the parent branch.
          if (self%factorMassGrowthConsolidate > 0.0d0) then
             ! Seek down through the branch until which find a node which either has a sibling (so can't be consolidated), or
             ! which has a mass sufficiently different from that of the starting node. Count how many such nodes we find.
             nodeChild   => nodeParent
             basic       => nodeChild %basic()
             basicParent => nodeParent%basic()
             countNodes  =  0_c_size_t          
             do while (                                                              &
                  &         associated(nodeChild%firstChild)                         &
                  &    .and.                                                         &
                  &    .not.associated(nodeChild%sibling   )                         &
                  &    .and.                                                         &
                  &      basic      %mass()*(1.0d0+self%factorMassGrowthConsolidate) &
                  &     >                                                            &
                  &      basicParent%mass()                                          &
                  &   )
                nodeChild  => nodeChild%firstChild
                basic      => nodeChild%basic     ()
                countNodes =  countNodes+1_c_size_t
             end do
             ! If we have found nodes that can be consolidated, remove the intervening nodes.
             if (countNodes > 1_c_size_t) then
                nodeChild => nodeParent%firstChild
                do while (countNodes > 1_c_size_t)
                   nodeGrandchild => nodeChild%firstChild
                   call nodeChild%destroy()
                   deallocate(nodeChild)
                   countNodes =  countNodes-1_c_size_t
                   nodeChild  => nodeGrandchild
                end do
                nodeParent%firstChild => nodeChild
                nodeChild %parent     => nodeParent
                do while (associated(nodeChild%sibling))
                   nodeChild        => nodeChild %sibling
                   nodeChild%parent => nodeParent
                end do
             end if
          end if
       end if
       nodeWork => nodeNext
    end do
    return
  end subroutine cleanSubsampleStubsNodeTreeInitialize
    
