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
  Implements a node operator class that determines if a node is on the most massive branch of its tree.
  !!}
  
  !![
  <nodeOperator name="nodeOperatorBranchMostMassive">
    <description>
      A node operator class that determines if a node is on the most massive branch of its tree. Intended to be paired with the
      \refClass{nodePropertyExtractorBranchMostMassive} class to extract this meta-data for output.
    </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorBranchMostMassive
     !!{
     A node operator class that tracks the maximum host halo mass which a node has occupied.
     !!}
     private
     integer :: isMostMassiveBranchID
   contains
     procedure :: nodeTreeInitialize => branchMostMassiveNodeTreeInitialize
  end type nodeOperatorBranchMostMassive
  
  interface nodeOperatorBranchMostMassive
     !!{
     Constructors for the \refClass{nodeOperatorBranchMostMassive} node operator class.
     !!}
     module procedure branchMostMassiveConstructorParameters
     module procedure branchMostMassiveConstructorInternal
  end interface nodeOperatorBranchMostMassive
  
contains

  function branchMostMassiveConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorBranchMostMassive} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodeOperatorBranchMostMassive)                :: self
    type(inputParameters              ), intent(inout) :: parameters
    
    self=nodeOperatorBranchMostMassive()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function branchMostMassiveConstructorParameters

  function branchMostMassiveConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorBranchMostMassive} node operator class.
    !!}
    implicit none
    type(nodeOperatorBranchMostMassive) :: self
    
    !![
    <addMetaProperty component="basic" name="isMostMassiveBranch" type="integer" id="self%isMostMassiveBranchID" isCreator="yes"/>
    !!]
    return
  end function branchMostMassiveConstructorInternal

  subroutine branchMostMassiveNodeTreeInitialize(self,node)
    !!{
    Determine if this node is on the most massive branch of its tree.
    !!}
    use :: Galacticus_Nodes   , only : nodeComponentBasic
    use :: Merger_Tree_Walkers, only : mergerTreeWalkerIsolatedNodes
    implicit none
    class(nodeOperatorBranchMostMassive), intent(inout), target  :: self
    type (treeNode                     ), intent(inout), target  :: node
    type (treeNode                     )               , pointer :: nodeOther
    class(nodeComponentBasic           )               , pointer :: basic     , basicOther
    type (mergerTreeWalkerIsolatedNodes)                         :: treeWalker

    basic => node %basic()
    call basic%integerRank0MetaPropertySet(self%isMostMassiveBranchID,1)
    treeWalker=mergerTreeWalkerIsolatedNodes(node%hostTree)
    do while (treeWalker%next(nodeOther))
       basicOther => nodeOther%basic()
       if (basicOther%time() <= basic%time() .and. basicOther%mass() > basic%mass()) then
          ! A more massive halo exists at a time earlier than the current halo. Therefore, the current halo is not on the most
          ! massive branch.
          call basic%integerRank0MetaPropertySet(self%isMostMassiveBranchID,0)
          return
       end if
    end do
    return
  end subroutine branchMostMassiveNodeTreeInitialize
