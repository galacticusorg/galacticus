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
Implements a node operator class that shifts node constrained branch status at node promotion.
!!}

  !![
  <nodeOperator name="nodeOperatorConstrainedBranch">
   <description>A node operator class that shifts node constrained branch status at node promotion.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorConstrainedBranch
     !!{
     A node operator class that shifts node constrained branch status at node promotion.
     !!}
     private
     integer :: isConstrainedID
   contains
     procedure :: nodePromote => constrainedBranchNodePromote
  end type nodeOperatorConstrainedBranch

  interface nodeOperatorConstrainedBranch
     !!{
     Constructors for the \refClass{nodeOperatorConstrainedBranch} node operator class.
     !!}
     module procedure constrainedBranchConstructorParameters
     module procedure constrainedBranchConstructorInternal
  end interface nodeOperatorConstrainedBranch

contains

  function constrainedBranchConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorConstrainedBranch} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodeOperatorConstrainedBranch)                :: self
    type(inputParameters              ), intent(inout) :: parameters
    
    self=nodeOperatorConstrainedBranch()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function constrainedBranchConstructorParameters

  function constrainedBranchConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorConstrainedBranch} node operator class.
    !!}
    implicit none
    type(nodeOperatorConstrainedBranch) :: self

    ! Use a matched "meta-property" as we used in the constrained tree build controller. This allows us to recover the stored
    ! "isConstrained" state of each node that was written by that object.
    !![
    <addMetaProperty component="basic" name="isConstrained" type="integer" id="self%isConstrainedID" isCreator="no"/>
    !!]
    return
  end function constrainedBranchConstructorInternal

  subroutine constrainedBranchNodePromote(self,node)
    !!{
    Act on node promotion.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(nodeOperatorConstrainedBranch), intent(inout) :: self
    type (treeNode                     ), intent(inout) :: node
    class(nodeComponentBasic           ), pointer       :: basic, basicParent

    basic       => node       %basic()
    basicparent => node%parent%basic()
    call basic%integerRank0MetaPropertySet(                                                               &
         &                                                                         self%isConstrainedID , &
         &                                 basicParent%integerRank0MetaPropertyGet(self%isConstrainedID)  &
         &                                )
    return
  end subroutine constrainedBranchNodePromote
