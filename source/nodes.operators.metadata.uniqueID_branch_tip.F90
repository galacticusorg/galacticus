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
Implements a node operator class that records the index of the branch tip for each node.
!!}

  !![
  <nodeOperator name="nodeOperatorUniqueIDBranchTip">
   <description>A node operator class that records the index of the branch tip for each node.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorUniqueIDBranchTip
     !!{
     A node operator class that records the index of the branch tip for each node.
     !!}
     private
     integer :: uniqueIDBranchTipID
   contains
     procedure :: nodeTreeInitialize => uniqueIDBranchTipNodeTreeInitialize
  end type nodeOperatorUniqueIDBranchTip

  interface nodeOperatorUniqueIDBranchTip
     !!{
     Constructors for the \refClass{nodeOperatorUniqueIDBranchTip} node operator class.
     !!}
     module procedure uniqueIDBranchTipConstructorParameters
     module procedure uniqueIDBranchTipConstructorInternal
  end interface nodeOperatorUniqueIDBranchTip

contains

  function uniqueIDBranchTipConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorUniqueIDBranchTip} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodeOperatorUniqueIDBranchTip)                :: self
    type(inputParameters              ), intent(inout) :: parameters
    
    self=nodeOperatorUniqueIDBranchTip()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function uniqueIDBranchTipConstructorParameters

  function uniqueIDBranchTipConstructorInternal() result(self)
    !!{
    Constructor for the \refClass{nodeOperatorUniqueIDBranchTip} node operator class which takes a parameter set as input.
    !!}
    implicit none
    type(nodeOperatorUniqueIDBranchTip) :: self
    
    !![
    <addMetaProperty component="basic" name="nodeUniqueIDBranchTip" type="longInteger" id="self%uniqueIDBranchTipID" isCreator="yes"/>
    !!]
    return
  end function uniqueIDBranchTipConstructorInternal

  subroutine uniqueIDBranchTipNodeTreeInitialize(self,node)
    !!{
    Initialize node branch tip indices.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(nodeOperatorUniqueIDBranchTip), intent(inout), target  :: self
    type (treeNode                     ), intent(inout), target  :: node
    type (treeNode                     )               , pointer :: nodeWork
    class(nodeComponentBasic           )               , pointer :: basic

    if (associated(node%firstChild)) return
    nodeWork => node
    do while (associated(nodeWork))
       basic => nodeWork%basic()
       call basic%longIntegerRank0MetaPropertySet(self%uniqueIDBranchTipID,node%uniqueID())
       if (nodeWork%isPrimaryProgenitor()) then
          nodeWork => nodeWork%parent
       else
          nodeWork => null()
       end if
    end do
    return
  end subroutine uniqueIDBranchTipNodeTreeInitialize
