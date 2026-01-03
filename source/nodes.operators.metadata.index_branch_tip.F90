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
Implements a node operator class that records the index of the branch tip for each node.
!!}

  !![
  <nodeOperator name="nodeOperatorIndexBranchTip">
   <description>A node operator class that records the index of the branch tip for each node.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorIndexBranchTip
     !!{
     A node operator class that records the index of the branch tip for each node.
     !!}
     private
     integer :: indexBranchTipID
   contains
     procedure :: nodeTreeInitialize => indexBranchTipNodeTreeInitialize
  end type nodeOperatorIndexBranchTip

  interface nodeOperatorIndexBranchTip
     !!{
     Constructors for the \refClass{nodeOperatorIndexBranchTip} node operator class.
     !!}
     module procedure indexBranchTipConstructorParameters
     module procedure indexBranchTipConstructorInternal
  end interface nodeOperatorIndexBranchTip

contains

  function indexBranchTipConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorIndexBranchTip} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodeOperatorIndexBranchTip)                :: self
    type(inputParameters           ), intent(inout) :: parameters
    
    self=nodeOperatorIndexBranchTip()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function indexBranchTipConstructorParameters

  function indexBranchTipConstructorInternal() result(self)
    !!{
    Constructor for the \refClass{nodeOperatorIndexBranchTip} node operator class which takes a parameter set as input.
    !!}
    implicit none
    type(nodeOperatorIndexBranchTip) :: self
    
    !![
    <addMetaProperty component="basic" name="nodeIndexBranchTip" type="longInteger" id="self%indexBranchTipID" isCreator="yes"/>
    !!]
    return
  end function indexBranchTipConstructorInternal

  subroutine indexBranchTipNodeTreeInitialize(self,node)
    !!{
    Initialize node branch tip indices.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(nodeOperatorIndexBranchTip), intent(inout), target  :: self
    type (treeNode                  ), intent(inout), target  :: node
    type (treeNode                  )               , pointer :: nodeWork
    class(nodeComponentBasic        )               , pointer :: basic

    if (associated(node%firstChild)) return
    nodeWork => node
    do while (associated(nodeWork))
       basic => nodeWork%basic()
       call basic%longIntegerRank0MetaPropertySet(self%indexBranchTipID,node%index())
       if (nodeWork%isPrimaryProgenitor()) then
          nodeWork => nodeWork%parent
       else
          nodeWork => null()
       end if
    end do
    return
  end subroutine indexBranchTipNodeTreeInitialize
