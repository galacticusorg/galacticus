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
  Implements a node operator class that collects and stores the mass accretion history of each node.
  !!}

  !![
  <nodeOperator name="nodeOperatorMassAccretionHistory">
   <description>A node operator class that collects and stores the mass accretion history of each node.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorMassAccretionHistory
     !!{
     A node operator class that collects and stores the mass accretion history of each node.
     !!}
     private
     integer :: massAccretionHistoryTimeID, massAccretionHistoryMassID
   contains
     procedure :: nodeInitialize => massAccretionHistoryNodeInitialize
  end type nodeOperatorMassAccretionHistory
  
  interface nodeOperatorMassAccretionHistory
     !!{
     Constructors for the \refClass{nodeOperatorMassAccretionHistory} node operator class.
     !!}
     module procedure massAccretionHistoryConstructorParameters
     module procedure massAccretionHistoryConstructorInternal
  end interface nodeOperatorMassAccretionHistory
  
contains

  function massAccretionHistoryConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorMassAccretionHistory} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodeOperatorMassAccretionHistory)                :: self
    type(inputParameters                 ), intent(inout) :: parameters
    
    self=nodeOperatorMassAccretionHistory()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function massAccretionHistoryConstructorParameters

  function massAccretionHistoryConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorMassAccretionHistory} node operator class.
    !!}
    use :: Galacticus_Nodes, only : defaultBasicComponent
    implicit none
    type(nodeOperatorMassAccretionHistory) :: self
   
    !![
    <addMetaProperty component="basic" name="massAccretionHistoryTime" id="self%massAccretionHistoryTimeID" rank="1" isCreator="yes"/>
    <addMetaProperty component="basic" name="massAccretionHistoryMass" id="self%massAccretionHistoryMassID" rank="1" isCreator="yes"/>
    !!]
    return
  end function massAccretionHistoryConstructorInternal

  subroutine massAccretionHistoryNodeInitialize(self,node)
    !!{
    Record the mass accretion history of the node.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class           (nodeOperatorMassAccretionHistory), intent(inout), target      :: self
    type            (treeNode                        ), intent(inout), target      :: node
    double precision                                  , dimension(:) , allocatable :: times     , masses
    type            (treeNode                        )               , pointer     :: nodeWork
    class           (nodeComponentBasic              )               , pointer     :: basic
    integer                                                                        :: countNodes
    
    if (associated(node%firstChild)) return
    nodeWork   => node
    countNodes =  1
    do while (nodeWork%isPrimaryProgenitor())
       countNodes =  countNodes       +1
       nodeWork   => nodeWork  %parent
    end do
    allocate(times (countNodes))
    allocate(masses(countNodes))
    nodeWork               => node
    countNodes             =  1
    basic                  => nodeWork%basic()
    times     (countNodes) =  basic   %time ()
    masses    (countNodes) =  basic   %mass ()
    do while (nodeWork%isPrimaryProgenitor())
       countNodes         =  countNodes         +1
       nodeWork           => nodeWork  %parent
       basic              => nodeWork  %basic ()
       times (countNodes) =  basic     %time  ()
       masses(countNodes) =  basic     %mass  ()
    end do
    basic => node%basic()
    call basic%floatRank1MetaPropertySet(self%massAccretionHistoryTimeID,times )
    call basic%floatRank1MetaPropertySet(self%massAccretionHistoryMassID,masses)
    return
  end subroutine massAccretionHistoryNodeInitialize
  
