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
  Implements a node operator class that tracks the maximum progenitor halo mass of a node.
  !!}
  
  !![
  <nodeOperator name="nodeOperatorMassProgenitorMaximum">
    <description>
      A node operator class that tracks the maximum progenitor halo mass of a node. Intended to be paired with the
      \refClass{nodePropertyExtractorMassProgenitorMaximum} class to extract those ages for output.
    </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorMassProgenitorMaximum
     !!{
     A node operator class that tracks the maximum progenitor halo mass which of a node.
     !!}
     private
     integer :: massProgenitorMaximumID
   contains
     procedure :: nodeInitialize => massProgenitorMaximumNodeInitialize
     procedure :: nodePromote    => massProgenitorMaximumNodePromote
  end type nodeOperatorMassProgenitorMaximum
  
  interface nodeOperatorMassProgenitorMaximum
     !!{
     Constructors for the \refClass{nodeOperatorMassProgenitorMaximum} node operator class.
     !!}
     module procedure massProgenitorMaximumConstructorParameters
     module procedure massProgenitorMaximumConstructorInternal
  end interface nodeOperatorMassProgenitorMaximum
  
contains

  function massProgenitorMaximumConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorMassProgenitorMaximum} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodeOperatorMassProgenitorMaximum)                :: self
    type(inputParameters                  ), intent(inout) :: parameters
    
    self=nodeOperatorMassProgenitorMaximum()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function massProgenitorMaximumConstructorParameters

  function massProgenitorMaximumConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorMassProgenitorMaximum} node operator class.
    !!}
    implicit none
    type(nodeOperatorMassProgenitorMaximum) :: self
    
    !![
    <addMetaProperty component="basic" name="massProgenitorMaximum" id="self%massProgenitorMaximumID" isEvolvable="no" isCreator="yes"/>
    !!]
    return
  end function massProgenitorMaximumConstructorInternal

  subroutine massProgenitorMaximumNodeInitialize(self,node)
    !!{
    Initialize the maximum progenitor mass of this node.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class           (nodeOperatorMassProgenitorMaximum), intent(inout), target  :: self
    type            (treeNode                         ), intent(inout), target  :: node
    type            (treeNode                         )               , pointer :: nodeProgenitor
    class           (nodeComponentBasic               )               , pointer :: basic         , basicProgenitor
    double precision                                                            :: massMaximum
    
    ! Find the maximum progenitor mass.
    nodeProgenitor => node
    massMaximum    =  0.0d0
    do while (associated(nodeProgenitor))
       basicProgenitor => nodeProgenitor%basic()
       massMaximum     =  max(massMaximum,basicProgenitor%mass())
       nodeProgenitor  => nodeProgenitor%firstChild
    end do
    basic => node%basic()
    call basic%floatRank0MetaPropertySet(self%massProgenitorMaximumID,massMaximum)
    return
  end subroutine massProgenitorMaximumNodeInitialize

  subroutine massProgenitorMaximumNodePromote(self,node)
    !!{
    Update the maximum host mass of this node as a result of node promotion.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(nodeOperatorMassProgenitorMaximum), intent(inout) :: self
    type (treeNode                         ), intent(inout) :: node
    class(nodeComponentBasic               ), pointer       :: basicParent, basic
    
    ! Adjust the maximum mass to that of the parent node.
    basic       => node       %basic()
    basicParent => node%parent%basic()
    call basic%floatRank0MetaPropertySet(                                                                     &
         &                                                                     self%massProgenitorMaximumID , &
         &                               basicParent%floatRank0MetaPropertyGet(self%massProgenitorMaximumID)  &
         &                              )
    return
  end subroutine massProgenitorMaximumNodePromote
