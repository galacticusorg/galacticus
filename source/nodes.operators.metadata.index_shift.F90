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
Implements a node operator class that shifts node indices at node promotion.
!!}

  !![
  <nodeOperator name="nodeOperatorIndexShift">
   <description>A node operator class that shifts node indices at node promotion.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorIndexShift
     !!{
     A node operator class that shifts node indices at node promotion.
     !!}
     private
   contains
     procedure :: nodePromote => indexShiftNodePromote
  end type nodeOperatorIndexShift

  interface nodeOperatorIndexShift
     !!{
     Constructors for the \refClass{nodeOperatorIndexShift} node operator class.
     !!}
     module procedure indexShiftConstructorParameters
  end interface nodeOperatorIndexShift

contains

  function indexShiftConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorIndexShift} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodeOperatorIndexShift)                :: self
    type(inputParameters       ), intent(inout) :: parameters
    
    self=nodeOperatorIndexShift()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function indexShiftConstructorParameters

  subroutine indexShiftNodePromote(self,node)
    !!{
    Act on node promotion.
    !!}
    implicit none
    class(nodeOperatorIndexShift), intent(inout) :: self
    type (treeNode              ), intent(inout) :: node
    !$GLC attributes unused :: self

    ! Shift the index from our node to the parent node.
    call node%parent%indexSet(node%index())
    return
  end subroutine indexShiftNodePromote
