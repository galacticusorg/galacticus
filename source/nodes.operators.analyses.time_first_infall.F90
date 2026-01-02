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
  Implements a node operator class that tracks the time of first infall for a node.
  !!}
  
  !![
  <nodeOperator name="nodeOperatorTimeFirstInfall">
    <description>
      A node operator class that tracks the time of first infall for a node. Intended to be paired with the
      \refClass{nodePropertyExtractorTimeFirstInfall} class to extract these times for output.
    </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorTimeFirstInfall
     !!{
     A node operator class that tracks the time of first infall for a node.
     !!}
     private
     integer :: timeFirstInfallID
   contains
     procedure :: nodeInitialize => timeFirstInfallNodeInitialize
  end type nodeOperatorTimeFirstInfall
  
  interface nodeOperatorTimeFirstInfall
     !!{
     Constructors for the \refClass{nodeOperatorTimeFirstInfall} node operator class.
     !!}
     module procedure timeFirstInfallConstructorParameters
     module procedure timeFirstInfallConstructorInternal
  end interface nodeOperatorTimeFirstInfall
  
contains

  function timeFirstInfallConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorTimeFirstInfall} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodeOperatorTimeFirstInfall)                :: self
    type(inputParameters            ), intent(inout) :: parameters
    
    self=nodeOperatorTimeFirstInfall()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function timeFirstInfallConstructorParameters

  function timeFirstInfallConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorTimeFirstInfall} node operator class.
    !!}
    implicit none
    type(nodeOperatorTimeFirstInfall) :: self
    
    !![
    <addMetaProperty component="basic" name="timeFirstInfall" id="self%timeFirstInfallID" isEvolvable="no" isCreator="yes"/>
    !!]
    return
  end function timeFirstInfallConstructorInternal
  
  subroutine timeFirstInfallNodeInitialize(self,node)
    !!{
    Initialize the time of first infall for this node.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class           (nodeOperatorTimeFirstInfall), intent(inout), target  :: self
    type            (treeNode                   ), intent(inout), target  :: node
    type            (treeNode                   )               , pointer :: nodeDescendant
    class           (nodeComponentBasic         )               , pointer :: basic
    double precision                                                      :: timeFirstInfall

    timeFirstInfall =  -1.0d0
    nodeDescendant  => node
    do while (nodeDescendant%isPrimaryProgenitor())
       nodeDescendant => nodeDescendant%parent
    end do
    if (associated(nodeDescendant%parent)) then
       basic           => nodeDescendant%parent%basic()
       timeFirstInfall =  basic                %time ()
    end if
    basic => node%basic()
    call basic%floatRank0MetaPropertySet(self%timeFirstInfallID,timeFirstInfall)
    return
  end subroutine timeFirstInfallNodeInitialize
