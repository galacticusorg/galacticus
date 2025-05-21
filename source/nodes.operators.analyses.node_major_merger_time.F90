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
  Implements a node operator class that computes the time of the most recent major merger between nodes.
  !!}

  !![
  <nodeOperator name="nodeOperatorNodeMajorMergerTime">
   <description>A node operator class that computes the time of the most recent major merger between nodes.</description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorNodeMajorMergerTime
     !!{
     A node operator class that computes the time of the most recent major merger between nodes.
     !!}
     private
     integer          :: nodeMajorMergerTimeID
     double precision :: fractionMassMajorMerger
   contains
     procedure :: nodeInitialize => nodeMajorMergerTimeNodeInitialize
     procedure :: nodesMerge     => nodeMajorMergerTimeNodesMerge
     procedure :: nodePromote    => nodeMajorMergerTimeNodePromote
  end type nodeOperatorNodeMajorMergerTime
  
  interface nodeOperatorNodeMajorMergerTime
     !!{
     Constructors for the \refClass{nodeOperatorNodeMajorMergerTime} node operator class.
     !!}
     module procedure nodeMajorMergerTimeConstructorParameters
     module procedure nodeMajorMergerTimeConstructorInternal
  end interface nodeOperatorNodeMajorMergerTime
  
contains

  function nodeMajorMergerTimeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodeOperatorNodeMajorMergerTime} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (nodeOperatorNodeMajorMergerTime)                :: self
    type            (inputParameters                ), intent(inout) :: parameters
    double precision                                                 :: fractionMassMajorMerger

    !![
    <inputParameter>
      <name>fractionMassMajorMerger</name>
      <defaultValue>0.25d0</defaultValue>
      <description>The mass ratio ($M_2/M_1$ where $M_2 &lt; M_1$) of merging halos above which the merger should be considered to be ``major''.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=nodeOperatorNodeMajorMergerTime(fractionMassMajorMerger)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function nodeMajorMergerTimeConstructorParameters

  function nodeMajorMergerTimeConstructorInternal(fractionMassMajorMerger) result(self)
    !!{
    Internal constructor for the \refClass{nodeOperatorNodeMajorMergerTime} node operator class.
    !!}
    implicit none
    type            (nodeOperatorNodeMajorMergerTime)                :: self
    double precision                                 , intent(in   ) :: fractionMassMajorMerger
    !![
    <constructorAssign variables="fractionMassMajorMerger"/>
    !!]
    
    !![
    <addMetaProperty component="basic" name="nodeMajorMergerTime" id="self%nodeMajorMergerTimeID" isEvolvable="no" isCreator="yes"/>
    !!]
    return
  end function nodeMajorMergerTimeConstructorInternal

  subroutine nodeMajorMergerTimeNodeInitialize(self,node)
    !!{
    Initialize nodeMajorMergerTime level data.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(nodeOperatorNodeMajorMergerTime), intent(inout), target  :: self
    type (treeNode                       ), intent(inout), target  :: node
    class(nodeComponentBasic             )               , pointer :: basic

    basic => node%basic()
    call basic%floatRank0MetaPropertySet(self%nodeMajorMergerTimeID,-1.0d0)
    return
  end subroutine nodeMajorMergerTimeNodeInitialize

  subroutine nodeMajorMergerTimeNodesMerge(self,node)
    !!{
    Record node major merger times.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(nodeOperatorNodeMajorMergerTime), intent(inout) :: self
    type (treeNode                       ), intent(inout) :: node
    class(nodeComponentBasic             ), pointer       :: basic, basicParent

    basic       => node       %basic()
    basicParent => node%parent%basic()
    if (basic%mass() >= self%fractionMassMajorMerger*basicParent%mass()) &
         &  call basicParent%floatRank0MetaPropertySet(self%nodeMajorMergerTimeID,basic%time())
    return
  end subroutine nodeMajorMergerTimeNodesMerge
 
  subroutine nodeMajorMergerTimeNodePromote(self,node)
    !!{
    Promote node major merger times.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(nodeOperatorNodeMajorMergerTime), intent(inout) :: self
    type (treeNode                       ), intent(inout) :: node
    class(nodeComponentBasic             ), pointer       :: basic, basicParent
    
    basic       => node       %basic()
    basicParent => node%parent%basic()
    if (basicParent%floatRank0MetaPropertyGet(self%nodeMajorMergerTimeID) > basic%floatRank0MetaPropertyGet(self%nodeMajorMergerTimeID)) &
         & call basic%floatRank0MetaPropertySet(self%nodeMajorMergerTimeID,basicParent%floatRank0MetaPropertyGet(self%nodeMajorMergerTimeID))
    return
  end subroutine nodeMajorMergerTimeNodePromote
