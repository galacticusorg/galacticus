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
Implements a merger tree filter which passes if the tree matches the given index.
!!}

  use :: Kind_Numbers, only : kind_int8
  
  !![
  <mergerTreeFilter name="mergerTreeFilterTreeIndex">
   <description>A merger tree filter which passes if the tree matches the given index.</description>
  </mergerTreeFilter>
  !!]
  type, extends(mergerTreeFilterClass) :: mergerTreeFilterTreeIndex
     !!{
     A merger tree filter class which passes if the tree matches the given index.
     !!}
     private
     integer(kind_int8) :: mergerTreeIndex
   contains
     procedure :: passes => treeIndexPasses
  end type mergerTreeFilterTreeIndex

  interface mergerTreeFilterTreeIndex
     !!{
     Constructors for the \refClass{mergerTreeFilterTreeIndex} merger tree filter class.
     !!}
     module procedure treeIndexConstructorParameters
     module procedure treeIndexConstructorInternal
  end interface mergerTreeFilterTreeIndex

contains
  
  function treeIndexConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerTreeFilterTreeIndex} merger tree filter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (mergerTreeFilterTreeIndex)                :: self
    type   (inputParameters          ), intent(inout) :: parameters
    integer(kind_int8                )                :: mergerTreeIndex

    !![
    <inputParameter>
      <name>mergerTreeIndex</name>
      <source>parameters</source>
      <description>The index of the merger tree to pass.</description>
    </inputParameter>
    !!]
    self=mergerTreeFilterTreeIndex(mergerTreeIndex)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function treeIndexConstructorParameters

  function treeIndexConstructorInternal(mergerTreeIndex) result(self)
    !!{
    Internal constructor for the \refClass{mergerTreeFilterTreeIndex} merger tree filter class.
    !!}
    implicit none
    type   (mergerTreeFilterTreeIndex)                :: self
    integer(kind_int8                ), intent(in   ) :: mergerTreeIndex
    !![
    <constructorAssign variables="mergerTreeIndex"/>
    !!]

    return
  end function treeIndexConstructorInternal

  logical function treeIndexPasses(self,tree) result(passes)
    !!{
    Implement a merger tree filter which passes if the index matches a target index.
    !!}
    implicit none
    class(mergerTreeFilterTreeIndex), intent(inout) :: self
    type (mergerTree               ), intent(in   ) :: tree
    
    passes=tree%index == self%mergerTreeIndex
    return
  end function treeIndexPasses
