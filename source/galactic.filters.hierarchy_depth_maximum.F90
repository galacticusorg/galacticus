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
Implements a filter which passes only halos below a specified hierarchy depth.
!!}

  !![
  <galacticFilter name="galacticFilterHierarchyDepthMaximum">
   <description>A filter which passes only isolated halos.</description>
  </galacticFilter>
  !!]
  type, extends(galacticFilterClass) :: galacticFilterHierarchyDepthMaximum
     !!{
     A galactic filter class which passes only isolated halos.
     !!}
     private
     integer :: depthHierarchyLargest, nodeHierarchyLevelMaximumID
   contains
     procedure :: passes => hierarchyDepthMaximumPasses
  end type galacticFilterHierarchyDepthMaximum

  interface galacticFilterHierarchyDepthMaximum
     !!{
     Constructors for the ``hierarchyDepthMaximum'' galactic filter class.
     !!}
     module procedure hierarchyDepthMaximumConstructorParameters
     module procedure hierarchyDepthMaximumConstructorInternal
  end interface galacticFilterHierarchyDepthMaximum

contains

  function hierarchyDepthMaximumConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``hierarchyDepthMaximum'' galactic filter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type   (galacticFilterHierarchyDepthMaximum)                :: self
    type   (inputParameters                    ), intent(inout) :: parameters
    integer                                                     :: depthHierarchyLargest

    !![
    <inputParameter>
      <name>depthHierarchyLargest</name>
      <source>parameters</source>
      <description>The largest value of hierarchy maximum depth to pass.</description>
    </inputParameter>
    !!]
    self=galacticFilterHierarchyDepthMaximum(depthHierarchyLargest)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function hierarchyDepthMaximumConstructorParameters

  function hierarchyDepthMaximumConstructorInternal(depthHierarchyLargest) result(self)
    !!{
    Internal constructor for the ``hierarchyDepthMaximum'' galactic filter class.
    !!}
    implicit none
    type   (galacticFilterHierarchyDepthMaximum)                :: self
    integer                                     , intent(in   ) :: depthHierarchyLargest
    !![
    <constructorAssign variables="depthHierarchyLargest"/>
    !!]

    !![
    <addMetaProperty component="basic" name="nodeHierarchyLevelMaximum" type="integer" id="self%nodeHierarchyLevelMaximumID" isCreator="no"/>
    !!]
    return
  end function hierarchyDepthMaximumConstructorInternal

  logical function hierarchyDepthMaximumPasses(self,node) result(passes)
    !!{
    Implement a galactic filter which passes only isolated halos.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(galacticFilterHierarchyDepthMaximum), intent(inout)          :: self
    type (treeNode                           ), intent(inout), target  :: node
    class(nodeComponentBasic                 )               , pointer :: basic

    basic  => node %basic                      (                                )
    passes =  basic%integerRank0MetaPropertyGet(self%nodeHierarchyLevelMaximumID) <= self%depthHierarchyLargest
    return
  end function hierarchyDepthMaximumPasses
