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
  Provides a node merger class implementing a multi level hierarchy.
  !!}

  !![
  <mergerTreeNodeMerger name="mergerTreeNodeMergerMultiLevelHierarchy">
   <description>A node merger class implementing a multi level hierarchy.</description>
  </mergerTreeNodeMerger>
  !!]
  type, extends(mergerTreeNodeMergerClass) :: mergerTreeNodeMergerMultiLevelHierarchy
     !!{
     Implementation of the multi-level hierarchy node merger class.
     !!}
     private
   contains
     procedure :: process  => multiLevelHierarchyProcess
  end type mergerTreeNodeMergerMultiLevelHierarchy

  interface mergerTreeNodeMergerMultiLevelHierarchy
     !!{
     Constructors for the \refClass{mergerTreeNodeMergerMultiLevelHierarchy} node merger class.
     !!}
     module procedure multiLevelHierarchyConstructorParameters
  end interface mergerTreeNodeMergerMultiLevelHierarchy

contains

  function multiLevelHierarchyConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerTreeNodeMergerMultiLevelHierarchy} node merger class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(mergerTreeNodeMergerMultiLevelHierarchy)                :: self
    type(inputParameters                        ), intent(inout) :: parameters

    self=mergerTreeNodeMergerMultiLevelHierarchy()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function multiLevelHierarchyConstructorParameters

  subroutine multiLevelHierarchyProcess(self,node)
    !!{
    Processes a node merging event, utilizing a multi level substructure hierarchy.
    !!}
    use :: Display            , only : displayGreen  , displayReset
    use :: Error              , only : Error_Report
    use :: Galacticus_Nodes   , only : treeNode
    use :: ISO_Varying_String , only : varying_string, operator(//), assignment(=)
    use :: String_Handling    , only : operator(//)
    implicit none
    class(mergerTreeNodeMergerMultiLevelHierarchy), intent(inout)          :: self
    type (treeNode                               ), intent(inout), target  :: node
    type (treeNode                               )               , pointer :: nodeChild    , nodeParent, &
         &                                                                    nodeSatellite
    type (varying_string                         )                         :: message
    !$GLC attributes unused :: self

    ! Get the parent node.
    nodeParent => node      %parent
    ! Uncouple node from the children of its parent.
    nodeChild  => nodeParent%firstChild
    if (.not.associated(nodeChild%sibling)) then
       message='attempting to make node '
       message=message//node%index()//' a satellite, but it is the primary progenitor'//char(10)
       message=message//'this can happen if branch jumps are allowed and the tree is postprocessed to remove nodes'//char(10)
       message=message//displayGreen()//'HELP:'//displayReset()//' to resolve this issue, either switch off postprocessing of the tree, or prevent'//char(10)
       message=message//'branch jumps by setting [mergerTreeReadAllowBranchJumps]=false'
       call Error_Report(message//{introspection:location})
    end if
    do while (.not.associated(nodeChild%sibling,node))
       nodeChild => nodeChild%sibling
    end do
    nodeChild%sibling => node%sibling
    ! Unset the sibling pointer for this node.
    node %sibling => null()
    ! Add it to the list of satellite nodes associated with its parent.
    if (associated(nodeParent%firstSatellite)) then
       nodeSatellite                => nodeParent%lastSatellite()
       nodeSatellite%sibling        => node
    else
       nodeParent   %firstSatellite => node
    end if
    return
  end subroutine multiLevelHierarchyProcess
