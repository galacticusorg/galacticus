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
  Provides a node merger class implementing a single level hierarchy.
  !!}

  !![
  <mergerTreeNodeMerger name="mergerTreeNodeMergerSingleLevelHierarchy">
   <description>
    A merger tree node merger class which maintains a single level hierarchy of substructure, i.e. it tracks only
    substructures, not sub-substructures or deeper levels. When a \gls{node} first becomes a satellite it is appended to the
    list of satellites associated with its host halo. If the \gls{node} contains its own satellites they will be detached from
    the \gls{node} and appended to the list of satellites of the new host (and assigned new merging times).
   </description>
  </mergerTreeNodeMerger>
  !!]
  type, extends(mergerTreeNodeMergerClass) :: mergerTreeNodeMergerSingleLevelHierarchy
     !!{
     Implementation of the standard merger tree evolver.
     !!}
     private
   contains
     procedure :: process  => singleLevelHierarchyProcess
  end type mergerTreeNodeMergerSingleLevelHierarchy

  interface mergerTreeNodeMergerSingleLevelHierarchy
     !!{
     Constructors for the \refClass{mergerTreeNodeMergerSingleLevelHierarchy} merger tree evolver.
     !!}
     module procedure singleLevelHierarchyConstructorParameters
  end interface mergerTreeNodeMergerSingleLevelHierarchy

contains

  function singleLevelHierarchyConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerTreeNodeMergerSingleLevelHierarchy} merger tree evolver class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(mergerTreeNodeMergerSingleLevelHierarchy)                :: self
    type(inputParameters                         ), intent(inout) :: parameters

    self=mergerTreeNodeMergerSingleLevelHierarchy()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function singleLevelHierarchyConstructorParameters

  subroutine singleLevelHierarchyProcess(self,node)
    !!{
    Processes a node merging event, utilizing a single level substructure hierarchy.
    !!}
    use :: Display            , only : displayGreen              , displayReset
    use :: Error              , only : Error_Report
    use :: Galacticus_Nodes   , only : treeNode
    use :: ISO_Varying_String , only : varying_string            , operator(//), assignment(=)
    use :: Satellite_Promotion, only : Satellite_Move_To_New_Host
    use :: String_Handling    , only : operator(//)
    implicit none
    class(mergerTreeNodeMergerSingleLevelHierarchy), intent(inout)          :: self
    type (treeNode                                ), intent(inout), target  :: node
    type (treeNode                                )               , pointer :: nodeChild    , nodeParent, &
         &                                                                     nodeSatellite
    type (varying_string                          )                         :: message
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
    ! Move any of its own satellites to become satellites of the parent and set their parent node pointers appropriately.
    do while (associated(node%firstSatellite))
       ! Move the satellite to the new parent.
       nodeSatellite => node%firstSatellite
       call Satellite_Move_To_New_Host(nodeSatellite,nodeParent)
    end do
    return
  end subroutine singleLevelHierarchyProcess
