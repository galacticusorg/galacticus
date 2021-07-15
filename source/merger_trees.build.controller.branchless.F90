!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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
Contains a module which implements a merger tree build controller class which builds branchless trees.
!!}

  !![
  <mergerTreeBuildController name="mergerTreeBuildControllerBranchless">
   <description>A merger tree build controller class which builds branchless trees.</description>
  </mergerTreeBuildController>
  !!]
  type, extends(mergerTreeBuildControllerClass) :: mergerTreeBuildControllerBranchless
     !!{     
     A merger tree build controller class which builds branchless trees. Specifically, whenever a branching occurs, only the
     primary progenitor is followed. This results in a tree which consists of the main branch, along with stubs of any branches
     off of the main branch, but those stubs do not grow full branches of their own.
     !!}
     private
  contains
     procedure :: control => controlBranchless
  end type mergerTreeBuildControllerBranchless

  interface mergerTreeBuildControllerBranchless
     !!{
     Constructors for the ``branchless'' merger tree build controller class.
     !!}
     module procedure branchlessConstructorParameters
  end interface mergerTreeBuildControllerBranchless

contains

  function branchlessConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``branchless'' merger tree build controller class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(mergerTreeBuildControllerBranchless)                :: self
    type(inputParameters                    ), intent(inout) :: parameters
  
    self=mergerTreeBuildControllerBranchless()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function branchlessConstructorParameters

  logical function controlBranchless(self,node,treeWalker_)
    !!{
    Skip side branches of a tree under construction.
    !!}
    implicit none
    class(mergerTreeBuildControllerBranchless), intent(inout)          :: self    
    type (treeNode                           ), intent(inout), pointer :: node
    class(mergerTreeWalkerClass              ), intent(inout)          :: treeWalker_
    !$GLC attributes unused :: self

    controlBranchless=.true.
    ! Move to the next node in the tree while such exists, and the current node is on a side branch.
    do while (controlBranchless.and.associated(node%parent).and..not.node%isPrimaryProgenitor())
       controlBranchless=treeWalker_%next(node)
    end do
    return
  end function controlBranchless
