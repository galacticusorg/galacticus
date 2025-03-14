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
Implements a filter which passes only constrained branch halos.
!!}

  !![
  <galacticFilter name="galacticFilterConstrainedBranch">
   <description>A filter which passes only constrained branch halos.</description>
  </galacticFilter>
  !!]
  type, extends(galacticFilterClass) :: galacticFilterConstrainedBranch
     !!{
     A galactic filter class which passes only constrained branch halos.
     !!}
     private
     integer :: isConstrainedID
   contains
     procedure :: passes => constrainedBranchPasses
  end type galacticFilterConstrainedBranch

  interface galacticFilterConstrainedBranch
     !!{
     Constructors for the ``constrainedBranch'' galactic filter class.
     !!}
     module procedure constrainedBranchConstructorParameters
     module procedure constrainedBranchConstructorInternal
  end interface galacticFilterConstrainedBranch

contains

  function constrainedBranchConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``constrainedBranch'' galactic filter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(galacticFilterConstrainedBranch)                :: self
    type(inputParameters                ), intent(inout) :: parameters

    self=galacticFilterConstrainedBranch()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function constrainedBranchConstructorParameters

  function constrainedBranchConstructorInternal() result(self)
    !!{
    Internal constructor for the ``constrainedBranch'' galactic filter class.
    !!}
    implicit none
    type(galacticFilterConstrainedBranch) :: self

    ! Use a matched "meta-property" as we used in the constrained tree build controller. This allows us to recover the stored
    ! "isConstrained" state of each node that was written by that object.
    !![
    <addMetaProperty component="basic" name="isConstrained" type="integer" id="self%isConstrainedID" isCreator="no"/>
    !!]
    return
  end function constrainedBranchConstructorInternal

  logical function constrainedBranchPasses(self,node) result(passes)
    !!{
    Implement a galactic filter which passes only constrained branch halos.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(galacticFilterConstrainedBranch), intent(inout)          :: self
    type (treeNode                       ), intent(inout), target  :: node
    class(nodeComponentBasic             )               , pointer :: basic

    basic  => node %basic                      (                    )
    passes =  basic%integerRank0MetaPropertyGet(self%isConstrainedID) == 1
    return
  end function constrainedBranchPasses
