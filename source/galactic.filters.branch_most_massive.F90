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
Implements a filter which passes only the most massive branch halo.
!!}

  !![
  <galacticFilter name="galacticFilterBranchMostMassive">
   <description>A filter which passes only the most massive branch halo.</description>
  </galacticFilter>
  !!]
  type, extends(galacticFilterClass) :: galacticFilterBranchMostMassive
     !!{
     A galactic filter class which passes only the most massive branch halo.
     !!}
     private
     integer :: isMostMassiveBranchID
   contains
     procedure :: passes => branchMostMassivePasses
  end type galacticFilterBranchMostMassive

  interface galacticFilterBranchMostMassive
     !!{
     Constructors for the \refClass{galacticFilterBranchMostMassive} galactic filter class.
     !!}
     module procedure branchMostMassiveConstructorParameters
     module procedure branchMostMassiveConstructorInternal
  end interface galacticFilterBranchMostMassive

contains

  function branchMostMassiveConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{galacticFilterBranchMostMassive} galactic filter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(galacticFilterBranchMostMassive)                :: self
    type(inputParameters                ), intent(inout) :: parameters

    self=galacticFilterBranchMostMassive()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function branchMostMassiveConstructorParameters

  function branchMostMassiveConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{galacticFilterBranchMostMassive} galactic filter class.
    !!}
    implicit none
    type(galacticFilterBranchMostMassive) :: self

    !![
    <addMetaProperty component="basic" name="isMostMassiveBranch" type="integer" id="self%isMostMassiveBranchID" isCreator="no"/>
    !!]
    return
  end function branchMostMassiveConstructorInternal

  logical function branchMostMassivePasses(self,node)
    !!{
    Implement a galactic filter which passes only main branch halos.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(galacticFilterBranchMostMassive), intent(inout)          :: self
    type (treeNode                       ), intent(inout), target  :: node
    class(nodeComponentBasic             )               , pointer :: basic

    basic                   => node %basic                      (                          )
    branchMostMassivePasses =  basic%integerRank0MetaPropertyGet(self%isMostMassiveBranchID) == 1
    return
  end function branchMostMassivePasses
