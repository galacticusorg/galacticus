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
Implements a filter which passes only main branch halos.
!!}

  !![
  <galacticFilter name="galacticFilterMainBranch">
   <description>A filter which passes only main branch halos.</description>
  </galacticFilter>
  !!]
  type, extends(galacticFilterClass) :: galacticFilterMainBranch
     !!{
     A galactic filter class which passes only main branch halos.
     !!}
     private
   contains
     procedure :: passes => mainBranchPasses
  end type galacticFilterMainBranch

  interface galacticFilterMainBranch
     !!{
     Constructors for the \refClass{galacticFilterMainBranch} galactic filter class.
     !!}
     module procedure mainBranchConstructorParameters
  end interface galacticFilterMainBranch

contains

  function mainBranchConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{galacticFilterMainBranch} galactic filter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(galacticFilterMainBranch)                :: self
    type(inputParameters         ), intent(inout) :: parameters

    self=galacticFilterMainBranch()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function mainBranchConstructorParameters

  logical function mainBranchPasses(self,node)
    !!{
    Implement a galactic filter which passes only main branch halos.
    !!}
    implicit none
    class(galacticFilterMainBranch), intent(inout)         :: self
    type (treeNode                ), intent(inout), target :: node
    !$GLC attributes unused :: self

    mainBranchPasses=node%isOnMainBranch()
    return
  end function mainBranchPasses
