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
Implements a filter which passes only halos that have always been isolated.
!!}

  !![
  <galacticFilter name="galacticFilterHaloAlwaysIsolated">
   <description>A filter which passes only halos that have always been isolated.</description>
  </galacticFilter>
  !!]
  type, extends(galacticFilterClass) :: galacticFilterHaloAlwaysIsolated
     !!{
     A galactic filter class which passes only halos that have always been isolated.
     !!}
     private
     integer :: nodeHierarchyLevelMaximumID
   contains
     procedure :: passes => haloAlwaysIsolatedPasses
  end type galacticFilterHaloAlwaysIsolated

  interface galacticFilterHaloAlwaysIsolated
     !!{
     Constructors for the {\normalfont \ttfamily haloAlwaysIsolated} galactic filter class.
     !!}
     module procedure haloAlwaysIsolatedConstructorParameters
     module procedure haloAlwaysIsolatedConstructorInternal
  end interface galacticFilterHaloAlwaysIsolated

contains

  function haloAlwaysIsolatedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily haloAlwaysIsolated} galactic filter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(galacticFilterHaloAlwaysIsolated)                :: self
    type(inputParameters                 ), intent(inout) :: parameters

    self=galacticFilterHaloAlwaysIsolated()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function haloAlwaysIsolatedConstructorParameters

  function haloAlwaysIsolatedConstructorInternal() result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily haloAlwaysIsolated} galactic filter class.
    !!}
    implicit none
    type(galacticFilterHaloAlwaysIsolated) :: self

    !![
    <addMetaProperty component="basic" name="nodeHierarchyLevelMaximum" type="integer" id="self%nodeHierarchyLevelMaximumID"/>
    !!]
    return
  end function haloAlwaysIsolatedConstructorInternal

  logical function haloAlwaysIsolatedPasses(self,node)
    !!{
    Implement a galactic filter which passes only isolated halos.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(galacticFilterHaloAlwaysIsolated), intent(inout)         :: self
    type (treeNode                        ), intent(inout), target :: node
    class(nodeComponentBasic              ), pointer               :: basic

    basic                    => node %basic                      (                                )
    haloAlwaysIsolatedPasses =  basic%integerRank0MetaPropertyGet(self%nodeHierarchyLevelMaximumID) == 0
    return
  end function haloAlwaysIsolatedPasses
