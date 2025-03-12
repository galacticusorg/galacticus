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
Implements a filter which passes only non-isolated halos.
!!}

  !![
  <galacticFilter name="galacticFilterHaloNotIsolated">
   <description>A filter which passes only non-isolated halos.</description>
  </galacticFilter>
  !!]
  type, extends(galacticFilterClass) :: galacticFilterHaloNotIsolated
     !!{
     A galactic filter class which passes only non-isolated halos.
     !!}
     private
   contains
     procedure :: passes => haloNotIsolatedPasses
  end type galacticFilterHaloNotIsolated

  interface galacticFilterHaloNotIsolated
     !!{
     Constructors for the ``haloNotIsolated'' galactic filter class.
     !!}
     module procedure haloNotIsolatedConstructorParameters
  end interface galacticFilterHaloNotIsolated

contains

  function haloNotIsolatedConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``haloNotIsolated'' galactic filter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(galacticFilterHaloNotIsolated)                :: self
    type(inputParameters              ), intent(inout) :: parameters

    self=galacticFilterHaloNotIsolated()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function haloNotIsolatedConstructorParameters

  logical function haloNotIsolatedPasses(self,node)
    !!{
    Implement a galactic filter which passes only isolated halos.
    !!}
    implicit none
    class(galacticFilterHaloNotIsolated), intent(inout)         :: self
    type (treeNode                     ), intent(inout), target :: node
    !$GLC attributes unused :: self

    haloNotIsolatedPasses=node%isSatellite()
    return
  end function haloNotIsolatedPasses
