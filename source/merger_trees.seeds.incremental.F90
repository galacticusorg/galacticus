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
Implements a merger tree random number seed in which the seed increases incrementally with tree index.
!!}

  !![
  <mergerTreeSeeds name="mergerTreeSeedsIncremental">
   <description>A merger tree random number seed in which the seed increases incrementally with tree index.</description>
  </mergerTreeSeeds>
  !!]
  type, extends(mergerTreeSeedsClass) :: mergerTreeSeedsIncremental
     !!{
     A merger tree random number seed in which the seed increases incrementally with tree index.
     !!}
     private
   contains
     procedure :: set => incrementalSet
  end type mergerTreeSeedsIncremental

  interface mergerTreeSeedsIncremental
     !!{
     Constructors for the {\normalfont \ttfamily incremental} merger tree seed class.
     !!}
     module procedure incrementalConstructorParameters
  end interface mergerTreeSeedsIncremental

contains

  function incrementalConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily incremental} merger tree seed class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(mergerTreeSeedsIncremental)                :: self
    type(inputParameters           ), intent(inout) :: parameters

    self=mergerTreeSeedsIncremental()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function incrementalConstructorParameters

  subroutine incrementalSet(self,tree)
    !!{
    Set the random number seed in the given {\normalfont \ttfamily tree}.
    !!}
    implicit none
    class(mergerTreeSeedsIncremental), intent(inout) :: self
    type (mergerTree                ), intent(inout) :: tree
    !$GLC attributes unused :: self

    ! Set the seed to the tree index, offset by the original seed.
    call tree%randomNumberGenerator_%seedSet(seed=tree%index,offset=.true.)
    return
  end subroutine incrementalSet
