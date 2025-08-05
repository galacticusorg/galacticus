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

!+    Contributions to this file made by:  Markus Haider.

  !!{
  Implements a satellite merging timescale class in which merging timescales are always infinite.
  !!}

  !![
  <satelliteMergingTimescales name="satelliteMergingTimescalesInfinite">
   <description>
    A satellite merging timescale class which always gives an infinite timescale for merging (technically, it returns a value
    close to the largest representable double precision floating point number which should be sufficiently close to infinity
    for practical purposes).
   </description>
  </satelliteMergingTimescales>
  !!]
  type, extends(satelliteMergingTimescalesClass) :: satelliteMergingTimescalesInfinite
     !!{
     A class implementing satellite merging timescales that are always infinite.
     !!}
     private
   contains
     procedure :: timeUntilMerging => infiniteTimeUntilMerging
  end type satelliteMergingTimescalesInfinite

  interface satelliteMergingTimescalesInfinite
     !!{
     Constructors for the \refClass{satelliteMergingTimescalesInfinite} satellite merging timescale class.
     !!}
     module procedure infiniteConstructorParameters
  end interface satelliteMergingTimescalesInfinite

contains

  function infiniteConstructorParameters(parameters) result(self)
    !!{
    A constructor for the {\normalfont \ttfamily infinite} satellite merging timescale class which builds the object from a
    parameter set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(satelliteMergingTimescalesInfinite)                :: self
    type(inputParameters                   ), intent(inout) :: parameters

    self=satelliteMergingTimescalesInfinite()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function infiniteConstructorParameters

  double precision function infiniteTimeUntilMerging(self,node,orbit)
    !!{
    Return a infinite timescale for satellite merging.
    !!}
    implicit none
    class(satelliteMergingTimescalesInfinite), intent(inout) :: self
    type (treeNode                          ), intent(inout) :: node
    type (keplerOrbit                       ), intent(inout) :: orbit
    !$GLC attributes unused :: self, node, orbit

    infiniteTimeUntilMerging=satelliteMergeTimeInfinite
    return
  end function infiniteTimeUntilMerging

