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
  Implements a satellite merging timescale class in which merging timescales are always zero.
  !!}

  !![
  <satelliteMergingTimescales name="satelliteMergingTimescalesZero">
   <description>
    A satellite merging timescale class which always gives a zero timescale for merging.
   </description>
  </satelliteMergingTimescales>
  !!]
  type, extends(satelliteMergingTimescalesClass) :: satelliteMergingTimescalesZero
     !!{
     A class implementing satellite merging timescales that are always zero.
     !!}
     private
   contains
     procedure :: timeUntilMerging => zeroTimeUntilMerging
  end type satelliteMergingTimescalesZero

  interface satelliteMergingTimescalesZero
     !!{
     Constructors for the \refClass{satelliteMergingTimescalesZero} satellite merging timescale class.
     !!}
     module procedure zeroConstructorParameters
  end interface satelliteMergingTimescalesZero

contains

  function zeroConstructorParameters(parameters) result(self)
    !!{
    A constructor for the {\normalfont \ttfamily zero} satellite merging timescale class which builds the object from a
    parameter set.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(satelliteMergingTimescalesZero)                :: self
    type(inputParameters               ), intent(inout) :: parameters

    self=satelliteMergingTimescalesZero()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function zeroConstructorParameters

  double precision function zeroTimeUntilMerging(self,node,orbit)
    !!{
    Return a zero timescale for satellite merging.
    !!}
    implicit none
    class(satelliteMergingTimescalesZero), intent(inout) :: self
    type (treeNode                      ), intent(inout) :: node
    type (keplerOrbit                   ), intent(inout) :: orbit
    !$GLC attributes unused :: self, node, orbit

    zeroTimeUntilMerging=0.0d0
    return
  end function zeroTimeUntilMerging
