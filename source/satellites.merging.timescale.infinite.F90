!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

  !% Implements a satellite merging timescale class in which merging timescales are always infinite.

  !# <satelliteMergingTimescales name="satelliteMergingTimescalesInfinite">
  !#  <description>Returns an infinite timescale for merging.</description>
  !# </satelliteMergingTimescales>
  type, extends(satelliteMergingTimescalesClass) :: satelliteMergingTimescalesInfinite
     !% A class implementing satellite merging timescales that are always infinite.
     private
   contains
     procedure :: timeUntilMerging => infiniteTimeUntilMerging
  end type satelliteMergingTimescalesInfinite

  interface satelliteMergingTimescalesInfinite
     !% Constructors for the {\normalfont \ttfamily infinite} satellite merging timescale class.
     module procedure infiniteConstructorParameters
  end interface satelliteMergingTimescalesInfinite

contains

  function infiniteConstructorParameters(parameters) result(self)
    !% A constructor for the {\normalfont \ttfamily infinite} satellite merging timescale class which builds the object from a
    !% parameter set.
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(satelliteMergingTimescalesInfinite)                :: self
    type(inputParameters                   ), intent(inout) :: parameters
    !GCC$ attributes unused :: parameters

    self=satelliteMergingTimescalesInfinite()
    return
  end function infiniteConstructorParameters

  double precision function infiniteTimeUntilMerging(self,node,orbit)
    !% Return a infinite timescale for satellite merging.
    implicit none
    class(satelliteMergingTimescalesInfinite), intent(inout) :: self
    type (treeNode                          ), intent(inout) :: node
    type (keplerOrbit                       ), intent(inout) :: orbit
    !GCC$ attributes unused :: self, node, orbit

    infiniteTimeUntilMerging=satelliteMergeTimeInfinite
    return
  end function infiniteTimeUntilMerging

