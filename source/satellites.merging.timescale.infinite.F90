!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

  !% Implements calculations of satellite merging times that are always infinite.
 
  !# <satelliteMergingTimescales name="satelliteMergingTimescalesInfinite">
  !#  <description>Returns an infinite timescale for merging.</description>
  !# </satelliteMergingTimescales>

  type, extends(satelliteMergingTimescalesClass) :: satelliteMergingTimescalesInfinite
     !% A class implementing satellite merging timescales that are always infinite.
     private
   contains
     final     ::                     infiniteDestructor
     procedure :: timeUntilMerging => infiniteTimeUntilMerging
  end type satelliteMergingTimescalesInfinite

  interface satelliteMergingTimescalesInfinite
     !% Constructors for the infinite merging timescale class.
     module procedure infiniteDefaultConstructor
  end interface satelliteMergingTimescalesInfinite

contains

  function infiniteDefaultConstructor()
    !% Default constructor for the infinite merging timescale class.
    use Galacticus_Display
    use Input_Parameters
    implicit none
    type(satelliteMergingTimescalesInfinite) :: infiniteDefaultConstructor

    return
  end function infiniteDefaultConstructor

  elemental subroutine infiniteDestructor(self)
    !% Default constructor for the infinite merging timescale class.
    implicit none
    type(satelliteMergingTimescalesInfinite), intent(inout) :: self

    ! Nothing to do.
    return
  end subroutine infiniteDestructor

  double precision function infiniteTimeUntilMerging(self,thisNode,thisOrbit)
    !% Return a zero timescale for satellite merging.
    use Galacticus_Nodes
    use Kepler_Orbits
    implicit none
    class           (satelliteMergingTimescalesInfinite), intent(inout)          :: self
    type            (treeNode                          ), intent(inout), pointer :: thisNode
    type            (keplerOrbit                       ), intent(inout)          :: thisOrbit
    double precision                                    , parameter              :: timeInfinite=1.0d30 ! Effective infinite time.

    infiniteTimeUntilMerging=timeInfinite
    return
  end function infiniteTimeUntilMerging

