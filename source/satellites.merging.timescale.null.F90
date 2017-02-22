!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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

  !% Implements calculations of satellite merging times that are always zero.
 
  !# <satelliteMergingTimescales name="satelliteMergingTimescalesNull">
  !#  <description>Returns a zero timescale for merging.</description>
  !# </satelliteMergingTimescales>

  type, extends(satelliteMergingTimescalesClass) :: satelliteMergingTimescalesNull
     !% A class implementing satellite merging timescales that are always zero.
     private
   contains
     final     ::                     nullDestructor
     procedure :: timeUntilMerging => nullTimeUntilMerging
  end type satelliteMergingTimescalesNull

contains

  elemental subroutine nullDestructor(self)
    !% Default constructor for the nul merging timescale class.
    implicit none
    type(satelliteMergingTimescalesNull), intent(inout) :: self
    !GCC$ attributes unused :: self
    
    ! Nothing to do.
    return
  end subroutine nullDestructor

  double precision function nullTimeUntilMerging(self,node,orbit)
    !% Return a zero timescale for satellite merging.
    use Galacticus_Nodes
    use Kepler_Orbits
    implicit none
    class(satelliteMergingTimescalesNull), intent(inout) :: self
    type (treeNode                      ), intent(inout) :: node
    type (keplerOrbit                   ), intent(inout) :: orbit
    !GCC$ attributes unused :: self, node, orbit
    
    nullTimeUntilMerging=0.0d0
    return
  end function nullTimeUntilMerging
