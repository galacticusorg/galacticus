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

  !% Implements calculations of satellite merging times using preset values.

  !# <satelliteMergingTimescales name="satelliteMergingTimescalesPreset" />

  type, extends(satelliteMergingTimescalesClass) :: satelliteMergingTimescalesPreset
     !% A class implementing preset satellite merging timescales.
     private
   contains
     final     ::                     presetDestructor
     procedure :: timeUntilMerging => presetTimeUntilMerging
  end type satelliteMergingTimescalesPreset

  interface satelliteMergingTimescalesPreset
     !% Constructors for the preset merging timescale class.
     module procedure presetDefaultConstructor
  end interface satelliteMergingTimescalesPreset

contains

  function presetDefaultConstructor()
    !% Default constructor for the preset merging timescale class.
    use Galacticus_Display
    use Input_Parameters
    implicit none
    type(satelliteMergingTimescalesPreset) :: presetDefaultConstructor

    return
  end function presetDefaultConstructor

  elemental subroutine presetDestructor(self)
    !% Default constructor for the preset merging timescale class.
    implicit none
    type(satelliteMergingTimescalesPreset), intent(inout) :: self

    ! Nothing to do.
    return
  end subroutine presetDestructor

  double precision function presetTimeUntilMerging(self,thisNode,thisOrbit)
    !% Return the timescale for merging satellites using the preset value.
    use Galacticus_Nodes
    use Kepler_Orbits
    implicit none
    class(satelliteMergingTimescalesPreset), intent(inout)          :: self
    type (treeNode                        ), intent(inout), pointer :: thisNode
    type (keplerOrbit                     ), intent(inout)          :: thisOrbit
    class(nodeComponentSatellite          )               , pointer :: thisSatellite

    ! Simply return the current time until merging as, by definition, this has been preset if this method is being used.
    thisSatellite          => thisNode     %satellite()
    presetTimeUntilMerging =  thisSatellite%mergeTime()
    return
  end function presetTimeUntilMerging
