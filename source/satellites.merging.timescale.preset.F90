!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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

  !% Implements a satellite merging timescale class which uses preset values for the timescale.

  !# <satelliteMergingTimescales name="satelliteMergingTimescalesPreset">
  !#  <description>This class assumes that merging times have been preset for every node (or, at least, every node which becomes a satellite). It therefore simply returns the preset merging time.</description>
  !# </satelliteMergingTimescales>
  type, extends(satelliteMergingTimescalesClass) :: satelliteMergingTimescalesPreset
     !% A class implementing preset satellite merging timescales.
     private
   contains
     procedure :: timeUntilMerging => presetTimeUntilMerging
  end type satelliteMergingTimescalesPreset

  interface satelliteMergingTimescalesPreset
     !% Constructors for the {\normalfont \ttfamily preset} satellite merging timescale class.
     module procedure presetConstructorParameters
  end interface satelliteMergingTimescalesPreset

contains

  function presetConstructorParameters(parameters) result(self)
    !% A constructor for the {\normalfont \ttfamily preset} satellite merging timescale class which builds the object from a
    !% parameter set.
    use Input_Parameters
    implicit none
    type(satelliteMergingTimescalesPreset)                :: self
    type(inputParameters                 ), intent(inout) :: parameters
    !GCC$ attributes unused :: parameters

    self=satelliteMergingTimescalesPreset()
    return
  end function presetConstructorParameters

  double precision function presetTimeUntilMerging(self,node,orbit)
    !% Return the timescale for merging satellites using the preset value.
    use Galacticus_Nodes, only : nodeComponentSatellite
    use Kepler_Orbits
    implicit none
    class(satelliteMergingTimescalesPreset), intent(inout) :: self
    type (treeNode                        ), intent(inout) :: node
    type (keplerOrbit                     ), intent(inout) :: orbit
    class(nodeComponentSatellite          ), pointer       :: satellite
    !GCC$ attributes unused :: self, orbit
    
    ! Simply return the current time until merging as, by definition, this has been preset if this method is being used.
    satellite              => node     %satellite()
    presetTimeUntilMerging =  satellite%mergeTime()
    return
  end function presetTimeUntilMerging
