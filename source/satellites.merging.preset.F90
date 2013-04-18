!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements calculations of satellite merging times using preset values.

module Satellite_Merging_Times_Preset
  !% Implements calculations of satellite merging times using preset values.
  implicit none
  private
  public :: Satellite_Time_Until_Merging_Preset_Initialize

contains

  !# <satelliteMergingMethod>
  !#  <unitName>Satellite_Time_Until_Merging_Preset_Initialize</unitName>
  !# </satelliteMergingMethod>
  subroutine Satellite_Time_Until_Merging_Preset_Initialize(satelliteMergingMethod,Satellite_Time_Until_Merging)
    !% Determine if this method is to be used and set pointer appropriately.
    use ISO_Varying_String
    implicit none
    type(varying_string),                 intent(in)    :: satelliteMergingMethod
    procedure(Satellite_Time_Until_Merging_Preset), pointer, intent(inout) :: Satellite_Time_Until_Merging

    if (satelliteMergingMethod == 'preset') Satellite_Time_Until_Merging => Satellite_Time_Until_Merging_Preset
    return
  end subroutine Satellite_Time_Until_Merging_Preset_Initialize

  double precision function Satellite_Time_Until_Merging_Preset(thisNode,thisOrbit)
    !% Return the timescale for merging satellites using the preset value.
    use Galacticus_Nodes
    use Kepler_Orbits
    implicit none
    type (treeNode              ), pointer, intent(inout) :: thisNode
    type (keplerOrbit),                 intent(inout) :: thisOrbit
    class(nodeComponentSatellite), pointer                :: thisSatelliteComponent

    ! Simply return the current time until merging as, by definition, this has been preset if this method is being used.
    thisSatelliteComponent => thisNode%satellite()
    Satellite_Time_Until_Merging_Preset=thisSatelliteComponent%mergeTime()
    return
  end function Satellite_Time_Until_Merging_Preset

end module Satellite_Merging_Times_Preset
