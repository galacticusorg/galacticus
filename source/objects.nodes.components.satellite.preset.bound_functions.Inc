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

!% Contains custom functions for the prest satellite component.

double precision function Node_Component_Satellite_Preset_Merge_Time(self)
  !% Return the {\tt mergeTime} property of the {\tt satellitePreset} component class.
  implicit none
  class(nodeComponentSatellitePreset), intent(inout) :: self
  type (treeNode                    ), pointer       :: selfNode
  class(nodeComponentBasic          ), pointer       :: selfBasicComponent

  selfNode => self%host()
  if (selfNode%isSatellite()) then
     selfBasicComponent => selfNode%basic()
     Node_Component_Satellite_Preset_Merge_Time=max(self%timeOfMergingData%value-selfBasicComponent%time(),0.0d0)
  else
     Node_Component_Satellite_Preset_Merge_Time=-1.0d0
  end if
  return
end function Node_Component_Satellite_Preset_Merge_Time

subroutine Node_Component_Satellite_Preset_Merge_Time_Set(self,setValue)
  !% Sets the {\tt mergeTime} property of the {\tt satellitePreset} component class.
  implicit none
  class           (nodeComponentSatellitePreset)         , intent(inout) :: self
  double precision                                       , intent(in   ) :: setValue
  type            (treeNode                    ), pointer                :: selfNode
  class           (nodeComponentBasic          ), pointer                :: selfBasicComponent

  selfNode           => self    %host ()
  selfBasicComponent => selfNode%basic()
  self%timeOfMergingData%value=setValue+selfBasicComponent%time()
  return
end subroutine Node_Component_Satellite_Preset_Merge_Time_Set

double precision function SatellitePresetMergeBoundMass(self)
  !% Return the satellite bound mass at the current time.
  use Histories
  use FGSL
  use Numerical_Interpolation
  implicit none
  class  (nodeComponentSatellitePreset), intent(inout) :: self
  class  (nodeComponentBasic          ), pointer       :: selfBasicComponent
  type   (treeNode                    ), pointer       :: selfNode
  integer                                              :: iTime
  logical                                              :: interpolationReset
  type   (history                     )                :: boundMassHistory
  type   (fgsl_interp_accel           )                :: interpolationAccelerator

  ! It does, so get the preset bound mass history.
  boundMassHistory=self%boundMassHistory()
  ! Check if the history exists.
  if (boundMassHistory%exists()) then
     ! It does, so find the preset time closest to the present one and return the mass at that time.
     selfNode           => self    %host ()
     selfBasicComponent => selfNode%basic()
     interpolationReset=.true.
     iTime=Interpolate_Locate(size(boundMassHistory%time),boundMassHistory%time,interpolationAccelerator&
          &,selfBasicComponent%time(),reset=interpolationReset,closest=.true.)
     SatellitePresetMergeBoundMass=boundMassHistory%data(iTime,1)
     call Interpolate_Done(interpolationAccelerator=interpolationAccelerator,reset=interpolationReset)
  else
     ! It does not, so return zero.
     SatellitePresetMergeBoundMass=0.0d0
  end if
  call boundMassHistory%destroy()
  return
end function SatellitePresetMergeBoundMass
