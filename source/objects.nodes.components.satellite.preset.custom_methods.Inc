!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


!% Contains custom functions for the prest satellite component.

double precision function Node_Component_Satellite_Preset_Time_Of_Merging(self)
  !% Return the {\tt timeOfMerging} property of the {\tt satellitePreset} component class.
  implicit none
  class(nodeComponentSatellitePreset), intent(inout) :: self
  type (treeNode                    ), pointer       :: selfNode

  selfNode => self%host()
  if (selfNode%isSatellite()) then
     Node_Component_Satellite_Preset_Time_Of_Merging=self%timeOfMergingData%value
  else
     Node_Component_Satellite_Preset_Time_Of_Merging=-1.0d0
  end if
  return
end function Node_Component_Satellite_Preset_Time_Of_Merging

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
  class(nodeComponentSatellitePreset), intent(inout) :: self
  double precision                   , intent(in   ) :: setValue
  type (treeNode                    ), pointer       :: selfNode
  class(nodeComponentBasic          ), pointer       :: selfBasicComponent

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
  class(nodeComponentSatellitePreset), intent(inout) :: self
  class(nodeComponentBasic          ), pointer       :: selfBasicComponent
  type (treeNode                    ), pointer       :: selfNode
  integer                                            :: iTime
  logical                                            :: interpolationReset
  type (history                     )                :: boundMassHistory
  type (fgsl_interp_accel           )                :: interpolationAccelerator
  
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
