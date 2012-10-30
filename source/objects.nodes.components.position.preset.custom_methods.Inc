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


!% Contains custom functions for the preset position component.

function PositionPresetPosition(self)
  !% Return the position of the node.
  use Numerical_Interpolation
  use FGSL
  implicit none
  double precision                  , dimension(:), allocatable  :: PositionPresetPosition
  class(nodeComponentPositionPreset), intent(inout)              :: self
  type (treeNode                   ), pointer                    :: selfNode
  class(nodeComponentBasic         ), pointer                    :: selfBasicComponent
  integer                           , parameter                  :: historyBegin=1,historyEnd=3
  integer                                                        :: iTime
  logical                                                        :: usingHistory,interpolationReset
  type (fgsl_interp_accel)                                       :: interpolationAccelerator
  type (history)                                                 :: positionHistory
  
  ! Check if this node has a position history attached to it.
  usingHistory   =.false.
  positionHistory=self%positionHistory()
  if (positionHistory%exists()) then
     selfNode           => self    %host ()
     selfBasicComponent => selfNode%basic()
     if (positionHistory%time(1) <= selfBasicComponent%time()) then
        interpolationReset=.true.
        iTime=Interpolate_Locate(size(positionHistory%time) ,positionHistory%time,interpolationAccelerator &
             &,selfBasicComponent%time(),reset=interpolationReset,closest=.true.)
        PositionPresetPosition=positionHistory%data(iTime,historyBegin:historyEnd)
        call Interpolate_Done(interpolationAccelerator=interpolationAccelerator,reset=interpolationReset)  
        usingHistory=.true.
     end if
  end if
  if (.not.usingHistory) PositionPresetPosition=self%positionData%value
  return
end function PositionPresetPosition

function PositionPresetVelocity(self)
  !% Return the velocity of the node.
  use Numerical_Interpolation
  use FGSL
  implicit none
  double precision                  , dimension(:), allocatable  :: PositionPresetVelocity
  class(nodeComponentPositionPreset), intent(inout)              :: self
  type (treeNode                   ), pointer                    :: selfNode
  class(nodeComponentBasic         ), pointer                    :: selfBasicComponent
  integer                           , parameter                  :: historyBegin=4,historyEnd=6
  integer                                                        :: iTime
  logical                                                        :: usingHistory,interpolationReset
  type (fgsl_interp_accel)                                       :: interpolationAccelerator
  type (history)                                                 :: positionHistory
  
  ! Check if this node has a position history attached to it.
  usingHistory   =.false.
  positionHistory=self%positionHistory()
  if (positionHistory%exists()) then
     selfNode           => self    %host ()
     selfBasicComponent => selfNode%basic()
     if (positionHistory%time(1) <= selfBasicComponent%time()) then
        interpolationReset=.true.
        iTime=Interpolate_Locate(size(positionHistory%time) ,positionHistory%time,interpolationAccelerator &
             &,selfBasicComponent%time(),reset=interpolationReset,closest=.true.)
        PositionPresetVelocity=positionHistory%data(iTime,historyBegin:historyEnd)
        call Interpolate_Done(interpolationAccelerator=interpolationAccelerator,reset=interpolationReset)  
        usingHistory=.true.
     end if
  end if
  if (.not.usingHistory) PositionPresetVelocity=self%velocityData%value
  return
end function PositionPresetVelocity
