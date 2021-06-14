!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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

!% Contains custom functions for the preset$+$orphans position component.

function PositionPresetOrphansPosition(self)
  !% Return the position of the node.
  implicit none
  double precision                                    , allocatable  , dimension(:) :: PositionPresetOrphansPosition
  class           (nodeComponentPositionPresetOrphans), intent(inout)               :: self
  type            (treeNode                          ), pointer                     :: node
  class           (nodeComponentBasic                ), pointer                     :: basic
  type            (history                           )                              :: positionHistory
  logical                                                                           :: positionIsKnown

  ! Check if the node has a known, preset position.
  positionIsKnown =  .false.
  node            => self%host()
  if (node%isSatellite()) then
     positionHistory=self%positionHistory()
     if (positionHistory%exists()) then
        basic => node%basic()
        positionIsKnown=positionHistory%time(size(positionHistory%time)) >= basic%time()
     else
        positionIsKnown=.false.
     end if
  else
     positionIsKnown=.true.
  end if
  ! If the node has no known position (i.e. is an orphan), check if it was previously assigned a position at this time. If so,
  ! return that position, otherwise, assign one now. If position is known, simply call the position method of our parent class.
  if (positionIsKnown) then
     PositionPresetOrphansPosition=self%nodeComponentPositionPreset%position      ()
  else
     PositionPresetOrphansPosition=self                            %positionOrphan()
  end if
  return
end function PositionPresetOrphansPosition

function PositionPresetOrphansVelocity(self)
  !% Return the velocity of the node.
  implicit none
  double precision                                    , allocatable  , dimension(:) :: PositionPresetOrphansVelocity
  class           (nodeComponentPositionPresetOrphans), intent(inout)               :: self
  type            (treeNode                          ), pointer                     :: node
  class           (nodeComponentBasic                ), pointer                     :: basic
  type            (history                           )                              :: positionHistory
  logical                                                                           :: positionIsKnown

  ! Check if the node has a known, preset position.
  positionIsKnown =  .false.
  node            => self%host()
  if (node%isSatellite()) then
     positionHistory=self%positionHistory()
     if (positionHistory%exists()) then
        basic => node%basic()
        positionIsKnown=positionHistory%time(size(positionHistory%time)) >= basic%time()
     else
        positionIsKnown=.false.
     end if
  else
     positionIsKnown=.true.
  end if
  ! If the node has no known position (i.e. is an orphan), check if it was previously assigned a position at this time. If so,
  ! return that position, otherwise, assign one now. If position is known, simply call the velocity method of our parent class.
  if (positionIsKnown) then
     PositionPresetOrphansVelocity=self%nodeComponentPositionPreset%velocity      ()
  else
     PositionPresetOrphansVelocity=self                            %velocityOrphan()
  end if
  return
end function PositionPresetOrphansVelocity
