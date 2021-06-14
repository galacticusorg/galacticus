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

!% Contains custom functions for the preset position component.

function PositionPresetPosition(self)
  !% Return the position of the node.
  use, intrinsic :: ISO_C_Binding          , only : c_size_t
  use            :: Numerical_Interpolation, only : interpolator
  implicit none
  double precision                             , allocatable  , dimension(:) :: PositionPresetPosition
  class           (nodeComponentPositionPreset), intent(inout)               :: self
  type            (treeNode                   ), pointer                     :: selfNode
  class           (nodeComponentBasic         ), pointer                     :: selfBasicComponent
  integer                                      , parameter                   :: historyBegin            =1, historyEnd=3
  integer         (c_size_t                   )                              :: iTime
  logical                                                                    :: usingHistory
  type            (interpolator               )                              :: interpolator_
  type            (history                    )                              :: positionHistory

  ! Check if this node has a position history attached to it.
  usingHistory   =.false.
  positionHistory=self%positionHistory()
  if (positionHistory%exists()) then
     selfNode           => self    %host ()
     selfBasicComponent => selfNode%basic()
     if (positionHistory%time(1) <= selfBasicComponent%time()) then
        interpolator_=interpolator(positionHistory%time)
        iTime=interpolator_%locate(selfBasicComponent%time(),closest=.true.)
        PositionPresetPosition=positionHistory%data(iTime,historyBegin:historyEnd)
        usingHistory=.true.
     end if
  end if
  if (.not.usingHistory) PositionPresetPosition=self%positionData
  return
end function PositionPresetPosition

function PositionPresetVelocity(self)
  !% Return the velocity of the node.
  use, intrinsic :: ISO_C_Binding          , only : c_size_t
  use            :: Numerical_Interpolation, only : interpolator
  implicit none
  double precision                             , allocatable  , dimension(:) :: PositionPresetVelocity
  class           (nodeComponentPositionPreset), intent(inout)               :: self
  type            (treeNode                   ), pointer                     :: selfNode
  class           (nodeComponentBasic         ), pointer                     :: selfBasicComponent
  integer                                      , parameter                   :: historyBegin            =4, historyEnd  =6
  integer         (c_size_t                   )                              :: iTime
  logical                                                                    :: usingHistory
  type            (interpolator               )                              :: interpolator_
  type            (history                    )                              :: positionHistory

  ! Check if this node has a position history attached to it.
  usingHistory   =.false.
  positionHistory=self%positionHistory()
  if (positionHistory%exists()) then
     selfNode           => self    %host ()
     selfBasicComponent => selfNode%basic()
     if (positionHistory%time(1) <= selfBasicComponent%time()) then
        interpolator_=interpolator(positionHistory%time)
        iTime=interpolator_%locate(selfBasicComponent%time(),closest=.true.)
        PositionPresetVelocity=positionHistory%data(iTime,historyBegin:historyEnd)
        usingHistory=.true.
     end if
  end if
  if (.not.usingHistory) PositionPresetVelocity=self%velocityData
  return
end function PositionPresetVelocity
