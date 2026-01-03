!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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

!!{
Contains custom functions for the preset satellite component.
!!}

double precision function SatellitePresetMergeBoundMass(self)
  !!{
  Return the satellite bound mass at the current time.
  !!}
  use            :: Histories              , only : history
  use, intrinsic :: ISO_C_Binding          , only : c_size_t
  use            :: Numerical_Interpolation, only : interpolator
  implicit none
  class  (nodeComponentSatellitePreset), intent(inout) :: self
  class  (nodeComponentBasic          ), pointer       :: selfBasicComponent
  type   (treeNode                    ), pointer       :: selfNode
  integer(c_size_t                    )                :: iTime
  type   (history                     )                :: boundMassHistory
  type   (interpolator                )                :: interpolator_

  ! Get the hosting node and the basic component.
  selfNode           => self    %host ()
  selfBasicComponent => selfNode%basic()
  ! Determine if this is a satellite node or not.
  if (selfNode%isSatellite()) then
     ! Node is a satellite, so get the preset bound mass history.
     boundMassHistory=self%boundMassHistory()
     ! Check if the history exists.
     if (boundMassHistory%exists()) then
        ! It does, so find the preset time closest to the present one and return the mass at that time.
        interpolator_=interpolator(boundMassHistory%time)
        iTime=interpolator_%locate(selfBasicComponent%time(),closest=.true.)
        SatellitePresetMergeBoundMass=boundMassHistory%data(iTime,1)
     else
        ! It does not, so return the basic mass.
        SatellitePresetMergeBoundMass=selfBasicComponent%mass()
     end if
  else
     ! Node is not a satellite, return the basic mass of the halo.
     SatellitePresetMergeBoundMass=selfBasicComponent%mass()
  end if
  call boundMassHistory%destroy()
  return
end function SatellitePresetMergeBoundMass

function SatellitePresetNodeIndex(self)
  !!{
  Return the satellite node index.
  !!}
  use            :: Histories              , only : longIntegerHistory
  use, intrinsic :: ISO_C_Binding          , only : c_size_t
  use            :: Kind_Numbers           , only : kind_int8
  use            :: Numerical_Interpolation, only : interpolator
  implicit none
  integer(kind=kind_int8              )                :: SatellitePresetNodeIndex
  class  (nodeComponentSatellitePreset), intent(inout) :: self
  class  (nodeComponentBasic          ), pointer       :: selfBasicComponent
  type   (treeNode                    ), pointer       :: selfNode
  integer(c_size_t                    )                :: iTime
  type   (longIntegerHistory          )                :: nodeIndexHistory
  type   (interpolator                )                :: interpolator_

  ! Get the node.
  selfNode         => self%host            ()
  ! Get the preset history.
  nodeIndexHistory =  self%nodeIndexHistory()
  ! For isolated nodes, or nodes with no nodeIndex history, return the current node index.
  if (.not.selfNode%isSatellite().or..not.nodeIndexHistory%exists()) then
     ! Return the standard node index.
     SatellitePresetNodeIndex=selfNode%index()
  else
     ! Node is a satellite with a node index history, so find the preset time closest to the
     ! present one and return the index at that time.
     selfBasicComponent => selfNode%basic()
     interpolator_=interpolator(nodeIndexHistory%time)
     iTime=interpolator_%locate(selfBasicComponent%time(),closest=.true.)
     SatellitePresetNodeIndex=nodeIndexHistory%data(iTime,1)
  end if
  call nodeIndexHistory%destroy()
  return
end function SatellitePresetNodeIndex
