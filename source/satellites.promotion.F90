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
Contains a module which handles events where a satellite is moved to a new host halo.
!!}

module Satellite_Promotion
  !!{
  Handles events where a satellite is moved to a new host halo.
  !!}
  implicit none
  private
  public :: Satellite_Move_To_New_Host

contains

  subroutine Satellite_Move_To_New_Host(satelliteNode,newHostNode)
    !!{
    Move {\normalfont \ttfamily satelliteNode} to be a satellite of {\normalfont \ttfamily newHostNode}.
    !!}
    use :: Display           , only : displayMessage    , displayVerbosity, verbosityLevelInfo
    use :: Galacticus_Nodes  , only : nodeComponentBasic, treeNode
    use :: ISO_Varying_String, only : assignment(=)     , operator(//)    , varying_string
    use :: String_Handling   , only : operator(//)
    implicit none
    type     (treeNode          ), intent(inout), target  :: satelliteNode    , newHostNode
    type     (treeNode          )               , pointer :: lastSatelliteNode
    class    (nodeComponentBasic), pointer                :: basic
    type     (varying_string    )                         :: message
    character(len=12            )                         :: label

    ! Report if necessary.
    if (displayVerbosity() >= verbosityLevelInfo) then
       basic => satelliteNode%basic()
       write (label,'(f12.6)') basic%time()
       message='Satellite node ['
       message=message//satelliteNode%index()//'] is being promoted to new host node ['//newHostNode%index()//'] at time '//trim(label)//' Gyr'
       call displayMessage(message)
    end if
    ! Allow arbitrary functions to act prior to the host change event.
    !![
    <eventHook name="satellitePreHostChange">
      <import>
	<module name="Galacticus_Nodes" symbols="treeNode"/>
      </import>
      <interface>
	type(treeNode), intent(inout), target :: satelliteNode, newHostNode
      </interface>
     <callWith>satelliteNode,newHostNode</callWith>
    </eventHook>
    !!]
    ! First remove from its current host.
    call satelliteNode%removeFromHost()
    ! Find attachment point for new host.
    if (associated(newHostNode%firstSatellite)) then
       lastSatelliteNode          => newHostNode%lastSatellite()
       lastSatelliteNode%sibling  => satelliteNode
    else
       newHostNode%firstSatellite => satelliteNode
    end if
    ! Set parent and sibling pointers.
    satelliteNode%parent  => newHostNode
    satelliteNode%sibling => null()
    ! Allow arbitrary functions to process the host change event.
    !![
    <eventHook name="satelliteHostChange">
      <import>
	<module name="Galacticus_Nodes" symbols="treeNode"/>
      </import>
      <interface>
	type(treeNode), intent(inout), target :: satelliteNode
      </interface>
     <callWith>satelliteNode</callWith>
    </eventHook>
    !!]
    return
  end subroutine Satellite_Move_To_New_Host

end module Satellite_Promotion
