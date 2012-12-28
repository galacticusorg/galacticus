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

!% Contains a module which handles events where a satellite is moved to a new host halo.

module Satellite_Promotion
  !% Handles events where a satellite is moved to a new host halo.
  implicit none
  private
  public :: Satellite_Move_To_New_Host
  
contains

  subroutine Satellite_Move_To_New_Host(satelliteNode,newHostNode)
    !% Move {\tt satelliteNode} to be a satellite of {\tt newHostNode}.
    use Galacticus_Nodes
    !# <include directive="satelliteHostChangeTask" type="moduleUse">
    include 'satellites.structures.host_change.moduleUse.inc'
    !# </include>
    implicit none
    type(treeNode), pointer, intent(inout) :: satelliteNode,newHostNode
    type(treeNode), pointer                :: lastSatelliteNode

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

    ! Allow arbitrary routines to process the host change event.
    !# <include directive="satelliteHostChangeTask" type="functionCall" functionType="void">
    !#  <functionArgs>satelliteNode</functionArgs>
    include 'satellites.structures.host_change.inc'
    !# </include>

    return
  end subroutine Satellite_Move_To_New_Host

end module Satellite_Promotion
