!% Contains a module which handles events where a satellite is moved to a new host halo.

module Satellite_Promotion
  !% Handles events where a satellite is moved to a new host halo.
  private
  public :: Satellite_Move_To_New_Host
  
contains

  subroutine Satellite_Move_To_New_Host(satelliteNode,newHostNode)
    !% Move {\tt satelliteNode} to be a satellite of {\tt newHostNode}.
    use Tree_Nodes
    !# <include directive="satelliteHostChangeTask" type="moduleUse">
    include 'satellites.structures.host_change.moduleUse.inc'
    !# </include>
    implicit none
    type(treeNode), pointer, intent(inout) :: satelliteNode,newHostNode
    type(treeNode), pointer                :: lastSatelliteNode

    ! First remove from its current host.
    call satelliteNode%removeFromHost()
    ! Find attachment point for new host.
    if (associated(newHostNode%satelliteNode)) then
       call newHostNode%lastSatellite(lastSatelliteNode)
       lastSatelliteNode%siblingNode => satelliteNode
    else
       newHostNode%satelliteNode => satelliteNode
    end if
    ! Set parent and sibling pointers.
    satelliteNode%parentNode  => newHostNode
    satelliteNode%siblingNode => null()

    ! Allow arbitrary routines to process the host change event.
    !# <include directive="satelliteHostChangeTask" type="code" action="subroutine">
    !#  <subroutineArgs>satelliteNode</subroutineArgs>
    include 'satellites.structures.host_change.inc'
    !# </include>

    return
  end subroutine Satellite_Move_To_New_Host

end module Satellite_Promotion
