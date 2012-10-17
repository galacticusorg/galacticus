!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements calculations of dark matter halo formation times.

module Dark_Matter_Halo_Formation_Times
  !% Implements calculations of dark matter halo formation times.
  implicit none
  private
  public :: Dark_Matter_Halo_Formation_Time

contains

  double precision function Dark_Matter_Halo_Formation_Time(thisNode,formationMassFraction)
    !% Returns the time at which the main branch progenitor of {\tt thisNode} first had a mass equal to {\tt
    !% formationMassFraction} of the current mass.
    use Tree_Nodes
    use Dark_Matter_Halo_Mass_Accretion_Histories
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, intent(in)             :: formationMassFraction
    type(treeNode),                  pointer :: workNode,formationNode
    double precision                         :: timeNode,massNode

    timeNode=Tree_Node_Time(thisNode)
    massNode=Tree_Node_Mass(thisNode)

    workNode => thisNode
    do while (associated(workNode))
       formationNode => workNode
       if (Tree_Node_Mass(workNode) <= formationMassFraction*massNode) exit
       workNode => workNode%childNode
    end do
    if (.not.associated(workNode)) then
       ! Find the formation time based on the mass accretion history.
       Dark_Matter_Halo_Formation_Time=Dark_Matter_Halo_Mass_Accretion_Time(formationNode,formationMassFraction*massNode)
    else
       ! Interpolate to get the exact time of formation.
       Dark_Matter_Halo_Formation_Time=                                      Tree_Node_Time(workNode)  &
            &                          +(Tree_Node_Time(workNode%parentNode)-Tree_Node_Time(workNode)) &
            &                          *(formationMassFraction*massNode     -Tree_Node_Mass(workNode)) &
            &                          /(Tree_Node_Mass(workNode%parentNode)-Tree_Node_Mass(workNode))
    end if

    return
  end function Dark_Matter_Halo_Formation_Time

end module Dark_Matter_Halo_Formation_Times
