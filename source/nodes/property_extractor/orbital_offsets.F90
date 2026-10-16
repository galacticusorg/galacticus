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

!+    Contributions to this file made by: Claude.

!!{RST
Contains a module which computes the position and velocity of a node relative to the top-level halo of its sub-halo hierarchy.
!!}

module Node_Orbital_Offsets
  !!{RST
  Computes the position and velocity of a node relative to the top-level halo of its sub-halo hierarchy.
  !!}
  implicit none
  private
  public :: Node_Orbital_Offset

contains

  subroutine Node_Orbital_Offset(node,time,position,velocity)
    !!{RST
    Compute the position and/or velocity of ``node`` relative to the top-level halo of its sub-halo hierarchy, by walking up
    through all host halos of the node and accumulating the position and velocity offsets of each from its host. If a host halo
    which is not a sub-halo is reached, and that halo is not on the main branch, the walk continues from the host halo, at the
    given ``time``, on the branch with which it will eventually merge.
    !!}
    use :: Galacticus_Nodes    , only : nodeComponentBasic, nodeComponentSatellite, treeNode
    use :: Numerical_Comparison, only : Values_Agree
    implicit none
    type            (treeNode              ), intent(inout), target                 :: node
    double precision                        , intent(in   )                         :: time
    double precision                        , intent(  out), dimension(3), optional :: position , velocity
    type            (treeNode              ), pointer                               :: nodeWork
    class           (nodeComponentBasic    ), pointer                               :: basic
    class           (nodeComponentSatellite), pointer                               :: satellite

    if (present(position)) position=0.0d0
    if (present(velocity)) velocity=0.0d0
    nodeWork => node
    do while (associated(nodeWork))
       satellite => nodeWork%satellite()
       if (present(position)) position=+position+satellite%position()
       if (present(velocity)) velocity=+velocity+satellite%velocity()
       if (nodeWork%isSatellite()) then
          ! Current node is a satellite, simply move to its parent.
          nodeWork => nodeWork%parent
       else
          ! Node is a host halo.
          if (nodeWork%isOnMainBranch()) then
             ! We are on the main branch - so we are done.
             nodeWork => null()
          else
             ! We are not on the main branch - find the host halo at this time on the branch with which we will merge.
             do while (nodeWork%isPrimaryProgenitor())
                nodeWork => nodeWork%parent
             end do
             nodeWork => nodeWork%parent
             basic    => nodeWork%basic ()
             do while (associated(nodeWork%firstChild).and..not.Values_Agree(basic%time(),time,relTol=1.0d-6))
                nodeWork => nodeWork%firstChild
                basic    => nodeWork%basic     ()
             end do
          end if
       end if
    end do
    return
  end subroutine Node_Orbital_Offset

end module Node_Orbital_Offsets
