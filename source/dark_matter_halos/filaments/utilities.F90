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

!+    Contributions to this file made by: Andrew Benson, Claude.

!!{RST
Contains a module which provides the axes of the cosmic web filaments feeding dark matter halos.
!!}

module Dark_Matter_Halo_Filaments
  !!{RST
  Provides the axes of the cosmic web filaments feeding dark matter halos.

  The axes of a halo are stored in a rank-1 ``float`` meta-property of its ``basic`` component, as the three direction cosines
  of each axis in turn, so that a halo fed by :math:`N` filaments stores :math:`3N` values. Each filament is a ray leaving the
  halo, so that a single filament passing through the halo is stored as two opposite axes. The halo has no orientation of its
  own, so the axes are expressed in the same arbitrary frame as the angles, :math:`(\theta,\phi)`, of the virial orbits of its
  satellites.

  The axes are shared by every node on a branch of the merger tree, so that a halo is fed along the same filaments throughout
  its history. They are drawn the first time that any node on the branch requires them, and stored on every node of the
  branch, from its earliest progenitor up to the node of which it is a secondary progenitor. They are therefore the same
  whether a satellite's host is identified with a node before or after the promotion of its primary progenitor, and the same
  however the nodes of the tree are ordered when they are first required. This matters because virial orbits are assigned at
  several points: during the construction of a tree, by tree operators, and when nodes are initialized.

  At present each branch is assigned a single filament passing through the halo, i.e. two opposite axes, with the filament
  oriented isotropically at random.
  !!}
  implicit none
  private
  public :: Dark_Matter_Halo_Filament_Axes

contains

  function Dark_Matter_Halo_Filament_Axes(node,filamentAxesID) result(axes)
    !!{RST
    Return the axes, as an array of shape :math:`(3,N)`, of the :math:`N` filaments feeding the halo of ``node``, drawing them
    for the branch of ``node`` if they do not yet exist. ``filamentAxesID`` is the identifier of the ``filamentAxes``
    meta-property of the ``basic`` component.
    !!}
    use :: Galacticus_Nodes        , only : nodeComponentBasic, treeNode
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision                    , allocatable  , dimension(:,:) :: axes
    type            (treeNode          ), intent(inout), target         :: node
    integer                             , intent(in   )                 :: filamentAxesID
    type            (treeNode          ), pointer                       :: nodeBranch
    class           (nodeComponentBasic), pointer                       :: basic
    double precision                    , allocatable  , dimension(:  ) :: directionCosines
    double precision                                   , dimension(3  ) :: axis
    double precision                                                    :: cosineTheta     , phi

    basic            => node %basic                    (              )
    directionCosines =  basic%floatRank1MetaPropertyGet(filamentAxesID)
    if (size(directionCosines) == 0) then
       ! No axes exist for this node. Find the earliest progenitor on its branch, which holds the axes of the branch if any
       ! node of the branch has been assigned them.
       nodeBranch => node
       do while (associated(nodeBranch%firstChild))
          nodeBranch => nodeBranch%firstChild
       end do
       basic            => nodeBranch%basic                    (              )
       directionCosines =  basic     %floatRank1MetaPropertyGet(filamentAxesID)
       if (size(directionCosines) == 0) then
          ! The branch has no axes, so draw a single filament passing through the halo, oriented isotropically.
          cosineTheta     =+2.0d0   *nodeBranch%hostTree%randomNumberGenerator_%uniformSample()-1.0d0
          phi             =+2.0d0*Pi*nodeBranch%hostTree%randomNumberGenerator_%uniformSample()
          axis            =[                                     &
               &            sqrt(1.0d0-cosineTheta**2)*cos(phi), &
               &            sqrt(1.0d0-cosineTheta**2)*sin(phi), &
               &                       cosineTheta               &
               &           ]
          directionCosines=[+axis,-axis]
       end if
       ! Store the axes on every node of the branch which lacks them, from its earliest progenitor upward.
       do while (associated(nodeBranch))
          basic => nodeBranch%basic()
          if (size(basic%floatRank1MetaPropertyGet(filamentAxesID)) == 0) &
               & call basic%floatRank1MetaPropertySet(filamentAxesID,directionCosines)
          if (.not.nodeBranch%isPrimaryProgenitor()) exit
          nodeBranch => nodeBranch%parent
       end do
    end if
    axes=reshape(directionCosines,[3,size(directionCosines)/3])
    return
  end function Dark_Matter_Halo_Filament_Axes

end module Dark_Matter_Halo_Filaments
