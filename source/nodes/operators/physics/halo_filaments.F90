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
  Implements a node operator which assigns each halo the axes of the cosmic web filaments which feed it.
  !!}

  !![
  <nodeOperator name="nodeOperatorHaloFilaments" docformat="rst">
   <description>
   A node operator which assigns each halo the axes of the cosmic web filaments along which it accretes, for use by
   :galacticus-class:`virialOrbitFilamentary`. The axes are stored in a meta-property, ``filamentAxes``, of the ``basic``
   component.

   Each branch of the merger tree is assigned a single filament passing through the halo---that is, a pair of opposite
   axes---oriented isotropically at random. Filaments persist: every node on a branch shares the axes of the branch, so that a
   halo accretes along the same filament throughout its history. A single filament is the limit of the number of filaments
   crossing a sphere of small radius about a density peak in the peak theory calculations of :cite:t:`codis_connectivity_2018`, while
   the number crossing a sphere of :math:`1.5 R_\mathrm{200}` in simulations is typically two for groups, rising to three to
   four for clusters :cite:p:`gouin_shape_2021`. A single filament is used for now because the mass dependence of the anisotropy
   of accretion assumed by :galacticus-class:`virialOrbitFilamentary` is calibrated for it: with the number of filaments
   increasing with mass as measured by :cite:t:`gouin_shape_2021`, the fraction of satellites accreted along filaments required
   to match the measured anisotropy is no longer a simple function of mass.

   Axes are assigned to all nodes when the tree is initialized. They are also assigned on demand, by
   :galacticus-class:`virialOrbitFilamentary`, to any branch which does not yet have them when an orbit is required---for
   example, for orbits assigned while a tree is being constructed---so that the axes do not depend on the order in which these
   events occur.
   </description>
  </nodeOperator>
  !!]
  type, extends(nodeOperatorClass) :: nodeOperatorHaloFilaments
     !!{RST
     A node operator which assigns each halo the axes of the cosmic web filaments which feed it.
     !!}
     private
     integer :: filamentAxesID
   contains
     procedure :: nodeTreeInitialize => haloFilamentsNodeTreeInitialize
  end type nodeOperatorHaloFilaments

  interface nodeOperatorHaloFilaments
     !!{RST
     Constructors for the :galacticus-class:`nodeOperatorHaloFilaments` node operator class.
     !!}
     module procedure haloFilamentsConstructorParameters
     module procedure haloFilamentsConstructorInternal
  end interface nodeOperatorHaloFilaments

contains

  function haloFilamentsConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`nodeOperatorHaloFilaments` node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(nodeOperatorHaloFilaments)                :: self
    type(inputParameters          ), intent(inout) :: parameters

    self=nodeOperatorHaloFilaments()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function haloFilamentsConstructorParameters

  function haloFilamentsConstructorInternal() result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`nodeOperatorHaloFilaments` node operator class.
    !!}
    implicit none
    type(nodeOperatorHaloFilaments) :: self

    !![
    <addMetaProperty component="basic" name="filamentAxes" id="self%filamentAxesID" rank="1" isEvolvable="no" isCreator="yes"/>
    !!]
    return
  end function haloFilamentsConstructorInternal

  subroutine haloFilamentsNodeTreeInitialize(self,node)
    !!{RST
    Assign filament axes to the branch of a node, if it does not already have them.
    !!}
    use :: Dark_Matter_Halo_Filaments, only : Dark_Matter_Halo_Filament_Axes
    implicit none
    class           (nodeOperatorHaloFilaments), intent(inout), target         :: self
    type            (treeNode                 ), intent(inout), target         :: node
    double precision                           , allocatable  , dimension(:,:) :: axes

    axes=Dark_Matter_Halo_Filament_Axes(node,self%filamentAxesID)
    return
  end subroutine haloFilamentsNodeTreeInitialize
