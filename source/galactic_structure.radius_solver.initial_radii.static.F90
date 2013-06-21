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

!% Contains a module which implements calculations of initial radius in which the dark matter halo is assumed to be static.

module Galactic_Structure_Initial_Radii_Static
  !% Implements calculations of initial radius in which the dark matter halois assumed to be static
  use Galacticus_Nodes
  implicit none
  private
  public :: Galactic_Structure_Initial_Radii_Static_Initialize

contains

  !# <galacticStructureRadiusSolverInitialRadiusMethod>
  !#  <unitName>Galactic_Structure_Initial_Radii_Static_Initialize</unitName>
  !# </galacticStructureRadiusSolverInitialRadiusMethod>
  subroutine Galactic_Structure_Initial_Radii_Static_Initialize(galacticStructureRadiusSolverInitialRadiusMethod&
       &,Galactic_Structure_Radius_Initial_Get,Galactic_Structure_Radius_Initial_Derivative_Get)
    !% Initializes the ``static'' initial radii module.
    use Input_Parameters
    use ISO_Varying_String
    implicit none
    type     (varying_string                                     ), intent(in   )          :: galacticStructureRadiusSolverInitialRadiusMethod
    procedure(Galactic_Structure_Radius_Initial_Static           ), intent(inout), pointer :: Galactic_Structure_Radius_Initial_Get
    procedure(Galactic_Structure_Radius_Initial_Derivative_Static), intent(inout), pointer :: Galactic_Structure_Radius_Initial_Derivative_Get

    if (galacticStructureRadiusSolverInitialRadiusMethod == 'static') then
       Galactic_Structure_Radius_Initial_Get            => Galactic_Structure_Radius_Initial_Static
       Galactic_Structure_Radius_Initial_Derivative_Get => Galactic_Structure_Radius_Initial_Derivative_Static
    end if
    return
  end subroutine Galactic_Structure_Initial_Radii_Static_Initialize

  double precision function Galactic_Structure_Radius_Initial_Static(thisNode,radius)
    !% Compute the initial radius in the dark matter halo assuming the halo is static.
    implicit none
    type            (treeNode), intent(inout), pointer :: thisNode
    double precision          , intent(in   )          :: radius

    Galactic_Structure_Radius_Initial_Static=radius
    return
  end function Galactic_Structure_Radius_Initial_Static

  double precision function Galactic_Structure_Radius_Initial_Derivative_Static(thisNode,radius)
    !% Compute the derivative of the initial radius in the dark matter halo assuming the halo is static.
    implicit none
    type            (treeNode), intent(inout), pointer :: thisNode
    double precision          , intent(in   )          :: radius

    Galactic_Structure_Radius_Initial_Derivative_Static=1.0d0
    return
  end function Galactic_Structure_Radius_Initial_Derivative_Static

end module Galactic_Structure_Initial_Radii_Static
