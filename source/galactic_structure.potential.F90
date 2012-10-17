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

!+    Contributions to this file made by:  Stéphane Mangeon, Andrew Benson.

!% Contains a module which implements calculations of gravitationl potential.

module Galactic_Structure_Potentials
  !% Implements calculations of the gravitational potential.
  use, intrinsic :: ISO_C_Binding
  use ISO_Varying_String
  use Tree_Nodes
  use Galactic_Structure_Options
  implicit none
  private
  public :: Galactic_Structure_Potential

contains

  double precision function Galactic_Structure_Potential(thisNode,radius,componentType,massType)
    !% Solve for the gravitational potential at a given radius. Assumes the galactic structure has already been computed.
    use Galacticus_Error
    use Input_Parameters
    use Dark_Matter_Profiles
    use Numerical_Constants_Physical
    use Numerical_Constants_Prefixes
    !# <include directive="potentialTask" type="moduleUse">
    include 'galactic_structure.potential.tasks.modules.inc'
    !# </include>
    implicit none
    type(treeNode),   intent(inout), pointer  :: thisNode
    integer,          intent(in),    optional :: componentType,massType
    double precision, intent(in)              :: radius
    integer                                   :: componentTypeActual,massTypeActual
    double precision                          :: componentPotential

    ! Determine which component type to use.
    if (present(componentType)) then
       componentTypeActual=componentType
    else
       componentTypeActual=componentTypeAll
    end if
    ! Determine which mass type to use.
    if (present(massType)) then
       massTypeActual     =massType
    else
       massTypeActual     =massTypeAll
    end if

    ! Initialize to zero potential.
    Galactic_Structure_Potential=0.0d0
    
    ! Call routines to supply the potential for all components.
    !# <include directive="potentialTask" type="code" action="subroutine">
    !#  <subroutineArgs>thisNode,radius,componentTypeActual,massTypeActual,componentPotential</subroutineArgs>
    !#  <subroutineAction>Galactic_Structure_Potential=Galactic_Structure_Potential+componentPotential</subroutineAction>
    include 'galactic_structure.potential.tasks.inc'
    !# </include>

    return
  end function Galactic_Structure_Potential
  
end module Galactic_Structure_Potentials
    
    
    
    








