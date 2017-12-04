!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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

!% Contains a module which implements calculations of sizes of galactic components (or more general components).

module Galactic_Structure_Radii
  !% Implements calculations of sizes of galactic components (or more general components).
  use ISO_Varying_String
  use Galacticus_Nodes
  implicit none
  private
  public :: Galactic_Structure_Radii_Solve, Galactic_Structure_Radii_Revert

  ! Flag to indicate if this module has been initialized.
  logical                                             :: galacticStructureRadiusSolverInitialized=.false.

  ! Name of cooling rate available method used.
  type     (varying_string                 )          :: galacticStructureRadiusSolverMethod

  ! Pointer to the function that actually does the calculation.
  procedure(Galactic_Structure_Radii_Solve ), pointer :: Galactic_Structure_Radii_Solve_Do  => null()
  procedure(Galactic_Structure_Radii_Revert), pointer :: Galactic_Structure_Radii_Revert_Do => null()

contains

  subroutine Galacticus_Structure_Radii_Initialize()
    !% Initialize galactic structure radius solver.
    use Galacticus_Error
    use Input_Parameters
    !# <include directive="galacticStructureRadiusSolverMethod" type="moduleUse">
    include 'galactic_structure.radius_solver.modules.inc'
    !# </include>
    implicit none

    if (.not.galacticStructureRadiusSolverInitialized) then
       !$omp critical(Galactic_Structure_Radii_Initialization)
       if (.not.galacticStructureRadiusSolverInitialized) then
          ! Get the galactic structure radii solver method parameter.
          !# <inputParameter>
          !#   <name>galacticStructureRadiusSolverMethod</name>
          !#   <cardinality>1</cardinality>
          !#   <defaultValue>var_str('adiabatic')</defaultValue>
          !#   <description>Selects the method to be used for solving for galactic structure.</description>
          !#   <source>globalParameters</source>
          !#   <type>string</type>
          !# </inputParameter>
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="galacticStructureRadiusSolverMethod" type="functionCall" functionType="void">
          !#  <functionArgs>galacticStructureRadiusSolverMethod,Galactic_Structure_Radii_Solve_Do,Galactic_Structure_Radii_Revert_Do</functionArgs>
          include 'galactic_structure.radius_solver.inc'
          !# </include>
          if (.not.(associated(Galactic_Structure_Radii_Solve_Do).and.associated(Galactic_Structure_Radii_Revert_Do))) &
               & call Galacticus_Error_Report('method '//char(galacticStructureRadiusSolverMethod)//' is unrecognized'//{introspection:location})
          galacticStructureRadiusSolverInitialized=.true.
       end if
       !$omp end critical(Galactic_Structure_Radii_Initialization)
    end if
    return
  end subroutine Galacticus_Structure_Radii_Initialize
  
  !# <preDerivativeTask>
  !#  <unitName>Galactic_Structure_Radii_Solve</unitName>
  !# </preDerivativeTask>
  subroutine Galactic_Structure_Radii_Solve(node)
    !% Solve for the radii of galactic components in {\normalfont \ttfamily node}.
    implicit none
    type(treeNode), intent(inout), target :: node

    call Galacticus_Structure_Radii_Initialize(    )
    call Galactic_Structure_Radii_Solve_Do    (node)
    return
  end subroutine Galactic_Structure_Radii_Solve
  
  subroutine Galactic_Structure_Radii_Revert(node)
    !% Revert the radii of galactic components in {\normalfont \ttfamily node} (if necessary to ensure that the structure solver
    !% will give the same result when called consequtively).
    implicit none
    type(treeNode), intent(inout), target :: node

    call Galacticus_Structure_Radii_Initialize(    )
    call Galactic_Structure_Radii_Revert_Do   (node)
    return
  end subroutine Galactic_Structure_Radii_Revert

end module Galactic_Structure_Radii
