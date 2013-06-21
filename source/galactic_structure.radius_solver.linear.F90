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

!% Contains a module which implements a ``linear'' galactic radii solver (no adiabatic contraction
!% and no self-gravity of baryons, and size simply scales in proportion to specific angular momentum).

module Galactic_Structure_Radii_Linear
  !% Implements a ``linear'' galactic radii solver (no adiabatic contraction and no self-gravity
  !% of baryons, and size simply scales in proportion to specific angular momentum).
  use Galacticus_Nodes
  use Galactic_Structure_Radius_Solver_Procedures
  implicit none
  private
  public :: Galactic_Structure_Radii_Linear_Initialize

  ! Module variables used to communicate current state of radius solver.
  type(treeNode), pointer :: haloNode 
  !$omp threadprivate(haloNode)
contains

  !# <galacticStructureRadiusSolverMethod>
  !#  <unitName>Galactic_Structure_Radii_Linear_Initialize</unitName>
  !# </galacticStructureRadiusSolverMethod>
  subroutine Galactic_Structure_Radii_Linear_Initialize(galacticStructureRadiusSolverMethod,Galactic_Structure_Radii_Solve_Do)
    !% Initializes the ``linear'' galactic radii solver module.
    use ISO_Varying_String
    implicit none
    type     (varying_string                       ), intent(in   )          :: galacticStructureRadiusSolverMethod 
    procedure(Galactic_Structure_Radii_Solve_Linear), intent(inout), pointer :: Galactic_Structure_Radii_Solve_Do   
    
    if (galacticStructureRadiusSolverMethod == 'linear') Galactic_Structure_Radii_Solve_Do => Galactic_Structure_Radii_Solve_Linear
    return
  end subroutine Galactic_Structure_Radii_Linear_Initialize

  subroutine Galactic_Structure_Radii_Solve_Linear(thisNode)
    !% Find the radii of galactic components in {\tt thisNode} using the ``linear'' method.
    use Galacticus_Error
    include 'galactic_structure.radius_solver.tasks.modules.inc'
    include 'galactic_structure.radius_solver.plausible.modules.inc'
    implicit none
    type            (treeNode              ), intent(inout), pointer :: thisNode                                              
    procedure       (Structure_Get_Template)               , pointer :: Radius_Get             =>null(), Velocity_Get=>null() 
    procedure       (Structure_Set_Template)               , pointer :: Radius_Set             =>null(), Velocity_Set=>null() 
    logical                                                          :: componentActive                                       
    double precision                                                 :: specificAngularMomentum                               
    
    ! Check that the galaxy is physical plausible. In this linear solver, we don't act on this.
    thisNode%isPhysicallyPlausible=.true.
    include 'galactic_structure.radius_solver.plausible.inc'
    include 'galactic_structure.radius_solver.tasks.inc'

    return
  end subroutine Galactic_Structure_Radii_Solve_Linear
  
  subroutine Solve_For_Radius(thisNode,specificAngularMomentum,Radius_Get,Radius_Set,Velocity_Get,Velocity_Set)
    !% Solve for the equilibrium radius of the given component.
    use Dark_Matter_Halo_Scales
    implicit none
    type            (treeNode              ), intent(inout), pointer :: thisNode                              
    double precision                        , intent(in   )          :: specificAngularMomentum               
    procedure       (Structure_Get_Template), intent(in   ), pointer :: Radius_Get             , Velocity_Get 
    procedure       (Structure_Set_Template), intent(in   ), pointer :: Radius_Set             , Velocity_Set 
    double precision                                                 :: radius                 , velocity     
    
    ! Return immediately if the specific angular momentum is zero.
    if (specificAngularMomentum <= 0.0d0) return
    
    ! Find the radius of the component, assuming radius scales linearly with angular momentum.
    velocity=Dark_Matter_Halo_Virial_Velocity(thisNode)
    radius  =specificAngularMomentum/velocity

    ! Set the component size to new radius and velocity.
    call Radius_Set  (thisNode,radius  )
    call Velocity_Set(thisNode,velocity)
    return
  end subroutine Solve_For_Radius

end module Galactic_Structure_Radii_Linear
