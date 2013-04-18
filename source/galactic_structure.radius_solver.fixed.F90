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

!% Contains a module which implements a ``fixed'' galactic radii solver in which sizes are always equal to
!% the halo virial radius multiplied by its spin parameter and a multiplicative constant.

module Galactic_Structure_Radii_Fixed
  !% Implements a ``fixed'' galactic radii solver in which sizes are always equal to           
  !% the halo virial radius multiplied by its spin parameter and a multiplicative constant.
  use Galactic_Structure_Radius_Solver_Procedures
  implicit none
  private
  public :: Galactic_Structure_Radii_Fixed_Initialize

  ! The ratio of galaxy size to the product of spin parameter and virial radius.
  double precision :: galacticStructureRadiiFixedFactor

contains

  !# <galacticStructureRadiusSolverMethod>
  !#  <unitName>Galactic_Structure_Radii_Fixed_Initialize</unitName>
  !# </galacticStructureRadiusSolverMethod>
  subroutine Galactic_Structure_Radii_Fixed_Initialize(galacticStructureRadiusSolverMethod,Galactic_Structure_Radii_Solve_Do)
    !% Initializes the ``fixed'' galactic radii solver module.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string),          intent(in)    :: galacticStructureRadiusSolverMethod
    procedure(Galactic_Structure_Radii_Solve_Fixed),          pointer, intent(inout) :: Galactic_Structure_Radii_Solve_Do
    
    if (galacticStructureRadiusSolverMethod == 'fixed') then
       Galactic_Structure_Radii_Solve_Do => Galactic_Structure_Radii_Solve_Fixed
       !@ <inputParameter>
       !@   <name>galacticStructureRadiiFixedFactor</name>
       !@   <defaultValue>$\sqrt{1/2}$ \citep{mo_formation_1998}</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The ratio of galaxy radius to $\lambda r_{\rm vir}$ in the ``fixed'' galactic structure radius solver algorithm.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('galacticStructureRadiiFixedFactor',galacticStructureRadiiFixedFactor,defaultValue=sqrt(0.5d0))
    end if
    return
  end subroutine Galactic_Structure_Radii_Fixed_Initialize

  subroutine Galactic_Structure_Radii_Solve_Fixed(thisNode)
    !% Find the radii of galactic components in {\tt thisNode} using the ``fixed'' method.
    use Galacticus_Error
    use Galacticus_Nodes
    include 'galactic_structure.radius_solver.tasks.modules.inc'
    include 'galactic_structure.radius_solver.plausible.modules.inc'
    implicit none
    type(treeNode),                    intent(inout), pointer :: thisNode
    procedure(Structure_Get_Template),                pointer :: Radius_Get => null(), Velocity_Get => null()
    procedure(Structure_Set_Template),                pointer :: Radius_Set => null(), Velocity_Set => null()
    logical                                                   :: componentActive
    double precision                                          :: specificAngularMomentum

    ! Check that the galaxy is physical plausible. In this fixed solver, we don't act on this.
    thisNode%isPhysicallyPlausible=.true.
    include 'galactic_structure.radius_solver.plausible.inc'

    ! Solve for radii.
    include 'galactic_structure.radius_solver.tasks.inc'

    return
  end subroutine Galactic_Structure_Radii_Solve_Fixed
  
  subroutine Solve_For_Radius(thisNode,specificAngularMomentum,Radius_Get,Radius_Set,Velocity_Get,Velocity_Set)
    !% Solve for the equilibrium radius of the given component.
    use Dark_Matter_Halo_Scales
    implicit none
    type(treeNode),                    pointer, intent(inout) :: thisNode
    double precision,                           intent(in)    :: specificAngularMomentum
    procedure(Structure_Get_Template), pointer, intent(in)    :: Radius_Get, Velocity_Get
    procedure(Structure_Set_Template), pointer, intent(in)    :: Radius_Set, Velocity_Set
    class    (nodeComponentSpin     ), pointer                :: thisSpinComponent
    double precision                                          :: radius,velocity

    ! Return immediately if the specific angular momentum is zero.
    if (specificAngularMomentum <= 0.0d0) return
    
    ! Find the radius of the component, assuming radius scales fixedly with angular momentum.
    thisSpinComponent => thisNode%spin()
    velocity=Dark_Matter_Halo_Virial_Velocity(thisNode)
    radius  =Dark_Matter_Halo_Virial_Radius  (thisNode)*thisSpinComponent%spin()*galacticStructureRadiiFixedFactor

    ! Set the component size to new radius and velocity.
    call Radius_Set  (thisNode,radius  )
    call Velocity_Set(thisNode,velocity)

    return
  end subroutine Solve_For_Radius

end module Galactic_Structure_Radii_Fixed
