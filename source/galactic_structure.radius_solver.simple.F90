!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements a simple galactic radii solver (no adiabatic contraction and no self-gravity of baryons).

module Galactic_Structure_Radii_Simple
  !% Implements a simple galactic radii solver (no adiabatic contraction and no self-gravity of baryons).
  use Galactic_Structure_Radius_Solver_Procedures
  implicit none
  private
  public :: Galactic_Structure_Radii_Simple_Initialize

  ! Module variables used to communicate current state of radius solver.
  type   (treeNode), pointer :: haloNode
  !$omp threadprivate(haloNode)
  ! Options controlling the solver.
  logical                    :: simpleRadiusSolverUseFormationHalo

contains

  !# <galacticStructureRadiusSolverMethod>
  !#  <unitName>Galactic_Structure_Radii_Simple_Initialize</unitName>
  !# </galacticStructureRadiusSolverMethod>
  subroutine Galactic_Structure_Radii_Simple_Initialize(galacticStructureRadiusSolverMethod,Galactic_Structure_Radii_Solve_Do)
    !% Initializes the ``simple'' galactic radii solver module.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type     (varying_string                       ), intent(in   )          :: galacticStructureRadiusSolverMethod
    procedure(Galactic_Structure_Radii_Solve_Simple), intent(inout), pointer :: Galactic_Structure_Radii_Solve_Do

    if (galacticStructureRadiusSolverMethod == 'simple') then
       Galactic_Structure_Radii_Solve_Do => Galactic_Structure_Radii_Solve_Simple
       !@ <inputParameter>
       !@   <name>simpleRadiusSolverUseFormationHalo</name>
       !@   <defaultValue>false</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Specifies whether or not the ``formation halo'' should be used when solving for the radii of galaxies.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('simpleRadiusSolverUseFormationHalo',simpleRadiusSolverUseFormationHalo,defaultValue=.false.)
    end if
    return
  end subroutine Galactic_Structure_Radii_Simple_Initialize

  subroutine Galactic_Structure_Radii_Solve_Simple(thisNode)
    !% Find the radii of galactic components in {\tt thisNode} using the ``simple'' method.
    use Galacticus_Error
    !# <include directive="radiusSolverTask" type="moduleUse">
    include 'galactic_structure.radius_solver.tasks.modules.inc'
    !# </include>
    !# <include directive="radiusSolverPlausibility" type="moduleUse">
    include 'galactic_structure.radius_solver.plausible.modules.inc'
    !# </include>
    implicit none
    type            (treeNode              ), intent(inout), pointer       :: thisNode
    logical                                 , parameter                    :: specificAngularMomentumRequired=.true.
    procedure       (Structure_Get_Template)               , pointer, save :: Radius_Get             =>null(), Velocity_Get=>null()
    procedure       (Structure_Set_Template)               , pointer, save :: Radius_Set             =>null(), Velocity_Set=>null()
    !$omp threadprivate(Radius_Get,Radius_Set,Velocity_Get,Velocity_Set)
    logical                                                                :: componentActive
    double precision                                                       :: specificAngularMomentum

    ! Check that the galaxy is physical plausible. In this simple solver, we don't act on this.
    thisNode%isPhysicallyPlausible=.true.
    !# <include directive="radiusSolverPlausibility" type="functionCall" functionType="void">
    !#  <functionArgs>thisNode,thisNode%isPhysicallyPlausible</functionArgs>
    include 'galactic_structure.radius_solver.plausible.inc'
    !# </include>

    ! Determine which node to use for halo properties.
    if (simpleRadiusSolverUseFormationHalo) then
       if (.not.associated(thisNode%formationNode)) call Galacticus_Error_Report('Galactic_Structure_Radii_Solve_Simple','no formation node exists')
       haloNode => thisNode%formationNode
    else
       haloNode => thisNode
    end if

    !# <include directive="radiusSolverTask" type="functionCall" functionType="void">
    !#  <functionArgs>thisNode,componentActive,specificAngularMomentumRequired,specificAngularMomentum,Radius_Get,Radius_Set,Velocity_Get,Velocity_Set</functionArgs>
    !#  <onReturn>if (componentActive) call Solve_For_Radius(thisNode,specificAngularMomentum,Radius_Get,Radius_Set,Velocity_Get,Velocity_Set)</onReturn>
    include 'galactic_structure.radius_solver.tasks.inc'
    !# </include>

    return
  end subroutine Galactic_Structure_Radii_Solve_Simple

  subroutine Solve_For_Radius(thisNode,specificAngularMomentum,Radius_Get,Radius_Set,Velocity_Get,Velocity_Set)
    !% Solve for the equilibrium radius of the given component.
    use Dark_Matter_Profiles
    implicit none
    type            (treeNode              ), intent(inout), pointer :: thisNode
    double precision                        , intent(in   )          :: specificAngularMomentum
    procedure       (Structure_Get_Template), intent(in   ), pointer :: Radius_Get             , Velocity_Get
    procedure       (Structure_Set_Template), intent(in   ), pointer :: Radius_Set             , Velocity_Set
    class           (darkMatterProfileClass)               , pointer :: darkMatterProfile_
    double precision                                                 :: radius                 , velocity

    ! Return immediately if the specific angular momentum is zero.
    if (specificAngularMomentum <= 0.0d0) return

    ! Get required objects.
    darkMatterProfile_ => darkMatterProfile()
    ! Find the radius in the dark matter profile with the required specific angular momentum
    radius=darkMatterProfile_%radiusFromSpecificAngularMomentum(haloNode,specificAngularMomentum)

    ! Find the velocity at this radius.
    velocity=darkMatterProfile_%circularVelocity(haloNode,radius)

    ! Set the component size to new radius and velocity.
    call Radius_Set  (thisNode,radius  )
    call Velocity_Set(thisNode,velocity)
    return
  end subroutine Solve_For_Radius

end module Galactic_Structure_Radii_Simple
