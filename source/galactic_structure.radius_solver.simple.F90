!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@caltech.edu>
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
!!
!!
!!    COPYRIGHT 2010. The Jet Propulsion Laboratory/California Institute of Technology
!!
!!    The California Institute of Technology shall allow RECIPIENT to use and
!!    distribute this software subject to the terms of the included license
!!    agreement with the understanding that:
!!
!!    THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE CALIFORNIA
!!    INSTITUTE OF TECHNOLOGY (CALTECH). THE SOFTWARE IS PROVIDED "AS-IS" TO
!!    THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY WARRANTIES OF
!!    PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A PARTICULAR USE OR
!!    PURPOSE (AS SET FORTH IN UNITED STATES UCC §2312-§2313) OR FOR ANY
!!    PURPOSE WHATSOEVER, FOR THE SOFTWARE AND RELATED MATERIALS, HOWEVER
!!    USED.
!!
!!    IN NO EVENT SHALL CALTECH BE LIABLE FOR ANY DAMAGES AND/OR COSTS,
!!    INCLUDING, BUT NOT LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF
!!    ANY KIND, INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST
!!    PROFITS, REGARDLESS OF WHETHER CALTECH BE ADVISED, HAVE REASON TO KNOW,
!!    OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.
!!
!!    RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF THE
!!    SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY CALTECH FOR
!!    ALL THIRD-PARTY CLAIMS RESULTING FROM THE ACTIONS OF RECIPIENT IN THE
!!    USE OF THE SOFTWARE.
!!
!!    In addition, RECIPIENT also agrees that Caltech is under no obligation
!!    to provide technical support for the Software.
!!
!!    Finally, Caltech places no restrictions on RECIPIENT's use, preparation
!!    of Derivative Works, public display or redistribution of the Software
!!    other than those specified in the included license and the requirement
!!    that all copies of the Software released be marked with the language
!!    provided in this notice.
!!
!!    This software is separately available under negotiable license terms
!!    from:
!!    California Institute of Technology
!!    Office of Technology Transfer
!!    1200 E. California Blvd.
!!    Pasadena, California 91125
!!    http://www.ott.caltech.edu


!% Contains a module which implements a simple galactic radii solver (no adiabatic contraction and no self-gravity of baryons).

module Galactic_Structure_Radii_Simple
  !% Implements a simple galactic radii solver (no adiabatic contraction and no self-gravity of baryons).
  use Tree_Nodes
  use Galactic_Structure_Radius_Solver_Procedures
  implicit none
  private
  public :: Galactic_Structure_Radii_Simple_Initialize

  ! Module variables used to communicate current state of radius solver.
  type(treeNode), pointer :: haloNode
  !$omp threadprivate(haloNode)

  ! Options controlling the solver.
  logical                 :: simpleRadiusSolverUseFormationHalo

contains

  !# <galacticStructureRadiusSolverMethod>
  !#  <unitName>Galactic_Structure_Radii_Simple_Initialize</unitName>
  !# </galacticStructureRadiusSolverMethod>
  subroutine Galactic_Structure_Radii_Simple_Initialize(galacticStructureRadiusSolverMethod,Galactic_Structure_Radii_Solve_Do)
    !% Initializes the ``simple'' galactic radii solver module.
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string),          intent(in)    :: galacticStructureRadiusSolverMethod
    procedure(),          pointer, intent(inout) :: Galactic_Structure_Radii_Solve_Do
    
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
    use Tree_Nodes
    !# <include directive="radiusSolverTask" type="moduleUse">
    include 'galactic_structure.radius_solver.tasks.modules.inc'
    !# </include>
    !# <include directive="radiusSolverPlausibility" type="moduleUse">
    include 'galactic_structure.radius_solver.plausible.modules.inc'
    !# </include>
    implicit none
    type(treeNode),                    intent(inout), pointer :: thisNode
    procedure(Structure_Get_Template), save,          pointer :: Radius_Get => null(), Velocity_Get => null()
    procedure(Structure_Set_Template), save,          pointer :: Radius_Set => null(), Velocity_Set => null()
    !$omp threadprivate(Radius_Get,Radius_Set,Velocity_Get,Velocity_Set)
    logical                                                   :: componentActive
    double precision                                          :: specificAngularMomentum

    ! Check that the galaxy is physical plausible. In this simple solver, we don't act on this.
    thisNode%isPhysicallyPlausible=.true.
    !# <include directive="radiusSolverPlausibility" type="code" action="subroutine">
    !#  <subroutineArgs>thisNode,thisNode%isPhysicallyPlausible</subroutineArgs>
    include 'galactic_structure.radius_solver.plausible.inc'
    !# </include>

    ! Determine which node to use for halo properties.
    if (simpleRadiusSolverUseFormationHalo) then
       if (.not.associated(thisNode%formationNode)) call Galacticus_Error_Report('Galactic_Structure_Radii_Solve_Simple','no formation node exists')
       haloNode => thisNode%formationNode
    else
       haloNode => thisNode
    end if

    !# <include directive="radiusSolverTask" type="code" action="subroutine">
    !#  <subroutineArgs>thisNode,componentActive,specificAngularMomentum,Radius_Get,Radius_Set,Velocity_Get,Velocity_Set</subroutineArgs>
    !#  <subroutineAction>if (componentActive) call Solve_For_Radius(thisNode,specificAngularMomentum,Radius_Get,Radius_Set,Velocity_Get,Velocity_Set)</subroutineAction>
    include 'galactic_structure.radius_solver.tasks.inc'
    !# </include>

    return
  end subroutine Galactic_Structure_Radii_Solve_Simple
  
  subroutine Solve_For_Radius(thisNode,specificAngularMomentum,Radius_Get,Radius_Set,Velocity_Get,Velocity_Set)
    !% Solve for the equilibrium radius of the given component.
    use Dark_Matter_Profiles
    implicit none
    type(treeNode),                    pointer, intent(inout) :: thisNode
    double precision,                           intent(in)    :: specificAngularMomentum
    procedure(Structure_Get_Template), pointer, intent(in)    :: Radius_Get, Velocity_Get
    procedure(Structure_Set_Template), pointer, intent(in)    :: Radius_Set, Velocity_Set
    double precision                                          :: radius,velocity
    
    ! Find the radius in the dark matter profile with the required specific angular momentum
    radius=Dark_Matter_Profile_Radius_from_Specific_Angular_Momentum(haloNode,specificAngularMomentum)

    ! Find the velocity at this radius.
    velocity=Dark_Matter_Profile_Circular_Velocity(haloNode,radius)

    ! Set the component size to new radius and velocity.
    call Radius_Set  (thisNode,radius  )
    call Velocity_Set(thisNode,velocity)
    return
  end subroutine Solve_For_Radius

end module Galactic_Structure_Radii_Simple
