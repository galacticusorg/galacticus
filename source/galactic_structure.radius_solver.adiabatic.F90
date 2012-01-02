!! Copyright 2009, 2010, 2011 Andrew Benson <abenson@caltech.edu>
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


!% Contains a module which implements an adiabatic contraction galactic radii solver (including self-gravity of baryons) using the
!% adiabatic contraction algorithm of \cite{gnedin_response_2004}.

module Galactic_Structure_Radii_Adiabatic
  !% Implements an adiabatic contraction galactic radii solver (including self-gravity of baryons) using the adiabatic
  !% contraction algorithm of \cite{gnedin_response_2004}.
  use Tree_Nodes
  use Galactic_Structure_Radius_Solver_Procedures
  implicit none
  private
  public :: Galactic_Structure_Radii_Adiabatic_Initialize

  ! Parameters of the adiabatic contraction algorithm.
  double precision        :: adiabaticContractionGnedinA,adiabaticContractionGnedinOmega&
       &,inverseAdiabaticContractionGnedinOmegaMinusOne

  ! Parameter controlling the accuracy of the solutions sought.
  double precision        :: adiabaticContractionSolutionTolerance

  ! Module variables used to communicate current state of radius solver.
  integer                 :: iterationCount,activeComponentCount
  double precision        :: fitMeasure,haloFraction
  type(treeNode), pointer :: haloNode
  !$omp threadprivate(iterationCount,activeComponentCount,fitMeasure,haloFraction,haloNode)

  ! Options controlling the solver.
  logical                 :: adiabaticContractionIncludeBaryonGravity,adiabaticContractionUseFormationHalo

contains

  !# <galacticStructureRadiusSolverMethod>
  !#  <unitName>Galactic_Structure_Radii_Adiabatic_Initialize</unitName>
  !# </galacticStructureRadiusSolverMethod>
  subroutine Galactic_Structure_Radii_Adiabatic_Initialize(galacticStructureRadiusSolverMethod,Galactic_Structure_Radii_Solve_Do)
    !% Initializes the ``adiabatic'' galactic radii solver module.
    use Input_Parameters
    use ISO_Varying_String
    implicit none
    type(varying_string),          intent(in)    :: galacticStructureRadiusSolverMethod
    procedure(),          pointer, intent(inout) :: Galactic_Structure_Radii_Solve_Do
    
    if (galacticStructureRadiusSolverMethod == 'adiabatic') then
       Galactic_Structure_Radii_Solve_Do => Galactic_Structure_Radii_Solve_Adiabatic
       ! Get parameters of the model.
       !@ <inputParameter>
       !@   <name>adiabaticContractionGnedinA</name>
       !@   <defaultValue>0.8 (\citealt{gustafsson_baryonic_2006}; from their Fig. 9, strong feedback case)</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The parameter $A$ appearing in the \cite{gnedin_response_2004} adiabatic contraction algorithm.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('adiabaticContractionGnedinA'    ,adiabaticContractionGnedinA    ,defaultValue=0.80d0)
       !@ <inputParameter>
       !@   <name>adiabaticContractionGnedinOmega</name>
       !@   <defaultValue>0.77 (\citealt{gustafsson_baryonic_2006}; from their Fig. 9, strong feedback case)</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The parameter $\omega$ appearing in the \cite{gnedin_response_2004} adiabatic contraction algorithm.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('adiabaticContractionGnedinOmega',adiabaticContractionGnedinOmega,defaultValue=0.77d0)
       !@ <inputParameter>
       !@   <name>adiabaticContractionIncludeBaryonGravity</name>
       !@   <defaultValue>true</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Specifies whether or not gravity from baryons is included when solving for sizes of galactic components in adiabatically contracted dark matter halos.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('adiabaticContractionIncludeBaryonGravity',adiabaticContractionIncludeBaryonGravity,defaultValue=.true.)
       !@ <inputParameter>
       !@   <name>adiabaticContractionUseFormationHalo</name>
       !@   <defaultValue>false</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Specifies whether or not the ``formation halo'' should be used when solving for the radii of galaxies.
       !@   </description>
       !@   <type>boolean</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('adiabaticContractionUseFormationHalo',adiabaticContractionUseFormationHalo,defaultValue=.false.)
       !@ <inputParameter>
       !@   <name>adiabaticContractionSolutionTolerance</name>
       !@   <defaultValue></defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     Maximum allowed mean fractional error in the radii of all components when seeking equilibrium solutions for galactic structure.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('adiabaticContractionSolutionTolerance',adiabaticContractionSolutionTolerance,defaultValue=1.0d-2)
       ! Store the exponent that we actually use in the equations.
       inverseAdiabaticContractionGnedinOmegaMinusOne=1.0d0/adiabaticContractionGnedinOmega-1.0d0
    end if
    return
  end subroutine Galactic_Structure_Radii_Adiabatic_Initialize

  subroutine Galactic_Structure_Radii_Solve_Adiabatic(thisNode)
    !% Find the radii of galactic components in {\tt thisNode} using the ``adiabatic'' method.
    use Tree_Nodes
    use Cosmological_Parameters
    use Galacticus_Error
    use Galactic_Structure_Enclosed_Masses
    use Galactic_Structure_Options
    include 'galactic_structure.radius_solver.tasks.modules.inc'
    implicit none
    type(treeNode),                    intent(inout), pointer :: thisNode
    integer,                           parameter              :: iterationMaximum=100
    procedure(Structure_Get_Template),                pointer :: Radius_Get => null(), Velocity_Get => null()
    procedure(Structure_Set_Template),                pointer :: Radius_Set => null(), Velocity_Set => null()
    logical                                                   :: componentActive,galaxyIsPhysicallyPlausible
    double precision                                          :: specificAngularMomentum

    ! Check that the galaxy is physical plausible. If not, do not try to solve for its structure.
    galaxyIsPhysicallyPlausible=.true.
    include 'galactic_structure.radius_solver.plausible.inc'
    if (galaxyIsPhysicallyPlausible) then
       ! Initialize the solver state.
       iterationCount=0
       fitMeasure    =2.0d0*adiabaticContractionSolutionTolerance

       ! Determine which node to use for halo properties.
       if (adiabaticContractionUseFormationHalo) then
          if (.not.associated(thisNode%formationNode)) call Galacticus_Error_Report('Galactic_Structure_Radii_Solve_Adiabatic','no formation node exists')
          haloNode => thisNode%formationNode
       else
          haloNode => thisNode
       end if

       ! Compute fraction of mass distribution as the halo. Truncate this to zero: we can get negative values if the ODE solver is
       ! exploring regimes of high baryonic mass, and this would cause problems.
       haloFraction=max(1.0d0-Galactic_Structure_Enclosed_Mass(thisNode,massType=massTypeGalactic)/Tree_Node_Mass(haloNode),0.0d0)
       
       ! Begin iteration to find a converged solution.
       do while (iterationCount <= 2 .or. ( fitMeasure > adiabaticContractionSolutionTolerance .and. iterationCount < iterationMaximum ) )
          iterationCount      =iterationCount+1
          activeComponentCount=0
          if (iterationCount > 1) fitMeasure=0.0d0
          include 'galactic_structure.radius_solver.tasks.inc'
          ! Check that we have some active components.
          if (activeComponentCount == 0) then
             fitMeasure=0.0d0
             exit
          else
             ! Normalize the fit measure by the number of active components.
             fitMeasure=fitMeasure/dble(activeComponentCount)
          end if
       end do
       ! Check that we found a converged solution.
       if (fitMeasure > adiabaticContractionSolutionTolerance) call Galacticus_Error_Report('Galactic_Structure_Radii_Solve_Adiabatic','failed to find converged solution')

    end if
    return
  end subroutine Galactic_Structure_Radii_Solve_Adiabatic

  subroutine Solve_For_Radius(thisNode,specificAngularMomentum,Radius_Get,Radius_Set,Velocity_Get,Velocity_Set)
    !% Solve for the equilibrium radius of the given component.
    use Dark_Matter_Profiles
    use Dark_Matter_Halo_Scales
    use Numerical_Constants_Physical
    use Galactic_Structure_Rotation_Curves
    use Galactic_Structure_Enclosed_Masses
    use Galactic_Structure_Options
    implicit none
    type(treeNode),                    pointer, intent(inout) :: thisNode
    double precision,                           intent(in)    :: specificAngularMomentum
    procedure(Structure_Get_Template), pointer, intent(in)    :: Radius_Get, Velocity_Get
    procedure(Structure_Set_Template), pointer, intent(in)    :: Radius_Set, Velocity_Set
    double precision                                          :: radius,velocity,virialRadius,angularMomentumC&
         &,angularMomentumCPrimed ,radiusInitial,haloMassInitial,darkMatterMassFinal,darkMatterVelocitySquared&
         &,baryonicVelocitySquared,radiusNew ,specificAngularMomentumPrimed

    ! Count the number of active comonents.
    activeComponentCount=activeComponentCount+1

    if (iterationCount == 1 .or. haloFraction <= 0.0d0) then
       ! On first iteration, see if we have a previous radius set for this component.
       radius=Radius_Get(thisNode)

       if (radius <= 0.0d0) then
          ! No previous radius was set, so make a simple estimate of sizes of all components ignoring adiabatic contraction and self-gravity.

          ! Find the radius in the dark matter profile with the required specific angular momentum
          radius=Dark_Matter_Profile_Radius_from_Specific_Angular_Momentum(haloNode,specificAngularMomentum)
          
          ! Find the velocity at this radius.
          velocity=Dark_Matter_Profile_Circular_Velocity(haloNode,radius)
       else
          ! A previous radius was set, so use it, and the previous circular velocity, as the initial guess.
          velocity=Velocity_Get(thisNode)
       end if

    else
       ! On subsequent iterations do the full calculation providing component has non-zero specific angular momentum.

       if (specificAngularMomentum <= 0.0d0) return
       
       ! Get current radius of the component.
       radius=Radius_Get(thisNode)

       ! Get the virial radius of the node.
       virialRadius=Dark_Matter_Halo_Virial_Radius(haloNode)

       ! Compute the angular momentum parameter, c.
       angularMomentumC=specificAngularMomentum**2/gravitationalConstantGalacticus

       ! Compute the primed angular momentum parameter, c'. To do this, we need r':
       !   r' = r_vir (r/r_vir/A)^(1/omega),
       ! from which we get c' as:
       !  c' = c (r'/r_vir)^(1-omega) / A.
       ! This reduces to:
       !  c' = (c/A) (r/r_vir/A)^(1/omega-1)
       angularMomentumCPrimed=(angularMomentumC/adiabaticContractionGnedinA)*((radius/virialRadius/adiabaticContractionGnedinA)&
            &**inverseAdiabaticContractionGnedinOmegaMinusOne)

       ! Solve for radius in halo with correct pseudo-specific angular momentum.
       specificAngularMomentumPrimed=dsqrt(angularMomentumCPrimed*gravitationalConstantGalacticus)
       radiusInitial=Dark_Matter_Profile_Radius_from_Specific_Angular_Momentum(haloNode,specificAngularMomentumPrimed)

       ! Compute mass within that radius.
       haloMassInitial=Dark_Matter_Profile_Enclosed_Mass(haloNode,radiusInitial)

       ! Compute dark matter mass within final radius.
       darkMatterMassFinal=haloMassInitial*haloFraction

       ! Compute dark matter contribution to rotation curve.
       darkMatterVelocitySquared=gravitationalConstantGalacticus*darkMatterMassFinal/radius

       ! Compute baryonic contribution to rotation curve.
       if (adiabaticContractionIncludeBaryonGravity) then
          baryonicVelocitySquared=Galactic_Structure_Rotation_Curve(thisNode,radius,massType=massTypeGalactic)**2
       else
          baryonicVelocitySquared=0.0d0
       end if

       ! Compute new estimate of velocity.
       velocity=dsqrt(darkMatterVelocitySquared+baryonicVelocitySquared)

       ! Compute new estimate of radius.
       if (radius > 0.0d0) then
          radiusNew=sqrt(specificAngularMomentum/velocity*radius)
       else
          radiusNew=specificAngularMomentum/velocity
       endif
       ! Compute a fit measure.
       if (radius > 0.0d0 .and. radiusNew > 0.0d0) fitMeasure=fitMeasure+dabs(dlog(radiusNew/radius))

       ! Set radius to new radius.
       radius=radiusNew

    end if

    ! Set the component size to new radius and velocity.
    call Radius_Set  (thisNode,radius  )
    call Velocity_Set(thisNode,velocity)
 
    return
  end subroutine Solve_For_Radius

end module Galactic_Structure_Radii_Adiabatic
