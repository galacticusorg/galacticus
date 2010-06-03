!! Copyright 2009, 2010, Andrew Benson <abenson@caltech.edu>
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






!% Contains a module which implements an adiabatic contraction galactic radii solver (including self-gravity of baryons) using the
!% adiabatic contraction algorithm of \cite{gnedin_response_2004}.

module Galactic_Structure_Radii_Adiabatic
  !% Implements an adiabatic contraction galactic radii solver (including self-gravity of baryons) using the adiabatic
  !% contraction algorithm of \cite{gnedin_response_2004}.
  use Tree_Nodes
  use Galactic_Structure_Radius_Solver_Procedures
  private
  public :: Galactic_Structure_Radii_Adiabatic_Initialize

  ! Parameters of the adiabatic contraction algorithm.
  double precision :: adiabaticContractionGnedinA,adiabaticContractionGnedinOmega

  ! Module variables used to communicate current state of radius solver.
  integer          :: iterationCount,activeComponentCount
  double precision :: fitMeasure,haloFraction
  !$omp threadprivate(iterationCount,activeComponentCount,fitMeasure,haloFraction)

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
       !@ </inputParameter>
       call Get_Input_Parameter('adiabaticContractionGnedinA'    ,adiabaticContractionGnedinA    ,defaultValue=0.80d0)
       !@ <inputParameter>
       !@   <name>adiabaticContractionGnedinOmega</name>
       !@   <defaultValue>0.77 (\citealt{gustafsson_baryonic_2006}; from their Fig. 9, strong feedback case)</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The parameter $\omega$ appearing in the \cite{gnedin_response_2004} adiabatic contraction algorithm.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('adiabaticContractionGnedinOmega',adiabaticContractionGnedinOmega,defaultValue=0.77d0)
    end if
    return
  end subroutine Galactic_Structure_Radii_Adiabatic_Initialize

  subroutine Galactic_Structure_Radii_Solve_Adiabatic(thisNode)
    !% Find the radii of galactic components in {\tt thisNode} using the ``adiabatic'' method.
    use Tree_Nodes
    use Tree_Node_Methods
    use Cosmological_Parameters
    use Galacticus_Error
    use Galactic_Structure_Enclosed_Masses
    use Galactic_Structure_Options
    include 'galactic_structure.radius_solver.tasks.modules.inc'
    implicit none
    type(treeNode),          intent(inout), pointer :: thisNode
    integer,                 parameter              :: iterationMaximum=100
    double precision,        parameter              :: fitMeasureAcceptable=1.0d-2
    logical                                         :: componentActive,galaxyIsPhysicallyPlausible
    double precision                                :: specificAngularMomentum

    ! Check that the galaxy is physical plausible. If not, do not try to solve for its structure.
    include 'galactic_structure.radius_solver.plausible.inc'
    if (galaxyIsPhysicallyPlausible) then

       ! Initialize the solver state.
       iterationCount=0
       fitMeasure    =2.0d0*fitMeasureAcceptable
       
       ! Compute fraction of mass distribution as the halo. Truncate this to zero: we can get negative values if the ODE solver is
       ! exploring regimes of high baryonic mass, and this would cause problems.
       haloFraction=max(1.0d0-Galactic_Structure_Enclosed_Mass(thisNode,massType=massTypeGalactic)/Tree_Node_Mass(thisNode),0.0d0)
       
       ! Begin iteration to find a converged solution.
       do while (iterationCount <= 2 .or. ( fitMeasure > fitMeasureAcceptable .and. iterationCount < iterationMaximum ) )
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
       if (fitMeasure > fitMeasureAcceptable) call Galacticus_Error_Report('Galactic_Structure_Radii_Solve_Adiabatic','failed to find converged solution')

    end if
    return
  end subroutine Galactic_Structure_Radii_Solve_Adiabatic
  
  subroutine Solve_For_Radius(thisNode,specificAngularMomentum)
    !% Solve for the equilibrium radius of the given component.
    use Dark_Matter_Profiles
    use Dark_Matter_Halo_Scales
    use Numerical_Constants_Physical
    use Galactic_Structure_Rotation_Curves
    use Galactic_Structure_Enclosed_Masses
    use Galactic_Structure_Options
    implicit none
    type(treeNode),   pointer, intent(inout) :: thisNode
    double precision,          intent(in)    :: specificAngularMomentum
    double precision                         :: radius,velocity,virialRadius,radiusPrimed,angularMomentumC,angularMomentumCPrimed&
         &,radiusInitial,haloMassInitial,darkMatterMassFinal,darkMatterVelocitySquared,baryonicVelocitySquared,radiusNew&
         &,specificAngularMomentumPrimed

    ! Count the number of active comonents.
    activeComponentCount=activeComponentCount+1

    if (iterationCount == 1) then
       ! On first iteration, make a simple estimate of sizes of all components ignoring adiabatic contraction and self-gravity.

       ! Find the radius in the dark matter profile with the required specific angular momentum
       radius=Dark_Matter_Profile_Radius_from_Specific_Angular_Momentum(thisNode,specificAngularMomentum)
       
       ! Find the velocity at this radius.
       velocity=Dark_Matter_Profile_Circular_Velocity(thisNode,radius)
       
    else
       ! On subsequent iterations do the full calculation providing component has non-zero specific angular momentum.

       if (specificAngularMomentum <= 0.0d0) return
       
       ! Get current radius of the component.
       radius=Radius_Get(thisNode)

       ! Get the virial radius of the node.
       virialRadius=Dark_Matter_Halo_Virial_Radius(thisNode)

       ! Compute the primed radius.
       radiusPrimed=virialRadius*((radius/virialRadius/adiabaticContractionGnedinA)**(1.0d0/adiabaticContractionGnedinOmega))

       ! Compute the angular momentum parameter, c.
       angularMomentumC=specificAngularMomentum**2/gravitationalConstantGalacticus

       ! Compute the primed angular momentum parameter, c'.
       angularMomentumCPrimed=angularMomentumC*((radiusPrimed/virialRadius)**(1.0d0-adiabaticContractionGnedinOmega))&
            &/adiabaticContractionGnedinA

       ! Solve for radius in halo with correct pseudo-specific angular momentum.
       specificAngularMomentumPrimed=dsqrt(angularMomentumCPrimed*gravitationalConstantGalacticus)
       radiusInitial=Dark_Matter_Profile_Radius_from_Specific_Angular_Momentum(thisNode,specificAngularMomentumPrimed)

       ! Compute mass within that radius.
       haloMassInitial=Dark_Matter_Profile_Enclosed_Mass(thisNode,radiusInitial)

       ! Compute dark matter mass within final radius.
       darkMatterMassFinal=haloMassInitial*haloFraction

       ! Compute dark matter contribution to rotation curve.
       darkMatterVelocitySquared=gravitationalConstantGalacticus*darkMatterMassFinal/radius

       ! Compute baryonic contribution to rotation curve.
       baryonicVelocitySquared=Galactic_Structure_Rotation_Curve(thisNode,radius,massType=massTypeGalactic)**2

       ! Compute new estimate of velocity.
       velocity=dsqrt(darkMatterVelocitySquared+baryonicVelocitySquared)

       ! Compute new estimate of radius.
       radiusNew=specificAngularMomentum/velocity

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
