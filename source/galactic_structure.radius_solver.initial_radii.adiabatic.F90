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

!% Contains a module which implements calculations of initial radius in the dark matter halo using the adiabatic contraction
!% algorithm of \cite{gnedin_response_2004}.

module Galactic_Structure_Initial_Radii_Adiabatic
  !% Implements calculations of initial radius in the dark matter halo using the adiabatic contraction algorithm of
  !% \cite{gnedin_response_2004}.
  use Galacticus_Nodes
  use Galactic_Structure_Options
  implicit none
  private
  public :: Galactic_Structure_Initial_Radii_Adiabatic_Initialize

  ! Parameters of the adiabatic contraction algorithm.
  double precision                      :: adiabaticContractionGnedinA,adiabaticContractionGnedinOmega

  ! Module scope quantities used in solving the initial radius root function.
  integer                   , parameter :: componentType=componentTypeAll,massType=massTypeBaryonic,weightBy=weightByMass,weightIndex=weightIndexNull
  logical                   , parameter :: haloLoaded=.false.
  double precision                      :: radiusShared,baryonicFinalTerm,baryonicFinalTermDerivative,darkMatterFraction,initialMassFraction,radiusFinal,radiusFinalMean,virialRadius,radiusFinalMeanSelfDerivative,radiusInitialMeanSelfDerivative,radiusInitial
  type            (treeNode), pointer   :: activeNode
  !$omp threadprivate(radiusShared,baryonicFinalTerm,baryonicFinalTermDerivative,darkMatterFraction,initialMassFraction,radiusFinal,radiusFinalMean,virialRadius,radiusFinalMeanSelfDerivative,radiusInitialMeanSelfDerivative,radiusInitial,activeNode)

contains

  !# <galacticStructureRadiusSolverInitialRadiusMethod>
  !#  <unitName>Galactic_Structure_Initial_Radii_Adiabatic_Initialize</unitName>
  !# </galacticStructureRadiusSolverInitialRadiusMethod>
  subroutine Galactic_Structure_Initial_Radii_Adiabatic_Initialize(galacticStructureRadiusSolverInitialRadiusMethod,Galactic_Structure_Radius_Initial_Get,Galactic_Structure_Radius_Initial_Derivative_Get)
    !% Initializes the ``adiabatic'' initial radii module.
    use Input_Parameters
    use ISO_Varying_String
    implicit none
    type     (varying_string  ),          intent(in   ) :: galacticStructureRadiusSolverInitialRadiusMethod
    procedure(double precision), pointer, intent(inout) :: Galactic_Structure_Radius_Initial_Get&
         &,Galactic_Structure_Radius_Initial_Derivative_Get
    
    if (galacticStructureRadiusSolverInitialRadiusMethod == 'adiabatic') then
       Galactic_Structure_Radius_Initial_Get            => Galactic_Structure_Radius_Initial_Adiabatic
       Galactic_Structure_Radius_Initial_Derivative_Get => Galactic_Structure_Radius_Initial_Derivative_Adiabatic
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
    end if
    return
  end subroutine Galactic_Structure_Initial_Radii_Adiabatic_Initialize

  subroutine Galactic_Structure_Radii_Initial_Adiabatic_Compute_Factors(thisNode,radius,computeGradientFactors)
    !% Compute various factors needed when solving for the initial radius in the dark matter halo using the adiabatic contraction
    !% algorithm of \cite{gnedin_response_2004}.
    use Dark_Matter_Profiles
    use Dark_Matter_Halo_Scales
    use Cosmological_Parameters
    use Galactic_Structure_Options
    use Numerical_Constants_Physical
    !# <include directive="rotationCurveTask" name="radiusSolverRotationCurveTask" type="moduleUse">
    !# <exclude>Dark_Matter_Profile_Structure_Tasks</exclude>
    include 'galactic_structure.radius_solver.initial_radii.adiabatic.rotation_curve.tasks.modules.inc'
    !# </include>
    !# <include directive="rotationCurveGradientTask" name="radiusSolverRotationCurveGradientTask" type="moduleUse">
    !# <exclude>Dark_Matter_Profile_Structure_Tasks</exclude>
    include 'galactic_structure.radius_solver.initial_radii.adiabatic.rotation_curve_gradient.tasks.modules.inc'
    !# </include>
    !# <include directive="enclosedMassTask" name="radiusSolverEnclosedMassTask" type="moduleUse">
    !# <exclude>Dark_Matter_Profile_Structure_Tasks</exclude>
    include 'galactic_structure.radius_solver.initial_radii.adiabatic.enclosed_mass.tasks.modules.inc'
    !# </include>
    implicit none
    type            (treeNode                         ), pointer, intent(inout) :: thisNode
    double precision                                   ,          intent(in   ) :: radius
    logical                                            ,          intent(in   ) :: computeGradientFactors    
    type            (treeNode                         ), pointer                :: currentNode
    class           (nodeComponentBasic               ), pointer                :: thisBasic
    double precision                                   , parameter              :: toleranceAbsolute=0.0d0,toleranceRelative=1.0d-3
    procedure       (Component_Enclosed_Mass          ), pointer                :: componentEnclosedMass
    procedure       (Component_Rotation_Curve         ), pointer                :: componentRotationCurve
    procedure       (Component_Rotation_Curve_Gradient), pointer                :: componentRotationCurveGradient
    double precision                                                            :: rotationCurveSquared &
         &,rotationCurveSquaredGradient ,componentVelocitySquaredGradient,baryonicMassTotal,baryonicMassSelfTotal,componentMass&
         &,componentVelocity

    ! Get the virial radius of the node.
    virialRadius=Dark_Matter_Halo_Virial_Radius(thisNode)
    ! Store the final radius and its orbit-averaged mean.
    radiusFinal    =                                     radius
    radiusFinalMean=Adiabatic_Solver_Mean_Orbital_Radius(radius)
    radiusShared   =radiusFinalMean
    ! Compute the baryonic contribution to the rotation curve.
    currentNode            => thisNode
    componentRotationCurve => Component_Rotation_Curve
    rotationCurveSquared=currentNode%mapDouble0(componentRotationCurve,reductionSummation)
    !# <include directive="rotationCurveTask" name="radiusSolverRotationCurveTask" type="functionCall" functionType="function" returnParameter="componentVelocity">
    !#  <exclude>Dark_Matter_Profile_Rotation_Curve_Task</exclude>
    !#  <functionArgs>currentNode,radiusFinalMean,componentType,massType,haloLoaded</functionArgs>
    !#  <onReturn>rotationCurveSquared=rotationCurveSquared+componentVelocity**2</onReturn>
    include 'galactic_structure.radius_solver.initial_radii.adiabatic.rotation_curve.tasks.inc'
    !# </include>
    baryonicFinalTerm=rotationCurveSquared*radiusFinalMean*radiusFinal/gravitationalConstantGalacticus
    ! Compute the baryonic contribution to the rotation curve.
    if (computeGradientFactors) then
       componentRotationCurveGradient => Component_Rotation_Curve_Gradient
       rotationCurveSquaredGradient   =currentNode%mapDouble0(componentRotationCurveGradient,reductionSummation)
       !# <include directive="rotationCurveGradientTask" name="radiusSolverRotationCurveGradientTask" type="functionCall" functionType="function" returnParameter="componentVelocitySquaredGradient">
       !#  <exclude>Dark_Matter_Profile_Rotation_Curve_Gradient_Task</exclude>
       !#  <functionArgs>thisNode,radiusFinalMean,componentType,massType,haloLoaded</functionArgs>
       !#  <onReturn>rotationCurveSquaredGradient=rotationCurveSquaredGradient+componentVelocitySquaredGradient</onReturn>
       include 'galactic_structure.radius_solver.initial_radii.adiabatic.rotation_curve_gradient.tasks.inc'
       !# </include>
       baryonicFinalTermDerivative=rotationCurveSquaredGradient*Adiabatic_Solver_Mean_Orbital_Radius_Derivative(radiusFinal)*radiusFinalMean*radiusFinal/gravitationalConstantGalacticus
    end if
    ! Compute the initial baryonic contribution from this halo, and any satellites.
    baryonicMassTotal=0.0d0
    currentNode => thisNode
    do while (associated(currentNode))
       componentEnclosedMass => Component_Enclosed_Mass
       baryonicMassTotal=baryonicMassTotal+currentNode%mapDouble0(componentEnclosedMass,reductionSummation)
       !# <include directive="enclosedMassTask" name="radiusSolverEnclosedMassTask" type="functionCall" functionType="function" returnParameter="componentMass">
       !#  <exclude>Dark_Matter_Profile_Enclosed_Mass_Task</exclude>
       !#  <functionArgs>currentNode,virialRadius,componentType,massType,weightBy,weightIndex,haloLoaded</functionArgs>
       !#  <onReturn>baryonicMassTotal=baryonicMassTotal+componentMass</onReturn>
       include 'galactic_structure.radius_solver.initial_radii.adiabatic.enclosed_mass.tasks.inc'
       !# </include>
       if (associated(currentNode,thisNode)) then
          baryonicMassSelfTotal=baryonicMassTotal
          do while (associated(currentNode%firstSatellite))
             currentNode => currentNode%firstSatellite
          end do
          if (associated(currentNode,thisNode)) currentNode => null()
       else
          if (associated(currentNode%sibling)) then
             currentNode => currentNode%sibling
             do while (associated(currentNode%firstSatellite))
                currentNode => currentNode%firstSatellite
             end do
          else
             currentNode => currentNode%parent
             if (associated(currentNode,thisNode)) currentNode => null()
          end if
       end if
    end do
    ! Limit masses to physical values.
    baryonicMassSelfTotal=max(baryonicMassSelfTotal,0.0d0)
    baryonicMassTotal    =max(baryonicMassTotal    ,0.0d0)
    ! Compute the dark matter fraction.
    thisBasic => thisNode%basic()
    darkMatterFraction =(Omega_Matter()-Omega_B())/Omega_Matter()+(baryonicMassTotal-baryonicMassSelfTotal)/thisBasic%mass()
    ! Compute the initial mass fraction.
    initialMassFraction=(Omega_Matter()-Omega_B())/Omega_Matter()+                   baryonicMassSelfTotal /thisBasic%mass()
    ! Store the current node.
    activeNode => thisNode
    return
  end subroutine Galactic_Structure_Radii_Initial_Adiabatic_Compute_Factors

  double precision function Galactic_Structure_Radius_Initial_Adiabatic(thisNode,radius)
    !% Compute the initial radius in the dark matter halo using the adiabatic contraction algorithm of
    !% \cite{gnedin_response_2004}.
    use Root_Finder
    use FGSL
    use Galacticus_Error
    use, intrinsic :: ISO_C_Binding
    implicit none
    type            (treeNode         ), pointer, intent(inout) :: thisNode
    double precision                   ,          intent(in   ) :: radius
    type            (fgsl_function    ), save                   :: rootFunction
    type            (fgsl_root_fsolver), save                   :: rootFunctionSolver
    !$omp threadprivate(rootFunction,rootFunctionSolver)
    double precision                   , parameter              :: toleranceAbsolute=0.0d0,toleranceRelative=1.0d-3
    double precision                                            :: radiusMinimum,radiusMaximum
    type            (c_ptr            )                         :: parameterPointer

    ! Compute the various factors needed by this calculation.
    call Galactic_Structure_Radii_Initial_Adiabatic_Compute_Factors(thisNode,radius,computeGradientFactors=.false.)
    ! Return radius unchanged if larger than the virial radius.
    if (radius >= virialRadius) then
       Galactic_Structure_Radius_Initial_Adiabatic=radius
       return
    end if
    ! Choose suitable minimum and maximum radii.
    radiusMinimum=radius
    radiusMaximum=virialRadius
    ! Decrease minimum radius as necessary
    do while (Galactic_Structure_Radius_Initial_Adiabatic_Solver(radiusMinimum,parameterPointer) > 0.0d0)
       radiusMinimum=0.5d0*radiusMinimum
    end do
    ! Check that solution is within bounds.
    if (Galactic_Structure_Radius_Initial_Adiabatic_Solver(radiusMaximum,parameterPointer) < 0.0d0) then
       Galactic_Structure_Radius_Initial_Adiabatic=virialRadius
       return
    end if
    ! Find the solution for initial radius.
    Galactic_Structure_Radius_Initial_Adiabatic=                         &
         & Root_Find(                                                    &
         &           radiusMinimum                                     , &
         &           radiusMaximum                                     , &
         &           Galactic_Structure_Radius_Initial_Adiabatic_Solver, &
         &           parameterPointer                                  , &
         &           rootFunction                                      , &
         &           rootFunctionSolver                                , &
         &           toleranceAbsolute                                 , &
         &           toleranceRelative                                   &
         &          )
    return
  end function Galactic_Structure_Radius_Initial_Adiabatic

  double precision function Galactic_Structure_Radius_Initial_Derivative_Adiabatic(thisNode,radius)
    !% Compute the derivative of the initial radius in the dark matter halo using the adiabatic contraction algorithm of
    !% \cite{gnedin_response_2004}.
    use Root_Finder
    use FGSL
    use Galacticus_Error
    use, intrinsic :: ISO_C_Binding
    implicit none
    type            (treeNode         ), pointer, intent(inout) :: thisNode
    double precision                   ,          intent(in   ) :: radius
    type            (fgsl_function    ), save                   :: rootFunction
    type            (fgsl_root_fsolver), save                   :: rootFunctionSolver
    !$omp threadprivate(rootFunction,rootFunctionSolver)
    double precision                   , parameter              :: toleranceAbsolute=0.0d0,toleranceRelative=1.0d-3
    double precision                                            :: radiusDerivativeMinimum,radiusDerivativeMaximum
    type            (c_ptr            )                         :: parameterPointer

    ! Compute the various factors needed by this calculation.
    call Galactic_Structure_Radii_Initial_Adiabatic_Compute_Factors(thisNode,radius,computeGradientFactors=.true.)
    ! Return unit derivative if radius is larger than the virial radius.
    if (radius >= virialRadius) then
       Galactic_Structure_Radius_Initial_Derivative_Adiabatic=1.0d0
       return
    end if
    ! Compute initial radius, and derivatives of initial and final mean radii.
    radiusFinal                    =                                                         radius
    radiusInitial                  =Galactic_Structure_Radius_Initial_Adiabatic    (thisNode,radius       )
    radiusInitialMeanSelfDerivative=Adiabatic_Solver_Mean_Orbital_Radius_Derivative(         radiusInitial)
    radiusFinalMeanSelfDerivative  =Adiabatic_Solver_Mean_Orbital_Radius_Derivative(         radius       )    
    ! Choose suitable minimum and maximum radii.
    radiusDerivativeMinimum=1.0d0
    radiusDerivativeMaximum=1.0d0
    ! Decrease minimum derivative as necessary
    do while (Galactic_Structure_Radius_Initial_Derivative_Adiabatic_Solver(radiusDerivativeMinimum,parameterPointer) > 0.0d0)
       radiusDerivativeMinimum=0.5d0*radiusDerivativeMinimum
       ! Catch cases where the gradient falls below zero (which should never happen).
       if (radiusDerivativeMinimum <= 0.0d0) call Galacticus_Error_Report('Galactic_Structure_Radius_Initial_Derivative_Adiabatic','gradient has reached zero [should not happen]')
    end do
    ! Decrease maximum derivative as necessary
    do while (Galactic_Structure_Radius_Initial_Derivative_Adiabatic_Solver(radiusDerivativeMaximum,parameterPointer) < 0.0d0)
       radiusDerivativeMaximum=2.0d0*radiusDerivativeMaximum
    end do
    ! Find the solution for initial radius.
    Galactic_Structure_Radius_Initial_Derivative_Adiabatic=                         &
         & Root_Find(                                                               &
         &           radiusDerivativeMinimum                                      , &
         &           radiusDerivativeMaximum                                      , &
         &           Galactic_Structure_Radius_Initial_Derivative_Adiabatic_Solver, &
         &           parameterPointer                                             , &
         &           rootFunction                                                 , &
         &           rootFunctionSolver                                           , &
         &           toleranceAbsolute                                            , &
         &           toleranceRelative                                              &
         &          )
    return
  end function Galactic_Structure_Radius_Initial_Derivative_Adiabatic

  function Galactic_Structure_Radius_Initial_Adiabatic_Solver(radiusInitial,parameterPointer) bind(c)
    !% Root function used in finding the initial radius in the dark matter halo when solving for adiabatic contraction.
    use, intrinsic :: ISO_C_Binding
    use Dark_Matter_Profiles
    implicit none
    real            (c_double), value :: radiusInitial
    type            (c_ptr   ), value :: parameterPointer
    real            (c_double)        :: Galactic_Structure_Radius_Initial_Adiabatic_Solver
    double precision                  :: darkMatterMassInitial,radiusInitialMean

    ! Find the initial mean orbital radius.
    radiusInitialMean    =Adiabatic_Solver_Mean_Orbital_Radius(           radiusInitial    )
    ! Get the mass of dark matter inside the initial radius.
    darkMatterMassInitial=Dark_Matter_Profile_Enclosed_Mass   (activeNode,radiusInitialMean)
    ! Compute the root function.
    Galactic_Structure_Radius_Initial_Adiabatic_Solver= &
         & darkMatterMassInitial                        &
         & *(                                           &
         &    initialMassFraction*radiusInitial         &
         &   -darkMatterFraction *radiusFinal           &
         &  )                                           &
         & -baryonicFinalTerm
    return
  end function Galactic_Structure_Radius_Initial_Adiabatic_Solver

  function Galactic_Structure_Radius_Initial_Derivative_Adiabatic_Solver(radiusInitialDerivative,parameterPointer) bind(c)
    !% Root function used in finding the derivative of the initial radius in the dark matter halo when solving for adiabatic
    !% contraction.
    use, intrinsic :: ISO_C_Binding
    use Dark_Matter_Profiles
    use Numerical_Constants_Math
    implicit none
    real            (c_double), value :: radiusInitialDerivative
    type            (c_ptr   ), value :: parameterPointer
    real            (c_double)        :: Galactic_Structure_Radius_Initial_Derivative_Adiabatic_Solver
    double precision                  :: darkMatterMassInitial,darkMatterDensityInitial,radiusInitialMean

    ! Find the initial mean orbital radius.
    radiusInitialMean       =Adiabatic_Solver_Mean_Orbital_Radius(           radiusInitial    )
    ! Get the mass of dark matter inside the initial radius.
    darkMatterMassInitial   =Dark_Matter_Profile_Enclosed_Mass   (activeNode,radiusInitialMean)
    ! Get the mass of dark matter inside the initial radius.
    darkMatterDensityInitial=Dark_Matter_Profile_Density         (activeNode,radiusInitialMean)
    ! Compute the root function.
    Galactic_Structure_Radius_Initial_Derivative_Adiabatic_Solver= &
         & darkMatterMassInitial                                   &
         & *(                                                      &
         &    initialMassFraction*radiusInitialDerivative          &
         &   -darkMatterFraction                                   &
         &  )                                                      &
         & +(                                                      &
         &    initialMassFraction*radiusInitial                    &
         &   -darkMatterFraction *radiusFinal                      &
         &  )                                                      &
         &  *4.0d0                                                 &
         &  *Pi                                                    &
         &  *radiusInitialMean**2                                  &
         &  *darkMatterDensityInitial                              &
         &  *radiusInitialMeanSelfDerivative                       &
         &  *radiusInitialDerivative                               &
         &  -baryonicFinalTerm                                     &
         &  *(                                                     &
         &     1.0d0                        /radiusFinal           &
         &    +radiusFinalMeanSelfDerivative/radiusFinalMean       &
         &   )                                                     &
         &  -baryonicFinalTermDerivative
    return
  end function Galactic_Structure_Radius_Initial_Derivative_Adiabatic_Solver

  double precision function Adiabatic_Solver_Mean_Orbital_Radius(radius)
    !% Returns the orbit averaged radius for dark matter corresponding the given {\tt radius} using the model of
    !% \cite{gnedin_response_2004}.
    implicit none
    double precision, intent(in) :: radius
    
    Adiabatic_Solver_Mean_Orbital_Radius=                                                        &
         &                                adiabaticContractionGnedinA                            &
         &                               *virialRadius                                           &
         &                               *(radius/virialRadius)**adiabaticContractionGnedinOmega
    return
  end function Adiabatic_Solver_Mean_Orbital_Radius

  double precision function Adiabatic_Solver_Mean_Orbital_Radius_Derivative(radius)
    !% Returns the derivative of the orbit averaged radius for dark matter corresponding the given {\tt radius} using the model of
    !% \cite{gnedin_response_2004}.
    implicit none
    double precision, intent(in) :: radius
    
    Adiabatic_Solver_Mean_Orbital_Radius_Derivative=                                             &
         &                        adiabaticContractionGnedinA                                    &
         &                       *(radius/virialRadius)**(adiabaticContractionGnedinOmega-1.0d0)
    return
  end function Adiabatic_Solver_Mean_Orbital_Radius_Derivative

  double precision function Component_Enclosed_Mass(component)
    !% Unary function returning the enclosed mass in a component. Suitable for mapping over components. Ignores the dark matter
    !% profile.
    implicit none
    class  (nodeComponent), intent(inout) :: component
 
    select type (component)
    class is (nodeComponentDarkMatterProfile)
       Component_Enclosed_Mass=0.0d0
    class default
       Component_Enclosed_Mass=component%enclosedMass(virialRadius,componentType ,massType,weightBy,weightIndex,haloLoaded)
    end select
    return
  end function Component_Enclosed_Mass
    
  double precision function Component_Rotation_Curve(component)
    !% Unary function returning the squared rotation curve in a component. Suitable for mapping over components.
    use Galacticus_Nodes
    implicit none
    class(nodeComponent), intent(inout) :: component
 
    select type (component)
    class is (nodeComponentDarkMatterProfile)
       Component_Rotation_Curve=0.0d0
    class default
       Component_Rotation_Curve=component%rotationCurve(radiusShared,componentType,massType,haloLoaded)**2
    end select
    return
  end function Component_Rotation_Curve
    
  double precision function Component_Rotation_Curve_Gradient(component)
    !% Unary function returning the squared rotation curve gradient in a component. Suitable for mapping over components.
    use Galacticus_Nodes
    implicit none
    class(nodeComponent), intent(inout) :: component
 
    select type (component)
    class is (nodeComponentDarkMatterProfile)
       Component_Rotation_Curve_Gradient=0.0d0
    class default
       Component_Rotation_Curve_Gradient=component%rotationCurveGradient(radiusShared,componentType,massType,haloLoaded)
    end select
    return
  end function Component_Rotation_Curve_Gradient

end module Galactic_Structure_Initial_Radii_Adiabatic
