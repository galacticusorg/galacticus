!% Contains a module which implements NFW halo profiles.

module Dark_Matter_Profiles_NFW
  !% Implements NFW halo profiles.
  use Tree_Nodes
  use FGSL
  private
  public :: Dark_Matter_Profile_NFW_Initialize, Dark_Matter_Profiles_NFW_State_Store, Dark_Matter_Profiles_NFW_State_Retrieve

  ! Minimum and maximum concentrations to tabulate.
  double precision                            :: concentrationMinimum   = 1.0d0
  double precision                            :: concentrationMaximum   =20.0d0
  ! Minimum and maximum radii to tabulate.
  double precision                            :: radiusMinimum          = 1.0d-3
  double precision                            :: radiusMaximum          = 1.0d+2
  double precision                            :: specificAngularMomentumMinimum
  double precision                            :: specificAngularMomentumMaximum
  ! Number of points per decade of concentration in NFW tabulations.
  integer,          parameter                 :: nfwTablePointsPerDecade       =100
  integer,          parameter                 :: nfwInverseTablePointsPerDecade=100
  ! Tables of NFW properties.
  logical                                     :: nfwTableInitialized=.false.,nfwInverseTableInitialized=.false.
  integer                                     :: nfwTableNumberPoints,nfwInverseTableNumberPoints
  double precision, allocatable, dimension(:) :: nfwConcentration,nfwEnergy,nfwRotationNormalization,nfwRadius&
       &,nfwSpecificAngularMomentum

  ! Interpolator variables.
  type(fgsl_interp)                           :: interpolationObject      ,interpolationInverseObject
  type(fgsl_interp_accel)                     :: interpolationAccelerator ,interpolationInverseAccelerator
  logical                                     :: interpolationReset=.true.,interpolationInverseReset=.true.

  ! Module variables used in integrations.
  double precision                            :: concentrationParameter

contains

  !# <darkMatterProfileMethod>
  !#  <unitName>Dark_Matter_Profile_NFW_Initialize</unitName>
  !# </darkMatterProfileMethod>
  subroutine Dark_Matter_Profile_NFW_Initialize(darkMatterProfileMethod,Dark_Matter_Profile_Energy_Get&
       &,Dark_Matter_Profile_Energy_Growth_Rate_Get,Dark_Matter_Profile_Rotation_Normalization_Get &
       &,Dark_Matter_Profile_Radius_from_Specific_Angular_Momentum_Get,Dark_Matter_Profile_Circular_Velocity_Get&
       &,Dark_Matter_Profile_Potential_Get,Dark_Matter_Profile_Enclosed_Mass_Get)
    !% Initializes the ``NFW'' halo profile module.
    use ISO_Varying_String
    use Tree_Node_Methods
    use Galacticus_Error
    implicit none
    type(varying_string),          intent(in)    :: darkMatterProfileMethod
    procedure(),          pointer, intent(inout) :: Dark_Matter_Profile_Energy_Get,Dark_Matter_Profile_Energy_Growth_Rate_Get&
         &,Dark_Matter_Profile_Rotation_Normalization_Get,Dark_Matter_Profile_Radius_from_Specific_Angular_Momentum_Get&
         &,Dark_Matter_Profile_Circular_Velocity_Get,Dark_Matter_Profile_Potential_Get,Dark_Matter_Profile_Enclosed_Mass_Get
    
    if (darkMatterProfileMethod == 'NFW') then
       Dark_Matter_Profile_Energy_Get                                => Dark_Matter_Profile_Energy_NFW
       Dark_Matter_Profile_Energy_Growth_Rate_Get                    => Dark_Matter_Profile_Energy_Growth_Rate_NFW
       Dark_Matter_Profile_Rotation_Normalization_Get                => Dark_Matter_Profile_Rotation_Normalization_NFW
       Dark_Matter_Profile_Radius_from_Specific_Angular_Momentum_Get => Radius_from_Specific_Angular_Momentum_NFW
       Dark_Matter_Profile_Circular_Velocity_Get                     => Dark_Matter_Profile_Circular_Velocity_NFW
       Dark_Matter_Profile_Potential_Get                             => Dark_Matter_Profile_Potential_NFW
       Dark_Matter_Profile_Enclosed_Mass_Get                         => Dark_Matter_Profile_Enclosed_Mass_NFW       
       ! Ensure that the dark matter profile component supports a "scale" property. Since we've been called with a treeNode to
       ! process, it should have been initialized by now.
       if (.not.associated(Tree_Node_Dark_Matter_Profile_Scale)) call&
            & Galacticus_Error_Report('Dark_Matter_Profile_NFW_Initialize','NFW dark matter profile requires a dark matter&
            & profile component that supports the "scale" property')
       ! Initialize the tabulations.
       call Dark_Matter_Profile_NFW_Tabulate
       call Dark_Matter_Profile_NFW_Inverse_Angular_Momentum
    end if
    return
  end subroutine Dark_Matter_Profile_NFW_Initialize

  subroutine Dark_Matter_Profile_NFW_Tabulate(concentration)
    !% Tabulate properties of the NFW halo profile which must be computed numerically.
    use Numerical_Ranges
    use Memory_Management
    implicit none
    double precision, intent(in), optional :: concentration
    integer                                :: iConcentration
    logical                                :: retabulate

    !$omp critical (Dark_Matter_Profile_NFW_Tabulate)
    retabulate=.not.nfwTableInitialized
    if (present(concentration)) then
       if (concentration < concentrationMinimum) then
          concentrationMinimum=0.5d0*concentration
          retabulate=.true.
       end if
       if (concentration > concentrationMaximum) then
          concentrationMaximum=2.0d0*concentration
          retabulate=.true.
       end if
    end if
    if (retabulate) then
       ! Decide how many points to tabulate and allocate table arrays.
       nfwTableNumberPoints=int(dlog10(concentrationMaximum/concentrationMinimum)*dble(nfwTablePointsPerDecade))+1
       if (allocated(nfwConcentration)) then
          call Dealloc_Array(nfwConcentration        )
          call Dealloc_Array(nfwEnergy               )
          call Dealloc_Array(nfwRotationNormalization)
       end if
       call Alloc_Array(nfwConcentration        ,nfwTableNumberPoints,'nfwConcentration'        )
       call Alloc_Array(nfwEnergy               ,nfwTableNumberPoints,'nfwEnergy'               )
       call Alloc_Array(nfwRotationNormalization,nfwTableNumberPoints,'nfwRotationNormalization')
       ! Create a range of concentrations.
       nfwConcentration=Make_Range(concentrationMinimum,concentrationMaximum,nfwTableNumberPoints,rangeType=rangeTypeLogarithmic)
       ! Loop over concentrations and populate tables.
       do iConcentration=1,nfwTableNumberPoints
          nfwEnergy(iConcentration)               =NFW_Profile_Energy(nfwConcentration(iConcentration))
          nfwRotationNormalization(iConcentration)=nfwConcentration(iConcentration)&
               &/Angular_Momentum_NFW_Scale_Free(nfwConcentration(iConcentration))
       end do
       ! Ensure interpolations get reset.
       interpolationReset=.true.
       ! Specify that tabulation has been made.
       nfwTableInitialized=.true.
    end if
    !$omp end critical (Dark_Matter_Profile_NFW_Tabulate)
    return
  end subroutine Dark_Matter_Profile_NFW_Tabulate

  subroutine Dark_Matter_Profile_NFW_Inverse_Angular_Momentum(specificAngularMomentum)
    !% Tabulates the specific angular momentum vs. radius in an NFW profile for rapid inversion.
    use Numerical_Ranges
    use Memory_Management
    implicit none
    double precision, intent(in), optional :: specificAngularMomentum
    integer                                :: iRadius
    logical                                :: retabulate

    !$omp critical (Dark_Matter_Profile_NFW_Inverse_AM)
    retabulate=.not.nfwInverseTableInitialized
    if (present(specificAngularMomentum)) then
       specificAngularMomentumMinimum=Specific_Angular_Momentum_NFW_Scale_Free(radiusMinimum)
       specificAngularMomentumMaximum=Specific_Angular_Momentum_NFW_Scale_Free(radiusMaximum)
       do while (specificAngularMomentum < specificAngularMomentumMinimum)
          radiusMinimum=0.5d0*radiusMinimum
          specificAngularMomentumMinimum=Specific_Angular_Momentum_NFW_Scale_Free(radiusMinimum)
          retabulate=.true.
       end do
       do while (specificAngularMomentum > specificAngularMomentumMaximum)
          radiusMaximum=2.0d0*radiusMaximum
          specificAngularMomentumMaximum=Specific_Angular_Momentum_NFW_Scale_Free(radiusMaximum)
          retabulate=.true.
       end do
    end if
    if (retabulate) then
       ! Decide how many points to tabulate and allocate table arrays.
       nfwInverseTableNumberPoints=int(dlog10(radiusMaximum/radiusMinimum)*dble(nfwInverseTablePointsPerDecade))+1
       if (allocated(nfwRadius)) then
          call Dealloc_Array(nfwRadius                 )
          call Dealloc_Array(nfwSpecificAngularMomentum)
       end if
       call Alloc_Array(nfwRadius                 ,nfwInverseTableNumberPoints,'nfwRadius'                 )
       call Alloc_Array(nfwSpecificAngularMomentum,nfwInverseTableNumberPoints,'nfwSpecificAngularMomentum')
       ! Create a range of radii.
       nfwRadius=Make_Range(radiusMinimum,radiusMaximum,nfwInverseTableNumberPoints,rangeType=rangeTypeLogarithmic)
       ! Loop over radii and populate tables.
       do iRadius=1,nfwInverseTableNumberPoints
          nfwSpecificAngularMomentum(iRadius)=Specific_Angular_Momentum_NFW_Scale_Free(nfwRadius(iRadius))
       end do
       ! Ensure interpolations get reset.
       interpolationInverseReset=.true.
       ! Specify that tabulation has been made.
       nfwInverseTableInitialized=.true.
    end if
    !$omp end critical (Dark_Matter_Profile_NFW_Inverse_AM)
    return
  end subroutine Dark_Matter_Profile_NFW_Inverse_Angular_Momentum

  double precision function Dark_Matter_Profile_Enclosed_Mass_NFW(thisNode,radius)
    !% Returns the enclosed mass (in $M_\odot$) in the dark matter profile of {\tt thisNode} at the given {\tt radius} (given in
    !% units of Mpc).
    use Tree_Nodes
    use Tree_Node_Methods
    use Dark_Matter_Halo_Scales
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, intent(in)             :: radius
    double precision                         :: radiusOverScaleRadius,virialRadiusOverScaleRadius

    radiusOverScaleRadius      =radius                                  /Tree_Node_Dark_Matter_Profile_Scale(thisNode)
    virialRadiusOverScaleRadius=Dark_Matter_Halo_Virial_Radius(thisNode)/Tree_Node_Dark_Matter_Profile_Scale(thisNode)
    Dark_Matter_Profile_Enclosed_Mass_NFW=Enclosed_Mass_NFW_Scale_Free(radiusOverScaleRadius,virialRadiusOverScaleRadius)&
         &*Tree_Node_Mass(thisNode)
    return
  end function Dark_Matter_Profile_Enclosed_Mass_NFW

  double precision function Dark_Matter_Profile_Potential_NFW(thisNode,radius)
    !% Returns the potential (in (km/s)$^2$) in the dark matter profile of {\tt thisNode} at the given {\tt radius} (given in
    !% units of Mpc).
    use Tree_Nodes
    use Tree_Node_Methods
    use Dark_Matter_Halo_Scales
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, intent(in)             :: radius
    double precision                         :: radiusOverScaleRadius,virialRadiusOverScaleRadius

    radiusOverScaleRadius            =radius                                  /Tree_Node_Dark_Matter_Profile_Scale(thisNode)
    virialRadiusOverScaleRadius      =Dark_Matter_Halo_Virial_Radius(thisNode)/Tree_Node_Dark_Matter_Profile_Scale(thisNode)
    Dark_Matter_Profile_Potential_NFW=(-1.0d0-virialRadiusOverScaleRadius*(dlog(1.0d0+radiusOverScaleRadius)&
         &/radiusOverScaleRadius-dlog(1.0d0+virialRadiusOverScaleRadius)/virialRadiusOverScaleRadius)/(dlog(1.0d0&
         &+virialRadiusOverScaleRadius)-virialRadiusOverScaleRadius/(1.0d0+virialRadiusOverScaleRadius)))&
         &*Dark_Matter_Halo_Virial_Velocity(thisNode)**2
    return
  end function Dark_Matter_Profile_Potential_NFW
  
  double precision function Dark_Matter_Profile_Circular_Velocity_NFW(thisNode,radius)
    !% Returns the circular velocity (in km/s) in the dark matter profile of {\tt thisNode} at the given {\tt radius} (given in
    !% units of Mpc). For an NFW halo this is independent of radius and therefore equal to the virial velocity.
    use Tree_Nodes
    use Dark_Matter_Halo_Scales
    use Numerical_Constants_Physical
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, intent(in)             :: radius

    if (radius > 0.0d0) then
       Dark_Matter_Profile_Circular_Velocity_NFW=dsqrt(gravitationalConstantGalacticus&
            &*Dark_Matter_Profile_Enclosed_Mass_NFW(thisNode,radius)/radius)
    else
       Dark_Matter_Profile_Circular_Velocity_NFW=0.0d0
    end if
    return
  end function Dark_Matter_Profile_Circular_Velocity_NFW
  
  double precision function Radius_from_Specific_Angular_Momentum_NFW(thisNode,specificAngularMomentum)
    !% Returns the radius (in Mpc) in {\tt thisNode} at which a circular orbit has the given {\tt specificAngularMomentum} (given
    !% in units of km s$^{-1}$ Mpc). For an NFW halo, the circular velocity is constant (and therefore equal to the virial
    !% velocity). Therefore, $r = j/V_{\rm virial}$ where $j$(={\tt specificAngularMomentum}) is the specific angular momentum and
    !% $r$ the required radius.
    use Tree_Nodes
    use Tree_Node_Methods
    use Dark_Matter_Halo_Scales
    use Numerical_Interpolation
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision, intent(in)             :: specificAngularMomentum
    double precision                         :: scaleRadius,specificAngularMomentumScaleFree

    if (specificAngularMomentum <= 0.0d0) then
       Radius_from_Specific_Angular_Momentum_NFW=0.0d0
       return
    end if

    ! Get the scale radius.
    scaleRadius=Tree_Node_Dark_Matter_Profile_Scale(thisNode)
    
    ! Compute the specific angular momentum in scale free units (using the scale length for distances the sqrt(G M(r_scale) /
    ! r_scale) for velocities).
    specificAngularMomentumScaleFree=specificAngularMomentum/scaleRadius/Dark_Matter_Profile_Circular_Velocity_NFW(thisNode,scaleRadius)

    ! Ensure that the interpolations exist and extend sufficiently far.
    call Dark_Matter_Profile_NFW_Inverse_Angular_Momentum(specificAngularMomentumScaleFree)

    ! Interpolate to get the dimensionless radius at which this specific angular momentum is found.
    !$omp critical(NFW_Inverse_Interpolation)
    Radius_from_Specific_Angular_Momentum_NFW=Interpolate(nfwInverseTableNumberPoints,nfwSpecificAngularMomentum,nfwRadius&
         &,interpolationInverseObject,interpolationInverseAccelerator,specificAngularMomentumScaleFree,reset&
         &=interpolationInverseReset)
    !$omp end critical(NFW_Inverse_Interpolation)

    ! Convert to a physical radius.
    Radius_from_Specific_Angular_Momentum_NFW=Radius_from_Specific_Angular_Momentum_NFW*scaleRadius

    return
  end function Radius_from_Specific_Angular_Momentum_NFW
  
  double precision function Dark_Matter_Profile_Rotation_Normalization_NFW(thisNode)
    !% Return the normalization of the rotation velocity vs. specific angular momentum relation.
    use Tree_Nodes
    use Tree_Node_Methods
    use Dark_Matter_Halo_Scales
    use Numerical_Interpolation
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision                         :: concentration

    ! Find the concentration parameter of this halo.
    concentration=Dark_Matter_Halo_Virial_Radius(thisNode)/Tree_Node_Dark_Matter_Profile_Scale(thisNode)

    ! Ensure that the interpolations exist and extend sufficiently far.
    call Dark_Matter_Profile_NFW_Tabulate(concentration)

    ! Find the energy by interpolation.
    !$omp critical(NFW_Interpolation)
    Dark_Matter_Profile_Rotation_Normalization_NFW=Interpolate(nfwTableNumberPoints,nfwConcentration,nfwRotationNormalization&
         &,interpolationObject,interpolationAccelerator,concentration,reset=interpolationReset)&
         &/Dark_Matter_Halo_Virial_Radius(thisNode)
    !$omp end critical(NFW_Interpolation)
    return
  end function Dark_Matter_Profile_Rotation_Normalization_NFW
  
  double precision function Dark_Matter_Profile_Energy_NFW(thisNode)
    !% Return the energy of an NFW halo density profile.
    use Tree_Nodes
    use Tree_Node_Methods
    use Dark_Matter_Halo_Scales
    use Numerical_Interpolation
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision                         :: concentration

    ! Find the concentration parameter of this halo.
    concentration=Dark_Matter_Halo_Virial_Radius(thisNode)/Tree_Node_Dark_Matter_Profile_Scale(thisNode)

    ! Ensure that the interpolations exist and extend sufficiently far.
    call Dark_Matter_Profile_NFW_Tabulate(concentration)

    ! Find the energy by interpolation.
    !$omp critical(NFW_Interpolation)
    Dark_Matter_Profile_Energy_NFW=Interpolate(nfwTableNumberPoints,nfwConcentration,nfwEnergy,interpolationObject&
         &,interpolationAccelerator,concentration,reset=interpolationReset)*Tree_Node_Mass(thisNode)&
         &*Dark_Matter_Halo_Virial_Velocity(thisNode)**2
    !$omp end critical(NFW_Interpolation)
    return
  end function Dark_Matter_Profile_Energy_NFW
  
  double precision function Dark_Matter_Profile_Energy_Growth_Rate_NFW(thisNode)
    !% Return the rate of change of the energy of an NFW halo density profile.
    use Tree_Nodes
    use Tree_Node_Methods
    use Dark_Matter_Halo_Scales
    use Numerical_Interpolation
    implicit none
    type(treeNode),   intent(inout), pointer :: thisNode
    double precision                         :: concentration,energy,energyGradient

    ! Find the concentration parameter of this halo.
    concentration=Dark_Matter_Halo_Virial_Radius(thisNode)/Tree_Node_Dark_Matter_Profile_Scale(thisNode)

    ! Ensure that the interpolations exist and extend sufficiently far.
    call Dark_Matter_Profile_NFW_Tabulate(concentration)

    ! Find the energy gradient by interpolation.
    !$omp critical(NFW_Interpolation)
    energy=Interpolate                   (nfwTableNumberPoints,nfwConcentration,nfwEnergy,interpolationObject&
         &,interpolationAccelerator,concentration,reset=interpolationReset)
    energyGradient=Interpolate_Derivative(nfwTableNumberPoints,nfwConcentration,nfwEnergy,interpolationObject&
         &,interpolationAccelerator,concentration,reset=interpolationReset)
    !$omp end critical(NFW_Interpolation)

    Dark_Matter_Profile_Energy_Growth_Rate_NFW=Dark_Matter_Profile_Energy_NFW(thisNode)&
         &*(Tree_Node_Mass_Accretion_Rate(thisNode)/Tree_Node_Mass(thisNode)+2.0d0 &
         &*Dark_Matter_Halo_Virial_Velocity_Growth_Rate(thisNode)/Dark_Matter_Halo_Virial_Velocity(thisNode)+(energyGradient&
         &*concentration/energy)*(Dark_Matter_Halo_Virial_Radius_Growth_Rate(thisNode)/Dark_Matter_Halo_Virial_Radius(thisNode)&
         &-Tree_Node_Dark_Matter_Profile_Scale_Growth_Rate(thisNode)/Tree_Node_Dark_Matter_Profile_Scale(thisNode)))

    return
  end function Dark_Matter_Profile_Energy_Growth_Rate_NFW

  double precision function Angular_Momentum_NFW_Scale_Free(concentration)
    !% Returns the total angular momentum (in units of the virial mass times scale radius times [assumed constant] rotation speed)
    !% in an NFW dark matter profile with given {\tt concentration}. This is given by:
    !% \begin{equation}
    !% J = \left. \int_0^c 4 \pi x^3 \rho(x) \d x \right/ \int_0^c 4 \pi x^2 \rho(x) \d x,
    !% \end{equation}
    !% where $x$ is radius in units of the scale radius and $c$ is concentration. This can be evaluated to give
    !% \begin{equation}
    !% J = \left. \left[ 1 + c - 2 \ln (1+c) - {1 \over 1+c} \right] \right/ \left[ \ln(1+c)-{c\over 1+c} \right].
    !% \end{equation}
    use Numerical_Constants_Math
    implicit none
    double precision, intent(in) :: concentration

    Angular_Momentum_NFW_Scale_Free=(1.0d0+concentration-2.0d0*dlog(1.0d0+concentration)-1.0d0/(1.0d0+concentration)) &
         &/(dlog(1.0d0+concentration)-concentration/(1.0d0+concentration))
    return
  end function Angular_Momentum_NFW_Scale_Free

  double precision function Specific_Angular_Momentum_NFW_Scale_Free(radius)
    !% Returns the specific angular momentum, normalized to unit scale length and unit velocity at the scale radius, at position
    !% {\tt radius} (in units of the scale radius) in an NFW profile.
    implicit none
    double precision, intent(in) :: radius

    Specific_Angular_Momentum_NFW_Scale_Free=dsqrt(radius*Enclosed_Mass_NFW_Scale_Free(radius,1.0d0))
    return
  end function Specific_Angular_Momentum_NFW_Scale_Free

  double precision function Enclosed_Mass_NFW_Scale_Free(radius,concentration)
    !% Returns the enclosed mass (in units of the virial mass) in an NFW dark matter profile with given {\tt concentration} at the
    !% given {\tt radius} (given in units of the scale radius).
    implicit none
    double precision, intent(in) :: radius,concentration
    double precision, parameter  :: minimumRadiusForExactSolution=1.0d-7

    if (radius >= minimumRadiusForExactSolution) then
       Enclosed_Mass_NFW_Scale_Free=(dlog(1.0d0+radius)-radius/(1.0d0+radius))
    else
       Enclosed_Mass_NFW_Scale_Free=(radius**2)*(0.5d0+radius*(-2.0d0/3.0d0+radius*(0.75d0+radius*(-0.8d0))))
    end if
    Enclosed_Mass_NFW_Scale_Free=Enclosed_Mass_NFW_Scale_Free/(dlog(1.0d0+concentration)-concentration/(1.0d0 &
         &+concentration))
    return
  end function Enclosed_Mass_NFW_Scale_Free
  
  double precision function Density_NFW_Scale_Free(radius,concentration)
    !% Returns the density (in units such that the virial mass and scale length are unity) in an NFW dark matter profile with
    !% given {\tt concentration} at the given {\tt radius} (given in units of the scale radius).
    use Numerical_Constants_Math
    implicit none
    double precision, intent(in) :: radius,concentration

    Density_NFW_Scale_Free=1.0d0/(dlog(1.0d0+concentration)-concentration/(1.0d0+concentration))/radius/(1.0d0+radius)**2/4.0d0/Pi
    return
  end function Density_NFW_Scale_Free
  
  double precision function NFW_Profile_Energy(concentration)
    !% Computes the total energy of an NFW profile halo of given {\tt concentration} using the methods of
    !% \citeauthor{cole_hierarchical_2000}~(\citeyear{cole_hierarchical_2000}; their Appendix~A).
    use, intrinsic :: ISO_C_Binding
    use Numerical_Constants_Math
    use Numerical_Integration
    implicit none
    double precision,                intent(in) :: concentration
    type(c_ptr)                                 :: parameterPointer
    type(fgsl_function)                         :: integrandFunction
    type(fgsl_integration_workspace)            :: integrationWorkspace
    double precision                            :: radiusMinimum,radiusMaximum,potentialEnergyIntegral,potentialEnergy&
         &,kineticEnergyIntegral,kineticEnergy,jeansEquationIntegral

    ! Compute the potential energy.
    radiusMinimum=0.0d0
    radiusMaximum=concentration
    concentrationParameter=concentration
    potentialEnergyIntegral=Integrate(radiusMinimum,radiusMaximum,Potential_Energy_Integrand&
         &,parameterPointer,integrandFunction,integrationWorkspace,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-3)
    call Integrate_Done(integrandFunction,integrationWorkspace)
    potentialEnergy=-0.5d0*(1.0d0/concentration+potentialEnergyIntegral)

    ! Compute the velocity dispersion at the virial radius.
    radiusMinimum=concentration
    radiusMaximum=100.0d0*concentration
    concentrationParameter=concentration
    jeansEquationIntegral=Integrate(radiusMinimum,radiusMaximum,Jeans_Equation_Integrand&
         &,parameterPointer,integrandFunction,integrationWorkspace,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-3)
    call Integrate_Done(integrandFunction,integrationWorkspace)

    ! Compute the kinetic energy.
    radiusMinimum=0.0d0
    radiusMaximum=concentration
    concentrationParameter=concentration
    kineticEnergyIntegral=Integrate(radiusMinimum,radiusMaximum,Kinetic_Energy_Integrand&
         &,parameterPointer,integrandFunction,integrationWorkspace,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-3)
    call Integrate_Done(integrandFunction,integrationWorkspace)
    kineticEnergy=2.0d0*Pi*(jeansEquationIntegral*concentration**3+kineticEnergyIntegral)
    
    ! Compute the total energy.
    NFW_Profile_Energy=(potentialEnergy+kineticEnergy)*concentration

    return
  end function NFW_Profile_Energy
  
  function Potential_Energy_Integrand(radius,parameterPointer) bind(c)
    !% Integrand for NFW profile potential energy.
    use, intrinsic :: ISO_C_Binding
    implicit none
    real(c_double)          :: Potential_Energy_Integrand
    real(c_double), value   :: radius
    type(c_ptr),    value   :: parameterPointer
    
    Potential_Energy_Integrand=(Enclosed_Mass_NFW_Scale_Free(radius,concentrationParameter)/radius)**2
    return
  end function Potential_Energy_Integrand
  
  function Kinetic_Energy_Integrand(radius,parameterPointer) bind(c)
    !% Integrand for NFW profile kinetic energy.
    use, intrinsic :: ISO_C_Binding
    implicit none
    real(c_double)          :: Kinetic_Energy_Integrand
    real(c_double), value   :: radius
    type(c_ptr),    value   :: parameterPointer

    Kinetic_Energy_Integrand=Enclosed_Mass_NFW_Scale_Free(radius,concentrationParameter)*Density_NFW_Scale_Free(radius&
         &,concentrationParameter)*radius
    return
  end function Kinetic_Energy_Integrand

  function Jeans_Equation_Integrand(radius,parameterPointer) bind(c)
    !% Integrand for NFW profile Jeans equation.
    use, intrinsic :: ISO_C_Binding
    implicit none
    real(c_double)          :: Jeans_Equation_Integrand
    real(c_double), value   :: radius
    type(c_ptr),    value   :: parameterPointer
    
    Jeans_Equation_Integrand=Enclosed_Mass_NFW_Scale_Free(radius,concentrationParameter)*Density_NFW_Scale_Free(radius &
         &,concentrationParameter)/radius**2
    return
  end function Jeans_Equation_Integrand

  !# <galacticusStateStoreTask>
  !#  <unitName>Dark_Matter_Profiles_NFW_State_Store</unitName>
  !# </galacticusStateStoreTask>
  subroutine Dark_Matter_Profiles_NFW_State_Store(stateFile,fgslStateFile)
    !% Write the tablulation state to file.
    implicit none
    integer,         intent(in) :: stateFile
    type(fgsl_file), intent(in) :: fgslStateFile

    write (stateFile) concentrationMinimum,concentrationMaximum,radiusMinimum,radiusMaximum
    return
  end subroutine Dark_Matter_Profiles_NFW_State_Store
  
  !# <galacticusStateRetrieveTask>
  !#  <unitName>Dark_Matter_Profiles_NFW_State_Retrieve</unitName>
  !# </galacticusStateRetrieveTask>
  subroutine Dark_Matter_Profiles_NFW_State_Retrieve(stateFile,fgslStateFile)
    !% Retrieve the tabulation state from the file.
    implicit none
    integer,         intent(in) :: stateFile
    type(fgsl_file), intent(in) :: fgslStateFile
    double precision            :: tCosmological

    ! Read the minimum and maximum tabulated times.
    read (stateFile) concentrationMinimum,concentrationMaximum,radiusMinimum,radiusMaximum
    ! Retabulate.
    nfwTableInitialized       =.false.
    nfwInverseTableInitialized=.false.
    call Dark_Matter_Profile_NFW_Tabulate
    call Dark_Matter_Profile_NFW_Inverse_Angular_Momentum
    return
  end subroutine Dark_Matter_Profiles_NFW_State_Retrieve
  
end module Dark_Matter_Profiles_NFW
